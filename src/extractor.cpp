#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <list>
#include <map>
#include <algorithm>
#include <locale>
#include <getopt.h>
#include <stdlib.h>
#include <unordered_map>

#include "fileroutines.h"
#include "search.h"
#include "stats.h"
#include "seq.h"

/*! \brief Read adapter patterns from an input stream
 *
 *  \param[in]  kmers_f     an input stream to read the patterns from
 *  \param[out] patterns    the vector to which the patterns are written
 */
void build_patterns(std::ifstream & kmers_f, std::vector <std::pair <std::string, Node::Type> > & patterns)
{
    std::string tmp;
    while (!kmers_f.eof()) {
        std::getline(kmers_f, tmp);
        std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
        if (!tmp.empty()) {
            size_t tab = tmp.find('\t');
            if (tab == std::string::npos) {
                patterns.push_back(std::make_pair(tmp, Node::Type::adapter));
            } else {
                patterns.push_back(std::make_pair(tmp.substr(0, tab), Node::Type::adapter));
            }
        }
    }
    kmers_f.close();
}

/*! \brief Check a read against patterns
 *
 *  \param[in]  read        a read sequence
 *  \param[in]  root        a root of the trie structure used for string matching
 *  \param[in]  patterns    a vector of patterns
 *  \param[in]  errors      the number of resolved mismatches between a read and a pattern
 *  \return                 the read type
 */
ReadType check_read(std::string const & read, Node * root, std::vector <std::pair<std::string, Node::Type> > const & patterns, int errors)
{
    if (errors) {
        return (ReadType)search_inexact(read, root, patterns, errors);
    } else {
        return (ReadType)search_any(read, root);
    }
}

/*! \brief Filter single-end reads by patterns
 *
 *  \param[in]  reads_f     an input stream of read sequences
 *  \param[out] bad_f       an output stream to write filtered out reads to
 *  \param[out] stats       statistics on processed reads
 *  \param[out] root        a root of the trie structure used to perform string matching
 *  \param[in]  patterns    a vector of patterns for read filtration
 *  \param[in]  errors      the number of resolved mismatches between a read and a pattern
 */
void filter_single_reads(std::ifstream & reads_f, std::ofstream & bad_f, 
                         Stats & stats, Node * root, std::vector <std::pair<std::string, Node::Type> > const & patterns, int errors)
{
    Seq read;
    int processed = 0;

    while (read.read_seq(reads_f)) {
        ReadType type = check_read(read.seq, root, patterns, errors);
        stats.update(type);
        if (type != ReadType::ok) {
            read.write_seq(bad_f);
        }

        processed += 1;
        if (processed % 1000000 == 0) {
            std::cout << "Processed: " << processed << std::endl;
        }
    }
}

/*! \brief Filter paired-end reads by patterns
 *
 *  \param[in]  reads1_f    an input stream of paired-end read 1 sequences
 *  \param[in]  reads2_f    an input stream of paired-end read 2 sequences
 *  \param[out] bad1_f      an output stream to write filtered out paired-end read 1 sequences to
 *  \param[out] bad2_f      an output stream to write filtered out paired-end read 2 sequences to
 *  \param[out] se1_f       an output stream to write filtered first parts of paired-end reads
 *  \param[out] se2_f       an output stream to write filtered second parts of paired-end reads
 *  \param[out] stats1      statistics on first parts of processed reads
 *  \param[out] stats2      statistics on second parts of processed reads
 *  \param[in]  root        a root of the trie structure used to perform string matching
 *  \param[in]  patterns    a vector of patterns for read filtration
 *  \param[in]  errors      the number of resolved mismatches between a read and a pattern
 *
 *  \remark The streams \p se1_f (and \p se2_f) correspond to paired-end reads which second
 *  (or first) part was filtered but the other one was left.
 */
void filter_paired_reads(std::ifstream & reads1_f, std::ifstream & reads2_f,
                         std::ofstream & bad1_f, std::ofstream & bad2_f,
                         std::ofstream & se1_f, std::ofstream & se2_f,
                         Stats & stats1, Stats & stats2,
                         Node * root, std::vector <std::pair<std::string, Node::Type> > const & patterns,
                         int errors)
{
    Seq read1;
    Seq read2;
    int processed = 0;

    while (read1.read_seq(reads1_f) && read2.read_seq(reads2_f)) {
        ReadType type1 = check_read(read1.seq, root, patterns, errors);
        ReadType type2 = check_read(read2.seq, root, patterns, errors);
        if (type1 == ReadType::ok && type2 == ReadType::ok) {
            
        } else {
            stats1.update(type1, false);
            stats2.update(type2, false);
            if (type1 == ReadType::ok) {
                read1.write_seq(se1_f);
                read2.write_seq(bad2_f);
            } else if (type2 == ReadType::ok) { 
                read1.write_seq(bad1_f);
                read2.write_seq(se2_f);
            } else {
                read1.write_seq(bad1_f);
                read2.write_seq(bad2_f);
            }
        }

        processed += 1;
        if (processed % 1000000 == 0) {
            std::cout << "Processed: " << processed << std::endl;
        }
    }
}

/*! \brief Print program parameters */
void print_help() 
{
    std::cout << "./extractor [-i raw_data.fastq | -1 raw_data1.fastq -2 raw_data2.fastq] -o output_dir --fragments fragments.dat -e errors" << std::endl;
}

/*! \brief The main function of the **extract** tool. */
int main(int argc, char ** argv)
{
    Node root('0');
    std::vector <std::pair<std::string, Node::Type> > patterns;

    std::string kmers, reads, out_dir;
    std::string reads1, reads2;
    char rez = 0;
    int errors = 0;
    const struct option long_options[] = {
        {"fragments",required_argument,NULL,'a'},
        {NULL,0,NULL,0}
    };

    while ((rez = getopt_long(argc, argv, "1:2:l:p:a:i:o:e:", long_options, NULL)) != -1) {
        switch (rez) {
        case 'a':
            kmers = optarg;
            break;
        case 'i':
            reads = optarg;
            break;
        case '1':
            reads1 = optarg;
            break;
        case '2':
            reads2 = optarg;
            break;
        case 'o':
            out_dir = optarg;
            break;
        case 'e':
            errors = std::atoi(optarg);
            break;
        case '?':
            print_help();
            return -1;
        }
    }

    if (errors < 0 || errors > 2) {
        std::cout << "possible errors count are 0, 1, 2" << std::endl;
        return -1;
    }

    if (kmers.empty() || out_dir.empty() || (
            reads.empty() &&
            (reads1.empty() || reads2.empty()))) {
        print_help();
        return -1;
    }
    
    if (!verify_directory(out_dir)) {
        std::cout << "Output directory does not exist, failed to create" << std::endl;
    }

    std::ifstream kmers_f (kmers.c_str());
    if (!kmers_f.good()) {
        std::cout << "Cannot open kmers file" << std::endl;
        print_help();
        return -1;
    }

    std::cout << "Building patterns..." << std::endl;
    build_patterns(kmers_f, patterns);

    if (patterns.empty()) {
        std::cout << "patterns are empty" << std::endl;
        return -1;
    }

    std::cout << "Building trie..." << std::endl;
    build_trie(root, patterns, errors);
	add_failures(root);


    std::cout << "Iterate reads..." << std::endl;
    if (!reads.empty()) {
        std::string reads_base = basename(reads);
        std::ifstream reads_f (reads.c_str());
        std::ofstream bad_f((out_dir + "/" + reads_base + ".filtered.fastq").c_str(), std::ofstream::out);

        if (!reads_f.good()) {
            std::cout << "Cannot open reads file" << std::endl;
            print_help();
            return -1;
        }

        if (!bad_f.good()) {
            std::cout << "Cannot open output file" << std::endl;
            print_help();
            return -1;
        }

        Stats stats(reads);

        filter_single_reads(reads_f, bad_f, stats, &root, patterns, errors);

        std::cout << stats;

        bad_f.close();
        reads_f.close();
    } else {
        std::string reads1_base = basename(reads1);
        std::string reads2_base = basename(reads2);
        std::ifstream reads1_f(reads1.c_str());
        std::ifstream reads2_f(reads2.c_str());

        std::string file_name_se1 = out_dir + "/" + reads1_base + ".se.fastq";
        std::string file_name_se2 = out_dir + "/" + reads2_base + ".se.fastq";
        std::string file_name_bad1 = out_dir + "/" + reads1_base + ".filtered.fastq";
        std::string file_name_bad2 = out_dir + "/" + reads2_base + ".filtered.fastq";

        std::ofstream se1_f(file_name_se1.c_str(),
                            std::ofstream::out);
        std::ofstream se2_f(file_name_se2.c_str(),
                            std::ofstream::out);
        std::ofstream bad1_f(file_name_bad1.c_str(),
                            std::ofstream::out);
        std::ofstream bad2_f(file_name_bad2.c_str(),
                            std::ofstream::out);

        std::cout << "Created output files:" << std::endl;
        std::cout << "\t" << file_name_se1 << std::endl;
        std::cout << "\t" << file_name_se2 << std::endl;
        std::cout << "\t" << file_name_bad1 << std::endl;
        std::cout << "\t" << file_name_bad2 << std::endl;


        if (!reads1_f.good() || !reads2_f.good()) {
            std::cout << "reads file is bad" << std::endl;
            print_help();
            return -1;
        }

        if (!bad1_f.good() || !bad2_f.good() ||
                !se1_f.good() || !se2_f.good()) {
            std::cout << "out file is bad" << std::endl;
            print_help();
            return -1;
        }

        Stats stats1(reads1);
        Stats stats2(reads2);

        filter_paired_reads(reads1_f, reads2_f,
                            bad1_f, bad2_f, se1_f, se2_f,
                            stats1, stats2,
                            &root, patterns, errors);

        std::cout << stats1;
        std::cout << stats2;

        bad1_f.close();
        bad2_f.close();
        reads1_f.close();
        reads2_f.close();
    }
}
