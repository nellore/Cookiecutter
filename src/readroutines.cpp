#include "readroutines.h"
#include <cmath>
#include <unordered_map>

/*! \brief Given a read sequence, calculate its DUST score
 *
 *  \param[in]  read    a read sequence
 *  \param[in]  k       the DUST algorithm parameter
 *  \return             the DUST score
 *
 *  \remark For more information on the DUST score, please check the following
 *  paper:
 *  Morgulis, Aleksandr, E. Michael Gertz, Alejandro A. Schäffer, and Richa
 *  Agarwala. "A fast and symmetric DUST implementation to mask low-complexity
 *  DNA sequences." *Journal of Computational Biology* 13, no. 5 (2006): 1028-1040.
 */
double get_dust_score(std::string const & read, int k)
{
    std::unordered_map <int, int> counts;
    static std::unordered_map <char, int> hashes = {{'N', 1},
        {'A', 2},
        {'C', 3},
        {'G', 4},
        {'T', 5}};
    unsigned int hash = 0;
    unsigned int max_pow = pow(10, k - 1);
    for (auto it = read.begin(); it != read.end(); ++it) {
        char c = (*it > 96) ? (*it - 32) : *it;
        hash = hash * 10 + hashes[c];
        if (it - read.begin() >= k - 1) {
            ++counts[hash];
            hash = hash - (hash / max_pow) * max_pow;
        }
    }
    double score = 0;
    double total = 0;
    for (auto it = counts.begin(); it != counts.end(); ++it) {
        score += it->second * (it->second - 1) / 2;
        total += score;
    }
    //    std::cout << (total / (read.size() - k + 1)) << std::endl;
    return (total / (read.size() - k + 1));
}

double get_mean_quality(std::string const & qual) {
    int sum = 0;
    int val = 0;
    double size = 0;
    for (auto it = qual.begin(); it != qual.end(); ++it) {
        val = (int)*it;
        sum += val;
        size += 1;
    }
    return sum/(double)size;
}

/*! \brief Check a read against patterns
 *
 *  \param[in]  read        a read sequence
 *  \param[in]  root        a root of the trie structure used for string matching
 *  \param[in]  patterns    a vector of patterns
 *  \param[in]  length      the read length threshold
 *  \param[in]  dust_k      the DUST algorithm parameter
 *  \param[in]  dust_cutoff the DUST score threshold
 *  \param[in]  errors      the number of resolved mismatches between a read and
 *                          a pattern
 *  \return                 the read type
 */
ReadType check_read(std::string const & read, std::string const & qual, Node * root, std::vector <std::pair<std::string, Node::Type> > const & patterns,
                    unsigned int length, int dust_k, int dust_cutoff, int errors, int mean_quality)
{
    if (length && read.size() < length) {
        return ReadType::length;
    }

    if (mean_quality > 0 && get_mean_quality(qual) < mean_quality) {
        return ReadType::mean_quality;
    }

    if (dust_cutoff && get_dust_score(read, dust_k) > dust_cutoff) {
        return ReadType::dust;
    }
    
    if (errors) {
        return (ReadType)search_inexact(read, root, patterns, errors);
    } else {
        return (ReadType)search_any(read, root);
    }
}
