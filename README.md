Description
----------------------

A tool for filtering reads by given kmers set. Initially created to filter reads with primers and adapters from large datasets.

Requirements
----------------------

- make
- gcc 4.7 and higher


Usage for adapter removing
--------------------------

./rm_reads <-i raw_data.fastq | -1 raw_data1.fastq -2 raw_data2.fastq> --fragments adapters.dat [-o output_dir --polyG 13 --length 50 --dust_cutoff cutoff --dust_k k -filterN]

    -i              input file
    -1              first input file for paired reads
    -2              second input file for paired reads
    -o              output directory (current directory by default)
    --polyG, -p     length of polyG/polyC tails (13 by default)
    --length, -l    minimum length cutoff (50 by default)
    --fragments  file with adapter kmers
    --dust_k, -k    window size for dust filter (not used by default)
    --dust_cutoff, -c   cutoff by dust score (not used by default)
    --errors, -e    maximum error count in match, possible values - 0, 1, 2 (0 by default)
    --filterN, -N   allow filter by N's in reads

Usage for reads removing by list with kmers
-------------------------------------------

./remove <-i raw_data.fastq | -1 raw_data1.fastq -2 raw_data2.fastq> --fragments fragments.dat -o output_dir

    -i              input file
    -1              first input file for paired reads
    -2              second input file for paired reads
    -o              output directory (current directory by default)
    --fragments     file with kmers

Usage for reads extraction by list with kmers
-------------------------------------------

./extractor <-i raw_data.fastq | -1 raw_data1.fastq -2 raw_data2.fastq> --fragments fragments.dat -o output_dir

    -i              input file
    -1              first input file for paired reads
    -2              second input file for paired reads
    -o              output directory (current directory by default)
    --fragments     file with kmers

Usage for reads separating into two groups by list with kmers
-------------------------------------------

./separate <-i raw_data.fastq | -1 raw_data1.fastq -2 raw_data2.fastq> --fragments fragments.dat -o output_dir

    -i              input file
    -1              first input file for paired reads
    -2              second input file for paired reads
    -o              output directory (current directory by default)
    --fragments     file with kmers


Input files
--------------------

Tools takes files with reads in fastq format as input. You can also use paired end reads. In case if you are using paired end reads, please, make sure that all reads from first file have correct pairs in second file.

```sh
run_batch.py -1 fastqA_1.fastq,fastqB_1.fastq -2 fastqA_2.fastq,fastqB_2.fastq -c remove --fragments fragments.dat -o filtered
```

Or with single end data:

```sh
run_batch_se.py -i fastqA.fastq,fastqB.fastq -c remove --fragments fragments.dat -o filtered
```


Creating  a kmer list from fasta file
-------------------------------------


```sh
make_library.py -i input_fasta.fa -o fragments.dat
```


Output files
--------------------

Tool creates following files in output directory:


1) file with correct reads (rm_reads, remove, separate)

```
input_prefix.ok.fastq       
```

2) file with reads, containing adapter kmers, N's, polyG/polyC tails or filtered by dust  (rm_reads, extractor separate)filter. Reason why read was filtered is given in the read id.

```
input_prefix.fitered.fastq  
```

3) for paired reads only. File with correct reads which have incorrect pair (rm_reads, separate).

```
input_prefix.se.fastq       
```


Project page
--------------------

Previous version of project can be found at https://github.com/allivi/rm_reads

Last version of project can be found at https://github.com/ad3002/Cookiecutter