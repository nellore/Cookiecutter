# Cookiecutter: a tool for kmer-based read filtering and extraction

Cookiecutter is a computational tool for filtering reads by their 
matches to specified k-mers. Originally it was created to
filter reads with primers and adapters from large sets of sequencing 
data.

## Requirements

To compile and use Cookiecutter, the following tools must be 
installed.

- make;
- gcc 4.7 or higher;
- python 2.7.

Cookiecutter is designed for use on Linux/UNIX and OS X systems.

## Installation

The package should be compiled from its source code using the 
provided Makefile in the following way.

```
git clone http://github.com/ad3002/Cookiecutter.git
cd Cookiecutter/src
make
sudo make install
```

If you do not have root access, you can use Cookiecutter from the `src` 
directory or specify another installation directory using `PREFIX`:

```
PREFIX=/my/dir make install
```

To uninstall Cookiecutter, use `make uninstall`. If an installation 
directory was specified using `PREFIX`, then it should also be 
specified 
for uninstalling: `PREFIX=/my/dir make uninstall`.

## Usage

Cookiecutter contains a number of subroutines for various tasks:

- **remove** searches given k-mers in reads and outputs the reads 
without any matches to the k-mers;

- **rm_reads** is an extension of **remove** that additionally provides 
options to filter reads by the presence of (C)n/(G)n tracks or 
unknown nucleotides, read length or low sequence complexity and 
outputs both filtered and unfiltered reads;

- **extract** searches given k-mers in reads and outputs the reads that
 matched the k-mers;
 
- **separate** searches given k-mers in reads and outputs both matched 
and unmatched reads to two separate files.

These subroutines may be launched directly or using the wrapper script 
**cookiecutter**. We recommend to use the wrapper because it allows to
process multiple input files in parallel mode and provides a 
convenient command-line interface to the subroutines. Also one may 
create k-mer libraries from FASTA files using the **cookiecutter 
make_library** tool.

Below we give examples of Cookiecutter usage. To get more information 
about the program options, use the `-h` argument: `cookiecutter -h`. 
It can also be applied to a specific subroutine, for example 
`cookiecutter rm_reads -h`.

### Creating a library of k-mers

A library of k-mers is necessary for all Cookiecutter subroutines. It 
can be created from a FASTA file using `cookiecutter make_library`.
For example, the command

```
cookiecutter make_library -i adapters.fa -o adapters.txt -l 5
```

will create the file *adapters.txt* of k-mers of length 5 bp from the
FASTA file *adapters.fa*.

### Removing reads by k-mers

Let us have a library of k-mers *adapters.txt* created as described 
above and a FASTQ file of single-end reads *raw_data.fastq*, and we 
would like to remove all reads containing any k-mers from the library.
It can be done using **remove** in the following way.

```
cookiecutter remove -i raw_data.fastq -f adapters.txt -o filtered
```

The output FASTQ file *raw_data.ok.fastq* will be created in the 
directory specified by the `-o` argument. It will contain the reads 
that do not include any of the specified matches.

### Extracting reads by k-mers

Let us have the same data set as in the
subsection above, but now we are to extract the reads matching 
any of the specified k-mers. For that, one should use **extract** in 
the same way as **remove**:

```
cookiecutter remove -i raw_data.fastq -f adapters.txt -o filtered
```

### Advanced read filtration

Let us have two FASTQ files of paired-end reads *raw_data_1.fastq* 
and *raw_data_2.fastq*. In addition to the k-mer presence filter, we 
would also like to filter them by the following criteria: read length,
presence of (G)n or (C)n tracks, sequence complexity (DUST) and
unknown nucleotides within a read. The **rm_reads** tool was designed
for such filtration.
 
```
cookiecutter rm_reads -1 raw_data_1.fastq -2 raw_data_2.fastq
    -f adapters.txt -o output_dir --polygc 13 --length 50
    --dust --filterN
```

Since we specified a pair of FASTQ files, the output files will also 
be paired. Read pairs are maintained if both paired-end read parts 
passed the filtration. If one part of a read passed the filtration but 
another 
failed it, then the passed part will be output to the file which 
name ends with *.se.fastq*.

### Read separation

Let us have the same paired-end FASTQ files *raw_data_1.fastq* and
*raw_data_2.fastq* as in the subsection above. We would like to 
separate reads matching the k-mer library from reads that do not 
match it. We will use the **separate** tool.

```
cookiecutter separate -1 raw_data_1.fastq -2 raw_data_2.fastq
    -f adapters.txt -o output_dir
```

### Processing multiple input files

Cookiecutter supports processing multiple input files (or pairs 
of input FASTQ files for paired-end reads) in parallel mode. For 
that, one should specify multiple input files in the arguments `-1`, 
`-2` or `-i` (see examples below).

```
cookiecutter remove -1 reads_a_1.fastq reads_b_1.fastq
    -2 reads_a_2.fastq reads_b_2.fastq -f adapters.txt
    -o output_dir
```

```
cookiecutter extract -i reads_a.fastq reads_b.fastq
    -f adapters.txt -o output_dir
```

Also one may specify multiple input FASTA files for the k-mer library
making tool.

```
cookiecutter make_library -i input_1.fa input_2.fa -o library.txt -l 5
```

## Repository and feedback

The latest version of Cookiecutter is publicly available at 
its [GitHub repository](https://github.com/ad3002/Cookiecutter). If 
you find any bugs or have any suggestions how to improve the tool, 
please find free to post issues at the repository. The earliest 
version of Cookiecutter can also be found at
[GitHub](https://github.com/allivi/rm_reads).
