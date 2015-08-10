#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# created: 05.06.2015
# author: Aleksey Komissarov
# contact: ad3002@gmail.com

import argparse
import os
import sys

if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Reproduce published analysis with small data '
                    'subset or with user defined fastq files.')
    parser.add_argument('-1', '--fastq1', help='Left fastq file',
                        required=False,
                        default="../demo/ERR194147_subset_1.fastq")
    parser.add_argument('-2', '--fastq2', help='Right fastq file',
                        required=False,
                        default="../demo/ERR194147_subset_2.fastq")
    args = vars(parser.parse_args())

    fastq_file1 = args["fastq1"]
    fastq_file2 = args["fastq2"]

    data = {
        "fastq1": fastq_file1,
        "fastq2": fastq_file2,
        "output_dir_1a": "../demo/temp_results_remove",
        "output_dir_1b": "../demo/temp_results_remove_dust",
        "output_dir_1c": "../demo/temp_results_rRNA",
        "output_dir_1d": "../demo/temp_results_mtDNA",
        "output_dir_1e": "../demo/temp_results_satDNA",
        "transc_fastq": "../demo/SRR100173_1.fastq",
    }

    command = "cookiecutter remove -1 %(fastq1)s -2 %(fastq2)s -o %(" \
              "output_dir_1a)s --fragments ../data/illumina.dat" % data
    print "Running analysis 1a (technical sequences)..."
    print command
    os.system(command)

    command = "cookiecutter rm_reads -1 %(fastq1)s -2 %(fastq2)s -o " \
              "%(output_dir_1b)s --polygc 13 --length 50 --fragments " \
              "../data/illumina.dat --dust_cutoff 3 --dust_k 4" % data
    print "Running analysis 1b (technical sequences and DUST filter)..."
    print command
    os.system(command)

    command = "cookiecutter remove -i %(transc_fastq)s -o " \
              "%(output_dir_1c)s --fragments ../data/rdna.dat" % data
    print "Running analysis 1c (rRNA removing from transcriptome " \
          "data)..."
    print command
    os.system(command)

    command = "cookiecutter extractor -1 %(fastq1)s -2 %(fastq2)s -o " \
              "%(output_dir_1d)s --fragments ../data/mtdna.dat" % data
    print "Running analysis 1d (mtDNA extracting)..."
    print command
    os.system(command)

    command = "cookiecutter separate -1 %(fastq1)s -2 %(fastq2)s -o " \
              "%(output_dir_1e)s --fragments ../data/alpha.dat" % data
    print "Running analysis 1e (alpha satDNA removing)..."
    print command
    os.system(command)
