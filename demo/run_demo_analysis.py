#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# created: 05.06.2015
# author: Aleksey Komissarov
# contact: ad3002@gmail.com

import argparse
import os
import subprocess
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

    if not os.path.isdir("../src"):
        print "Please, run demo script from demo directory."
        sys.exit(2)

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

    try:
        commands = [
            ["make", "clean"],
            ["make"],
        ]
        os.chdir("../src")
        for command in commands:
            print command
            subprocess.check_call(command, shell=True)
    except subprocess.CalledProcessError:
        sys.stderr.write('failed to compile\n')

    command_names = dict(a='removing technical sequences',
                         b='removing technical sequences and '
                           'applying DUST filter',
                         c='rRNA removing from transcriptome data',
                         d='extracting mtDNA',
                         e='removing alpha satDNA')

    command_launches = dict(
        a='../src/remove -1 %(fastq1)s -2 %(fastq2)s -o %('
          'output_dir_1a)s --fragments ../data/illumina.dat',
        b='../src/rm_reads -1 %(fastq1)s -2 %(fastq2)s -o %('
          'output_dir_1b)s --polyG 13 --length 50 --fragments '
          '../data/illumina.dat --dust_cutoff 3 --dust_k 4',
        c='../src/remove -i %(transc_fastq)s -o %(output_dir_1c)s '
          '--fragments ../data/rdna.dat',
        d='../src/extractor -1 %(fastq1)s -2 %(fastq2)s -o %('
          'output_dir_1d)s --fragments ../data/mtdna.dat',
        e='../src/separate -1 %(fastq1)s -2 %(fastq2)s -o %('
          'output_dir_1e)s --fragments ../data/alpha.dat'
    )

    for label in sorted(command_launches.iterkeys()):
        command = command_launches[label]
        print 'Running analysis 1{} ({})'.format(label,
                                                 command_names[label])
        print 'Command 1{}: {}'.format(label, command)
        try:
            subprocess.check_call(command % data, shell=True)
            print 'Command 1{}: completed!'.format(label)
        except subprocess.CalledProcessError:
            print 'Command 1{}: failed!'.format(label)
