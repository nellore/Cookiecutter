#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian, Aleksey Komissarov
# gaik (dot) tamazian (at) gmail (dot) com

import argparse
import logging
import re
import subprocess
import time
from collections import defaultdict

logging.basicConfig()
logger = logging.getLogger(__name__)


class ParallelLauncher(object):
    """
    Run a program with the specified parameters in the parallel mode.
    """
    def __init__(self, program, input_files, args=None, threads=1):
        """
        Create a launcher object.

        :param program: a program to be launched
        :param input_files: a list of input files for parallel runs
        :param threads: the number of parallel threads
        :param args: a dictionary of arguments to be passed to the runs
        :type: str
        :type input_files: list
        :type threads: int
        :type args: dict
        """
        self.__program = program
        self.__input_files = input_files
        self.__args = args
        self.__threads = threads

        # form a list of commands from the specified parameters.
        self.__cmd_list = [self.__program + ' ' +
                           map(lambda x: ' '.join(map(str, x)),
                               files) +
                           map(lambda x: ' '.join(map(str, x)),
                               files.iteritems())
                           for files in self.__input_files]

    def launch(self):
        """
        Perform the launch.
        """
        running = []
        while self.__cmd_list:
            while len(running) > self.__threads:
                logger.info('checking %d launched processes',
                            len(running))
                for i, process in enumerate(running):
                    returncode = process.poll()
                    if returncode is not None:
                        if returncode == 0:
                            logger.info('process succeeded (%d '
                                        'processes remain)',
                                        len(self.__cmd_list))
                        else:
                            logger.info('process failed with the error '
                                        'code %d (%d processes remain)',
                                        process.returncode,
                                        len(self.__cmd_list))
                        running[i] = None
                        running = [x for x in running if x is not None]
                        break
                time.sleep(1)
            curr_cmd = self.__cmd_list.pop()
            running.append(subprocess.Popen(curr_cmd, shell=True))
            logger.info('new command launched: %s', curr_cmd)


class Extractor(ParallelLauncher):
    """
    Launch the extractor tool.
    """
    def __init__(self, files, fragments, output, threads):
        super(Extractor, self).__init__(
            'extractor', files,
            dict(zip(('-f', '-o'), (fragments, output))),
            threads
        )


class Remove(ParallelLauncher):
    """
    Launch the remove tool.
    """
    def __init__(self, files, fragments, output, threads):
        super(Remove, self).__init__(
            'remove', files,
            dict(zip(('-f', '-o'), (fragments, output))),
            threads
        )


class RmReads(ParallelLauncher):
    """
    Launch the rn_reads tool.
    """
    def __init__(self, files, fragments, output, polygc, length,
                 dust, dust_k, dust_cutoff, filter_n, threads):
        super(RmReads, self).__init__(
            'rm_reads', files,
            dict(zip('-f', '-o', '-p', '-l', '-d', '-k', '-c', '-N'),
                 (fragments, output, polygc, length, dust, dust_k,
                  dust_cutoff, filter_n)),
            threads
        )


class Separate(ParallelLauncher):
    """
    Launch the separate tool.
    """
    def __init__(self, files, fragments, output, threads):
        super(Separate, self).__init__(
            'separate', files,
            dict(zip(('-f', '-o'), (fragments, output))),
            threads
        )


def get_revcomp(seq):
    """
    Given a nucleotide sequence, return its reverse complement.

    :param seq: a nucleotide sequence
    :type seq: str
    :return: a reverse complement of the specified sequence
    :rtype: str
    """
    c = dict(zip('ATCGNatcgn[]', 'TAGCNtagcn]['))
    return ''.join(c.get(nucleotide, '') for nucleotide in reversed(
        seq))


def sc_iter_fasta_brute(file_name):
    """
    Iterate over a FASTA file.

    :param file_name: a name of a FASTA file
    :type file_name: str
    :return: a tuple of a sequence and its header
    :rtype: tuple
    """
    seq_header = None
    seq = []
    with open(file_name) as file_handler:
        data = file_handler.readlines()
        for line in data:
            if line.startswith(">"):
                if seq:
                    yield (seq_header, "".join(seq))
                seq_header = line
                seq = []
                continue
            seq.append(line)
        if seq or seq_header:
            yield (seq_header, "".join(seq))


def create_kmer_file(fasta_name, output_name, kmer_length):
    """
    Write the list of all kmers of the specified length from the
    specified input file to the given output file.

    :param fasta_name: a name of a FASTA file from which sequences
        k-mers will be obtained
    :param output_name: a name of an output file where k-mers will be
        written to
    :param kmer_length: the length of obtained k-mers
    :type fasta_name: str
    :type output_name: str
    :type kmer_length: int
    """
    kmers = defaultdict(int)
    for header, sequence in sc_iter_fasta_brute(fasta_name):
        sequence = sequence.upper()
        sequence = re.sub("\s+", "", sequence, re.S | re.M)
        for i in xrange(0, len(sequence)-kmer_length+1):
            kmer = sequence[i:i+kmer_length]
            kmers[kmer] += 1
            rkmer = get_revcomp(kmer)
            kmers[rkmer] += 1
    with open(output_name, "w") as fh:
        for kmer in kmers:
            fh.write("%s\t%s\n" % (kmer, kmers[kmer]))


def cookiecutter():
    """
    Wrapper around tools of the Cookiecutter package.
    """
    parser = argparse.ArgumentParser(
        description='Cookiecutter: a kmer-based read filtration and '
                    'extraction tool.')
    subparsers = parser.add_subparsers(dest='command')

    # Parser for the extractor tool

    extractor_parser = subparsers.add_parser(
        'extractor',
        description='Extracts reads according to a given list of '
                    'kmers and outputs only the reads that matched '
                    'the list.'
    )

    extractor_required = extractor_parser.add_argument_group(
        'required arguments')

    extractor_io = extractor_required.add_mutually_exclusive_group(
        required=True
    )
    extractor_io.add_argument('-i', '--input',
                              help='a FASTQ file of single-end reads')
    extractor_io.add_argument('-1', '--fastq1',
                              help='a FASTQ file of the first '
                                   'paired-end reads')

    extractor_parser.add_argument('-2', '--fastq2',
                                  help='a FASTQ file of the second '
                                       'paired-end reads')
    extractor_required.add_argument('-f', '--fragments', required=True,
                                    help='a file of k-mers')
    extractor_required.add_argument('-o', '--output', required=True,
                                    help='a directory for output files')

    # Parser for the remove tool.

    remove_parser = subparsers.add_parser(
        'remove',
        description='Removes reads according to a given list of kmers '
                    'and outputs only reads without any matches to '
                    'the provided kmer list.'
    )

    remove_required = remove_parser.add_argument_group(
        'required arguments')

    remove_io = remove_required.add_mutually_exclusive_group(
        required=True
    )
    remove_io.add_argument('-i', '--input',
                           help='a FASTQ file of single-end reads')
    remove_io.add_argument('-1', '--fastq1',
                           help='a FASTQ file of the first '
                                'paired-end reads')

    remove_parser.add_argument('-2', '--fastq2',
                               help='a FASTQ file of the second '
                                    'paired-end reads')
    remove_required.add_argument('-f', '--fragments', required=True,
                                 help='a file of k-mers')
    remove_required.add_argument('-o', '--output', required=True,
                                 help='a directory for output files')

    # Parser for the rm_reads tool.

    rm_reads_parser = subparsers.add_parser(
        'rm_reads',
        description='The rm_reads tool is an extended version of '
                    'remove '
                    'enhanced with the DUST filter, removing reads '
                    'containing (G)n- and (C)n-tracks and unknown '
                    'nucleotides and filtering reads by their length; '
                    'also its output includes both filtered and '
                    'unfiltered reads.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    rm_reads_required = rm_reads_parser.add_argument_group(
        'required arguments')

    rm_reads_io = rm_reads_required.add_mutually_exclusive_group(
        required=True
    )
    rm_reads_io.add_argument('-i', '--input',
                             default=argparse.SUPPRESS,
                             help='a FASTQ file of single-end reads')
    rm_reads_io.add_argument('-1', '--fastq1',
                             default=argparse.SUPPRESS,
                             help='a FASTQ file of the first '
                                  'paired-end reads')

    rm_reads_parser.add_argument('-2', '--fastq2',
                                 default=argparse.SUPPRESS,
                                 help='a FASTQ file of the second '
                                      'paired-end reads')
    rm_reads_parser.add_argument('-p', '--polygc', type=int,
                                 default=13,
                                 help='the polyG/polyC sequence '
                                      'length cutoff')
    rm_reads_parser.add_argument('-l', '--length', type=int,
                                 default=50,
                                 help='the read length cutoff')
    rm_reads_parser.add_argument('-d', '--dust', action='store_true',
                                 help='use the DUST filter')
    rm_reads_parser.add_argument('-c', '--dust_cutoff', type=int,
                                 default=2,
                                 help='the score cutoff for the DUST '
                                      'filter')
    rm_reads_parser.add_argument('-k', '--dust_k', type=int,
                                 default=4,
                                 help='the window size for the DUST '
                                      'filter')
    rm_reads_parser.add_argument('-N', '--filterN',
                                 action='store_true',
                                 help='filter reads by the presence '
                                      'of Ns')

    rm_reads_required.add_argument('-f', '--fragments', required=True,
                                   default=argparse.SUPPRESS,
                                   help='a file of k-mers')
    rm_reads_required.add_argument('-o', '--output', required=True,
                                   default=argparse.SUPPRESS,
                                   help='a directory for output files')

    # Parser for the separate tool.

    separate_parser = subparsers.add_parser(
        'separate',
        description='Outputs both matched and not matched '
                    'reads in separate files.'
    )

    separate_required = separate_parser.add_argument_group(
        'required arguments')

    separate_io = separate_required.add_mutually_exclusive_group(
        required=True
    )
    separate_io.add_argument('-i', '--input',
                             help='a FASTQ file of single-end reads')
    separate_io.add_argument('-1', '--fastq1',
                             help='a FASTQ file of the first '
                                  'paired-end reads')

    separate_parser.add_argument('-2', '--fastq2',
                                 help='a FASTQ file of the second '
                                      'paired-end reads')
    separate_required.add_argument('-f', '--fragments', required=True,
                                   help='a file of k-mers')
    separate_required.add_argument('-o', '--output', required=True,
                                   help='a directory for output files')

    # Parser for the make_library tool.

    make_library_parser = subparsers.add_parser(
        'make_library',
        description='Create a library of k-mers from the specified '
                    'FASTA file.'
    )

    make_library_required = make_library_parser.add_argument_group(
        'required_arguments')
    make_library_required.add_argument('-i', '--input',
                                       help='a FASTA file',
                                       required=True)
    make_library_required.add_argument('-o', '--output',
                                       help='an output file of '
                                            'k-mers',
                                       required=True)
    make_library_required.add_argument('-l', '--length', type=int,
                                       help='the length of generated '
                                            'k-mers',
                                       required=True)

    args = parser.parse_args()

    if args.command == 'make_library':
        create_kmer_file(args.input, args.output, args.length)


if __name__ == '__main__':
    cookiecutter()
