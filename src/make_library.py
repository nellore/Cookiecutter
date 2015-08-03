#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# created: 10.10.2015
# author: Aleksey Komissarov
# contact: ad3002@gmail.com

import argparse
import re
from collections import defaultdict


def get_revcomp(seq):
    """Return reverse complementary sequence.

    >>> get_revcomp('AT CG')
    'CGAT'

    """
    c = dict(zip('ATCGNatcgn[]', 'TAGCNtagcn]['))
    return ''.join(c.get(nucleotide, '') for nucleotide in reversed(
        seq))


def sc_iter_fasta_brute(file_name):
    """ Iter over fasta file."""
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


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Create a filf of kmers from a FASTA file')
    parser.add_argument('-i', '--input', help='FASTA input',
                        required=True)
    parser.add_argument('-o', '--output', help='kmer output',
                        required=True)
    parser.add_argument('-k', '--k', help='length of kmers',
                        required=True)
    args = vars(parser.parse_args())
    fasta_file = args["input"]
    output = args["output"]
    k = int(args["k"])

    kmers = defaultdict(int)
    for header, sequence in sc_iter_fasta_brute(fasta_file):
        sequence = sequence.upper()
        sequence = re.sub("\s+", "", sequence, re.S | re.M)
        for i in xrange(0, len(sequence)-k+1):
            kmer = sequence[i:i+k]
            kmers[kmer] += 1
            rkmer = get_revcomp(kmer)
            kmers[rkmer] += 1
    with open(output, "w") as fh:
        for kmer in kmers:
            fh.write("%s\t%s\n" % (kmer, kmers[kmer]))
