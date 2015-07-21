#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 10.10.2015
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com

import sys
import argparse
import re
from collections import defaultdict

def get_revcomp(sequence):
    '''Return reverse complementary sequence.

    >>> complementary('AT CG')
    'CGAT'

    '''
    c = dict(zip('ATCGNatcgn[]', 'TAGCNtagcn]['))
    return ''.join(c.get(nucleotide, '') for nucleotide in reversed(sequence))


def sc_iter_fasta_brute(file_name):
    """ Iter over fasta file."""
    header = None
    seq = []
    with open(file_name) as fh:
        data = fh.readlines()
        for line in data:
            if line.startswith(">"):
                if seq:
                    sequence = "".join(seq)
                    yield (header, sequence)
                header = line
                seq = []
                continue
            seq.append(line)
        if seq or header:
            sequence = "".join(seq)
            yield (header, sequence)


if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='Fasta file to kmers file.')
	parser.add_argument('-i','--input', help='Fasta input', required=True)
	parser.add_argument('-o','--output', help='Kmer output', required=True)
	parser.add_argument('-k','--k', help='K lengths', required=True)
	args = vars(parser.parse_args())
	fasta_file = args["input"]
	output = args["output"]
	k = int(args["k"])

	kmers = defaultdict(int)
	for header, sequence in sc_iter_fasta_brute(fasta_file):
		sequence = sequence.upper()
		sequence = re.sub("\s+", "", sequence, re.S|re.M)
		for i in xrange(0, len(sequence)-k+1):
			kmer = sequence[i:i+k]
			kmers[kmer] += 1
			rkmer = get_revcomp(kmer)
			kmers[rkmer] += 1
	with open(output, "w") as fh:
		for kmer in kmers:
			fh.write("%s\t%s\n" % (kmer, kmers[kmer]))


