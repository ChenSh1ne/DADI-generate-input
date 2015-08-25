#! /usr/bin/python

import argparse
import os
import vcf
from pyfaidx import Fasta



parser = argparse.ArgumentParser()
#parser.add_argument('-N', help="Number of cores", type=int, default=6)
parser.add_argument('-ref',\
    help="Fasta file with reference genome", \
    type=str,\
    default='ref_genome.fa')
parser.add_argument('-out_ref',\
    help="Fasta file with outgroup reference genome",\
    type=str,\
    default='out_genome.fa')
parser.add_argument('-in_vcf',\
    help="Vcf file with multiple samples",\
    type=str,\
    default='snv.vcf')
parser.add_argument('-synt',\
    help="File with synteny fragments beetwen reference and outgroup genome",\
	type=str,\
	default='query_list.txt')

args = parser.parse_args()

ref_scaffolds = Fasta(args.ref)
out_ref_scaffolds = Fasta(args.out_ref)

in_vcf = vcf.Reader(open(args.in_vcf, 'r'))

for record in in_vcf:
    #print len(record.samples)
    print record.num_unknown