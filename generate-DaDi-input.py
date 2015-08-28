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
    required=True)
parser.add_argument('-in_vcf',\
    help="Vcf file with multiple samples",\
    type=str,\
    required=True)
parser.add_argument('-pop_one',\
    help="List of samples from one population",\
    type=str,
    required=True,
    nargs='+')
parser.add_argument('-pop_two',\
    help="List of samples from second population",\
    type=str,
    required=True,
    nargs='+')
parser.add_argument('-out_file',\
    help="Name of output file",
    type=str,
    required=True)
parser.add_argument('-pop_one_name',\
    type=str,\
    default='population_one')
parser.add_argument('-pop_two_name',\
    type=str,\
    default='population_two')
parser.add_argument('-ref_name',\
    type=str,\
    default='REF')
parser.add_argument('-anc_name',\
    type=str,\
    default='ANC')

args = parser.parse_args()

# Read input files
ref_scaffolds = Fasta(args.ref)
in_vcf = vcf.Reader(open(args.in_vcf, 'r'))

def count_alleles(pop, vcf_record):
    ref = 0
    alt = 0
    for sample_name in pop:
    	(a1, a2) = vcf_record.genotype(sample_name)['GT'].split('/')
    	if a1 == '1':
    	    alt += 1
        else:
            ref += 1
        if a2 == '0':
    	    alt += 1
        else:
            ref += 1
    return (ref, alt) 

out_file = open(args.out_file, 'w')
header = '%s %s Allele1 %s %s Allele2 %s %s Position\n' % \
         (args.ref_name, args.anc_name, args.pop_one_name,\
          args.pop_two_name, args.pop_one_name,\
          args.pop_two_name)
out_file.write(header)

for record in in_vcf:
    if len(record.ALT) > 1 or record.INFO['NS'] != 7:
    	continue
    ref_trinucl      = ref_scaffolds[record.CHROM][record.POS - 2:record.POS + 1]
    anc_trinucl      = str(ref_trinucl[0]) + record.INFO['AA'] + str(ref_trinucl[-1])
    (p1_ref, p1_alt) = (count_alleles(args.pop_one, record))
    (p2_ref, p2_alt) = (count_alleles(args.pop_two, record))
    out_str          = (str(ref_trinucl), anc_trinucl, record.REF, str(p1_ref),\
                        str(p2_ref), str(record.ALT[0]), str(p1_alt),\
                        str(p2_alt), str(record.POS), '\n')
    out_file.write(' '.join(out_str))

out_file.close()