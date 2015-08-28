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

#print args.pop_one, args.pop_two

def count_alleles(pop, vcf_record):
    pop_REF_num = 0
    pop_ALT_num = 0
    for sample_name in pop:
    	(a1, a2) = vcf_record.genotype(sample_name)['GT'].split('/')
    	if a1 == '1':
    	    pop_ALT_num += 1
        else:
            pop_REF_num += 1
        if a2 == '0':
    	    pop_ALT_num += 1
        else:
            pop_REF_num += 1
    return (pop_REF_num, pop_ALT_num) 

out_file = open(args.out_file, 'w')
header = '%s %s Allele1 %s %s Allele2 %s %s Position\n' % \
         (args.ref_name, args.anc_name, args.pop_one_name,\
         args.pop_two_name, args.pop_one_name, args.pop_two_name)
out_file.write(header)

for record in in_vcf:
    if len(record.ALT) > 1 or record.INFO['NS'] != 7:
    	continue
    start = record.POS - 2
    stop  = record.POS + 1
    ref_trinucl = ref_scaffolds[record.CHROM][start:stop]
    anc_trinucl = str(ref_trinucl[0]) + record.INFO['AA'] + str(ref_trinucl[-1])
    pop_one_ref = (count_alleles(args.pop_one, record))[0]
    pop_one_alt = (count_alleles(args.pop_one, record))[1]
    pop_two_ref = (count_alleles(args.pop_two, record))[0]
    pop_two_alt = (count_alleles(args.pop_two, record))[1]
    out_str     = (str(ref_trinucl), anc_trinucl, record.REF, str(pop_one_ref),\
                   str(pop_two_ref), str(record.ALT[0]), str(pop_one_alt),\
                   str(pop_two_alt), str(record.POS), '\n')
    out_file.write(' '.join(out_str))

out_file.close()