#!/usr/bin/env python

from dev_make_otus import ReadPair, CombinedReads
import sys
import glob

sample_folder = sys.argv[1]

def main(sample_folder):
	
	if sample_folder[-1] != '/':
		sample_folder += '/'

	forwards = glob.glob(sample_folder + '*R1*.fastq')
	paired = [(x, x.replace('_R1_', '_R2_')) for x in forwards]
	
	read_pairs = [ReadPair(*x) for x in paired]

main(sample_folder)
	
	
	


