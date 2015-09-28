#!/usr/bin/env python

from make_otus import ReadPair, CombinedReads
import sys
import glob

#Avoid .pyc for now
sys.dont_write_bytecode = True

sample_folder = sys.argv[1]


def main(sample_folder):
	
	if sample_folder[-1] != '/':
		sample_folder += '/'

	forwards = glob.glob(sample_folder + '*_R1_*.*q')
	paired = [(x, x.replace('_R1_', '_R2_')) for x in forwards]

	read_pairs = [ReadPair(*x) for x in paired]
	combined = CombinedReads()
	
	for sample in read_pairs:
		sample.cleanup()

	combined.cleanup()
		
	

if __name__ == "__main__":
	main(sample_folder)
