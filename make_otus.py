#!/usr/bin/env python

import sys
import magic
import subprocess
import os

####################################################

fwd_read = sys.argv[1]
rev_read = sys.argv[2]
fwd_n_bases = 17
rev_n_bases = 21

#Forward read adapter seqeunce
F_ADAPTER = 'ATCTCGTATGCCGTCTTCTGCTTG'
#Reverse read adapter sequence
R_ADAPTER = 'TCTCGGTGGTCGCCGTATCATT'

####################################################



def determine_file_type(file):
	"""Check that the files are not compressed"""
	return 'ASCII text' in magic.from_file(file) 


def trim(infile, n_bases):
    """Trim away primer sequences (n number of bases at the beginning of the sequence)"""

    outfile = os.path.splitext(infile)[0] + "_trim." + 'fastq'

    with open(outfile,'w') as o:
	    with open(infile, 'r') as f:
	        flag = 0            #flag next line as sequence
	        for line in f.readlines():
	            if flag == 0:
	                if "@HWI" in line or line == "+\n":
	                    o.write(line)
	                    flag = 1 
	            else:
	                o.write(line[n_bases:])
	                flag = 0 
	
    return outfile


def remove_adapters(fwd, rev):
	"""Remove reads with adapters"""

	fwd_out = os.path.splitext(fwd)[0] + '_adp.fastq'
	rev_out = os.path.splitext(rev)[0] + '_adp.fastq'

	cmd = ['cutadapt', '-a', F_ADAPTER, '-A', R_ADAPTER, '-o', fwd_out, '-p', rev_out, '--discard-trimmed', fwd, rev] 
	p = subprocess.Popen(cmd)
	p.wait()

	return (fwd_out, rev_out)


def merge_reads(fwd, rev):
	"""Merge forward and reverse reads"""

	merged_file = os.path.splitext(fwd)[0] + '_merged.fastq'
	cmd = ['usearch', '--fastq_mergepairs', fwd, '-reverse', rev, '-fastqout', merged_file, '-fastq_merge_maxee', '1.0']
	p = subprocess.Popen(cmd)
	p.wait()
	return merged_file


def rename_headers(file):
	"""Add sample information to the fastq-headers"""
	pass


def dereplicate(file):
	"""Perform dereplication for usearch"""
	derep_out = os.path.splitext(file)[0] + '_derep.fastq'
	
	cmd = ['usearch', '-derep_fulllength', file, '-fastqout', derep_out, '-sizeout', '-minseqlength', '64']
	p = subprocess.Popen(cmd)
	p.wait()

	return derep_out


def sort_reads(file):
	"""Perform sortbysize and remove singletons"""
	sorted_out = os.file.splitext(file)[0] + '_sort.fastq'

	cmd = ['usearch', '-sortbysize', file, '-fastqout', sorted_out, '-minsize', '2']
	p = subprocess.Popen(cmd)
	p.wait()
	
	return sorted_out


def cluster_otus(file, radius=3):
	"""Perfor clustering with usearch"""

	identity = str(100 - radius)
	otus = 'otus_' + identity + '.fasta'
	up_out = 'uparse_out.csv'

	cmd = ['usearch', '-cluster_otus', file, '-otus', otus, '-uparseout', up_out, '-relabel', 'OTU_', '-sizein', '-sizeout']
	p = subprocess.Popen(cmd)
	p.wait()
	
	return otus



def main():

	if determine_file_type(fwd_read) and determine_file_type(rev_read):
		#run all steps
		pass

	else:
		raise ValueError('Non-readable file. Make sure input is uncompressed.')


if __name__ == '__main__':
	main()


