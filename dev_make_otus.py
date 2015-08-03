#!/usr/bin/env python

import os
import magic
import subprocess

class ReadPair(object):

	__all__ = set()
	
	# Forward primer length
	fwd_n_bases = 17
	# Reverse primer lenght
	rev_n_bases = 21
	# Forward read adapter seqeunce
	F_ADAPTER = 'ATCTCGTATGCCGTCTTCTGCTTG'
	# Reverse read adapter sequence
	R_ADAPTER = 'TCTCGGTGGTCGCCGTATCATT'
	

	def __init__(self, fwd_read, rev_read):

		self.__class__.__all__.add(self)

		self.fwd_read = fwd_read
		self.rev_read = rev_read
		
		self.fwd_adp, self.rev_adp = self.remove_adapters()
		self.merged = self.merge_reads(self.fwd_adp, self.rev_adp)	


	def trim_primers(self):
		pass


	def remove_adapters(self):
	
		fwd_adp_out = os.path.splitext(self.fwd_read)[0] + '_adp.fastq'
		rev_adp_out = os.path.splitext(self.rev_read)[0] + '_adp.fastq'	
		cmd = ['cutadapt',
				'-a', ReadPair.F_ADAPTER,
				'-A', ReadPair.R_ADAPTER,
				'-o', fwd_adp_out,
				'-p', rev_adp_out,
				'--discard-trimmed',
				self.fwd_read, self.rev_read]
		p = subprocess.Popen(cmd)
		p.wait()

		return (fwd_adp_out, rev_adp_out)
		
		
	def merge_reads(self, fwd, rev):
		"""Merge forward and reverse reads"""
		merged_file = os.path.splitext(fwd)[0] + '_merged.fastq' #TODO fix filename
		cmd = ['usearch', 
				'--fastq_mergepairs', fwd, 
				'-reverse', rev, 
				'-fastqout', merged_file, 
				'-fastq_merge_maxee', '1.0']

		p = subprocess.Popen(cmd)
		p.wait()
		return merged_file



class CombinedReads(object):
	"""Class for setting up the otus"""

	def __init__(self):
		self.fasta = self.combine_merged_reads()
		self.derep = self.dereplicate(self.fasta)
		self.sorted = self.sort_reads(self.derep)
		self.otus = self.cluster_otus(self.sorted)

	def combine_merged_reads(self):
		"""Get all the merged reads from the samples and concat them"""
		
		files = [x.merged for x in ReadPair.__all__]

		combined_fasta = 'all_seqs.fasta'

		with open('all_seqs.fasta', 'w') as o:
			for fname in files:
				with open(fname) as f:
					for line in f:
						o.write(line)
		
		return combined_fasta


	def dereplicate(self, file):
		"""Perform dereplication for usearch"""

		derep_out = os.path.splitext(file)[0] + '_derep.fastq'

		cmd = ['usearch', '-derep_fulllength', file, '-fastqout', derep_out, '-sizeout', '-minseqlength', '64']
		p = subprocess.Popen(cmd)
		p.wait()

		return derep_out


	def sort_reads(self, file):
		"""Perform sortbysize and remove singletons"""
		sorted_out = os.path.splitext(file)[0] + '_sort.fastq'

		cmd = ['usearch', '-sortbysize', file, '-fastqout', sorted_out, '-minsize', '2']
		p = subprocess.Popen(cmd)
		p.wait()

		return sorted_out

	def cluster_otus(self, file, radius=3):
		"""Perfor clustering with usearch"""

		identity = str(100 - radius)
		otus = 'otus_' + identity + '.fasta'
		up_out = 'uparse_out.csv'

		cmd = ['usearch', '-cluster_otus', file, '-otus', otus, '-uparseout', up_out, '-relabel', 'OTU_', '-sizein', '-sizeout']
		p = subprocess.Popen(cmd)
		p.wait()

		return otus

