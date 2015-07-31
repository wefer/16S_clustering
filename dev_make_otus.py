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
		#self.merged = self.merge_reads()	


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
		
		
	def merge_reads(self):
		"""Merge forward and reverse reads"""
		merged_file = os.path.splitext(self.fwd_adp_out)[0] + '_merged.fastq'
		cmd = ['usearch', 
				'--fastq_mergepairs', self.fwd_adp, 
				'-reverse', self.rev_adp, 
				'-fastqout', merged_file, 
				'-fastq_merge_maxee', '1.0']

		p = subprocess.Popen(cmd)
		p.wait()
		return merged_file
 		


class CombinedReads(object):

	def __init__(self):
		pass


	def combine_merged_reads(self):
		"""Get all the merged reads from the samples and concat them"""
		files = [x.merged for x in ReadPair.__all__]
		pass


	def dereplicate(self):
		pass


	def sort_reads(self):
		pass

	
	def cluster_otus(self):
		pass


