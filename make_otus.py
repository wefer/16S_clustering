#!/usr/bin/env python

import os
import subprocess
import parse_rdp

class ReadPair(object):
	"""Read pair forward and reverse"""

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

		self.sampleid = os.path.basename(fwd_read).split('_')[0]
		
		self.fwd_adp, self.rev_adp = self.remove_adapters()
		self.fwd_trim, self.rev_trim = self.trim_primers(self.fwd_adp, ReadPair.fwd_n_bases),  self.trim_primers(self.rev_adp, ReadPair.rev_n_bases)
		self.merged = self.merge_reads(self.fwd_trim, self.rev_trim)
		self.reheaded = self.add_barcode(self.merged, self.sampleid)


	def trim_primers(self, infile, n_bases):
		"""Trim away primer sequences (n number of bases at the beginning of the sequence)"""

		outfile = os.path.splitext(infile)[0] + "_trim." + 'fastq'

		with open(outfile,'w') as o:
			with open(infile, 'r') as f:
				flag = 0	#flag next line as sequence or quality string
				for line in f.readlines():
					if flag == 0:
						if "@HWI" in line or line == "+\n":
							o.write(line)
							flag = 1 
					else:
						o.write(line[n_bases:])
						flag = 0 
      
		return outfile


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
		merged_file = os.path.splitext(fwd)[0] + '_merged.fastq'
		cmd = ['usearch', 
				'--fastq_mergepairs', fwd, 
				'-reverse', rev, 
				'-fastqout', merged_file, 
				'-fastq_merge_maxee', '1.0']

		p = subprocess.Popen(cmd)
		p.wait()
		return merged_file

	
	def add_barcode(self, fastq, sampleid):
		"""
		Add sample ID in the fastq header
		[originalheader];barcodelabel=[sample_id]
		"""
		outfile = os.path.splitext(fastq)[0] + '_rehead.fastq' 

		with open(outfile, 'w') as o:
			with open(fastq, 'r') as f: 
				for line in f.readlines():
					if line[:4] == "@HWI":
						line = line.rstrip() + ";barcodelabel={}\n".format(sampleid)
					o.write(line)   

		return outfile
	

	def cleanup(self):
		"""Remove all intermediate files create"""
		
		files = [self.fwd_adp, self.rev_adp, self.fwd_trim, self.rev_trim, self.merged, self.reheaded]
		for file in files:
			try:
				os.remove(file)
			except OSError:
				pass


class CombinedReads(object):
	"""Class for setting up the otus"""

	def __init__(self, radius=1):
		self.radius = radius
		self.fasta = self.combine_merged_reads()
		self.derep = self.dereplicate(self.fasta)
		self.sorted = self.sort_reads(self.derep)
		self.otus = self.cluster_otus(self.sorted, self.radius)
		self.readmap = self.map_to_clusters(self.fasta, self.otus, self.radius)
		self.taxa = self.assign_taxa(self.otus)
		self.otu_table = self.convert_to_otu_table(self.readmap)
		self.parse_taxa()

	def combine_merged_reads(self):
		"""Get all the merged reads from the samples and concat them"""
		
		files = [x.reheaded for x in ReadPair.__all__]

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


	def cluster_otus(self, file, radius=1):
		"""Perform clustering with usearch"""

		identity = str(100 - radius)
		otus = 'otus_' + identity + '.fasta'
		up_out = 'uparse_out.csv'

		cmd = ['usearch', '-cluster_otus', file, '-otus', otus, '-uparseout', up_out, '-relabel', 'OTU_', '-sizein', '-sizeout']
		p = subprocess.Popen(cmd)
		p.wait()

		return otus

	
	def map_to_clusters(self, fasta, otus, radius):
		"""Map the merged reads to the cluster centroids"""
		
		id = 1 - (float(radius)/100)
		readmap = "readmap.uc"
		cmd = ['usearch',
				'-usearch_global', fasta, 
				'-db', otus, 
				'-strand', 'both', 
				'-id', str(id), 
				'-uc', readmap, 
				'-maxaccepts', '8', 
				'-maxrejects', '64', 
				'-top_hit_only']
		
		p = subprocess.Popen(cmd)
		p.wait()

		return readmap


	def convert_to_otu_table(self, readmap):
		"""Convert to a sensible format"""
		
		otu_table = 'otu_table.csv'	
		with open(otu_table, 'w') as o:	

			cmd = ['python', 'scripts/uc2otutab.py', self.readmap]
			p = subprocess.Popen(cmd, stdout = o)
			p.wait()

		return otu_table


	def assign_taxa(self, centroids):
		"""Assign taxomony to centroid sequences using RDP"""
		
		outfile = "otu_taxonomy.csv" 
		rdp_cls_binary = "/home/hugo.wefer/rdp_classifier_2.10.1/dist/classifier.jar" 
		
		cmd = ['java', '-jar', rdp_cls_binary, 'classify', centroids, '-o', outfile]
		
		p = subprocess.Popen(cmd)
		p.wait()
		
		return outfile
	

	def parse_taxa(self):
		"""Parse out the important parts of the rdp output"""

		outfile = 'parsed_taxa.csv'

		with open(outfile, 'w') as o:
			for item in parse_rdp.reformat(self.taxa).items():
				o.write('\t'.join(item))
				o.write('\n')

		return outfile


	def cleanup(self):
		"""Remove intermediate files"""

		files = [self.fasta, self.derep, self.sorted, self.taxa]
		for file in files:
			try:
				os.remove(file)
			except OSError:
				pass
		


