#!/usr/bin/env python

def reformat(infile):
	"""Reformat the output from RDP classifier"""

	taxa_dict = {}
	with open(infile, 'r') as f:
		for line in f.readlines():
			content = line.split('\t')
			otu_id = content[0]
			taxa = ';'.join(content[5::3])
			taxa_dict[otu_id] = taxa

	return taxa_dict


