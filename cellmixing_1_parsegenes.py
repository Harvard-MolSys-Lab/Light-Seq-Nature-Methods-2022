# Changed to count ALL genes using gene_id. sometimes two ENSEMBL Ids can have the same gene name.
# NL
import os
import glob
import csv
from collections import defaultdict

import pandas as pd
import pysam

from utils import setNumbersDF

# Seq files are not hosted on the github platform due to size limitations.
# This script is meant to be run on the .bam files for the mouse human mixing experimens (436, 438)
# Set the folder for the .bam files that have been deduped from UMI tools.
# Files should have *dedup.bam. Script also uses pysam which requires an indexed ".bam.bai" file.
# So run the samtools indexer in the folder as well.

PATH_BAMFILES = "436_438_outFiles_deduped/"
PATH_GENES = "results_cellmixing_genes/"
PATH_NUMBERS = "results_cellmixing_numbers/"

BAMFILES = sorted(glob.glob('%s*_dedup.bam' % PATH_BAMFILES))
DICTIONARY = "genes_dict_combined/genesDict_types_lengths_combined.csv"

B1 = "GTTAGG"
B2 = "TATGGA"

if os.path.isdir("results_cellmixing_genes") is False:
	os.mkdir("results_cellmixing_genes")

if os.path.isdir("results_cellmixing_numbers") is False:
	os.mkdir("results_cellmixing_numbers")


class Parser(object):
	"""
	Parse out BC1 or BC2 for sequences
	Also sets the dataframe to be ENSG or ENSMUSG depending on the input reads
	"""
	def __init__(self, file, READS):
		self.file = file
		self.filename = (file.split('/')[-1])
		self.fname_short = self.filename.split("_")[0]
		self.reads = READS
		self.parseBC()
		self.species = '_species_'
		
	def parseBC(self):
		self.bc1 = []
		self.bc2 = []
		for idx, read in enumerate(self.reads):
			bc = read.query_name.split("_")[1]
			if B1 in bc:
				self.bc1.append(idx)
			if B2 in bc:
				self.bc2.append(idx)
		
		return self.bc1, self.bc2
				
	def geneCount(self, bc_set):
		genecounts = defaultdict(lambda: 0)

		for _, read in enumerate(bc_set):
			gene_ID = self.reads[read].get_tag("XT")
			genecounts[gene_ID] += 1
		
		return genecounts
	
	# Determines whether you entered Mouse or Human Reads
	# Will also return a string for later use.
	def setDF(self):
		df = pd.read_csv(DICTIONARY, index_col=0)
		df["gene_count"] = 0
		m_index = list(map(lambda name: "ENSMUSG" in name, list(df.index))).index(True)  # Finds first instance of a mouse gene

		if "ENSG" in self.reads[0].get_tag("XT"):
			df = df[:m_index]
			self.species = "_Human_"
		elif "ENSMUSG" in self.reads[0].get_tag("XT"):
			df = df[m_index:]
			self.species = "_Mouse_"
			
		return df

	def geneCountAll(self, bc_set):
		genecounts = self.geneCount(bc_set)
		df = self.setDF()

		for item in genecounts:
			count = genecounts.get(item)
			df.loc[item, "gene_count"] = count

		return df


# Find first instance of mChr1
def indexHumanMouse(samfile):
	temp_index = samfile.references.index('mchr1')

	gindex_human = range(0, temp_index)
	gindex_mouse = range(temp_index, len(samfile.references))
	
	return gindex_human, gindex_mouse


def parseReads(gindex_human, gindex_mouse, samfile):
	# Parse out the reads for the entire bamfile
	# Each item is an alignment file
	HUMAN_READS = []
	for idx in gindex_human:
		iter = samfile.fetch(tid=idx, until_eof=True)
		for reads in iter:
			HUMAN_READS.append(reads)
			
	MOUSE_READS = []
	for idx in gindex_mouse:
		iter = samfile.fetch(tid=idx, until_eof=True)
		for reads in iter:
			MOUSE_READS.append(reads)
	
	return HUMAN_READS, MOUSE_READS


def main():
	for bamfile in BAMFILES:

		samfile = pysam.AlignmentFile(bamfile, "rb")
		gindex_human, gindex_mouse = indexHumanMouse(samfile)
		HUMAN_READS, MOUSE_READS = parseReads(gindex_human, gindex_mouse, samfile)

		HUMAN_MAPS = Parser(bamfile, HUMAN_READS)
		MOUSE_MAPS = Parser(bamfile, MOUSE_READS)

		# HUMAN_MAPS.geneCountAll(HUMAN_MAPS.bc1).to_csv(
		# 	PATH_GENES + HUMAN_MAPS.fname_short + HUMAN_MAPS.species + "bc1_AllCounts.csv"
		# 	)
		HUMAN_MAPS.geneCountAll(HUMAN_MAPS.bc2).to_csv(
			PATH_GENES + HUMAN_MAPS.fname_short + HUMAN_MAPS.species + "bc2_AllCounts.csv"
			)
		MOUSE_MAPS.geneCountAll(MOUSE_MAPS.bc1).to_csv(
			PATH_GENES + MOUSE_MAPS.fname_short + MOUSE_MAPS.species + "bc1_AllCounts.csv"
			)
		# MOUSE_MAPS.geneCountAll(MOUSE_MAPS.bc2).to_csv(
		# 	PATH_GENES + MOUSE_MAPS.fname_short + MOUSE_MAPS.species + "bc2_AllCounts.csv"
		# 	)

		fname_short = HUMAN_MAPS.fname_short  # Either maps name is fine
		numbers_df = setNumbersDF.setDF(HUMAN_MAPS, MOUSE_MAPS, samfile.mapped)

		numbers_df.to_csv(PATH_NUMBERS + fname_short + "_NUMBERS.csv")
		print("Finished with file: " + fname_short)


if __name__ == '__main__':
	main()
