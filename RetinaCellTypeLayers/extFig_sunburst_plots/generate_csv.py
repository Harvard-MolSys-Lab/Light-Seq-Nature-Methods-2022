# Parse out numbers from the outputted log files to generate the sunburst plot.

import re
import glob

import pandas as pd
import numpy as np


LOGS = 'logs/'


# Parse the log file from umi-tools extract barcode function
def parse_umitools_extract(file):

	regexes = {
		"total_reads": r"Input Reads: (\d+)",
		"barcoded_reads": r"regex matches read1: (\d+)",
		"no_barcode_reads": r"regex does not match read1: (\d+)",
	}
	
	parsed_data = {}
	with open(file) as f:
		for line in f:
			for k, r in regexes.items():
				r_search = re.search(r, line)
				if r_search:
					parsed_data[k] = float(r_search.group(1))
	
	return parsed_data


def parse_umitools_dedup(file):

	regexes = {
		"total_reads": r"Input Reads: (\d+)",
		"UMI": r"Number of reads out: (\d+)",
	}
	
	parsed_data = {}
	with open(file) as f:
		for line in f:
			for k, r in regexes.items():
				r_search = re.search(r, line)
				if r_search:
					parsed_data[k] = float(r_search.group(1))
	
	return parsed_data


# Parse Log file from STAR aligner 
# Regex strings and functions are modeled after the MultiQC parsing function.
# See https://github.com/ewels/MultiQC
# Added strings for unmapped numbers
def parse_star_report(file):
	
	regexes = {
		"total_reads": r"Number of input reads \|\s+(\d+)",
		"uniquely_mapped": r"Uniquely mapped reads number \|\s+(\d+)",
		"multimapped": r"Number of reads mapped to multiple loci \|\s+(\d+)",
		"multimapped_toomany": r"Number of reads mapped to too many loci \|\s+(\d+)",
		"unmapped_tooshort": r"Number of reads unmapped: too short \|\s+(\d+)",
		"unmapped_other": r"Number of reads unmapped: other \|\s+(\d+)",
		"unmapped_mismatches": r"Number of reads unmapped: too many mismatches \|\s+(\d+)",
	}
	
	parsed_data = {}
	with open(file) as f:
		for line in f:
			for k, r in regexes.items():
				r_search = re.search(r, line, re.MULTILINE)
				if r_search:
					parsed_data[k] = float(r_search.group(1))
	
	return parsed_data


# Parse the log file from featurecounts
def parse_featurecounts_summary(file):

	regexes = {
		"Assigned": r"Assigned\t(\d+)",
		"Unassigned_NoFeatures": r"Unassigned_NoFeatures\t(\d+)",
		"Unassigned_Ambiguity": r"Unassigned_Ambiguity\t(\d+)",
	}
	
	parsed_data = {}
	with open(file) as f:
		for line in f:
			for k, r in regexes.items():
				r_search = re.search(r, line,)
				if r_search:
					parsed_data[k] = float(r_search.group(1))
	
	return parsed_data


def parse_geneassigned_file(file):
	
	rRNA_dict = {
		'ENSMUSG00000119584.1': 'Rn18s',
		'ENSMUSG00000064337.1': 'mt-Rnr1',
	}
	
	gene_assigned = pd.read_csv(file, sep='\t', header=1)
	# Sort by highest number of counts (last column)
	gene_assigned.sort_values(
		gene_assigned.columns[-1], ascending=False, inplace=True
	)
	
		# Take the GeneID and counts of the top two entries
	gene1 = tuple(
		(
			rRNA_dict[gene_assigned.iloc[0].Geneid], #gene name
			gene_assigned.iloc[0][-1], #gene counts
		)
	)
	gene2 = tuple(
		(
			rRNA_dict[gene_assigned.iloc[1].Geneid], #gene name
			gene_assigned.iloc[1][-1], #gene counts
		)
	)
	
	return gene1, gene2


def file_name(path):
	fname = path.split('/')[-1].split('_')[0]
	return fname


def main():

	umi_extract_log = sorted(glob.glob('%s*extract.log' % LOGS))
	featurecounts_summary = sorted(glob.glob('%s*GeneAssigned.summary' % LOGS))
	gene_assigned_table = sorted(glob.glob('%s*GeneAssigned' % LOGS))
	STAR_report = sorted(glob.glob('%s*Log.final.out' % LOGS))
	umi_dedup_log = sorted(glob.glob('%s*Dedup.log' % LOGS))
	

	for idx in range(len(umi_extract_log)):

		umi_extract_parsed = parse_umitools_extract(
			umi_extract_log[idx]
		)
		featurecounts_parsed = parse_featurecounts_summary(
			featurecounts_summary[idx]
		)
		gene_assigned_parsed = parse_geneassigned_file(
			gene_assigned_table[idx]
		)
		star_parsed = parse_star_report(
			STAR_report[idx]
		)
		umi_dedup_parsed = parse_umitools_dedup(
			umi_dedup_log[idx]
		)

		### Set DF ###
		RNA_df = pd.DataFrame(
			index = np.arange(9),
			columns=[
				'Reads', 
				'Barcode', 
				'Mapping', 
				'Assignment', 
				'RNA type', 
				'UMI'
			]
		)
		###

		#Fill in the first 5 rows that just have one output, Fairly straightforward.

		RNA_df.at[0,'Reads'] = umi_extract_parsed['no_barcode_reads']
		RNA_df.at[0,'Barcode']  = 'No Barcode'

		RNA_df.at[1,'Reads'] = star_parsed['multimapped'] + star_parsed['multimapped_toomany']
		RNA_df.at[1,'Barcode'] = 'Has Barcode'
		RNA_df.at[1,'Mapping'] = 'Multimap'

		RNA_df.at[2,'Reads'] = star_parsed['unmapped_tooshort'] + star_parsed['unmapped_other'] + star_parsed['unmapped_mismatches']
		RNA_df.at[2,'Barcode'] = 'Has Barcode'
		RNA_df.at[2,'Mapping'] = 'Unmapped'

		RNA_df.at[3,'Reads'] = featurecounts_parsed['Unassigned_Ambiguity']
		RNA_df.at[3,'Assignment'] = 'Unassigned Ambiguity'
		RNA_df.at[3,'Barcode']  = 'Has Barcode'
		RNA_df.at[3,'Mapping']  = 'Unique Map'

		RNA_df.at[4,'Reads'] = featurecounts_parsed['Unassigned_NoFeatures']
		RNA_df.at[4,'Barcode']  = 'Has Barcode'
		RNA_df.at[4,'Mapping']  = 'Unique Map'
		RNA_df.at[4,'Assignment'] = 'Unassigned No Features'

		#Remaining counts parses out the assigned genes
		RNA_df.at[5,'Reads'] = umi_dedup_parsed['UMI']
		RNA_df.at[5,'Barcode']  = 'Has Barcode'
		RNA_df.at[5,'Mapping']  = 'Unique Map'
		RNA_df.at[5,'Assignment'] = 'Gene Assigned'
		RNA_df.at[5,'RNA type'] = 'RNA Transcripts'
		RNA_df.at[5,'UMI'] = 'UMI'

		# The sunburst plot requires the data to be structured with preceding "roots" for each "branch"
		# This slice of the layer was ultimately removed for the final figure
		RNA_df.at[6,'Reads'] = umi_dedup_parsed['total_reads'] - umi_dedup_parsed['UMI']
		RNA_df.at[6,'Barcode']  = 'Has Barcode'
		RNA_df.at[6,'Mapping']  = 'Unique Map'
		RNA_df.at[6,'Assignment'] = 'Gene Assigned'
		RNA_df.at[6,'RNA type'] = 'RNA Transcripts'
		RNA_df.at[6,'UMI'] = 'Duplicate UMIs (Placeholder)'

		RNA_df.at[7,'Reads'] = gene_assigned_parsed[0][1]
		RNA_df.at[7,'Barcode']  = 'Has Barcode'
		RNA_df.at[7,'Mapping']  = 'Unique Map'
		RNA_df.at[7,'Assignment'] = 'Gene Assigned'
		RNA_df.at[7,'RNA type'] = gene_assigned_parsed[0][0]

		RNA_df.at[8,'Reads'] = gene_assigned_parsed[1][1]
		RNA_df.at[8,'Barcode']  = 'Has Barcode'
		RNA_df.at[8,'Mapping']  = 'Unique Map'
		RNA_df.at[8,'Assignment'] = 'Gene Assigned'
		RNA_df.at[8,'RNA type'] = gene_assigned_parsed[1][0]

		#Show final DF
		file_prefix = file_name(umi_extract_log[idx])
		print(file_prefix)
		print(RNA_df)

		RNA_df.to_csv('%s_RNADist.csv' %file_prefix, index=False)

if __name__ == '__main__':
	main()

	