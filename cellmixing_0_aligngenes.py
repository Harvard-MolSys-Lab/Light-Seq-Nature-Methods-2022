import os
import glob

PATH_BAMFILES = "outFiles_genome_alignPerGene/"
FILE_GTF      = "/Users/ninning/o2Download/STAR_mergedIndex/gencode.v35andvM25.annotation.gff3"
BAMFILES      = sorted(glob.glob('%s*sortedByCoord.out.bam' %PATH_BAMFILES))
PATH_OUTFILE  = "436_438_outFiles_deduped/"
#UMI dedup options
OUTPUT_STATS  = False

def getFileNamesAll(files):
    filenames = [file.split('/')[-1].split('.bam')[0] for file in files]
    print (filenames)

def getFileName(bamfile):
	filename = (bamfile.split('/')[-1].split('.bam')[0])
	return filename

def featureCounts(bamfile):
	filename = getFileName(bamfile)
	os.system('featureCounts -a %s -o %s -R BAM %s' % \
			 (FILE_GTF, \
			  PATH_BAMFILES + filename + '.gene_assigned', \
			  bamfile))

def sortIndex(bamfile):
	fcountFile = bamfile + '.featureCounts.bam'
	filename = getFileName(bamfile)

	print("Sorting File: %s" % (fcountFile))
	os.system('samtools sort %s -o %s' % \
			 ( fcountFile, PATH_BAMFILES + filename + '_sorted.bam'))

	print("Indexing Sorted File: %s" % (PATH_BAMFILES + filename + '_sorted.bam'))
	os.system('samtools index %s' % \
			 ( PATH_BAMFILES + filename + '_sorted.bam'))

def umiDedup(bamfile):
	#Make sure you have the sorted.bam and index file present in the directory
	filename = getFileName(bamfile)

	if OUTPUT_STATS:
		umistr = 'umi_tools dedup --per-gene --gene-tag=XT --assigned-status-tag=XS -I %s --output-stats=%s -S %s --log %s' % \
				 (PATH_BAMFILES + filename + '_sorted.bam', \
				  PATH_OUTFILE + filename + '_dedupstats', \
				  PATH_OUTFILE + filename + '_dedup.bam',  \
				  PATH_OUTFILE + filename + '_deduplog.log')
		print('\n' + 'Running: ' + umistr + '\n')
		os.system(umistr)
	else:
		umistr = 'umi_tools dedup --per-gene --gene-tag=XT --assigned-status-tag=XS -I %s -S %s --log %s' % \
				 (PATH_BAMFILES + filename + '_sorted.bam', \
				  PATH_OUTFILE + filename + '_dedup.bam',  \
				  PATH_OUTFILE + filename + '_deduplog.log')
		print('\n' + 'Running: ' + umistr + '\n')
		os.system(umistr)  

def main():

	print('Getting file names: ')
	FILE_NAMES = getFileNamesAll(BAMFILES)

	#Run analysis one file at a time
	for bamfile in BAMFILES:
		featureCounts(bamfile)
		sortIndex(bamfile)
		umiDedup(bamfile)

if __name__ == '__main__':
	main()