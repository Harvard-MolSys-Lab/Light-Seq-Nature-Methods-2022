import glob
import os

filelist = glob.glob('*_Dedup.bam')

for file in filelist:
	print("Indexing File: %s" % file)
	os.system('samtools index %s' % file)
