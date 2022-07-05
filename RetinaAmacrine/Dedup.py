import glob
import os
import pprint
import pysam
import statistics
from collections import defaultdict

GENE_FILES = sorted(glob.glob('outFiles/*Aligned.sortedByCoord.out.bam'))
GFF_FILE = '../gencode.vM27.annotation.gff3'

OUT_DIR = 'outFiles/'

GENES_TO_FILTER = ['ENSMUSG00000119584.1', 'ENSMUSG00000064337.1']

def createSmallerGeneFile(filePrefix, geneFile):
  outGeneFile = '%s.bam' % filePrefix

  gF = pysam.AlignmentFile(geneFile)
  outGF = pysam.AlignmentFile(outGeneFile, 'w', template=gF)
  
  filteredCounts = defaultdict(int)

  for read in gF.fetch(until_eof=True):
    if read.has_tag('XT'):
      xtTag = read.get_tag('XT')

      if xtTag in GENES_TO_FILTER:
        filteredCounts[xtTag] += 1
      else:
        outGF.write(read)

  gF.close()
  outGF.close()

  print('Filtered out %s' % ', '.join(['%d of %s' % \
                                       (filteredCounts[gene], gene) \
                                       for gene in GENES_TO_FILTER]))

  return outGeneFile

def deduplicate(filePrefix, geneFile):
  assignedGeneFile = '%s_GeneAssigned' % filePrefix
  countedGeneFile = '%s.featureCounts.bam' % geneFile  
  sortedGeneFile = '%s_Sorted.bam' % filePrefix
  dedupGeneFile = '%s_Dedup.bam' % filePrefix
  dedupLogFile = '%s_Dedup.log' % filePrefix

  # Use feature-counts to count genes
  print('  Counting genes with featureCounts...')
  os.system('featureCounts -a %s -o %s -R BAM %s' % \
            (GFF_FILE, assignedGeneFile, geneFile))

  # Now, filter out reads from the two problem loci
  filteredCountsFile = createSmallerGeneFile(filePrefix, countedGeneFile)

  # De-duplicate
  print('  Deduping...')

  # Filter, sort, index the bam files then deduplicate
  # NOTE: the filtering is not random
  os.system('samtools sort %s -o %s' % (filteredCountsFile,
                                        sortedGeneFile))

  os.system('samtools index %s' % sortedGeneFile)

  os.system(('umi_tools dedup --per-gene --gene-tag=XT ' \
             '--assigned-status-tag=XS -I %s -S %s -L %s') \
             % (sortedGeneFile, dedupGeneFile,
                dedupLogFile))

def main():  
  # Deduplicate each file individually
  for geneFile in GENE_FILES:
    filePrefix = '%s%s' % (OUT_DIR, geneFile.split('/')[1].split('_R1_')[0])
    deduplicate(filePrefix, geneFile)

if __name__ == '__main__':
  main()
