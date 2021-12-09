from collections import defaultdict 
import glob
import gzip
import os
import pysam

from Bio import SeqIO
from multiprocessing import Pool
from skbio.alignment import StripedSmithWaterman

IN_FILES = sorted(glob.glob('inFiles/*_R1_*.fastq.gz'))
OUT_DIR = 'outFiles/'

# Copied from https://www.ncbi.nlm.nih.gov/nuccore/L29345.1?report=fasta
# on 5/23/2021. Took reverse complement.
GFP_SEQUENCE = 'ATACACTCCAGTAGCCTATTTAATAAGAAAATAGCCCCTATTAATAACATCAATAAATTAT' \
               'TCATAAAATTTTAATGAATCTATAAATATATAAATAAAGTCTCAGCCTGAATTTAACCAGG' \
               'AACCCTGAGAATTTAGTAATTGTTCGGACACTTTAGTGTCAATTGGAAGTCTGGACATTTA' \
               'TTTGTATAGTTCATCCATGCCATGTGTAATCCCAGCAGCTGTTACAAACTCAAGAAGGATC' \
               'ATGTGATCTCTCTTTTCGTTGGGATCTTTGGAAAGGGCAGATTGTGTGGACAGGTAATGGT' \
               'TGTCTGGTAAAAGGACAGGGCCATCGCCAATTGGAGTATTTTGTTGATAATGGTCTGCTAA' \
               'TTGAACGCTTCCATCTTTAATGTTGTGTCTAATTTTGAAGTTAACTTTGATTCCATTCTTT' \
               'GGTTTGTCTGCCATGATGTATACATTATGTGAGTTATAGTTGTATTCCATTTTGTGTCCAA' \
               'GAATGTTTCCATCTTCTTTAAAATCAATACCTTTTAACTCGATTCTATTAACAAGGGTATC' \
               'ACCTTCAAACTTGACTTCAGCACGTGTCTTGTAGTTCCCGTCATCTTTGTAAAATATAGTT' \
               'CTTTCCTGTACATAACCTTCGGGCATGGCACTCTTGAAAAAGTCATGCTGTTTCATATGAT' \
               'CTGGGTATCTTGAAAAGCATTGAACACCATAAGAGAAAGTAGTGACAAGTGTTGGCCATGG' \
               'AACAGGTAGCTTCCCAGTAGTGCAAATAAATTTAAGGGTAAGTTTTCCGTATGTTGCATCA' \
               'CCTTCACCCTCTCCACTGACAGAGAATTTTTGCCCATTAACATCGCCATCTAATTCAACAA' \
               'GAATTGGGACAACTCCAGTGAAAAGTTCTTCTCCTTTACTCATCTTTGTTATCTTTTATTC' \
               'GTGTGTA'

# Copied from SnapGene https://www.snapgene.com/resources/plasmid-files/?set=
# fluorescent_protein_genes_and_plasmids&plasmid=EGFP on 5/27/2021. Took reverse
# complement.
EGFP_SEQUENCE = 'TTACTTGTACAGCTCGTCCATGCCGAGAGTGATCCCGGCGGCGGTCACGAACTCCAGCAG' \
                'GACCATGTGATCGCGCTTCTCGTTGGGGTCTTTGCTCAGGGCGGACTGGGTGCTCAGGTA' \
                'GTGGTTGTCGGGCAGCAGCACGGGGCCGTCGCCGATGGGGGTGTTCTGCTGGTAGTGGTC' \
                'GGCGAGCTGCACGCTGCCGTCCTCGATGTTGTGGCGGATCTTGAAGTTCACCTTGATGCC' \
                'GTTCTTCTGCTTGTCGGCCATGATATAGACGTTGTGGCTGTTGTAGTTGTACTCCAGCTT' \
                'GTGCCCCAGGATGTTGCCGTCCTCCTTGAAGTCGATGCCCTTCAGCTCGATGCGGTTCAC' \
                'CAGGGTGTCGCCCTCGAACTTCACCTCGGCGCGGGTCTTGTAGTTGCCGTCGTCCTTGAA' \
                'GAAGATGGTGCGCTCCTGGACGTAGCCTTCGGGCATGGCGGACTTGAAGAAGTCGTGCTG' \
                'CTTCATGTGGTCGGGGTAGCGGCTGAAGCACTGCACGCCGTAGGTCAGGGTGGTCACGAG' \
                'GGTGGGCCAGGGCACGGGCAGCTTGCCGGTGGTGCAGATGAACTTCAGGGTCAGCTTGCC' \
                'GTAGGTGGCATCGCCCTCGCCCTCGCCGGACACGCTGAACTTGTGGCCGTTTACGTCGCC' \
                'GTCCAGCTCGACCAGGATGGGCACCACCCCGGTGAACAGCTCCTCGCCCTTGCTCACCAT'

# Copied from SnapGene https://www.snapgene.com/resources/plasmid-files/?set=
# fluorescent_protein_genes_and_plasmids&plasmid=Emerald_GFP on 11/8/2021. Took
# reverse complement.
EMGFP_SEQUENCE = 'TTACTTGTACAGCTCGTCCATGCCGAGAGTGATCCCGGCGGCGGTCACGAACTCCAGCA' \
                 'GGACCATGTGATCGCGCTTCTCGTTGGGGTCTTTGCTCAGGGCGGACTGGGTGCTCAGG' \
                 'TAGTGGTTGTCGGGCAGCAGCACGGGGCCGTCGCCGATGGGGGTGTTCTGCTGGTAGTG' \
                 'GTCGGCGAGCTGCACGCTGCCGTCCTCGATGTTGTGGCGGGTCTTGAAGTTCACCTTGA' \
                 'TGCCGTTCTTCTGCTTGTCGGCGGTGATATAGACCTTGTGGCTGTTGTAGTTGTACTCC' \
                 'AGCTTGTGCCCCAGGATGTTGCCGTCCTCCTTGAAGTCGATGCCCTTCAGCTCGATGCG' \
                 'GTTCACCAGGGTGTCGCCCTCGAACTTCACCTCGGCGCGGGTCTTGTAGTTGCCGTCGT' \
                 'CCTTGAAGAAGATGGTGCGCTCCTGGACGTAGCCTTCGGGCATGGCGGACTTGAAGAAG' \
                 'TCGTGCTGCTTCATGTGGTCGGGGTAGCGGGCGAAGCACTGCACGCCGTAGGTCAAGGT' \
                 'GGTCACGAGGGTGGGCCAGGGCACGGGCAGCTTGCCGGTGGTGCAGATGAACTTCAGGG' \
                 'TCAGCTTGCCGTAGGTGGCATCGCCCTCGCCCTCGCCGGACACGCTGAACTTGTGGCCG' \
                 'TTTACGTCGCCGTCCAGCTCGACCAGGATGGGCACCACCCCGGTGAACAGCTCCTCGCC' \
                 'CTTGCTCACCAT'

def mapGFPs(filePrefix, trimmedR1File, sortedFile):
  with open('%s_results_compare.txt' % filePrefix, 'a') as resultsFile:

    # Load sorted file to see which sequences successfully mapped
    print('  Loading sorted file %s...' % sortedFile)

    samFile = pysam.AlignmentFile(sortedFile, "rb")

    samIter = samFile.fetch(until_eof=True)

    # Track which queries successfully mapped
    hits = set()

    for read in samIter:
      hits.add(read.query_name)

    print('Mapping to GFP sequences...')

    gfpSeqs = {'GFP': GFP_SEQUENCE,
              'eGFP': EGFP_SEQUENCE,
              'emGFP': EMGFP_SEQUENCE}
    gfpUMIs = {'GFP': [],
               'eGFP': [],
               'emGFP': []}
    gfpMaps = {'GFP': defaultdict(lambda: 0),
               'eGFP': defaultdict(lambda: 0),
               'emGFP': defaultdict(lambda: 0)}

    with gzip.open(trimmedR1File, 'rt') as trimHandle:
      for record in SeqIO.parse(trimHandle, 'fastq'):
        if not record.id in hits:
          
          query = StripedSmithWaterman(str(record.seq))
          for gSeq in gfpSeqs:
            alignment = query(gfpSeqs[gSeq])

            if alignment.optimal_alignment_score >= 40:
              extractedBarcode = record.id.split('_')[1]
              extractedUMI = record.id.split('_')[2]
              
              if not extractedUMI in gfpUMIs[gSeq]:
                gfpMaps[gSeq][extractedBarcode] += 1
                gfpUMIs[gSeq].append(extractedUMI)
                resultsFile.write('  Map to %s (score: %d, bar: %s): %s\n' % \
                                  (gSeq, alignment.optimal_alignment_score,
                                   extractedBarcode, str(record.seq)))

    for gSeq in gfpSeqs:
      for barcodeSeq in gfpMaps[gSeq]:
        resultsFile.write('  %s, barcode %s: %d maps to %s' % \
                          (filePrefix, barcodeSeq,
                           gfpMaps[gSeq][barcodeSeq], gSeq))

def main():

  parseFiles = []
  for inFile in IN_FILES:
    filePrefix = '%s%s' % (OUT_DIR, inFile.split('/')[1].split('_R1_')[0])
    trimmedR1File = '%s_R1_trimmed.fastq.gz' % filePrefix

    for genome in ['_Merged']: # Only filtered out merged genome maps
      sortedFile = '%s%s_Sorted.bam' % (filePrefix, genome)
      parseFiles.append(('%s%s' % (filePrefix, genome), trimmedR1File,
                                   sortedFile))


  with Pool() as pool:
    pool.starmap(mapGFPs, parseFiles)

if __name__ == '__main__':
  main()
