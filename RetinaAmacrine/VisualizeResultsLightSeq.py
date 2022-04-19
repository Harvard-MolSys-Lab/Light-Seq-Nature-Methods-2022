from collections import defaultdict 
import glob
import pandas as pd
import pprint
import pysam
import statistics

DEDUP_FILES = sorted(glob.glob('outFiles/*_Dedup.bam'))
GFF_FILE = '../gencode.vM27.annotation.gff3'
OUT_FILE_FREQ = 'TranscriptFrequencies.csv'
OUT_FILE = 'ReorderedLightSeq.csv'

def readGFF():
  genes = {}
  
  print('Reading GFF file...')
  with open(GFF_FILE) as gffHandle:
    for line in gffHandle:
      if line[0] != '#' and line.split('\t')[2] == 'gene':
        genes[line.split('gene_id=')[1].split(';')[0]] = line
  
  return genes

def main():
  filesToAnalyze = DEDUP_FILES

  genes = readGFF()
  allGenes = []
  allConditionBarcodes = []
  frequencies = defaultdict(lambda: defaultdict(int)) # default to 0 count

  # Create dictionary of gene frequencies
  for geneFile in filesToAnalyze:  
    print('Reading %s...' % geneFile)

    condition = geneFile.split('/')[-1].split('_Dedup')[0]
    samFile = pysam.AlignmentFile(geneFile, "rb")
    samIter = samFile.fetch(until_eof=True)

    for read in samIter:
      barcodeSeq = read.query_name.split('_')[1]
      geneName = read.get_tag('XT')

      condBar = '%s/%s' % (condition, barcodeSeq)

      frequencies[geneName][condBar] += 1

      if not condBar in allConditionBarcodes:
        allConditionBarcodes.append(condBar)

  # Now, create CSV
  allData = defaultdict(list)

  allData['Gene'] = frequencies.keys()

  print('Compiling data...')
  for gene in allData['Gene']:
    allData['Gene name'].append(genes[gene].split('gene_name=')[1] \
                                .split(';')[0])
    for condBar in allConditionBarcodes:
      allData[condBar].append(frequencies[gene][condBar])

  lightFreqs = pd.DataFrame.from_dict(allData)
  lightFreqs.to_csv(OUT_FILE_FREQ, index=False)

  # Now, re-order data according to layer, replicate information
  del lightFreqs['Gene']


  lightFreqs.rename(columns={'Gene name': 'Gene',
                             'TLS23A/GTTAGG': 'Th_1',
                             'TLS23B/GTTAGG': 'Th_2',
                             'TLS23D/GTTAGG': 'Th_3',
                             'TLS23E/GTTAGG': 'Th_4',
                             'TLS23F/GTTAGG': 'Th_5',
                             'TLS23A/AGGGTA': 'Am_1',
                             'TLS23B/AGGGTA': 'Am_2',
                             'TLS23D/AGGGTA': 'Am_3',
                             'TLS23E/AGGGTA': 'Am_4',
                             'TLS23F/AGGGTA': 'Am_5'},
                    inplace=True)

  lightFreqs = lightFreqs[['Gene', 'Th_1', 'Th_2', 'Th_3', 'Th_4',
                           'Th_5', 'Am_1', 'Am_2', 'Am_3',
                           'Am_4', 'Am_5']]

  lightFreqs.set_index('Gene', inplace=True)
  # Add together counts with duplicate indices
  lightFreqs = lightFreqs.groupby(lightFreqs.index).sum()

  lightFreqs.to_csv(OUT_FILE)

if __name__ == '__main__':
  main()
