from collections import defaultdict 
import glob
import pandas as pd
import pprint
import pysam
import statistics

DEDUP_FILES = sorted(glob.glob('outFiles/*_Dedup.bam'))
GFF_FILE = '../gencode.vM27.annotation.gff3'
OUT_FILE_FREQ = 'LightSeq/TranscriptFrequencies.csv'
OUT_FILE = 'LightSeq/ReorderedLightSeq.csv'

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
                             'LS8A/AGGGTA': 'RGC_1',
                             'LS8B/AGGGTA': 'RGC_2',
                             'LS8C/AGGGTA': 'RGC_3',
                             'LS8D/AGGGTA': 'RGC_4',
                             'LS8A/GTTAGG': 'BP_1',
                             'LS8B/GTTAGG': 'BP_2',
                             'LS8C/GTTAGG': 'BP_3',
                             'LS8D/GTTAGG': 'BP_4',
                             'LS8A/TATGGA': 'ONL_1',
                             'LS8B/TATGGA': 'ONL_2',
                             'LS8C/TATGGA': 'ONL_3',
                             'LS8D/TATGGA': 'ONL_4'},
                    inplace=True)

  lightFreqs = lightFreqs[['Gene', 'RGC_1', 'RGC_2', 'RGC_3', 'RGC_4',
                           'BP_1', 'BP_2', 'BP_3', 'BP_4',
                           'ONL_1', 'ONL_2', 'ONL_3', 'ONL_4']]

  lightFreqs.set_index('Gene', inplace=True)
  # Add together counts with duplicate indices
  lightFreqs = lightFreqs.groupby(lightFreqs.index).sum()

  lightFreqs.to_csv(OUT_FILE)

if __name__ == '__main__':
  main()
