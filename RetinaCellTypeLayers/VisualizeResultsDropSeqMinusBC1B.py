from collections import defaultdict 
import pandas as pd

DROP_FREQ = 'DropSeq/DS_matrix.csv'
DROP_CLUSTERS = 'DropSeq/DS_clusters.csv'
OUT_FILE = 'DropSeq/DropSeqFrequenciesByCellTypeMinusBC1B.csv'
OUT_COUNT_FILE = 'DropSeq/DropSeqCellCountsMinusBC1B.csv'
OUT_LAYER_FILE = 'DropSeq/DropSeqLayerFrequenciesMinusBC1B.csv'
OUT_LAYER_COUNTS = 'DropSeq/DropSeqLayerCellCountsMinusBC1B.csv'

CELL_TYPE_CLUSTERS = {'Rods': [20, 18], 'Cones': [21], 'MG': [2],
                      'BCs': [1, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15]}
                      
DROP_REPS = ['Bipolar1', 'Bipolar2', 'Bipolar3', 'Bipolar4', 'Bipolar5',
             'Bipolar6']
DROP_CONDS = [{'Rods': 0.972, 'Cones': 0.028},
              {'BCs': 0.722, 'MG': 0.278}]
DROP_NAMES = ['ONL', 'BP']

def main():
  dropFreqs = pd.read_csv(DROP_FREQ, index_col=0)
  dropClusters = pd.read_csv(DROP_CLUSTERS)

  dropClusters['Replicate'] = [cell.split('_')[0] \
                               for cell in dropClusters['Cell']]

  cellTypeFreqs = defaultdict(list)
  cellTypeFreqs['Gene'] = dropFreqs.index
  cellCounts = defaultdict(lambda: [0])

  for cellType in CELL_TYPE_CLUSTERS:
    # Get list of clusters for this cell type
    clusters = CELL_TYPE_CLUSTERS[cellType]
    clustered = dropClusters.loc[dropClusters['Cluster'].isin(clusters)]

    for dropRep in DROP_REPS:
      replicateClustered = clustered.loc[dropClusters['Replicate'] == dropRep]
      condStr = '%s/%s' % (dropRep, cellType)
      cellTypeFreqs[condStr] = dropFreqs[dropFreqs.columns.intersection( \
                                 replicateClustered['Cell'])].sum(axis = 1)
      cellCounts[condStr] = [len(replicateClustered['Cell'])]

  # Now, save the dataframe
  cellTypeFreqsDF = pd.DataFrame.from_dict(cellTypeFreqs) 
  cellTypeFreqsDF.to_csv(OUT_FILE, index='Gene', header=True)

  countsDF = pd.DataFrame.from_dict(cellCounts)
  countsDF.to_csv(OUT_COUNT_FILE, index=False)
  
  # Add together counts with duplicate indices, if any
  cellTypeFreqsDF = cellTypeFreqsDF.groupby(cellTypeFreqsDF.index).sum()

  simulatedData = defaultdict(list)
  simulatedCellCounts = defaultdict(lambda: [0])

  for geneInd, gene in enumerate(cellTypeFreqsDF.index):
    simulatedData['Gene'].append(gene)

    for condInd in range(len(DROP_CONDS)):
      for repInd, dropRep in enumerate(DROP_REPS):
        dropVal = 0
        cellCount = 0
        for cellType in DROP_CONDS[condInd]:
          condStr = '%s/%s' % (dropRep, cellType)
          dropVal += round(cellTypeFreqsDF[condStr][geneInd] * \
                           DROP_CONDS[condInd][cellType])
          cellCount += cellCounts['%s/%s' % (dropRep, cellType)][0]

        simulatedData['%s_%s' % (DROP_NAMES[condInd],
                                 repInd + 1)].append(dropVal)
        simulatedCellCounts['%s_%s' % (DROP_NAMES[condInd],
                                       repInd + 1)][0] = cellCount

  pd.DataFrame.from_dict(simulatedData).to_csv(OUT_LAYER_FILE, index=False)
  pd.DataFrame.from_dict(simulatedCellCounts).to_csv(OUT_LAYER_COUNTS,
                                                     index=False)


if __name__ == '__main__':
  main()
