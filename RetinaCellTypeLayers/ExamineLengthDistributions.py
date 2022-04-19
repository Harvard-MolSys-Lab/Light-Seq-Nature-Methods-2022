from collections import defaultdict
import collections
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

CONDITIONS = ['LS8A', 'LS8B', 'LS8C', 'LS8D', 'LS8E', 'LS8F', 'LS8H']
BARCODES = ['AGGGTA', 'GTTAGG', 'TATGGA']
IN_FILE_FREQ = 'LightSeq/TranscriptFrequencies.csv'
GENE_ASSIGNED_FILE = 'outFiles/%s_GeneAssigned'
BINS = [[0, 1475],
        [1475, 2130],
        [2130, 2730],
        [2730, 3339],
        [3339, 3991],
        [3991, 4724],
        [4724, 5688],
        [5688, 6986],
        [6986, 9226],
        [9226, 116825]]
ORDERED_BINS = ['0-1475', '1475-2130', '2130-2730',
                '2730-3339', '3339-3991', '3991-4724',
                '4724-5688', '5688-6986', '6986-9226',
                '9226-116825']
PLOT_DIR = 'Plots/'

def main():
  transcriptCounts = pd.read_csv(IN_FILE_FREQ, index_col=0)

  for cond in CONDITIONS:
    assignFile = pd.read_csv(GENE_ASSIGNED_FILE % cond, sep='\t',
                             skiprows=1, index_col=0)
    for bar in BARCODES:
      print('Condition %s/%s' % (cond, bar))

      lengths = []
      binCountVals = defaultdict(list)
      binRPKMVals = defaultdict(list)
      totalReads = transcriptCounts['%s/%s' % (cond, bar)].sum()
      print('total reads: %d' % totalReads)
      millionScalingFactor = totalReads / 1000000.0

      for gene in transcriptCounts.index:
        geneLength = assignFile['Length'][gene]
        tranCount = transcriptCounts['%s/%s' % (cond, bar)][gene]

        for _ in range(tranCount):
          lengths.append(geneLength)

        for i in range(len(BINS)):
          if BINS[i][0] <= geneLength < BINS[i][1] and tranCount > 0:
            binCountVals['%d-%d' % (BINS[i][0], BINS[i][1])].append(tranCount)
            binRPKMVals['%d-%d' % (BINS[i][0], BINS[i][1])].append( \
                tranCount / (geneLength / 1000.0) / millionScalingFactor)

      # Check that length of list matches number of transcripts counted
      print('  Length of lengths list: %d' % len(lengths))
      print('  Total counts for condition: %d' % \
              transcriptCounts['%s/%s' % (cond, bar)].sum())

      plt.figure()
      plt.hist(lengths, histtype='step', bins=50)
      plt.yscale('log')
      plt.xlabel('Transcript length (nt)')
      plt.ylabel('Read counts')
      plt.title('%s/%s' % (cond, bar))
      plt.savefig('%s%s_%s_Hist.svg' % (PLOT_DIR, cond, bar))
      plt.close()

      plt.figure()
      listOfListsCounts = [binCountVals[key] for key in ORDERED_BINS]
      ax = sns.boxplot(data=listOfListsCounts)
      ax.set_xticklabels(ORDERED_BINS)
      plt.xticks(rotation = 30)
      plt.xlabel('Transcript length (nt)')
      plt.ylabel('Read counts')
      plt.yscale('log')
      plt.tight_layout()
      plt.savefig('%s%s_%s_CountsPerBin.svg' % (PLOT_DIR, cond, bar))
      plt.close()

      plt.figure()
      listOfListsRPKM = [binRPKMVals[key] for key in ORDERED_BINS]
      ax = sns.boxplot(data=listOfListsRPKM)
      ax.set_xticklabels(ORDERED_BINS)
      plt.xticks(rotation = 30)
      plt.xlabel('Transcript length (nt)')
      plt.ylabel('RPKM')
      plt.title('%s/%s' % (cond, bar))
      plt.yscale('log')
      plt.tight_layout()
      plt.savefig('%s%s_%s_RPKMPerBin.svg' % (PLOT_DIR, cond, bar))
      plt.close()

if __name__ == '__main__':
  main()