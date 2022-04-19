import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import statistics

from collections import defaultdict 
from scipy import stats

DROP_FREQS = 'DropSeq/DropSeqLayerFrequencies.csv'
DROP_CELL_COUNTS = 'DropSeq/DropSeqLayerCellCounts.csv'
DROP_PADJ_FILE = 'DropSeq/BPvsONL.csv'
LIGHT_FREQS = 'LightSeq/ReorderedLightSeq.csv'
LIGHT_PADJ_FILE = 'LightSeq/BPvsONL.csv'
LIGHT_CELL_COUNTS = 'LightSeq/LightSeqCellCounts.csv'

GENE_ASSIGN_FILE = 'outFiles/LS8A_GeneAssigned'
IN_FILE_FREQ = 'LightSeq/TranscriptFrequencies.csv'

DROP_REPS = ['_1', '_2', '_3', '_4', '_5', '_6']
LIGHT_REPS = ['_1', '_2', '_3', '_4']

COLOR_MAP = matplotlib.cm.get_cmap('viridis_r')

BINS = [[0, 2130],
        [2130, 3339],
        [3339, 4724],
        [4724, 6986],
        [6986, 116825]]
ORDERED_BINS = ['0-2130', '2130-3339', '3339-4724',
                '4724-6986', '6986-116825']

def plotEnrichments(lightEnrich, dropEnrich, includedGenes, annotByGene, fignum,
                    geneLengths, tranCounts, padjVals, xlim=None, ylim=None):
  # Calculate means and replicate errors for enrichment values
  lightEnrichMeans = []
  lightEnrichErrs = []
  dropEnrichMeans = []
  dropEnrichErrs = []
  lengths = []

  for gene in includedGenes:
    lightMean = statistics.mean(lightEnrich[gene])
    lightErr = statistics.stdev(lightEnrich[gene])
    dropMean = statistics.mean(dropEnrich[gene])
    dropErr = statistics.stdev(dropEnrich[gene])

    lightEnrichMeans.append(lightMean)
    lightEnrichErrs.append(lightErr)
    dropEnrichMeans.append(dropMean)
    dropEnrichErrs.append(dropErr)

    tranNames = tranCounts.index[tranCounts['Gene name'] == gene].tolist()

    # This should not print anything, but double checks unique mapping
    if len(tranNames) != 1:
      print('Gene %s has multiple transcript names: %s' % (gene, tranNames))

    tranLength = geneLengths[tranNames[0]]
    lengths.append(tranLength)

  plt.figure()
  ax = plt.gca()

  logLengths = [math.log10(length) for length in lengths]
  norm = matplotlib.colors.Normalize(vmin=min(logLengths), vmax=max(logLengths))
  colors = [COLOR_MAP(norm(logLengths[i])) \
              for i in range(len(lightEnrichMeans))]

  plt.scatter(lightEnrichMeans, dropEnrichMeans, alpha=0.4, c=logLengths, s=0,
              cmap=COLOR_MAP, zorder=2)

  print('Minimum length: %d nt, Maximum length: %d nt' % \
          (min(lengths), max(lengths)))
  print('Minimum log length: %f, Maximum log length: %f' % \
          (min(logLengths), max(logLengths)))

  # Plot
  for i in range(len(lightEnrichMeans)):
    ax.errorbar(lightEnrichMeans[i], dropEnrichMeans[i],
                xerr=lightEnrichErrs[i], yerr=dropEnrichErrs[i], fmt='o',
                color=colors[i], ms=0, alpha=0.4, zorder=1)

  plt.xlabel('Light-Seq ONL - BCL, approx reads per cell')
  plt.ylabel('Drop-Seq ONL - BCL, approx reads per cell')

  if xlim is not None:
    plt.xlim(xlim)
  if ylim is not None:
    plt.ylim(ylim)

  plt.colorbar()

  for i in range(len(lightEnrichMeans)):
    if includedGenes[i] in annotByGene:
      if lightEnrichMeans[i] > 0:
        ax.scatter(lightEnrichMeans[i], dropEnrichMeans[i],
                   s=2, c='k', zorder=3)
        ax.annotate(includedGenes[i], (lightEnrichMeans[i],
                                       dropEnrichMeans[i]),
                    horizontalalignment='left', verticalalignment='bottom',
                    color='k', zorder=4)
      else:
        ax.scatter(lightEnrichMeans[i], dropEnrichMeans[i],
                   s=2, c='k', zorder=3)
        ax.annotate(includedGenes[i], (lightEnrichMeans[i],
                                       dropEnrichMeans[i]),
                    horizontalalignment='left', verticalalignment='bottom',
                    color='k', zorder=4)

  plt.grid(True)
  ax.set_aspect('equal')

  plt.savefig('Plots/%d.svg' % fignum)

def plotBipolarOnly(lightBCL, dropBCL, genesToPlot, geneLengths, tranCounts,
                    padjVals):
  # Calculate means and replicate errors for enrichment values
  lightMeans = []
  lightErrs = []
  dropMeans = []
  dropErrs = []

  ratios = []
  lengths = []

  print('Plotting %d BCL enriched genes' % len(genesToPlot))

  for gene in genesToPlot:
    lightMean = statistics.mean(lightBCL[gene])
    lightErr = statistics.stdev(lightBCL[gene])
    dropMean = statistics.mean(dropBCL[gene])
    dropErr = statistics.stdev(dropBCL[gene])

    lightMeans.append(lightMean)
    lightErrs.append(lightErr)
    dropMeans.append(dropMean)
    dropErrs.append(dropErr)

    ratios.append(dropMean / lightMean)

    tranNames = tranCounts.index[tranCounts['Gene name'] == gene].tolist()

    # This should not print anything, but double checks unique mapping
    if len(tranNames) != 1:
      print('Gene %s has multiple transcript names: %s' % (gene, tranNames))

    tranLength = geneLengths[tranNames[0]]
    lengths.append(tranLength)

  plt.figure()
  ax = plt.gca()

  logLengths = [math.log10(length) for length in lengths]
  norm = matplotlib.colors.Normalize(vmin=min(logLengths), vmax=max(logLengths))
  colors = [COLOR_MAP(norm(logLengths[i])) for i in range(len(lightMeans))]
  zVals = [1.0 / padjVals[gene] for gene in genesToPlot]

  plt.scatter(lightMeans, dropMeans, s=0, alpha=0.7, c=logLengths,
              cmap=COLOR_MAP)

  print('Minimum length: %d nt, Maximum length: %d nt' % \
          (min(lengths), max(lengths)))
  print('Minimum log length: %f, Maximum log length: %f' % \
          (min(logLengths), max(logLengths)))

  for i in range(len(lightMeans)):
    ax.errorbar(lightMeans[i], dropMeans[i], xerr=lightErrs[i],
                yerr=dropErrs[i], fmt='o', ms=0, alpha=0.7, color=colors[i],
                zorder=zVals[i])

  plt.xlabel('Estimated Light-Seq reads per BCL cell')
  plt.ylabel('Estimated Drop-Seq reads per BCL cell')

  xStart = 3e-3
  xEnd = 20
  yStart = 3e-3
  yEnd = 20

  ax.set_xscale('log')
  ax.set_yscale('log')
  plt.xlim([xStart, xEnd])
  plt.ylim([yStart, yEnd])
  plt.colorbar()

  ax.set_aspect('equal')

  plt.savefig('Plots/BCLonly.svg')

  print('%f +/- %f drop-seq vs. light-seq counts' % (statistics.mean(ratios),
                                                     statistics.stdev(ratios)))
  print('  median of %f' % statistics.median(ratios))
  r, p = stats.pearsonr(lightMeans, dropMeans)
  print('pearsonr = %f, p = %f' % (r, p))

  binLightSeq = defaultdict(list)
  binDropSeq = defaultdict(list)
  binRatios = defaultdict(list)

  for i in range(len(lightMeans)):
    for j in range(len(BINS)):
      if BINS[j][0] <= lengths[i] < BINS[j][1]:
        binStr = '%d-%d' % (BINS[j][0], BINS[j][1])
        binLightSeq[binStr].append(lightMeans[i])
        binDropSeq[binStr].append(dropMeans[i])
        binRatios[binStr].append(ratios[i])

  for bin in BINS:
    binStr = '%d-%d' % (bin[0], bin[1])
    lightMeansBin = binLightSeq[binStr]
    dropMeansBin = binDropSeq[binStr]
    ratiosBin = binRatios[binStr]
    r, p = stats.pearsonr(lightMeansBin, dropMeansBin)
    print('  For bin %s' % binStr)
    print('    pearsonr = %f, p = %f' % (r, p))
    print('    %f +/- %f drop-seq vs. light-seq counts' % \
             (statistics.mean(ratiosBin), statistics.stdev(ratiosBin)))
    print('      median of %f' % statistics.median(ratiosBin))

  plt.figure()
  listOfListsRatios = [binRatios[key] for key in ORDERED_BINS]
  ax = sns.boxplot(data=listOfListsRatios)
  ax.set_xticklabels(ORDERED_BINS)
  plt.xticks(rotation = 30)
  plt.xlabel('Transcript length (nt)')
  plt.ylabel('Drop-Seq vs. Light-Seq ratio')
  plt.yscale('log')
  plt.tight_layout()
  plt.savefig('Plots/BCL_Ratios.svg')
  plt.close()

def main():
  enrichLight = pd.read_csv(LIGHT_PADJ_FILE, index_col=0)
  enrichedLightGenes = enrichLight.loc[enrichLight['padj'] < 0.05].index
  lightCellCounts = pd.read_csv(LIGHT_CELL_COUNTS)

  enrichDrop = pd.read_csv(DROP_PADJ_FILE, index_col=0)
  enrichedDropGenes = enrichDrop.loc[enrichDrop['padj'] < 0.05].index
  dropCellCounts = pd.read_csv(DROP_CELL_COUNTS)

  enrichCounts = defaultdict(list)
  for gene in enrichLight.index.union(enrichDrop.index):
    lightStr = 'L not enriched'
    dropStr = 'D not enriched'

    if gene not in enrichLight.index:
      lightStr = 'L not present '
    elif enrichLight['padj'][gene] < 0.05:
      if enrichLight['log2FoldChange'][gene] > 0:
        lightStr = 'L BCL enriched'
      elif enrichLight['log2FoldChange'][gene] < 0:
        lightStr = 'L ONL enriched'
      else:
        print('something wrong')

    if gene not in enrichDrop.index:
      dropStr = 'D not present '
    elif enrichDrop['padj'][gene] < 0.05:
      if enrichDrop['log2FoldChange'][gene] > 0:
        dropStr = 'D BCL enriched'
      elif enrichDrop['log2FoldChange'][gene] < 0:
        dropStr = 'D ONL enriched'
      else:
        print('something wrong')

    enrichCounts['%s, %s' % (lightStr, dropStr)].append(gene)

  print('%d total genes in LS, %d enriched' % (len(enrichLight.index),
                                               len(enrichedLightGenes)))
  print('%d total genes in DS, %d enriched' % (len(enrichDrop.index),
                                               len(enrichedDropGenes)))

  print('\n'.join(['%s: %d' % (cond, len(enrichCounts[cond])) \
                   for cond in enrichCounts]))

  # Include only genes that were significantly enriched in both datasets
  allEnrichedGenes = list(set(enrichedLightGenes) & set(enrichedDropGenes))

  dropFreqs = pd.read_csv(DROP_FREQS, index_col='Gene')
  lightFreqs = pd.read_csv(LIGHT_FREQS, index_col='Gene')

  # Calculate subtraction-based enrichment per replicate
  includedGenesBCL = []
  lightBCL = defaultdict()
  dropBCL = defaultdict()

  includedGenes = []
  lightEnrich = defaultdict()
  dropEnrich = defaultdict()

  for gene in allEnrichedGenes:
    lightValsONL = []
    lightValsBCL = []

    for lightRep in LIGHT_REPS:
      condONL = '%s%s' % ('ONL', lightRep)
      lightValsONL.append(lightFreqs[condONL][gene] / \
                          lightCellCounts[condONL][0])

      condBCL = '%s%s' % ('BP', lightRep)
      lightValsBCL.append(lightFreqs[condBCL][gene] / \
                          lightCellCounts[condBCL][0])

    dropValsONL = []
    dropValsBCL = []

    for dropRep in DROP_REPS:
      condONL = '%s%s' % ('ONL', dropRep)
      dropValsONL.append(dropFreqs[condONL][gene] / dropCellCounts[condONL][0])

      condBCL = '%s%s' % ('BP', dropRep)
      dropValsBCL.append(dropFreqs[condBCL][gene] / dropCellCounts[condBCL][0])

    lightBCL[gene] = lightValsBCL
    dropBCL[gene] = dropValsBCL
    includedGenesBCL.append(gene)

    # Include all genes, even if only show up in one assay
    lightEnrich[gene] = list(np.array(lightValsONL) - np.array(lightValsBCL))
    dropEnrich[gene] = list(np.array(dropValsONL) - np.array(dropValsBCL))
    includedGenes.append(gene)      

  enrichedLightGenesBipolarOnly = enrichLight.loc[enrichLight['padj'] < 0.05] \
                                  .loc[enrichLight['log2FoldChange'] > 0].index
  enrichedDropGenesBipolarOnly = enrichDrop.loc[enrichDrop['padj'] < 0.05] \
                                 .loc[enrichDrop['log2FoldChange'] > 0].index
  genesToPlot = list(set(enrichedLightGenesBipolarOnly) & \
                     set(enrichedDropGenesBipolarOnly))
  assignFile = pd.read_csv(GENE_ASSIGN_FILE, sep='\t', skiprows=1, index_col=0)
  transcriptCounts = pd.read_csv(IN_FILE_FREQ, index_col=0)
  plotBipolarOnly(lightBCL, dropBCL, genesToPlot, assignFile['Length'],
                  transcriptCounts, enrichLight['padj'])

  markers1 = ['Rho', 'Sag', 'Gnat1', 'Gngt1', 'Pdc', 'Gnb1', 'Pde6g', 'Calm1',
              'Pcp2', 'Glul', 'Vsx2', 'Trpm1']
  plotEnrichments(lightEnrich, dropEnrich, includedGenes, markers1, 1,
                  assignFile['Length'], transcriptCounts, enrichLight['padj'])

  markers2 = markers1 + ['Vsx2', 'Cnga1', 'Pde6a', 'Crx', 'Prox1', 'Grm6',
                         'Dpysl3']
  plotEnrichments(lightEnrich, dropEnrich, includedGenes, markers2, 2,
                  assignFile['Length'], transcriptCounts, enrichLight['padj'],
                  xlim=[-1.5, 1.5], ylim=[-1.5, 1.5])

if __name__ == '__main__':
  main()
