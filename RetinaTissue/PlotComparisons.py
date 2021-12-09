from collections import defaultdict 
import pandas as pd
import statistics
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

DROP_FREQS = 'DropSeq/DropSeqLayerFrequencies.csv'
DROP_CELL_COUNTS = 'DropSeq/DropSeqLayerCellCounts.csv'
DROP_PADJ_FILE = 'DropSeq/BPvsONL.csv'
LIGHT_FREQS = 'LightSeq/ReorderedLightSeq.csv'
LIGHT_PADJ_FILE = 'LightSeq/BPvsONL.csv'
LIGHT_CELL_COUNTS = 'LightSeq/LightSeqCellCounts.csv'

DROP_REPS = ['_1', '_2', '_3', '_4', '_5', '_6']
LIGHT_REPS = ['_1', '_2', '_3', '_4']

MAIN_COLOR = np.array([1.0, 0.0, 1.0])

def plotEnrichments(lightEnrich, dropEnrich, includedGenes, includeErrors=False,
                    xlim=None, ylim=None, fignum=None, annotByGene=None):
  # Calculate means and replicate errors for enrichment values
  lightEnrichMeans = []
  lightEnrichErrs = []
  dropEnrichMeans = []
  dropEnrichErrs = []

  for gene in includedGenes:
    lightMean = statistics.mean(lightEnrich[gene])
    lightErr = statistics.stdev(lightEnrich[gene])
    dropMean = statistics.mean(dropEnrich[gene])
    dropErr = statistics.stdev(dropEnrich[gene])

    lightEnrichMeans.append(lightMean)
    lightEnrichErrs.append(lightErr)
    dropEnrichMeans.append(dropMean)
    dropEnrichErrs.append(dropErr)

  plt.figure()
  ax = plt.gca()

  sizes = [50] * len(includedGenes)

  # Plot
  if includeErrors:
    ax.errorbar(lightEnrichMeans, dropEnrichMeans, xerr=lightEnrichErrs,
                yerr=dropEnrichErrs, fmt='o', color='k', ms=0, alpha=0.15,
                zorder=1)
  ax.scatter(lightEnrichMeans, dropEnrichMeans, color='k', alpha=0.15, zorder=2)

  plt.xlabel('Light-Seq ONL - BCL, approx reads per cell')
  plt.ylabel('Drop-Seq ONL - BCL, approx reads per cell')

  if xlim is not None:
    plt.xlim(xlim)
  if ylim is not None:
    plt.ylim(ylim)

  if annotByGene is not None:
    for i in range(len(lightEnrichMeans)):
      if includedGenes[i] in annotByGene:
        if lightEnrichMeans[i] > 0:
          ax.scatter(lightEnrichMeans[i], dropEnrichMeans[i], color='m',
                     zorder=3)
          ax.annotate(includedGenes[i], (lightEnrichMeans[i],
                                         dropEnrichMeans[i]),
                      horizontalalignment='left', verticalalignment='bottom',
                      color='m', zorder=4)
        else:
          ax.scatter(lightEnrichMeans[i], dropEnrichMeans[i], color='c',
                     zorder=3)
          ax.annotate(includedGenes[i], (lightEnrichMeans[i],
                                         dropEnrichMeans[i]),
                      horizontalalignment='left', verticalalignment='bottom',
                      color='c', zorder=4)

  plt.grid(True)
  ax.set_aspect('equal')

  if fignum is not None:
    plt.savefig('Plots/%d.svg' % fignum)

def plotBipolarOnly(lightBCL, dropBCL, genesToPlot, annotList=None):
  # Calculate means and replicate errors for enrichment values
  lightMeans = []
  lightErrs = []
  dropMeans = []
  dropErrs = []

  ratios = []

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

  print('%f +/- %f drop-seq vs. light-seq counts' % (statistics.mean(ratios),
                                                     statistics.stdev(ratios)))
  print('  median of %f' % statistics.median(ratios))

  meanData = pd.DataFrame.from_dict({'Light-Seq': lightMeans,
                                     'Drop-Seq': dropMeans})

  plt.figure()
  ax = plt.gca()

  ax.errorbar(lightMeans, dropMeans, xerr=lightErrs, yerr=dropErrs, fmt='o',
              ms=0, alpha=0.5, color='c', zorder=1)
  ax.scatter(lightMeans, dropMeans, color='c', alpha=0.5, zorder=2)

  for i, gene in enumerate(genesToPlot):
    if annotList is None or gene in annotList:
      ax.annotate(gene, (lightMeans[i], dropMeans[i]),
                  horizontalalignment='left', verticalalignment='bottom')
      ax.scatter(lightMeans[i], dropMeans[i], color='k', zorder=3)

  plt.xlabel('Estimated Light-Seq reads per BCL cell')
  plt.ylabel('Estimated Drop-Seq reads per BCL cell')

  xStart = 3e-3
  xEnd = 20
  yStart = 3e-3
  yEnd = 20

  r, p = stats.pearsonr(meanData['Light-Seq'], meanData['Drop-Seq'])
  print('pearsonr = %f, p = %f' % (r, p))

  ax.set_xscale('log')
  ax.set_yscale('log')
  plt.xlim([xStart, xEnd])
  plt.ylim([yStart, yEnd])

  ax.set_aspect('equal')

  plt.savefig('Plots/BCLonly.svg')

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
        lightStr = 'L ONL enriched'
      elif enrichLight['log2FoldChange'][gene] < 0:
        lightStr = 'L BCL enriched'
      else:
        print('something wrong')

    if gene not in enrichDrop.index:
      dropStr = 'D not present '
    elif enrichDrop['padj'][gene] < 0.05:
      if enrichDrop['log2FoldChange'][gene] > 0:
        dropStr = 'D ONL enriched'
      elif enrichDrop['log2FoldChange'][gene] < 0:
        dropStr = 'D BCL enriched'
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
  markers1 = ['Slitrk2', 'Glul', 'Calm1', 'Vsx2', 'Trpm1', 'Pcp2', 'Gng13',
              'Prox1', 'Grm6', 'Rgs', 'Rgs16', 'Chtnap2']
  plotBipolarOnly(lightBCL, dropBCL, genesToPlot, annotList = markers1)

  markers2 = ['Rho', 'Sag', 'Gnat1', 'Gngt1', 'Pdc', 'Gnb1', 'Pde6g', 'Calm1',
              'Pcp2', 'Glul', 'Vsx2', 'Trpm1']
  plotEnrichments(lightEnrich, dropEnrich, includedGenes, includeErrors=True,
                  annotByGene=markers2, fignum=1)

  markers3 = markers2 + ['Vsx2', 'Cnga1', 'Pde6a', 'Crx', 'Prox1', 'Grm6',
                         'Dpysl3']
  plotEnrichments(lightEnrich, dropEnrich, includedGenes, includeErrors=True,
                  xlim=[-2, 2], ylim=[-2, 2], annotByGene=markers3, fignum=2)


if __name__ == '__main__':
  main()
