import glob
import os

IN_FILES = sorted(glob.glob('inFiles/*_R1_*.fastq.gz'))
OUT_DIR = 'outFiles/'

MERGED_INDEX_PATH = '../211102_GRCh38andGRCm39_Indexed'
HUMAN_INDEX_PATH = '../211101_GRCh38_Indexed'
MOUSE_INDEX_PATH = '../211101_GRChm39_Indexed'
STAR_PATH = '../STAR-2.7.9a/source/STAR'


def mapSeqs(trimmedR1File, fileR1Prefix, indexPath):
  print('  Mapping to merged genome...')
  os.system(('%s --runThreadN 12 --quantMode TranscriptomeSAM GeneCounts ' \
             'GeneCounts --genomeDir %s --readFilesIn %s ' \
             '--outFileNamePrefix %s --readFilesCommand zcat ' \
             '--outFilterMultimapNmax 1 ' \
             '--outSAMtype BAM SortedByCoordinate') % \
               (STAR_PATH, indexPath, trimmedR1File, fileR1Prefix))

def main():
  parseFilesMerged = []
  parseFilesHuman = []
  parseFilesMouse = []

  for inFile in IN_FILES:
    filePrefix = '%s%s' % (OUT_DIR, inFile.split('/')[1].split('_R1_')[0])
    trimmedR1File = '%s_R1_trimmed.fastq.gz' % filePrefix

    mergedR1Prefix = '%s_R1_Merged' % filePrefix
    humanR1Prefix = '%s_R1_Human' % filePrefix
    mouseR1Prefix = '%s_R1_Mouse' % filePrefix

    parseFilesMerged.append((trimmedR1File, mergedR1Prefix, MERGED_INDEX_PATH))
    parseFilesHuman.append((trimmedR1File, humanR1Prefix, HUMAN_INDEX_PATH))
    parseFilesMouse.append((trimmedR1File, mouseR1Prefix, MOUSE_INDEX_PATH))

  for pFiles in parseFilesMerged:
    mapSeqs(*pFiles)

  for pFiles in parseFilesHuman:
    mapSeqs(*pFiles)
    
  for pFiles in parseFilesMouse:
    mapSeqs(*pFiles)

if __name__ == '__main__':
  main()
