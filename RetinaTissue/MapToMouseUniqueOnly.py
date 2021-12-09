import glob
import os

IN_FILES = sorted(glob.glob('inFiles/*_R1.fastq.gz'))
OUT_DIR = 'outFiles/'

MOUSE_INDEX_PATH = '../211101_GRChm39_Indexed'
STAR_PATH = '../STAR-2.7.9a/source/STAR'

def mapMouseSeqs(trimmedR1File, mouseR1Prefix):
    # Map to mouse
    print('  Mapping to mouse genome...')
    os.system(('%s --runThreadN 12 --quantMode TranscriptomeSAM GeneCounts ' \
               'GeneCounts --genomeDir %s --readFilesIn %s ' \
               '--outFileNamePrefix %s --readFilesCommand zcat ' \
               '--outFilterMultimapNmax 1 ' \
               '--outSAMtype BAM SortedByCoordinate') % \
                 (STAR_PATH, MOUSE_INDEX_PATH, trimmedR1File, mouseR1Prefix))

def main():
  parseFiles = []
  for inFile in IN_FILES:
    filePrefix = '%s%s' % (OUT_DIR, inFile.split('/')[-1].split('_')[0])
    trimmedR1File = '%s_R1_trimmed.fastq.gz' % filePrefix
    mouseR1Prefix = '%s_R1_Mouse' % filePrefix
    parseFiles.append((trimmedR1File, mouseR1Prefix))

  for pFiles in parseFiles:
    mapMouseSeqs(*pFiles)

if __name__ == '__main__':
  main()
