import os
import glob

from utils import *

OUTFOLDER = 'simOut/'

### 
filelist = globBAMfiles(OUTFOLDER)

for file in filelist:
    sortedfile = sortBAM(file)
    indexBAM(sortedfile, type="sorted")
    umiDedup(sortedfile)

#### Index deduped bam files for analysis (if needed)

dedupfiles = globBAMfiles(OUTFOLDER, type="dedup")

for file in dedupfiles:
    indexBAM(file, type="dedup")