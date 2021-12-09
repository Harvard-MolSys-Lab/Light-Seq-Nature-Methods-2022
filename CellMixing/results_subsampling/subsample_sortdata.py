import os, sys
parentdir = os.path.dirname(os.path.realpath("."))
sys.path.append(parentdir)

from utils.utils import *

### 
filelist = globBAMfiles()

for file in filelist:
    sortedfile = sortBAM(file)
    indexBAM(sortedfile, type="sorted")
    umiDedup(sortedfile)

#### Index deduped bam files for analysis (if needed)

dedupfiles = globBAMfiles(type="dedup")

for file in dedupfiles:
    indexBAM(file, type="dedup")