import os, sys
parentdir = os.path.dirname(os.path.realpath("."))
sys.path.append(parentdir)
print(parentdir)

import pysam
import pandas as pd
import numpy as np
from utils.utils import *

####
bamfile = sorted(globBAMfiles(type="dedup"))
#print("Deduped and indexed bam files: ")
#print(bamfile)

umis = []
for file in bamfile:
    samfile = pysam.AlignmentFile(file, "rb")
    umis.append(samfile.mapped)

## get fraction values ##
temp = []
for file in bamfile:
    temp.append(file.split("_")[1])
    fracs = sorted(list(set(temp)))

umis_arr = np.array(umis)

#reshape array to match dataframe
d = np.transpose(umis_arr.reshape(12,5))

data = pd.DataFrame(d, columns = fracs)
print(data)

data.to_csv("subsampled.csv")