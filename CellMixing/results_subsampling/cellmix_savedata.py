import os
import glob

import pysam
import pandas as pd
import numpy as np


def globBAMfiles(path, **kwargs):
    if kwargs.get("type"):
        filelist = glob.glob(
            "%s*_%s.bam" % (path, kwargs.get('type'))
        )
    else:
        filelist = glob.glob(
            "%s*[!dedup][!sorted][!featurecounts].bam" % path
        )
    return filelist

OUTFOLDER = 'simOut/'

####
bamfile = sorted(globBAMfiles(OUTFOLDER, type="dedup"))

## get umis and fraction values ##
umis = []
fracs = []

for file in bamfile:

    samfile = pysam.AlignmentFile(file, "rb")
    umis.append(samfile.mapped)
    fracs.append(file.split("_")[1])

fracs = sorted(list(set(fracs)))
umis_arr = np.array(umis)

#reshape array to match dataframe
d = np.transpose(umis_arr.reshape(12,5))

data = pd.DataFrame(d, columns = fracs)
print(data)

data.to_csv("cellmix_subsampled.csv")
