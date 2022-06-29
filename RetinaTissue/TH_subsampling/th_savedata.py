import os, sys
import pysam
import pandas as pd
import numpy as np
from utils import *

TH_GENE     = "ENSMUSG00000000214"
CARTPT_GENE = "ENSMUSG00000021647"

OUTFOLDER = 'simOut/'

####
bamfile = sorted(globBAMfiles(OUTFOLDER, type="dedup"))


## get umis and fraction values ##
umis = []
fracs = []
th_gene = []
cart_gene = []

for file in bamfile:
    
    th_count = 0
    cart_count = 0

    samfile = pysam.AlignmentFile(file, "rb")
    umis.append(samfile.mapped)
    fracs.append(file.split("_")[1])

    for read in samfile.fetch():
        if TH_GENE in read.get_tag('XT'):
            th_count+=1
    th_gene.append(th_count)

    for read in samfile.fetch():
        if CARTPT_GENE in read.get_tag('XT'):
            cart_count+=1
    cart_gene.append(cart_count)

fracs = sorted(list(set(fracs)))
umis_arr = np.array(umis)
th_gene = np.array(th_gene)
cart_gene = np.array(cart_gene)

#reshape array to match dataframe
d = np.transpose(umis_arr.reshape(12,5))
dth = np.transpose(th_gene.reshape(12,5))
dcart = np.transpose(cart_gene.reshape(12,5))

data = pd.DataFrame(d, columns = fracs)
print(data)

data.to_csv("subsampled.csv")

data = pd.DataFrame(dth, columns = fracs)
data.to_csv("th_sub.csv")

data = pd.DataFrame(dcart, columns = fracs)
data.to_csv("cart_sub.csv")
