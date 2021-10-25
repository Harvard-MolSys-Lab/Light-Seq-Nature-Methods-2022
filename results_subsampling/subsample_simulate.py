import pysam
import numpy as np
import gzip
import re
import random

FRAC = np.array([0.005, 0.01, 0.025, 0.05, 0.075, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 0.9])
SIMS = 5

samfile = pysam.AlignmentFile("S436E_R1_MergedAligned.sortedByCoord.out_sorted.bam", "rb")

#Build name index for the reference file, gonna take up 16 gigs
print("Indexing the alignment file...")
nameindex = pysam.IndexedReads(samfile)
nameindex.build()

def qname(var):
    name = re.findall(
        "(?<=@).*",str(var)
    )\
    [0].split(" ")[0]
    
    return name

def subSim(frac, sim):
    subsamfile = pysam.AlignmentFile("subsampled_F%s_S%s.bam" % (str(frac), str(sim)), "wb", template=samfile)
    with gzip.open("S436E_headers.fastq.gz") as header:
        for idx,line in enumerate(header):
            if random.random()<=frac:
                try:
                    if nameindex.find(qname(line)):
                        iter = nameindex.find(qname(line))
                        subsamfile.write(next(iter))
                except KeyError:
                    pass

    subsamfile.close()

for frac in FRAC:
    print("Simulating FRAC: %s" %str(frac))
    for sim in range(SIMS):
        print("Simulating SIM: %s" %str(sim))
        subSim(frac, sim)