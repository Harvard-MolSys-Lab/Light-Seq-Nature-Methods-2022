# Subsampling one condition from the TH amacrine experiment
# Plot the Th and Cartpt gene as seperate graphs.

import os
import glob
import random

import numpy as np
import pysam

FRAC = np.array(
    [0.005, 0.01, 0.025, 0.05, 0.075, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 0.9]
)
SIMS = 5
OUTFOLDER = 'simOut/'

samfile = pysam.AlignmentFile("TLS23A_Sorted.bam", "rb")


def subsim(frac, sim):
    subsamfile = pysam.AlignmentFile(
        "%ssubsampled_F%s_S%s.bam" % (OUTFOLDER, str(frac), str(sim)), 
        "wb", 
        template=samfile
    )
    
    for read in samfile.fetch():
        if random.random()<=frac:
            subsamfile.write(read)
    
    subsamfile.close()


def main():
    for frac in FRAC:
        print("Simulating FRAC: %s" % str(frac))
        for sim in range(SIMS):
            print("Simulating SIM: %s" % str(sim))
            subsim(frac, sim)


if __name__ == '__main__':
    main()
