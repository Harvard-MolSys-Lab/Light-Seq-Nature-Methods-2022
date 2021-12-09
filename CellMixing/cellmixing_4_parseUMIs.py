import os
import glob

import pandas as pd
import pysam

PATH_BAMFILES = "outFiles/"
PATH_NUMBERS = "results_cellmixing_numbers/"

HUMAN = sorted(glob.glob('%s*Human_Dedup.bam' % PATH_BAMFILES))
MOUSE = sorted(glob.glob('%s*Mouse_Dedup.bam' % PATH_BAMFILES))

# Remove the neg controls
HUMAN = HUMAN[1:]
MOUSE = MOUSE[1:]


def parse_umis():
    samfiles = [pysam.AlignmentFile(file, "rb") for file in HUMAN+MOUSE]
    umis = [samfile.mapped for samfile in samfiles]

    return umis


def set_umi_df(umis):
    data = pd.DataFrame(umis, 
    columns = ["Umis"], 
    index=[
        "Human_1", 
        "Human_2", 
        "Human_3",
        "Human_deepseq", 
        "Mouse_1", 
        "Mouse_2", 
        "Mouse_3", 
        "Mouse_deepseq", 
        ]
    )

    return data


def main():
    umis = parse_umis()
    umi_df = set_umi_df(umis)

    umi_df.to_csv(PATH_NUMBERS + "Cell_umi_numbers.csv")


if __name__ == '__main__':
    main()
    