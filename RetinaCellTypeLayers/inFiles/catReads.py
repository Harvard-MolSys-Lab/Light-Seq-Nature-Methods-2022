import os
import glob

INDEX_READS_PATH = "../indexedReads/"
CONDITIONS = {'S1':'LS8A',
              'S2':'LS8B',
              'S3':'LS8C',
              'S4':'LS8D',
              'S5':'LS8E',
              'S6':'LS8F',
              'S7':'LS8H',
              }

indexed_reads = sorted(glob.glob(INDEX_READS_PATH + "*R1.fastq.gz"))
parsed_reads = sorted(glob.glob("LS8*R1.fastq.gz"))

for cond in list(CONDITIONS.keys()):

    cond_files = [file for file in indexed_reads if cond in file]
    
    LS_cond = CONDITIONS.get(cond)
    parsed_cond_files = [file for file in parsed_reads if LS_cond in file]
    
    catstring = "cat %s %s %s %s > %s_Merged_R1.fastq.gz" % (
        cond_files[0],
        cond_files[1],
        parsed_cond_files[0],
        parsed_cond_files[1],
        LS_cond
    )

    print(catstring)
    os.system(catstring)