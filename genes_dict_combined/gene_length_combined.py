# create a combined gene dictionary from the summary file from featurecounts

import csv

FILE = "436E_R1_MergedAligned.sortedByCoord.out.gene_assigned"  # Any of the gene_assigned files will be fine

print("Reading only gene id and gene length columns from: " + FILE)

genes = {}
with open(FILE) as file:
    for idx, lines in enumerate(file):
        if idx < 1:  # skip first line
            pass
        else:
            dupdate = {lines.split('\t')[0] : lines.split('\t')[5]}
            genes.update(dupdate)

with open("gene_length_combined.csv", "w") as f: 
    for key in genes.keys(): 
        f.write("%s,%s\n" % (key, genes[key]))
