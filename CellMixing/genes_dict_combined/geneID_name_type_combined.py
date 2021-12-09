import csv
import re

GFF_FILE = "gencode.v38andvM27.annotation.gff3"

#save gene_types
genes_types={}
with open(GFF_FILE) as gff:
    
    #first line of file
    dupdate = {"gene_id" : "gene_name" + ',' + "gene_type"}
    genes_types.update(dupdate)

    for line in gff:
        if line[0] != '#' and 'gene_id=' in line:
            gene_id = re.findall("(?<=gene_id=)[^;]*",line)[0]
            gene_name = re.findall("(?<=gene_name=)[^;]*",line)[0]
            gene_type = re.findall("(?<=gene_type=)[^;]*",line)[0]
            
            dupdate = {gene_id : gene_name + ',' + gene_type}
            genes_types.update(dupdate)

with open("genesDict_types_combined.csv", "w") as f: 
    for key in genes_types.keys(): 
        f.write("%s,%s\n" % (key , genes_types[key]))