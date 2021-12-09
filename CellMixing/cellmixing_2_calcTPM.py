import csv
import glob
import os

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt


# Note: Analyze 436 and 438 seperately as they use different polymerases and 
# starting material for amplificaiton, could affect analysis.
PATH_CSV = "results_cellmixing_genes/"
PATH_RESULTS = "results_cellmixing_geneExpression/"

# Take only the 436D, E, F files for analysis
HUMAN_436 = [
    PATH_CSV + "436D_Human_bc2_AllCounts.csv",
    PATH_CSV + "436E_Human_bc2_AllCounts.csv",
    PATH_CSV + "436F_Human_bc2_AllCounts.csv"
]
MOUSE_436 = [
    PATH_CSV + "436D_Mouse_bc1_AllCounts.csv",
    PATH_CSV + "436E_Mouse_bc1_AllCounts.csv",
    PATH_CSV + "436F_Mouse_bc1_AllCounts.csv"
]

if os.path.isdir("results_cellmixing_geneExpression") is False:
    os.mkdir("results_cellmixing_geneExpression")


class TPM:
    """
    Set of functions to calc TPM
    Transcripts per kilobase per million reads
    """
    def __init__(self, file):
        self.file = file
        self.filename = (file.split('/')[-1])
        self.fname_short = self.filename.split("_")[0]
        self.data = pd.read_csv(file)
        self.data.index.name = self.fname_short
        self.getSpecies()
        
    def getSpecies(self):
        if "ENSG" in self.data["gene_id"][0]:
            self.species = "_human"
        if "ENSMUSG" in self.data["gene_id"][0]:
            self.species = "_mouse"
    
    # Reads per kilobase
    def RPK(self):
        self.data["RPK"]=0
        rpk = self.data["gene_count"] / (self.data["gene_length"]/1000)
        self.data["RPK"] = rpk
        
    def TPM(self):
        self.RPK()
        RPKsum = self.data["RPK"].sum()
        
        self.data["TPM"] = 0
        tpm = (self.data["RPK"] / RPKsum) * 1e6
        self.data["TPM"] = tpm
    
    def log2TPM(self):
        self.TPM()
        self.data["log2(TPM+1)"] = 0
        
        logtpm = np.log2(self.data["TPM"] + 1)
        self.data["log2(TPM+1)"] = logtpm


# Take log2TPM values and calc the mean, std of the 3 sequencing replicates
def calcMeanTPM(df):
    df["meanTPM"] = 0
    df["meanTPM"] = df.loc[:, "436D":"436F"].mean(axis=1)
    df["stdTPM"] = 0
    df["stdTPM"] = df.loc[:, "436D":"436F"].std(axis=1)
    
    return df


def main():
    # Create a new dataframe to index the values from the three replicates
    # Keep the gene_id, name and type columns
    
    temp_df = pd.read_csv(HUMAN_436[0])
    human_436_tpm = temp_df.loc[:, "gene_id":"gene_type"]
    human_436_tpm.index.name = "Human_436_tpm"

    temp_df = pd.read_csv(MOUSE_436[0])
    mouse_436_tpm = temp_df.loc[:, "gene_id":"gene_type"]
    mouse_436_tpm.index.name = "Mouse_436_tpm"

    # Add new col names to the summary dataframe
    col_names = ["436D", "436E", "436F"]
    
    for idx,file in enumerate(HUMAN_436):
        human_geneExp = TPM(file)
        human_geneExp.log2TPM()

        human_geneExp.data.to_csv(
            PATH_RESULTS 
            + human_geneExp.fname_short 
            + human_geneExp.species 
            + ".csv"
            )

        print("Saving: " + human_geneExp.fname_short + human_geneExp.species)
        
        human_436_tpm[col_names[idx]] = 0
        human_436_tpm[col_names[idx]] = human_geneExp.data["log2(TPM+1)"]    
    
    for idx,file in enumerate(MOUSE_436):
        mouse_geneExp = TPM(file)
        mouse_geneExp.log2TPM()

        mouse_geneExp.data.to_csv(
            PATH_RESULTS 
            + mouse_geneExp.fname_short 
            + mouse_geneExp.species 
            + ".csv"
            )
        
        print("Saving: " + mouse_geneExp.fname_short + mouse_geneExp.species)
        
        mouse_436_tpm[col_names[idx]] = 0
        mouse_436_tpm[col_names[idx]] = mouse_geneExp.data["log2(TPM+1)"]

    human_436_tpm = calcMeanTPM(human_436_tpm)
    human_436_tpm.to_csv(PATH_RESULTS + "Human_436_meanTPM.csv")

    mouse_436_tpm = calcMeanTPM(mouse_436_tpm)
    mouse_436_tpm.to_csv(PATH_RESULTS + "Mouse_436_meanTPM.csv")


if __name__ == '__main__':
    main()

