#Plot and save the top 200 genes by relative expression level.
import pandas as pd
import seaborn as sns
import numpy as np
import glob
import os

import matplotlib.pyplot as plt

PATH = "results_cellmixing_geneExpression/"
TOP_GENES = 200
TPM_THRESH = 1 #Look at only detected genes for now. Shoulder point for the mean of TPM values is also ~1

def getFileName(file):
	fname = file.split("/")[-1].split("_")[0:2]
	fname = "_".join(fname)
	return fname

#Apply a cutoff to remove zero counts
def tpmCutoff(df):
	newdf = df[
		(df.loc[:, "436D":"436F"]>=TPM_THRESH).all(axis=1)
	]
	
	return newdf

#calculate pearson corr of only the three replicates 436 D, E, F
def corrMatrix(df, **kwargs):
	df_length = kwargs.get("top")
	print("CorrCoeff of: " + str(len(df.iloc[0:df_length])) + " items")
	
	corr = df.iloc[0:df_length] \
			  .loc[:, "436D":"436F"] \
			  .corr(method="pearson")
	
	return corr

def pairGridPlot(df, **kwargs):
	#slice only columns we want
	if kwargs.get("hue"):
		col = ["436D", "436E", "436F", "top_gene"]
	else:
		col = ["436D", "436E", "436F"]
		
	df_plot = df.loc[:, col]
	g = sns.PairGrid(df_plot, hue=kwargs.get("hue"))
	g.map_lower(sns.histplot, cbar=True, bins=50, hue=None)
	g.map_upper(sns.scatterplot, alpha=0.8, s=8)
	g.map_diag(sns.histplot, kde=False, hue=None, bins=50)

def main():

	files = [PATH + "Human_436_meanTPM.csv", PATH + "Mouse_436_meanTPM.csv"]

	for file in files:
		fname = getFileName(file)
		
		tpm_df = pd.read_csv(file, index_col=0)
		tpm_df = tpm_df.sort_values("meanTPM", ascending=False)
		tpm_df_1 = tpmCutoff(tpm_df)
		
		#Save files of selected top genes
		tpm_df_1.iloc[0:TOP_GENES].to_csv(
			PATH + fname + "_top%s.csv" %TOP_GENES
		)
		
		#Save corr matrices for both datasets
		corrMatrix(tpm_df_1, top=TOP_GENES).to_csv(
			PATH + fname + "_top%s_corr.csv" %TOP_GENES
		)
		corrMatrix(tpm_df_1).to_csv(
			PATH + fname + "_corr.csv"
		)
		
		#Mark the Top genes which will be highlighted in the scatterplot
		tpm_df_1["top_gene"]=0
		tpm_df_1.iloc[0:TOP_GENES].loc[:,"top_gene"]=1
		
		#Plot and Save scatterplot
		pairGridPlot(tpm_df_1, hue="top_gene")
		plt.savefig(PATH + fname + "_plot.pdf")

if __name__ == '__main__':
	main()
