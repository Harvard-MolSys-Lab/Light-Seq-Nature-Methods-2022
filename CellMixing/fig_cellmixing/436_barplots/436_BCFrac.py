import pandas as pd
import seaborn as sns
import numpy as np
import glob
import os
import csv

import matplotlib.pyplot as plt

PATH = "../../results_cellmixing_numbers/"

CSV = sorted(glob.glob(PATH+'*NUMBERS.csv'))

df_436 = []
for file in CSV[1:4]:
    df_436.append(pd.read_csv(file, index_col=0))

#Enter in the number for eGFP reads, order is experiment D, E, F
eGFP_bc2 = np.array([300, 811, 546])
eGFP_bc1 = np.array([21, 67, 41])

def calcFrac(data):
    h_bc2 = data.loc["Human_BC2"] / data.loc["Human_TotalReads"]
    m_bc1 = data.loc["Mouse_BC1"] / data.loc["Mouse_TotalReads"]
    
    return h_bc2, m_bc1

#construct a new dataframe for plotting scatterplots of fraction of correct reads
hm_df = pd.DataFrame(columns = ["Human BC2", "Mouse BC1", "Mapped Read"], index=range(9))
hm_df.index.name = "Fraction of Reads"

#construct the dataframe
for idx,df in enumerate(df_436):
    h_bc2, m_bc1 = calcFrac(df)

    hm_df.loc[idx, "Human BC2"] = h_bc2[0]
    hm_df.loc[idx, "Mouse BC1"] = 1-h_bc2[0]
    hm_df.loc[idx, "Mapped Read"] = "Human"
    
    hm_df.loc[idx+3, "Mouse BC1"] = m_bc1[0]
    hm_df.loc[idx+3, "Human BC2"] = 1-m_bc1[0]
    hm_df.loc[idx+3, "Mapped Read"] = "Mouse"
    
    hm_df.loc[idx+6, "Human BC2"] = eGFP_bc2[idx] / (eGFP_bc2[idx] + eGFP_bc1[idx])
    hm_df.loc[idx+6, "Mouse BC1"] = eGFP_bc1[idx] / (eGFP_bc2[idx] + eGFP_bc1[idx])
    hm_df.loc[idx+6, "Mapped Read"] = "eGFP"

#Save df
hm_df.to_csv("436_BCFraction.csv")

#plot and save figure
fig_dims = (4, 4)
fig, ax = plt.subplots(figsize=fig_dims)
sns.scatterplot(data=hm_df, x="Human BC2", y="Mouse BC1", hue="Mapped Read", alpha=0.5, s=300)
ax.set(ylim=(0,1))
ax.set(xlim=(0,1))
plt.savefig("436_BCFraction.pdf")