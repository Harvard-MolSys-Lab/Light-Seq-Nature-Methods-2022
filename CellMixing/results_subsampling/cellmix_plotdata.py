import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

REAL_DATA_X = 1
REAL_DATA_Y = 529028
FRAC = np.array(
    [0.005, 0.01, 0.025, 0.05, 0.075, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 0.9]
)

df = pd.read_csv('cellmix_subsampled.csv', index_col=0)

# Plotting
plt.errorbar(FRAC, df.mean(), yerr=df.std(), linestyle='-', marker='.', markersize=15)
plt.xlabel('Fraction of Reads')
plt.ylabel('UMIs')
plt.plot(REAL_DATA_X, REAL_DATA_Y, marker='.', markersize=15, markerfacecolor='m')
plt.savefig('cellmixing_subsampled.pdf')