# Script to generate the sunburst plots
# Note, plotly plots and saving the files requires seperate python packages that can interfere with LSEnv
# Thus this script was run with a new "LS-plots" environment

import plotly.express as px
import plotly.io as pio
import pandas as pd
import glob

CSV = sorted(glob.glob('*.csv'))

for file in CSV:
	
	df = pd.read_csv(file)
	file_prefix = file.split('_')[0]
	
	fig = px.sunburst(df, 
		path=['Barcode', 'Mapping', 'Assignment', 'RNA type','UMI'], 
		values='Reads',
		color='Mapping'
	)
	
	print("Plotting...%s" % file_prefix)
	pio.write_image(fig, "%s_sunburst.pdf" % file_prefix, width=1024, height=1024, scale=1)