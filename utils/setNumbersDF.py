import pandas as pd

def getNumbers(human_maps, mouse_maps, total_reads):
	df = [
		len(human_maps.bc1),
		len(human_maps.bc2),
		len(mouse_maps.bc1),
		len(mouse_maps.bc2),
		len(human_maps.bc1) + len(human_maps.bc2),
		len(mouse_maps.bc1) + len(mouse_maps.bc2),
		total_reads
		]

	return df

def setDF(human_maps, mouse_maps, total_reads):
	df = getNumbers(human_maps, mouse_maps, total_reads)

	data = pd.DataFrame(df, 
	columns = [human_maps.fname_short], 
	index=[
		"Human_BC1", 
		"Human_BC2", 
		"Mouse_BC1", 
		"Mouse_BC2", 
		"Human_TotalReads", 
		"Mouse_TotalReads", 
		"TotalReads"
		]
	)

	return data
