import skimage.io
from skimage.measure import regionprops, label
import glob
import os
import csv

import pandas as pd


PATH = "nikon_BF_ROI/"

ROI = sorted(glob.glob(PATH+"*ROI.tif"))

def parseWells(ROI):
    well = []
    for roi in ROI:
        fname = roi.split("/")[-1]
        well.append(fname.split("-")[0])
    
    well = sorted(set(well))
    return well

def calcArea(ROI):
    roi_area = []
    for roi in ROI:
        roi_temp = skimage.io.imread(roi)
        roi_temp = label(roi_temp)
        roi_temp = regionprops(roi_temp)
        roi_temp = roi_temp[0].area
        
        roi_area.append(roi_temp)
    
    return roi_area


wells = parseWells(ROI)
roi_area = calcArea(ROI)

#chunk the list into increments of 3 so they match with the table
roi_area_list = [roi_area[idx:idx+3] for idx in range(0,len(roi_area),3)]

#save items
df = pd.DataFrame(index=["BP","ONL","RGC"],columns=wells)

for idx, well in enumerate(wells):
    df[well] = roi_area_list[idx]

print("Current dataframe")
print(df)

df.to_csv("retina_areas.csv")
