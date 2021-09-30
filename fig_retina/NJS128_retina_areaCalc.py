import skimage.io
from skimage.measure import regionprops, label
import glob
import os
import csv

import pandas as pd


PATH = "NJS128_retina_HX7_nikon images/"

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

#save items
df = pd.DataFrame(index=["area1","area2","area3"],columns=wells)
df["A3"] = roi_area[0:3]
df["B2"] = roi_area[3:6]
df["B3"] = roi_area[6:9]
df["C3"] = roi_area[9:12]

print("Current dataframe")
print(df)

df.to_csv("NJS128_retina_areas.csv")
