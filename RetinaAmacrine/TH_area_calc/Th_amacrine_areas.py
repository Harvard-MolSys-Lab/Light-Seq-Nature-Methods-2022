import skimage.io
from skimage.measure import regionprops, label
import glob
import os
import csv
from collections import defaultdict


PATH = "Images/Barcode2/"
ROI = sorted(glob.glob(PATH+"*ROI.tif"))

COND = {
    'Slide2-A1':'Th_1',
    'Slide2-A2':'Th_2',
    'Slide2-C1':'Th_3',
    'Slide2-C2':'Th_4',
    'Slide1-C1':'Th_5',
    }


def getCond(ROI):
    cond_list = []
    for roi in ROI:
        key = next(
            (key for key in COND.keys() if key in roi), None
            )
        cond_list.append(COND.get(key))
    
    return cond_list


def calcArea(ROI):
    roi_area = []
    for roi in ROI:

        roi_img = regionprops(
            label(
                skimage.io.imread(roi)
                )
            )
        
        fov_area = []
        for roi in roi_img:
            fov_area.append(roi.area)
        
        roi_area.append(fov_area)
    
    return roi_area


roi_area = calcArea(ROI)
cond_list = getCond(ROI)

roi_dict = defaultdict(list)
for cond,area in list(zip(cond_list, roi_area)):
    roi_dict[cond].extend(area)

roi_dict.pop(None, None) #remove the None condition

with open('Th_amacrine_areas.csv', 'w') as f:
    writer = csv.writer(f, delimiter=',')
    
    for key in sorted(roi_dict.keys()):
        row = [key] + roi_dict.get(key)
        print(row)
        writer.writerow(row)
        
    #Sum up areas
    for key in sorted(roi_dict.keys()):
        row = [key] + [sum(roi_dict.get(key))]
        print(row)
        writer.writerow(row)
