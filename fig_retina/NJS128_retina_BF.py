import skimage.morphology
import skimage.io
import skimage.transform
from skimage.filters import threshold_otsu
import numpy as np

import matplotlib.pyplot as plt

PATH = "NJS128_retina_HX7_nikon images/"

#The Nikon camera is flipped compared to the ImgExpress camera, so flip the image to match them.
def rotate180(arr):
    arr180 = skimage.transform.rotate(arr, 180, preserve_range=True)
    return arr180

def normalizeBF(arr):
    thresh = threshold_otsu(arr)
    pxforeground = arr[arr>thresh]
    
    imgnorm = (arr-arr.min()) / (np.percentile(pxforeground, 99) - arr.min())
    
    return imgnorm

def normalizeBF95(arr):
    thresh = threshold_otsu(arr)
    pxforeground = arr[arr>thresh]

    imgnorm = (arr-arr.min()) / (np.percentile(pxforeground, 95) - arr.min())
    
    return imgnorm

#Figure uses BF from well B2
b2BF = skimage.io.imread(PATH + 'B2-BF.tif')
b2BF = normalizeBF95(b2BF)
b2BF = rotate180(b2BF)
b2BF = skimage.color.gray2rgb(b2BF)
plt.imshow(b2BF)
plt.axis("off")
plt.savefig('B2-1.png', dpi=1000, bbox_inches='tight', transparent='True', pad_inches=0)

b2BF = skimage.io.imread(PATH + 'B2-Area2-BF.tif')
b2BF = normalizeBF(b2BF)
b2BF = rotate180(b2BF)
b2BF = skimage.color.gray2rgb(b2BF)
plt.imshow(b2BF)
plt.axis("off")
plt.savefig('B2-2.png', dpi=1000, bbox_inches='tight', transparent='True', pad_inches=0)

b2BF = skimage.io.imread(PATH + 'B2-Area3-BF.tif')
b2BF = normalizeBF(b2BF)
b2BF = rotate180(b2BF)
b2BF = skimage.color.gray2rgb(b2BF)
plt.imshow(b2BF)
plt.axis("off")
plt.savefig('B2-3.png', dpi=1000, bbox_inches='tight', transparent='True', pad_inches=0)