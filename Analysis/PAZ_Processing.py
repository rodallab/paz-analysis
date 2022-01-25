'''
Assorted simple image processing functions
'''
import numpy as np
import tifffile as tif
from scipy.ndimage.filters import gaussian_filter as gblur
from skimage.filters import threshold_li as li
from scipy.ndimage.morphology import binary_fill_holes as fill_holes 
from scipy.ndimage.morphology import distance_transform_edt as edt


def segment(img, sigma):
    img_norm = img.astype(np.float32)
    for i, chan in enumerate(img_norm):
        img_norm[i]=chan/chan.mean()
    img_sum = img_norm.sum(axis=0)
    # Blur, threshold, and fill holes
    mask = gblur(img_sum, sigma)
    thresh = li(mask)
    mask = mask>thresh
    mask = fill_holes(mask)
    return mask

def sdt(mask):
    edm = edt(mask)
    sdm = edm -1
    #sdm[~mask]=-1
    return sdm

def norm_by_mask(img, mask=None, sigma=1):
    if type(mask)==type(None):
        mask = segment(img, sigma)
    fg_pix = img[:,mask>0]
    means = fg_pix.mean(axis=1)
    return np.moveaxis(img, 0, -1)/means