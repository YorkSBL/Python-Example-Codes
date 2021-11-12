#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EXearFingerprint.py

[IN PROGRESS]
Purpose: Attempt to create a Python-y version of EXhaircells2B.m to 
analyze the Anolis SEM image and relate back to tonotopic map

o Requires file ACsb5LhcMasked.jpg
o Borrow some image processing syntax from Hill (2015) ch.6.8.2
o Using scipy & Pillow packages to deal w/ image processing; installed via
  https://anaconda.org/anaconda/pillow)
 https://anaconda.org/anaconda/scipy
o Using code bits from here too:
  https://www.geeksforgeeks.org/working-images-python/
  https://stackoverflow.com/questions/14263050/segment-an-image-using-python-and-pil-to-calculate-centroid-and-rotations-of-mul

Created on Wed Jun 30 09:49:59 2021
@author: CB
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import pylab
from PIL import Image
import scipy
from scipy import ndimage
from matplotlib import image

# -------------------------------- 
root= '../Files/'         # root path to file
fileN = 'ACsb5LhcMasked.jpg'   # file name

# -------------------------------- 

# ==== bookeeping
fname = os.path.join(root,fileN)
img  = Image.open(fname)
imgA = image.imread(fname)  # imports as an array

"""
# ==== [borrowed bits are below to do segmentation]
im = np.where(img > 128, 0, 1)
label_im, num = ndimage.label(im)
slices = ndimage.find_objects(label_im)
centroids = ndimage.measurements.center_of_mass(im, label_im, xrange(1,num+1))

angles = []
for s in slices:
    height, width = label_im[s].shape
    opp = height - np.where(im[s][:,-1]==1)[0][-1] - 1
    adj = width - np.where(im[s][-1,:]==1)[0][0] - 1
    angles.append(np.degrees(np.arctan2(opp,adj)))
#print 'centers:', centroids
#print 'angles:', angles
"""

# ====
pylab.imshow(img)