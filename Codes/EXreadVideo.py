#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EXreadVideo.py

Purpose: Read in a video file and convert to array for sectioning/analysis

o Uses an Anolis Digimorph video as the source
o Reqs. the skvideo Python package (and does not use OpenCV); REF:
https://stackoverflow.com/questions/42163058/how-to-turn-a-video-into-numpy-array
o Note that the dimensions of the video file yield (I think):
    sag: 200 slices
    horiz: 192 slices
    coronal: 440 slices

NOTE: this code uses the package skvideo

Created on Thu Sep 30 14:33:56 2021
@author: C. Bergevin
"""

import skvideo.io  
import numpy as np
import matplotlib.pyplot as plt

# ==================================
# [User Params]
fileN= "AnolisSag.mp4"
slcCnum= 200   # slice # (i.e., indx) re CORONAL slice {250}
slcHnum= 157   # slice # (i.e., indx) re HORIZONTAL slice {157}
slcSnum= 107   # slice # (i.e., indx) re SAGITTAL slice {147}
# ==================================

# --- load in the data
D = skvideo.io.vread(fileN)  
# --- extract representative slices
sliceC= D[:,:,slcCnum,1]  # coronal slice
sliceH= D[:,slcHnum,:,1]  # horizontal slice
sliceS= D[slcSnum,:,:,1]  # sagittal slice
# --- visualize (kludge!)
#fig1, ax1 = plt.subplots(2,2)
#fig1= plt.imshow(sliceS)
#ax1.grid(False)
plt.close("all")
fig1 = plt.figure()
# -- coronal
fig1.add_subplot(3,1,1)
plt.imshow(sliceC,cmap='bone')
plt.title("Coronal")
ax = plt.gca()
ax.axes.xaxis.set_visible(False)
ax.axes.yaxis.set_visible(False)
plt.grid(False)
# -- horizontal
fig1.add_subplot(3,1,2)
plt.imshow(sliceH,cmap='viridis')
plt.title("Horizontal")
ax = plt.gca()
ax.axes.xaxis.set_visible(False)
ax.axes.yaxis.set_visible(False)
plt.grid(False)
# -- saggital
fig1.add_subplot(3,1,3)
plt.imshow(sliceS,cmap='jet')
plt.title("Sagittal")
ax = plt.gca()
ax.axes.xaxis.set_visible(False)
ax.axes.yaxis.set_visible(False)
plt.grid(False)
cbar = plt.colorbar()


# ===== save file?
if 1==0:
    np.save('Zena.npy',D)
    




