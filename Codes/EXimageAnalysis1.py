#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EXimageAnalysis1.py

Purpose: Read in an image file, dither (i.e., add some noise), and 
then use "linear filtering for [...] denoising" as per the example in
Kutz (2013) ch.14.2

o Useful refs:
https://www.geeksforgeeks.org/reading-images-in-python/
http://scipy-lectures.org/intro/scipy/auto_examples/solutions/plot_fft_image_denoise.html
https://scipython.com/book/chapter-6-numpy/examples/blurring-an-image-with-a-two-dimensional-fft/
https://stackoverflow.com/questions/59812522/fourier-transform-inverse-fourier-transform-python

NOTE: this code uses the package scikit-image

Created on Sun Oct 17 12:16:41 2021
@author: C. Bergevin
"""

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
import random
import math
from math import e
from skimage import color
from skimage import io
from scipy import fftpack
from matplotlib.colors import LogNorm

# ----------- [User Params]
fileN= './Images/catINbox.jpg'
#fileN= './Images/grating1.png'
noiseF= 0.1   # noise factor {0.1?}
bs= 0.0001   # Gaussian filter (SD) factor

# ----------- 
  
# ---- read in file and convery to grayscale
#img = mpimg.imread(fileN)
#img= img[:,:,0]
img = io.imread(fileN)
imgG = color.rgb2gray(img)
M= img.shape[0]  # numb. or rows
N= img.shape[1]  # numb. or columns
# ---- create noisy version
noise = np.random.normal(0,noiseF,imgG.shape)
imgN= imgG+ noise

# ==== deal w/ FFT
imgNspec = fftpack.fft2(imgN)
# --- now shift the zero-frequency component to spectrum's center 
imgNspec = np.fft.fftshift(imgNspec)
imgNspecM= np.abs(imgNspec)

# ==== deal w/ filtering (adapted from Matlab via Kutz ch.14.2)
# --- first create the "Gaussian filter"
x = np.linspace(0,N-1,N)
y = np.linspace(0,M-1,M)
xv,yv = np.meshgrid(x,y)
# --- now use those arrays to create the filter via a 2-D Gaussian
F= e**(-bs*(xv-np.round(N/2))**2 -bs*(yv-np.round(M/2))**2)
# --- now apply the filter (via multip. in spect. domain)
#filtM= imgNspecM*F  # apply filter
#filtT= filtM*e**(1j*np.angle(imgNspec))   # repackage (incl. phase)
#imgF= np.fft.ifft2(filtT)  # convert back to spatial domain
# --- alt. vers. (which seems to work!)
filtT= imgNspec*F  # apply filter
fftx = np.fft.ifftshift(filtT)
ffty = np.fft.ifft2(fftx)
imgF = np.abs(ffty)

# ========= VISUALIZE
plt.close("all")
# ---- display image
if 1==0:
    fig1= plt.figure()
    plt.imshow(imgN,cmap='bone')
    plt.grid(False)
    cbar = plt.colorbar()

# ---- display spectal bits
fig2,ax2 = plt.subplots(2,2)
# -- noisy img.
im1= ax2[0,0].imshow(imgN,cmap='bone')
ax2[0,0].set_title('Noisy original')
plt.colorbar(im1,ax=ax2[0,0])
ax2[0,0].axis('off')
# -- |FFT| of noisy img.
im2=ax2[0,1].imshow(imgNspecM,norm=LogNorm(vmin=5))
ax2[0,1].grid(False)
ax2[0,1].axes.xaxis.set_visible(False)
ax2[0,1].axes.yaxis.set_visible(False)
#ax2[0,1].set_title('|FFT| of noisy original (Shifted; Zero freq @ center)')
ax2[0,1].set_title('|FFT| of noisy original')
plt.colorbar(im2,ax=ax2[0,1])
# -- filter
im3=ax2[1,0].imshow(np.log(F))
ax2[1,0].axis('off')
ax2[1,0].set_title('Filter')
plt.colorbar(im3,ax=ax2[1,0])
# -- filtered img.
im4= ax2[1,1].imshow(imgF,cmap='bone')
ax2[1,1].set_title('Filt. image')
plt.colorbar(im4,ax=ax2[1,1])
ax2[1,1].axis('off')

