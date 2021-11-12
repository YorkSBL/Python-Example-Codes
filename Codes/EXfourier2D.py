#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EXfourier2D.py

[IN PROGRESS]
Purpose: 

o Useful refs:
+ EXimageAnalysis1.py
+ https://thepythoncodingbook.com/2021/08/30/2d-fourier-transform-in-python-and-fourier-synthesis-of-images/
+ https://www.djmannion.net/psych_programming/vision/draw_gratings/draw_gratings.html



Created on Mon Oct 18 12:44:57 2021
@author: CB
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy import fftpack
from matplotlib.colors import LogNorm

# ----------- [User Params]

noiseF= 0.1   # noise factor {0.1?}
bs= 0.0001   # Gaussian filter (SD) factor
M= 1000  #  num of pixels for (square) image length
ell= 10  # wavelength (i.e., # of pixels/period)
angle= -1*np.pi / 2.1  # angle (rads) re vertical [-pi,pi]
logM= 1   # boolean (1= plot mag on log scale)
# ----------- 
  
# ==== create grating
x = np.arange(-np.round(M/2), np.round(M/2), 1)
X, Y = np.meshgrid(x, x)
grating = np.sin(2*np.pi*(X*np.cos(angle) + Y*np.sin(angle)) /ell)
imgN= grating

# ==== deal w/ FFT
imgNspec = fftpack.fft2(imgN)
# --- now shift the zero-frequency component to spectrum's center 
imgNspec = np.fft.fftshift(imgNspec)
imgNspecM= np.abs(imgNspec)
imgNspecP= np.angle(imgNspec)



# ========= VISUALIZE
plt.close("all")
# ---- display image
if 1==0:
    fig1= plt.figure()
    plt.imshow(imgN,cmap='bone')
    plt.grid(False)
    cbar = plt.colorbar()

# ---- display spectal bits
fig2,ax2 = plt.subplots(3)
fig2.set_size_inches(3,6)
# -- noisy img.
im1= ax2[0].imshow(imgN,cmap='bone')
ax2[0].set_title('Grating')
plt.colorbar(im1,ax=ax2[0])
ax2[0].axis('off')
# -- |FFT| 
if logM==1:
    im2=ax2[1].imshow(imgNspecM,norm=LogNorm(vmin=5))
else:
    im2=ax2[1].imshow(imgNspecM)

ax2[1].grid(False)
ax2[1].axes.xaxis.set_visible(False)
ax2[1].axes.yaxis.set_visible(False)
ax2[1].set_title('|FFT| (cent. shift.)')
plt.colorbar(im2,ax=ax2[1])
# -- phase
im3=ax2[2].imshow(imgNspecP)
ax2[2].axis('off')
ax2[2].set_title('Phase')
plt.colorbar(im3,ax=ax2[2])

"""
# ---- display spectal bits
fig2,ax2 = plt.subplots(2,2)
# -- noisy img.
im1= ax2[0,0].imshow(imgN,cmap='bone')
ax2[0,0].set_title('Grating')
plt.colorbar(im1,ax=ax2[0,0])
ax2[0,0].axis('off')
# -- |FFT| 
im2=ax2[1,0].imshow(imgNspecM,norm=LogNorm(vmin=5))
ax2[1,0].grid(False)
ax2[1,0].axes.xaxis.set_visible(False)
ax2[1,0].axes.yaxis.set_visible(False)
ax2[1,0].set_title('|FFT| of noisy original (cent. shift.)')
plt.colorbar(im2,ax=ax2[1,0])
# -- phase
im3=ax2[1,1].imshow(imgNspecP)
ax2[1,1].axis('off')
ax2[1,1].set_title('Phase')
plt.colorbar(im3,ax=ax2[1,1])
"""


