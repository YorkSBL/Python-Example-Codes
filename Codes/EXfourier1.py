#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EXfourier1.py

Purpose: Learn Pythonn's syntax for 1-D FFTs using a noisy sinusoid 
as the waveform
    
o Lifting bits from my EXfftVSrfft.m
o Note the distiction creating the time arrays as follows:
    + Matlab: > t=[0:1/SR:(Npoints-1)/SR];
    + Python: >>> t = np.arange(0,Npoints/SR,1/SR)
o Be careful re:
    bad: noise= np.random.randn(Npoints,1)
    good: noise= np.random.randn(Npoints)

Created on Fri Jun 18 17:32:50 2021
@author: CB
"""

import numpy as np
import matplotlib.pyplot as plt
from numpy.fft import rfft

# -------------------------------- 
SR= 44100;         # sample rate [Hz]
Npoints= 8192;     # length of fft window (# of points) [should ideally be 2^N]
                   # [time window will be the same length]
f= 2510.0;         # Frequency (for waveforms w/ tones) [Hz] {2514}
noiseA= 0.1;     # amplitude factor for relative (flat-spectrum) noise
# -------------------------------- 

# ==== bookeeping
df = SR/Npoints;  
fQ= np.ceil(f/df)*df;   # quantized natural freq.
t = np.arange(0,Npoints/SR,1/SR)  # create an array of time points, Npoints long
dt= 1/SR;  # spacing of time steps
freq= np.arange(0,(Npoints+1)/2,1)    # create a freq. array (for FFT bin labeling)
freq= SR*freq/Npoints;
# ===== create base signal + noise
if (1==1):   
    base= np.cos(2*np.pi*f*t)   # non-quantized sinusoid
else :     
    base= np.cos(2*np.pi*fQ*t)  # quantized sinusoid
noise= np.random.randn(Npoints);      # Gaussian
signal= base+ noiseA*noise;
# ===== deal w/ FFT & whatnot
spec= rfft(signal)
specM = abs(spec)
specMdb = 20*np.log10(specM)
# --- visualize
fig1, ax1 = plt.subplots(2,1)
#fig1=plt.plot([])  # set up fig
ax1[0].plot(t,signal,'k.-',label='wf')
ax1[0].set_xlabel('Time [s]')  
ax1[0].set_ylabel('Signal [arb]') 
ax1[0].set_title('dealing w/ rfft in Python')
ax1[0].grid()
ax1[1].plot(freq/1000,specMdb,'r.-',label='X')
ax1[1].set_xlabel('Frequency [kHz]')  
ax1[1].set_ylabel('Magnitude [dB]') 
ax1[1].grid()
ax1[1].set_xlim([0, 2.5*f/1000])
fig1.tight_layout(pad=1.5)



