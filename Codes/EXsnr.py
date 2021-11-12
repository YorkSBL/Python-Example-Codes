#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EXsnr.py

Purpose: Bury a sinusoid in noise, determine the assoc. SNR,
and plot both the time waveform and spectral magnitude

Note: The noise doesn't quite have unit amplitude (there should
be some sort of normalization factor in there)

Created on Thu Nov 11 12:39:51 2021
@author: C. Bergevin
"""

import numpy as np
import matplotlib.pyplot as plt
from numpy.fft import rfft
from numpy.fft import irfft

# ========= [User Params]  ==================

fs= 1343   # sinusoidal freq [Hz]
As= 1.1      # amplitude of sinusoid
SR= 44100;         # sample rate [Hz]
Npts= 8192;     # length of fft window (# of points) 
#                     [should ideally be 2^N]
                                           
# ====================================
# ==== bookeeping
# --- create a freq. array (for FFT bin labeling)
freq= np.arange(0,(Npts+1)/2,1)    
freq= SR*freq/Npts
df = SR/Npts   # freq. bin width
t= np.linspace(0,(Npts-1)/SR,Npts)   # time array


# --- create noise (two types)
if 1==0:
    # flat-spectrum random phase noise in spectral domain
    nA= np.ones(len(freq))
    nP= (2*np.pi)*np.random.uniform(0,1,len(freq))
    noiseS= nA*np.exp(1j*nP)    # spectral rep. of noise
    noiseT= irfft(noiseS)   # convert to time domain
else:
    # noisy noise ;-)
    noiseT= np.random.normal(0,1,Npts)
        
signal= noiseT+ As*np.sin(fs*2*np.pi*t)
signalS= rfft(signal)  # compute spectrum

# ===== print vals. to screen
print('------')
print(f'SNR = {20*np.log10(As/1)} dB')
# Note: As coded, the noise has ~unit amplitude (hence /1)


# ===== visualize
plt.close("all")
# -----
fig3, ax3 = plt.subplots(2,1)
ax3[0].plot(t,signal,'k-',label='Sig.')
ax3[0].grid()
ax3[0].set_xlabel('Time [s]')  
ax3[0].set_ylabel('Signal [arb]') 
ax3[1].plot(freq/1000,20*np.log10(abs(signalS)),'k-',label='Unfilt.')
ax3[1].set_xlabel('Frequency [kHz]')  
ax3[1].set_ylabel('Magnitude [dB]')
ax3[1].grid()
ax3[1].legend()
ax3[1].set_xscale('log')
fig3.tight_layout(pad=1.5)