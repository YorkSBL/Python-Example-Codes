#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EXlowpass.py

[IN PROGRESS]
Plot the impulse response and transfer function for a low-pass
filter (SOL to PHYS 4030 F21 HW8; see deVries (1994) ch.6.4
for theory

TO DO
o Not sure how to scale the TF properly to get a "gain" function
of unity for low (i.e., "pass") freqs. Laplace transforms (Kludge
fixed at the moment)
        
REFS
+ deVries (1994) ch.6.4
+ https://en.wikipedia.org/wiki/Low-pass_filter

        
        
Created on Wed Nov 10 12:38:08 2021
@author: C. Bergevin
"""

import numpy as np
import matplotlib.pyplot as plt
#from numpy.fft import fft
from numpy.fft import rfft
from numpy.fft import irfft

# ========= [User Params]  ==================
R = 1000   # resistance [omhs]
C= 1.6*10**(-7)     # capacitance [F]


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
# ----
wC= 1/(R*C)  # cut-off freq. (rads)
fC= 1/(2*np.pi*R*C)  # cut-off freq. (Hz)
#IR= (1/wC)*np.exp(-t*wC)   # impulse response
IR= np.exp(-t*wC)   # impulse response
TF = rfft(IR)   # transfer function 
normA= 1/np.abs(TF[0])  # normaliz. factor for TF (KLUDGE!)
TF= normA*TF

# --- also express TF as a Laplace transform (ref not provided!)
s= -1j*2*np.pi*freq  # complex freq. (purely imag. here)
TFl= wC/(s+wC)

# ==== convolution analysis
# --- create a few sinusoidal waveforms
wf0= np.sin(freq[25]*2*np.pi*t)  # lowest (non-zero) freq
wfCF= np.sin(fC*2*np.pi*t)  # at cutoff freq
wf1= np.sin(2*fC*2*np.pi*t)  # at twice the cutoff freq
# --- convolve 2*CF waveform
convR1= np.convolve(wf1,IR)
convR1= convR1[0:Npts]   # toss redundant parts
convR1spec = rfft(convR1)
# --- similar for others
convR0= np.convolve(wf0,IR)
convR0spec= rfft(convR0[0:Npts])
convRcf= np.convolve(wfCF,IR)
convRcfSpec= rfft(convRcf[0:Npts])

# ==== create broadband noisy signal and "filter" via 
# convolution in spectral domain (i.e., digitally low-pass
# filter a broadband signal)
# --- create flat-spectrum random phase noise in spectral domain
nA= np.ones(len(freq))
nP= (2*np.pi)*np.random.uniform(0,1,len(freq))
noiseS= nA*np.exp(1j*nP)    # spectral rep. of noise
# --- filter (in spect. domain)
noiseSf= noiseS*TF   
# --- convert to time domain
noise= irfft(noiseS) 
noiseF= irfft(noiseSf)   


# ===== print vals. to screen
print('------')
print(f'Cutoff freq = {fC} Hz')

# =======================================
# ===== visualize
plt.close("all")

# --- IR timecourse
if 1==1:
    fig0 = plt.subplots()
    fig0= plt.plot(t,IR)
    fig0= plt.title('Impulse response (IR)') 
    fig0= plt.grid()
    fig0= plt.xlabel('Time [s]') 
    fig0= plt.ylabel('Impulse response: V_out [V]') 
    
# --- transfer fuction (spectral domain)    
if 1==1:
    fig1, ax1 = plt.subplots(2,1)
    ax1[0].plot(freq/1000,20*np.log10(abs(TF)),'k-',label='TF')
    ax1[0].plot(freq/1000,20*np.log10(abs(TFl)),'r--',label='Laplace')
    ax1[0].set_xlabel('Frequency [kHz]')  
    ax1[0].set_ylabel('Magnitude [dB]') 
    ax1[0].legend()
    ax1[0].axvline(fC/1000, color='g', linestyle='dashed', lw=1)
    ax1[0].grid()
    ax1[0].set_xscale('log')
    ax1[0].set_title('Transfer function')
    ax1[1].plot(freq/1000,(np.angle(TF))/(2*np.pi),'k-')
    ax1[1].plot(freq/1000,(np.angle(TFl))/(2*np.pi),'r--')
    ax1[1].axvline(fC/1000, color='g', linestyle='dashed', lw=1)
    ax1[1].set_xlabel('Frequency [kHz]')  
    ax1[1].set_ylabel('Phase [cycs]') 
    ax1[1].grid()
    ax1[1].set_xscale('log')
    fig1.tight_layout(pad=1.5)
 
# --- convolutions
if 1==1:
    fig2, ax2 = plt.subplots(2,1)
    ax2[0].plot(convR1,'k-',lw=2,label='2CF')
    ax2[0].plot(convR0[0:Npts],'b-',lw=1,label='CF')
    ax2[0].plot(convRcf[0:Npts],'r-',lw=1,label='Low')
    ax2[0].set_title('Convolution between sinusoid and IR') 
    ax2[0].grid()
    ax2[0].set_xlabel('Time shift [s]') 
    ax2[0].set_ylabel('Conv. [arb]')
    ax2[1].plot(freq/1000,20*np.log10(abs(convR1spec)),'k-',lw=2,label='2CF')
    ax2[1].plot(freq/1000,20*np.log10(abs(convR0spec)),'b-',lw=1,label='CF')
    ax2[1].plot(freq/1000,20*np.log10(abs(convRcfSpec)),'r-',lw=1,label='Low')
    ax2[1].set_xscale('log')
    ax2[1].grid()
    ax2[1].legend()
    ax2[1].set_xlabel('Frequency [kHz]')  
    ax2[1].set_ylabel('Magnitude [dB]') 
    fig2.tight_layout(pad=1.5)
    
# ---- filtered signal
if 1==1:
    fig3, ax3 = plt.subplots(2,1)
    ax3[0].plot(t,noise,'k-',label='Unfilt.')
    ax3[0].plot(t,noiseF,'r-',label='Filt.')
    ax3[0].grid()
    ax3[0].set_xlabel('Time [s]')  
    ax3[0].set_ylabel('Signal [arb]') 
    ax3[0].set_title('Low-pass filtering of (noisy) broadband signal')
    ax3[1].plot(freq/1000,20*np.log10(abs(noiseS)),'k-',label='Unfilt.')
    ax3[1].plot(freq/1000,20*np.log10(abs(noiseSf)),'r-',label='Filt.')
    ax3[1].set_xlabel('Frequency [kHz]')  
    ax3[1].set_ylabel('Magnitude [dB]')
    ax3[1].grid()
    ax3[1].legend()
    ax3[1].set_xscale('log')
    fig3.tight_layout(pad=1.5)




