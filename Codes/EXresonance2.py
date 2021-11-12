#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EXresonance2.py

[IN PROGRESS]
Purpose: Numerically integrate the DDHO using a hard-coded RK4
dx/dt = y
dy/dt = -gamma*y - (w^2)*x + A*sin(wd*t)

o Cycle through different DDHO drive freqs so to reconstruct the 
"resonance curve"
o built off of EXddho.py (formerly EXresonance.py)
o Note: Be careful about the relative val. of the params. (e.g., if w is 
large, than A needs to be large too) as otherwise the relative
force terms are small and transients become problematic in the
spectrum  
o Need to "find" a spectral peak for each wd. Using np.max, but can 
also use np.findpeaks. Ref. page is:
https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.find_peaks.html
o                                         
                                            

Created on 2021.09.29
@author:CB
"""

import math
import numpy as np
import random
import matplotlib.pyplot as plt
from numpy.fft import rfft
from scipy.signal import find_peaks


# ----------- [User Params]
# --- params
gamma = 100.0      # damping term
w = 1000     # natural freq (rads/s)
# --- sinusoidal drive
A = 100      # drive amplitude
wdNUM= 45     # # of drive freqs to run
wdR = w*np.linspace(0.65,1.35,wdNUM+1)    # drive freqs
# --- ICs
rIC = 0       # boolen to randomize ICs (set to 1 to turn on)
x0 = 0.0       # initial x (assuming rIC=0)
y0 = 0.0       # initial y  (assuming rIC=0)
# --- time & spec params
tE = 5        # integration time (assumes tS=0)
SR= 5*10**3     # "sample rate" (i.e., step-size = 1/SR) [Hz] {5*10^5}
Npoints= 2*8192  # of point for "steady-state" waveform to compute spectrum of
# -----------


# ----------- [define ODE system]
def f(r,t):
    x = r[0]
    y = r[1]
    fx = y
    fy = -gamma*y - (w**2)*x + A*np.sin(wd*t)
    #print(f't = {str(A*np.sin(wd*t))}')
    #fy = -(w**2)*x + A*np.sin(wd*t)  # Driven HO
    #fy = -(w**2)*x - mu*y + A*np.sin(wd*t)  # DDHO
    return np.array([fx,fy],float)


# ----------- [define RK4]
def rk4(r,t,h):
    k1 = h*f(r,t)
    k2 = h*f(r+0.5*k1,t+0.5*h)
    k3 = h*f(r+0.5*k2,t+0.5*h)
    k4 = h*f(r+k3,t+h)
    r += (k1+2*k2+2*k3+k4)/6
    return r

# --- randomize ICs?
if rIC==1:
    x0 = random.uniform(-1,1)
    y0 = random.uniform(-1,1)

# ===== bookkeeping
tS = 0         # start time
N = tE*SR     # num of time points
h = (tE-tS)/N   # step-size
tpoints = np.arange(tS,tE,h)   # time point array
r = [x0,y0]
maxA= []

# --- print a few vals to screen
#print(f'Natural freq: f = {str(w/2*np.pi)} Hz')
#print(f'Drive freq: fd = {str(wd/2*np.pi)} Hz')

# ==== loop to go through each drive freq.
for n in range(0,len(wdR)):
    wd= wdR[n]
    x = []
    y = [] 
    # ---- integration loop (for a given drive freq)
    indx = 0
    for t in tpoints:
        x = np.insert(x,indx,r[0])
        y = np.insert(y,indx,r[1])
        r = rk4(r,t,h)
        indx = indx+1
        
    # ---- deal w/ FFT & whatnot
    df = SR/Npoints
    freq= np.arange(0,(Npoints+1)/2,1)    # create a freq. array (for FFT bin labeling)
    freq= SR*freq/Npoints;
    signal = x[indx-Npoints:indx]  # grab last nPts
    spec= rfft(signal)
    specM = abs(spec)
    specMdb = 20*np.log10(specM)
    # --- extract peak mag. and store aways
    maxA= np.insert(maxA,n,np.max(specM))



# ===== visualize
plt.close("all")
# --- Fig.1: (last wd run) time waveforms & phase space
if 1==1:
    fig1, ax1 = plt.subplots(2,1)
    ax1[0].plot(tpoints,x,'b-',label='X')
    ax1[0].plot(tpoints[indx-Npoints:indx],signal,'r.',label='SS bit')
    ax1[0].set_xlabel('Time')  
    ax1[0].set_ylabel('Position x') 
    ax1[0].set_title('DDHO responses')
    ax1[0].grid()
    #ax1[0].legend(loc=1)
    # ~~ Fig.1B: phase space
    ax1[1].plot(x,y,'r-',label='X')
    ax1[1].set_xlabel('x')  
    ax1[1].set_ylabel('dx/dt') 
    ax1[1].grid()
    fig1.tight_layout(pad=1.5)

# --- Fig.2: (last wd run) spectrum
if 1==1:
    fig2 = plt.subplots()
    fig2= plt.plot(freq/(w/(2*np.pi)),specMdb,'r.-',label='X')
    fig2= plt.xlabel('Normaliz. Freq. [f/fo]')  
    fig2= plt.ylabel('Magnitude [dB]') 
    fig2= plt.grid()
    fig2= plt.xlim([0,3])
    
# --- Fig.3: (last wd run) time waveforms & phase space
if 1==1:
    fig3, ax3 = plt.subplots()
    fig3= plt.plot(wdR/w,maxA,'bo-')
    fig3= plt.xlabel('Normaliz. Freq. [f/fo]')  
    fig3= plt.ylabel('Steady-state amplitude') 
