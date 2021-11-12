#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EXvdp3.py

[IN PROGRESS; sinusoidal bit doesn't seem to work right...]
Purpose: Numerically integrate the van der Pol osc using a hard-coded 
RK4
dx/dt = y
dy/dt = mu*(1-x^2)*y + (w^2)*x + A*sin(wd*t)

o same as EXvdp2.py except adding in means to compute FFT of 
steady-state portion
o also weaves in bits/updates from EXvdpNoisy1.py (e.g., randomize
 the ICs)                                               

Created on 2021.07.13
@author:CB
"""

import math
import numpy as np
import random
import matplotlib.pyplot as plt
from numpy.fft import rfft




# ----------- [User Params]
# --- params
mu = 1.2      # (negative) damping term
w = 10     # natural freq (rads/s)
# --- sinusoidal drive
A = 10      #drive amplitude
wd = 0.55*w     # drive freq
# --- ICs
rIC = 1       # boolen to randomize ICs (set to 1 to turn on)
x0 = 1.0       # initial x (assuming rIC=0)
y0 = 0.0       # initial y  (assuming rIC=0)
# --- time & spec params
tE = 20        # integration time (assumes tS=0)
SR= 5*10**3     # "sample rate" (i.e., step-size = 1/SR) [Hz] {5*10^5}
Npoints= 5*8192  # of point for "steady-state" waveform to compute spectrum of
# -----------


# ----------- [define ODE system]
def f(r,t):
    x = r[0]
    y = r[1]
    fx = y
    fy = mu*(1-x**2)*y - (w**2)*x + A*np.sin(wd*t)
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
x = []
y = [] 
r = [x0,y0]

# --- print a few vals to screen
print(f'Natural freq: f = {str(w/2*np.pi)} Hz')
print(f'Drive freq: fd = {str(wd/2*np.pi)} Hz')

# --- integration loop
indx = 0
for t in tpoints:
    x = np.insert(x,indx,r[0])
    y = np.insert(y,indx,r[1])
    r = rk4(r,t,h)
    indx = indx+1

# ===== deal w/ FFT & whatnot
df = SR/Npoints
freq= np.arange(0,(Npoints+1)/2,1)    # create a freq. array (for FFT bin labeling)
freq= SR*freq/Npoints;
signal = x[indx-Npoints:indx]  # grab last nPts
spec= rfft(signal)
specM = abs(spec)
specMdb = 20*np.log10(specM)

# ===== visualize
# ~~ Fig.1A: time waveforms
fig1, ax1 = plt.subplots(2,1)
ax1[0].plot(tpoints,x,'b-',label='X')
ax1[0].plot(tpoints[indx-Npoints:indx],signal,'r.',label='SS bit')
ax1[0].set_xlabel('Time')  
ax1[0].set_ylabel('Position x') 
ax1[0].set_title('van der Pol System')
ax1[0].grid()
#ax1[0].legend(loc=1)
# ~~ Fig.1B: phase space
ax1[1].plot(x,y,'r-',label='X')
ax1[1].set_xlabel('x')  
ax1[1].set_ylabel('dx/dt') 
ax1[1].grid()
fig1.tight_layout(pad=1.5)
# ~~ Fig.2: Spectrum
fig2 = plt.subplots()
fig2= plt.plot(freq/1000,specMdb,'r.-',label='X')
fig2= plt.xlabel('Frequency [kHz]')  
fig2= plt.ylabel('Magnitude [dB]') 
fig2= plt.grid()
fig2= plt.xlim([0, 5.5*(w/(2*np.pi))/1000])
