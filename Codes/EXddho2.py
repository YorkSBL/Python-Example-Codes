#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EXddho2.py

[IN PROGRESS; not yet presently using the blackbox]
Purpose: Deal w/ the DDHO via blackbox scipy solver
dx/dt = y
dy/dt = -gamma*y - (w^2)*x + A*sin(wd*t)
where gamma=b/m, w= sqrt(k/m), and Ad= A/m

o modified from EXddho.py and 
o includes means to compute FFT of  steady-state portion
o Be careful about the relative val. of the params. (e.g., if w is 
large, than A needs to be large too) as otherwise the relative
force terms are small and transients become problematic in the
spectrum   
o Ref site:
https://pythonnumericalmethods.berkeley.edu/notebooks/chapter22.06-Python-ODE-Solvers.html                                       
 https://sam-dolan.staff.shef.ac.uk/mas212/notebooks/ODE_Example.html                                        

Created on 2021.09.29
@author:CB
"""

import math
import numpy as np
import random
import matplotlib.pyplot as plt
from numpy.fft import rfft
from scipy.signal import find_peaks
from scipy.integrate import solve_ivp


# ----------- [User Params]
# --- params

# --- original code
"""
gamma = 10.1      # damping term
w = 1000     # natural freq (rads/s)
# --- sinusoidal drive
A = 1      #drive amplitude
wd = 2.25*w     # drive freq/natural freq. (e.g., wd=1 drives at resonance)
"""
# --- new code
m= 1  # oscillator mass [g]
b= 0.01  # damping coeffic.
k= 5000  # stiffness [N/m]
# --- sinusoidal drive
Ad = 10      #drive amplitude
wd = 1.0005     # drive freq/natural freq. (e.g., wd=1 drives at resonance)


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
x = []
y = [] 
r = [x0,y0]

# ----

m= m/1000  # convert mass to kg
gamma = b/m      # damping term
w = np.sqrt(k/m)     # natural freq (rads/s)
A= Ad/m
wd = wd*w   # get wd to be an actual (not rel.) freq.


# --- print a few vals to screen
print(f'Natural freq: f = {str(w/(2*np.pi))} Hz')
print(f'Drive freq: fd = {str(wd/(2*np.pi))} Hz')





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
# --- extract peak mag.
maxA= np.max(specM)

# ===== visualize
plt.close("all")
# ~~ Fig.1A: time waveforms
fig1, ax1 = plt.subplots(1)
ax1.plot(tpoints,x,'b-',label='X')
ax1.plot(tpoints[indx-Npoints:indx],signal,'r.',label='SS bit')
"""
fig1, ax1 = plt.subplots(2,1)
ax1[0].plot(tpoints,x,'b-',label='X')
ax1[0].plot(tpoints[indx-Npoints:indx],signal,'r.',label='SS bit')
ax1[0].set_xlabel('Time')  
ax1[0].set_ylabel('Position x') 
ax1[0].set_title('Title')
ax1[0].grid()
#ax1[0].legend(loc=1)
# ~~ Fig.1B: phase space
ax1[1].plot(x,y,'r-',label='X')
ax1[1].set_xlabel('x')  
ax1[1].set_ylabel('dx/dt') 
ax1[1].grid()
fig1.tight_layout(pad=1.5)
"""
# ~~ Fig.2: Spectrum
fig2 = plt.subplots()
fig2= plt.plot(freq/(w/(2*np.pi)),specMdb,'r.-',label='X')
fig2= plt.xlabel('Normaliz. Freq. [f/fo]')  
fig2= plt.ylabel('Magnitude [dB]') 
fig2= plt.grid()
fig2= plt.xlim([0,3])
