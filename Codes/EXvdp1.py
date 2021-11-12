#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EXvdp1.py

Purpose: Numerically integrate the van der Pol osc using a hard-coded 
RK4
dx/dt = y
dy/dt = mu*(1-x^2)*y + (w^2)*x

o x and y vals passed along via array r
o this code marks my first attempt to use the "def" functionality
o Borrowing some syntax EXlotkaVolterra.py
o see also Newman (2012) exercise 8.6e (pg.352)

Created on Mon Jun 14 09:53:05 2021
@author:CB
"""

import math
import numpy as np
import random
import matplotlib.pyplot as plt



# ----------- [User Params]
# --- params
mu = 3.75      # 
w = 1.9     # natural freq
# --- ICs
x0 = 1.0       # initial x
y0 = 0.0       # initial y
# --- time params
tS = 0         # start time
tE = 20        # end time
N = 1000;      # num of time points to use
# -----------

# ----------- [define ODE system]
def f(r,t):
    x = r[0]
    y = r[1]
    fx = y
    fy = mu*(1-x**2)*y - (w**2)*x
    return np.array([fx,fy],float)

# --- bookkeeping
h = (tE-tS)/N   # step-size
tpoints = np.arange(tS,tE,h)   # time point array
x = []
y = [] 
#r = np.array([x0,y0],float)
r = [x0,y0]
indx = 0
# --- integration loop
for t in tpoints:
    x = np.insert(x,indx,r[0])
    y = np.insert(y,indx,r[1])
    k1 = h*f(r,t)
    k2 = h*f(r+0.5*k1,t+0.5*h)
    k3 = h*f(r+0.5*k2,t+0.5*h)
    k4 = h*f(r+k3,t+h)
    r += (k1+2*k2+2*k3+k4)/6
    indx = indx+1
  
    
# --- visualize
# time waveforms
fig1, ax1 = plt.subplots(2,1)
ax1[0].plot(tpoints,x,'b-',label='X')
#plt.xlim([-1,1]) 
ax1[0].set_xlabel('Time')  
ax1[0].set_ylabel('Position x') 
ax1[0].set_title('van der Pol System')
ax1[0].grid()
#ax1[0].legend(loc=1)
# phase space
ax1[1].plot(x,y,'r-',label='X')
ax1[1].set_xlabel('x')  
ax1[1].set_ylabel('dx/dt') 
ax1[1].grid()
fig1.tight_layout(pad=1.5)