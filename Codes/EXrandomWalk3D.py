#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EXrandomWalk3D.py
Created on Fri Sep 24 16:05:56 2021
@author: C. Bergevin

o computes an ensemble 3D random walk (as well as a corr. 1-D walks)
o the 3-D walk is computed in two different ways
 a. a "1-D" walk in 3-D (where for each step, the walker chooses
    two angles and then moves a dist. SS from the last step
 b. "modified" version where the walker take three steps: one along each
    x, y, and z (each of dist. SS) and thus a total dist. of 
    sqrt(x^2 + y^2 + z^2) = SS*sqrt(3)
--> all this culminates in the MSD plot for what the "slope" should 
    be (i.e., 2*D for a, but 6*D for b)
o Note that the diffusion const. is propto. SS^2 (actually, D=SS^2/tau)

"""

import numpy as np
import random
import matplotlib.pyplot as plt
#import sys
#from mpl_toolkits.mplot3d import Axes3D

# ===============================================
M= 300         # numb of walkers to simulate
N= 100         # numb of steps to simulate per walker
SS= 3   # step-size [m]
tau= 0.1   # time per step [s]
IC= [0,0,0]   # initial position of a walker
# ===============================================

Snum= np.linspace(0,N-1,num=N)
t= Snum*tau   # array of time vals.
D= SS**2/(2*tau)   # determine diffusion const.
# --- initialize arrays to hold the #s
x= np.empty((M,N))
y= np.empty((M,N))
z= np.empty((M,N))
r= np.empty((M,N))   # radial position
rM= np.empty((M,N))   # "modified" radial position
x1D= np.empty((M,N))
y1D= np.empty((M,N))
z1D= np.empty((M,N))
phi= np.empty((M,N))
theta= np.empty((M,N))
# ==== loop through for each walker
for m in range(0,M):
    # --- for the m'th walker, now go through steps 2:N
    for n in range(0,N):
        # --- determine random angles
        tmp1= random.uniform(0,1)
        tmp2= random.uniform(0,1)
        tmp3= random.uniform(0,1)
        phi[m,n]= 2*np.pi*tmp1
        theta[m,n]= 2*np.pi*tmp2
        # --- update re the last step
        if n==0:
            x[m,n]= IC[0]
            y[m,n]= IC[1]
            z[m,n]= IC[2]
            x1D[m,n]= IC[0]
            y1D[m,n]= IC[1]
            z1D[m,n]= IC[2]
        else:  
            x[m,n]= x[m,n-1]+ SS*np.cos(theta[m,n-1])* np.sin(phi[m,n-1])
            y[m,n]= y[m,n-1]+ SS*np.sin(theta[m,n-1])* np.sin(phi[m,n-1])
            z[m,n]= z[m,n-1]+ SS* np.cos(phi[m,n-1])      
            # --- deal w/ 1-D case(s)
            if tmp1<0.5:
                x1D[m,n]= x1D[m,n-1]- SS
            else:
                x1D[m,n]= x1D[m,n-1]+ SS
            if tmp2<0.5:
                y1D[m,n]= y1D[m,n-1]- SS
            else:
                y1D[m,n]= y1D[m,n-1]+ SS
            if tmp3<0.5:
                z1D[m,n]= z1D[m,n-1]- SS
            else:
                z1D[m,n]= z1D[m,n-1]+ SS
        # --- determine/store away radial position (kludge; assumes IC is origin)
        r[m,n]= np.sqrt(x[m,n]**2 + y[m,n]**2 + z[m,n]**2)
        # --- "modified" vers
        rM[m,n]= np.sqrt(x1D[m,n]**2 + y1D[m,n]**2 + z1D[m,n]**2)
    
# ==== calculate mean squared dist. as function of time
MSD= np.average(np.square(r),axis=0)
MSD1D= np.average(np.square(x1D),axis=0)
MSDm= np.average(np.square(rM),axis=0)

# ==== plot one of the 3-D walkers
if 1==1:
    fig1= plt.figure()
    ax1= fig1.add_subplot(111, projection='3d')
    ax1.plot(x[m,:],y[m,:],z[m,:],'b.-',lw=1)
    plt.show()

# ==== plot the MSD
fig0= plt.figure()
ax0= fig0.add_subplot()
ax0.plot(t,MSD,label="3-D")
ax0.plot(t,MSD1D,'r--',label="1-D")
ax0.plot(t,MSDm,'k-',label="modified (true) 3-D")
ax0= plt.xlabel('Time [s]')
ax0= plt.ylabel('MSD [m]')
fig0.legend()
fig0= plt.grid()
print(f'1-D slope = {str(2*D)}')
print(f'3-D slope = {str(6*D)}')


#sys.exit() 
