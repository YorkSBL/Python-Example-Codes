#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EXanscombe.py
Purpose: Code to reproduce Anscombe's quartet fits 
(as per https://en.wikipedia.org/wiki/Anscombe's_quartet#/media/File:Anscombe%27s_quartet_3.svg
o  [see also EXregressAnscombe.m]
o Linear regression done via Bevington (ch.6 eqns.6.13 and 6.23)
o calculation for r^2 (coefficient of determination) via:
    http://en.wikipedia.org/wiki/Coefficient_of_determination
Created on Thu Jun 10 10:03:20 2021
@author: CB
"""
#import math
import numpy as np
#import random
import matplotlib.pyplot as plt
# ==========================
# --- specify the associated data (here as a list)
A =[[10,    8,   13,    9,   11,   14,    6,    4,   12,    7,    5],
    [8.04, 6.95, 7.58, 8.81, 8.33, 9.96, 7.24, 4.26, 10.84, 4.82, 5.68],
   [10, 8, 13, 9, 11, 14, 6, 4,   12,    7,    5],
    [9.14, 8.14, 8.74, 8.77, 9.26, 8.1, 6.13, 3.1, 9.13, 7.26, 4.74],
   [10,    8,   13,    9,   11,   14,    6,    4,   12,    7,    5],
    [7.46, 6.77, 12.74, 7.11, 7.81, 8.84, 6.08, 5.39, 8.15, 6.42, 5.73],
    [8,    8,    8,    8.,    8,    8,    8,   19,    8,    8,    8],
    [6.58, 5.76, 7.71, 8.84, 8.47, 7.04, 5.25, 12.5, 5.56, 7.91, 6.89]]
# ==========================
# --- initialize arrays (via numpy) for key vals of each Anscombe pair)
m = np.array([])
b = np.array([])
sigma_m = np.array([])
sigma_b = np.array([])
Rsquared= np.array([])
# --- directly compute best fits for each set
for n in [0,1,2,3]:
    xT = A[2*n]   # x-vals
    yT = A[2*n+1] # y-vals
    N = len(xT)
    # --- compute slope and intercept via linear regression
    Delta = N*np.sum(np.square(xT))- np.square(np.sum(xT))
    bT = (1/Delta)*(np.sum(np.square(xT))*np.sum(yT)- \
                     np.sum(xT)*np.sum(np.multiply(xT,yT)))
    mT = (1/Delta)*(N*np.sum(np.multiply(xT,yT)) - \
                    np.sum(xT)*np.sum(yT))
    # --- compute associated uncertainties
    sigma_mT= np.sqrt(1/Delta*sum(np.square(xT)))
    sigma_bT= np.sqrt(N/Delta )  
    # --- also determine r^2
    meanY= np.sum(yT)/N;
    SStot= np.sum(np.square(yT-meanY))
    SSres= np.sum(np.square(yT-bT-np.multiply(mT,xT)))
    RsquaredT= 1- SSres/SStot
    # --- store away vals (for each Anscombe pair)
    m = np.insert(m,n,mT)
    b = np.insert(b,n,bT)
    Rsquared = np.insert(Rsquared,n,RsquaredT)
    
# --- print relevant vals to screen
for n in [0,1,2,3]:
    print(f'Set {str(n+1)}: y={str(m[n])}*x + {str(b[n])} (R^2 = {str(Rsquared[n])})')
# --- create model curves
x1= A[0]
ym1= b[0]+ np.multiply(m[0],x1)
x2= A[2]
ym2= b[1]+ np.multiply(m[1],x2)
x3= A[4]
ym3= b[2]+ np.multiply(m[2],x3)
x4= A[6]
ym4= b[3]+ np.multiply(m[3],x4)

# --- visualize
plt.close("all")
fig, axs = plt.subplots(2, 2)
axs[0,0].plot(A[0],A[1],'o')
axs[0,0].plot(x1,ym1,'r-')
axs[0,0].set_title('Set 1')
axs[0,0].grid()
axs[0,1].plot(A[2],A[3],'o')
axs[0,1].plot(x2,ym2,'r-')
axs[0,1].set_title('Set 2')
axs[0,1].grid()
axs[1,0].plot(A[4],A[5],'o')
axs[1,0].plot(x3,ym3,'r-')
axs[1,0].set_title('Set 3')
axs[1,0].grid()
#axs[1,0].xlabel('x') 
axs[1,1].plot(A[6],A[7],'o',label='data')
axs[1,1].plot(x4,ym4,'r-',label='model')
axs[1,1].set_title('Set 4')
axs[1,1].grid()
fig.tight_layout(pad=1.5)
axs[1,0].set(xlabel='x', ylabel='y')
legend = axs[1,1].legend(loc=4)
plt.show()

