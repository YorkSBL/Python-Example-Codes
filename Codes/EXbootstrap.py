#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EXbootstrap.py

[IN PROGRESS]
Purpose: Read in a data file and do some boostrapping analysis 
(this one deals w/ estimating Planck's const)     
--> This is a SOL to HW7 for PHYS 4030 F21 (see the assignment
 pdf for more background; 4030F21HW7.pdf]
                                              
 Notes:
o This requires the external file planckdata.xlsx  
o Code is parsed up into two sections: 1. pre-bootstrap and
subsequently 2. bootstrap. 
o For the pre-bootstrap, at present, the uncertainties 
assoc. w/ the regression are not factored in for estimating
h (it could/should though!)

REFS
+ https://www.scienceinschool.org/article/2014/planck/
+ https://www.cs.toronto.edu/~frossard/post/linear_regression/
+ https://en.wikipedia.org/wiki/Propagation_of_uncertainty

Created on Tue Nov  2 13:56:32 2021
@author: C. Bergevin
"""


import numpy as np
import random
import matplotlib.pyplot as plt
import pandas as pd


# ========= [User Params]  ==================
fileN= "./Other/planckdata.xlsx"
thresh= 0.01   # threshold to "trim" the raw #s
# --- bootstrap-related vals
#n= 20  
Nbs=10000   # numb. of bootstraps to do for each color data set

# ==== physical consts. (used in analysis)
ec=  1.6022* 10**(-19)  # unit electron charge [C]
c= 2.99792458* 10**(8)   # speed of light [m/s]
hEXACT= 6.62607015* 10**(-34)   # exact val. for h [J/Hz]
# --- wavelengths 
lR= 623* 10**(-9)  # red wavelengths [m]
lO= 586* 10**(-9)  # orange
lG= 567* 10**(-9)  # green
lB= 467* 10**(-9)  # blue


# ====================================

# ==== bookkeeping
# ---- read in data from file
dR = np.array(pd.read_excel(fileN,sheet_name="RED"))
dO = np.array(pd.read_excel(fileN,sheet_name="ORANGE"))
dG = np.array(pd.read_excel(fileN,sheet_name="GREEN"))
dB = np.array(pd.read_excel(fileN,sheet_name="BLUE"))
# ---- "trim" the data (via a thresholding)
dRt= dR[dR[:,1]>thresh]
dOt= dO[dO[:,1]>thresh]
dGt= dG[dG[:,1]>thresh]
dBt= dB[dB[:,1]>thresh]
# NOTE: the "x" data (i.e., voltage) for red is dRt[:,0] while
# the "y" data (i.e., current) is dRt[:,1] (similarly for orange et al)

# =======================================
# ======= 1. Pre-bootstrap: linear regress. on (all) trimmed data
# --- initialize arrays (via numpy) for key vals of each Anscombe pair)
m = np.array([])
b = np.array([])
sigma_m = np.array([])
sigma_b = np.array([])
Rsquared= np.array([])
# ---- loop thru for each color (and directly compute regression)
for n in [0,1,2,3]:
    if n==0:
        xT= dRt[:,0]
        yT= dRt[:,1]
    elif n==1:
        xT= dOt[:,0]
        yT= dOt[:,1]
    elif n==2:
        xT= dGt[:,0]
        yT= dGt[:,1]
    elif n==3:
        xT= dBt[:,0]
        yT= dBt[:,1]
    
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
    # --- store away vals (for each data set)
    m = np.insert(m,n,mT)
    b = np.insert(b,n,bT)
    Rsquared = np.insert(Rsquared,n,RsquaredT)
    
# --- create model curves (for plotting)
ymR= b[0]+ np.multiply(m[0],dRt[:,0])
ymO= b[1]+ np.multiply(m[1],dOt[:,0])
ymG= b[2]+ np.multiply(m[2],dGt[:,0])
ymB= b[3]+ np.multiply(m[3],dBt[:,0])
# --- find x-intercept (i.e., Vt)
vTR= -b[0]/m[0]
vTO= -b[1]/m[1]
vTG= -b[2]/m[2]
vTB= -b[3]/m[3]
# ==== put pieces together to estimate h (via linear regress of vT)
# --- repackage
xF1= [1/lR,1/lO,1/lG,1/lB]  # the "x"-data to fit
yF1= [vTR,vTO,vTG,vTB]  # the "y"-data to fit
N = len(xF1)
Delta = N*np.sum(np.square(xF1))- np.square(np.sum(xF1))
m1 = (1/Delta)*(N*np.sum(np.multiply(xF1,yF1)) - \
                    np.sum(xF1)*np.sum(yF1))
h1est= m1*ec/c  # ... and finally, our estimate for h

# ===== print vals. to screen
print('--- 1. Pre-Boostrap ---')
print(f'Exact value for h = {hEXACT} J/Hz')
print(f'Pre-bootstrap est. for h = {h1est} J/Hz')

# =======================================
# ======= 2. Bootstrap:
# --- initialize arrays for bootstrapped vT vals. [R O G B]
vTBS = np.array([])   # bootstrapped avg. vT 
vTBSstd = np.array([])  # assoc. uncertainty
vT_hist= np.zeros((4,Nbs))  # raw bootstrapped vT vals. (to histogram)
# ==== loop through each color
for mm in [0,1,2,3]:
    # --- load in appropriate trimmed data for a given color
    if mm==0:
        vals= dRt
    elif mm==1:
        vals=dOt
    elif mm== 2:
        vals=dGt
    elif mm== 3:
        vals=dBt

    # --- initialize arrays to store vals    
    mBS = np.array([])
    bBS = np.array([])
    vTbsN = np.array([])
    nbs= len(vals) 
    indx= np.arange(nbs)  # create array index 
    # === now loop through each bootstrap
    for nn in np.linspace(0,Nbs-1,Nbs):
        # === grab a resampled array 
        indxBS= np.random.choice(indx,replace=1,size=nbs)
        xT= vals[indxBS,0] # volt. vals
        yT= vals[indxBS,1] # curr. vals
        N = len(xT)
        # === calc. linear regress. 
        Delta = N*np.sum(np.square(xT))- np.square(np.sum(xT))
        bT = (1/Delta)*(np.sum(np.square(xT))*np.sum(yT)- \
                         np.sum(xT)*np.sum(np.multiply(xT,yT)))
        mT = (1/Delta)*(N*np.sum(np.multiply(xT,yT)) - \
                        np.sum(xT)*np.sum(yT))
         # --- store away vals (for each bootstrapped set)
        mBS = np.insert(mBS,int(nn),mT)
        bBS = np.insert(bBS,int(nn),bT)
        vTbsN= np.insert(vTbsN,int(nn),-bT/mT)
        
    # ==== final stats + bookkeeping
    # --- first calc. mean and std across the bootstrapped resamples
    mBSavg= np.mean(mBS)
    mBSstd= np.std(mBS)
    bBSavg= np.mean(bBS)
    bBSstd= np.std(bBS) 
    # ---- now use those to compute bootstrapped vT and uncertainty
    # (Note: need to use propagation of error here! easy enough to 
    # calculate directly or use generic formula for f=A/B)
    vTavg= -bBSavg/mBSavg
    vTstd= vTavg*np.sqrt((bBSstd/bBSavg)**2 + (mBSstd/mBSavg)**2)
    # --- store avgd. vals (comparable to yF1)
    vTBS = np.insert(vTBS,int(mm),vTavg)
    vTBSstd = np.insert(vTBSstd,int(mm),vTstd)
    # --- also store raw vT vals for each color (used to histogram)
    vT_hist[mm,:] = vTbsN
 
 # ==== now use bootstraped vals to determine h 
xF2= [1/lR,1/lO,1/lG,1/lB]  # the "x"-data to fit
yF2= vTBS  # the "y"-data to fit
N = len(xF2)
Delta = N*np.sum(np.square(xF2))- np.square(np.sum(xF2))
m2 = (1/Delta)*(N*np.sum(np.multiply(xF2,yF2)) - \
                    np.sum(xF2)*np.sum(yF2))
sigma_m2= np.sqrt(1/Delta*sum(np.square(xF2)))
h2= m2*ec/c  # ... and finally, our estimate for h
sigma_h2= h2*sigma_m2
    
    
# ===== print vals. to screen
print('--- 2. Boostrap ---')
print(f'Bootstrap est. for h = {h2} +/- {sigma_h2} J/Hz')

# =======================================
# ===== visualize
plt.close("all")

# --- plot data + linear fits to ALL (trimmed) data
if 1==1:
    fig, axs = plt.subplots(2, 2)
    axs[0,0].plot(dR[:,0],dR[:,1],'r.',alpha=0.1)
    axs[0,0].plot(dRt[:,0],dRt[:,1],'ro',alpha=0.3)
    axs[0,0].plot(dRt[:,0],ymR,'r-', linewidth=2)
    axs[0,0].set_title('Red') 
    axs[0,0].grid()
    axs[0,0].set_xlabel('Voltage [V]') 
    axs[0,0].set_ylabel('Current [A]') 
    axs[0,1].plot(dO[:,0],dO[:,1],'.', color='orange',alpha=0.2)
    axs[0,1].plot(dOt[:,0],dOt[:,1],'o',color='orange',alpha=0.3)
    axs[0,1].plot(dOt[:,0],ymO,'-',linewidth=2,color='orange')
    axs[0,1].set_title('Orange') 
    axs[0,1].grid()
    axs[0,1].set_xlabel('Voltage [V]') 
    axs[0,1].set_ylabel('Current [A]') 
    #axs[1,0].plot(dG[:,0],dG[:,1],'go')
    axs[1,0].plot(dG[:,0],dG[:,1],'.', color='green',alpha=0.2)
    axs[1,0].plot(dGt[:,0],dGt[:,1],'o',color='green',alpha=0.3)
    axs[1,0].plot(dGt[:,0],ymG,'-',linewidth=2,color='green')
    axs[1,0].set_title('Green') 
    axs[1,0].grid()
    axs[1,0].set_xlabel('Voltage [V]') 
    axs[1,0].set_ylabel('Current [A]')
    axs[1,1].plot(dB[:,0],dB[:,1],'.', color='blue',alpha=0.2)
    axs[1,1].plot(dBt[:,0],dBt[:,1],'o',color='blue',alpha=0.3)
    axs[1,1].plot(dBt[:,0],ymB,'-',linewidth=2,color='blue')
    axs[1,1].set_title('Blue') 
    axs[1,1].grid()
    axs[1,1].set_xlabel('Voltage [V]') 
    axs[1,1].set_ylabel('Current [A]')      
    fig.tight_layout(pad=1.5)


# --- histogram of bootstrapped vT vals for each color
if 1==1:
    fig, axs = plt.subplots(2, 2)
    axs[0,0].hist(vT_hist[0,:],density=1,bins=30,color = "red",lw=0)
    axs[0,0].axvline(vTBS[0], color='k', linestyle='dashed', lw=2)
    axs[0,0].set_title('vT: Red') 
    axs[0,0].set_xlabel('vT') 
    axs[0,0].set_ylabel('Count')
    axs[0,1].hist(vT_hist[1,:],density=1,bins=30,color = "orange",lw=0)
    axs[0,1].axvline(vTBS[1], color='k', linestyle='dashed', lw=2)
    axs[0,1].set_title('vT: Orange') 
    axs[0,1].set_xlabel('vT') 
    axs[0,1].set_ylabel('Count') 
    axs[1,0].hist(vT_hist[2,:],density=1,bins=30,color = "green",lw=0)
    axs[1,0].axvline(vTBS[2], color='k', linestyle='dashed', lw=2)
    axs[1,0].set_title('vT: Green') 
    axs[1,0].set_xlabel('vT') 
    axs[1,0].set_ylabel('Count') 
    axs[1,1].hist(vT_hist[3,:],density=1,bins=30,color = "blue",lw=0)
    axs[1,1].axvline(vTBS[3], color='k', linestyle='dashed', lw=2)
    axs[1,1].set_title('vT: Blue') 
    axs[1,1].set_xlabel('vT') 
    axs[1,1].set_ylabel('Count') 
    fig.tight_layout(pad=1.5)
