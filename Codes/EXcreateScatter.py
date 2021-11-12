#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EXcreateScatter.py
Created on Mon Sep 20 16:07:20 2021
@author: CB

[IN PROGRESS]

Purpose: Motivated by ch.3.8 from WS Cleveland's book Visualizing Data, 
alng w/ his figs.3.45-3.47, we attempt here to:
    o read in data extracted from fig.3.45 via deplot.m
    o simulate other (undertiminable) aspects so to recreate fig.3.46
    o do some basic analysis on such
    o create some associated plots (e.g., box&whisker)
    
NOTES: 
+ For each data point, there are 8x as many data referred to re
the Hersh study. 
+ This code is intended as SOL to HW4 for PHYS 4030 F21
+ Ref. for original study:
    A. H. Hersh. The Effect of Temperature upon the Heterozygotes in 
    the Bar Series of Drosophila. Journal of Experimental Zoology,
    39, 55-71,1924.
+ np.sum(dist)= 823 (i.e., total # of points as from the original 1924 study by Hersh
+ Useful ref re boxplots (& Python):
https://towardsdatascience.com/create-and-customize-boxplots-with-pythons-matplotlib-to-get-lots-of-insights-from-your-data-d561c9883643
+ 

                   
"""

import os
import sys
import numpy as np
#import random
import matplotlib.pyplot as plt

# ===============================================
root= '../SOL/'         # root path to file
fileN = 'extractedPoints.txt'   # file name (of extracted points)

sclSD= 0.75  # scale factor to tighten vert. spread for simulated data
jittD= 0.2   # spread factor for jittered vers.
# ===============================================

# ---- bookeeping re file stuff
fname = os.path.join(root,fileN)
data = np.loadtxt(fname)
Tm = [15, 17, 19, 21, 23, 25, 27, 29, 31] # list of meas. temps
# ---- repackage/label the #s
T = data[:,0]  # temperature [C]
facet = data[:,1]  # facet #
# round to correct for jitter introduced via deplot (that is, both 
# temp and facet # are integers)
T = np.round(T)   
facet= np.round(facet)
# ===== now fiddle somewhat to simulate the distrib. in Fig.3.46
dist= [94, 65, 78, 93, 77, 123, 134, 91,68]  # made up temp sample dist.
# --- first create simulated temp array (kludge?)
Ts= []
Ts= np.append(Ts,15*np.ones(dist[0]))
Ts= np.append(Ts,17*np.ones(dist[1]))
Ts= np.append(Ts,19*np.ones(dist[2]))
Ts= np.append(Ts,21*np.ones(dist[3]))
Ts= np.append(Ts,23*np.ones(dist[4]))
Ts= np.append(Ts,25*np.ones(dist[5]))
Ts= np.append(Ts,27*np.ones(dist[6]))
Ts= np.append(Ts,29*np.ones(dist[7]))
Ts= np.append(Ts,31*np.ones(dist[8]))
# --- now use the sampled data to constrain the simulated data; first 
# sequester diff. temps (NOTE: these are the facet # for a given
# temp as extracted via dplot, not the simulated facet #s, for that 
# see facetS## as created further below
data15= facet[np.where(T==15)]
data17= facet[np.where(T==17)]
data19= facet[np.where(T==19)]
data21= facet[np.where(T==21)]
data23= facet[np.where(T==23)]
data25= facet[np.where(T==25)]
data27= facet[np.where(T==27)]
data29= facet[np.where(T==29)]
data31= facet[np.where(T==31)]

# ==== finally create simulated facet array
facetS= []
#facetS= np.append(facetS,np.round(np.random.normal(5,1,dist[0])))
# --- T= 15 C
facetS15= np.round(np.random.normal(\
    np.mean(data15),sclSD*np.std(data15),dist[0]))
facetS= np.append(facetS,facetS15)
# --- T= 17 C
facetS17= np.round(np.random.normal(\
    np.mean(data17),sclSD*np.std(data17),dist[1]))
facetS= np.append(facetS,facetS17)
# --- T= 19 C
facetS19= np.round(np.random.normal(\
    np.mean(data19),sclSD*np.std(data19),dist[2]))
facetS= np.append(facetS,facetS19)
# --- T= 21 C
facetS21= np.round(np.random.normal(\
    np.mean(data21),sclSD*np.std(data21),dist[3]))
facetS= np.append(facetS,facetS21)
# --- T= 23 C
facetS23= np.round(np.random.normal(\
    np.mean(data23),sclSD*np.std(data23),dist[4]))
facetS= np.append(facetS,facetS23)
# --- T= 25 C
facetS25= np.round(np.random.normal(\
    np.mean(data25),sclSD*np.std(data25),dist[5]))
facetS= np.append(facetS,facetS25)
# --- T= 27 C
facetS27= np.round(np.random.normal(\
    np.mean(data27),sclSD*np.std(data27),dist[6]))
facetS= np.append(facetS,facetS27)
# --- T= 29 C
facetS29= np.round(np.random.normal(\
    np.mean(data29),sclSD*np.std(data29),dist[7]))
facetS= np.append(facetS,facetS29)
# --- T= 31 C
facetS31= np.round(np.random.normal(\
    np.mean(data31),sclSD*np.std(data31),dist[8]))
facetS= np.append(facetS,facetS31)    
# ==== add jitter to simulated data
TsJ= np.random.normal(Ts,jittD)
facetSJ= np.random.normal(facetS,jittD)

# ======== now do linear regression (using bits from EXanscombe.py)
# for both jittered and unjittered vers.
# --- initialize arrays (via numpy) for key vals 
m = np.array([])
b = np.array([])
sigma_m = np.array([])
sigma_b = np.array([])
Rsquared= np.array([])
# --- directly compute best fits for each set
for n in [0,1]:
    if n==0:  # unjittered vers.
        xT = Ts   # x-vals
        yT = facetS # y-vals
    else:   # jittered vers
        xT = TsJ   # x-vals
        yT = facetSJ # y-vals
        
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
    # --- store away vals 
    m = np.insert(m,n,mT)
    b = np.insert(b,n,bT)
    sigma_m = np.insert(sigma_m,n,sigma_mT)
    sigma_b = np.insert(sigma_b,n,sigma_bT)
    Rsquared = np.insert(Rsquared,n,RsquaredT)

# === determine means & variance for each temp (used to plot later)
mu15= np.average(facetS15)
sd15= np.std(facetS15)
mu17= np.average(facetS17)
sd17= np.std(facetS17)
mu19= np.average(facetS19)
sd19= np.std(facetS19)
mu21= np.average(facetS21)
sd21= np.std(facetS21)
mu23= np.average(facetS23)
sd23= np.std(facetS23)
mu25= np.average(facetS25)
sd25= np.std(facetS25)
mu27= np.average(facetS27)
sd27= np.std(facetS27)
mu29= np.average(facetS29)
sd29= np.std(facetS29)
mu31= np.average(facetS31)
sd31= np.std(facetS31)

# --- compile vals (kludge)
mu= [mu15,mu17,mu19,mu21,mu23,mu25,mu27,mu29,mu31]
sd= [sd15,sd17,sd19,sd21,sd23,sd25,sd27,sd29,sd31]


# =================================================
# [visualize and print out #s as needed]

# ===== recreate Fig.3.45
fig1, ax1 = plt.subplots()
ax1.plot(T,facet,'ko')
ax1.set_xlabel('Temperature [C]')  
ax1.set_ylabel('Facet number') 
ax1.grid()
ax1.set_title('Re-creation of Cleveland (1993) Fig.3.45') 
ax1.set_xlim([14.5, 31.5])
ax1.set_ylim([-17, 9]) 
# ===== recreate Fig.3.46
fig2, ax2 = plt.subplots()
ax2.scatter(TsJ,facetSJ,s=15, facecolors='none', edgecolors='k')
ax2.set_xlabel('Jittered temperature [C]')  
ax2.set_ylabel('Jittered facet number') 
ax2.grid()
ax2.set_title(os.path.join('Re-creation of Cleveland (1993) Fig.3.46')) 
ax2.set_xlim([14.5, 31.5])
ax2.set_ylim([-17, 9])  
# ===== recreate Fig.3.47
fig3, ax3 = plt.subplots()
# --- repackage as a list to feed into boxplot (Note: need to use
# np.array along w/ dtype=object so to avoid a warning)
dataBP= np.array([facetS15,facetS17,facetS19,facetS21,facetS23,
            facetS25,facetS27,facetS29,facetS31],dtype=object)
# ---
ax3.boxplot(dataBP, notch="True")
# --- relabel x-axis (likely kludge)
positions = [1, 2, 3, 4, 5, 6, 7, 8, 9]
labels = ["15","17","19","21","23","25","27","29","31"]
plt.xticks(positions, labels)  
ax3.set_xlabel('Temperature [C]')  
ax3.set_ylabel('Facet number') 
ax3.grid()

# --- print relevant vals to screen
#print(f'Linear fit to unjittered sim. data: y={str(m[0])}*x + {str(b[0])} (R^2 = {str(Rsquared[0])})')
#print(f'Linear fit to jittered sim. data: y={str(m[1])}*x + {str(b[1])} (R^2 = {str(Rsquared[1])})')

print(f'Linear fit to unjittered sim. data: y={str(m[0])}*x + {str(b[0])} (R^2 = {str(Rsquared[0])})')
print(f'Linear fit to jittered sim. data: y={str(m[1])}*x + {str(b[1])} (R^2 = {str(Rsquared[1])})')

# ===== plot means+ variance along w/ regression curve
fig4, ax4 = plt.subplots()
ax4.scatter(TsJ,facetSJ,s=1,marker='o',facecolors='none', \
            edgecolors='0.8')
# --- plot means & STD (as errorbars)
ax4.plot(Tm,mu,'r',marker='s',linestyle='None')
fig4=plt.errorbar(Tm,mu,sd,ecolor='r',linestyle='None')
 # --- plot regression curve
ax4.plot(Tm,(np.float64(Tm)*m[0]+b[0]))
ax4.set_xlabel('Temperature [C]')  
ax4.set_ylabel('Facet number')  
ax4.legend(['means & var (non-jitt)','regression (non-jitt)','jittered sim. data'])
fig4=plt.title('Regression re non-jittered data')
#sys.exit() 
    
    
    
    
    
    
    