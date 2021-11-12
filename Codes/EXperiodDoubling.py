#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EXperiodDoubling.py
[IN PROGRESS]

Purpose: Create/plot the period doubling cascade for the logistic map

o Using EXlogisticBIF.m as a template

Created on Mon Jun 14 13:35:50 2021
@author: CB
"""

import math
import numpy as np
import random
import matplotlib.pyplot as plt

# ----------- [User Params]
# --- params
range= [2,4]   # min and max values to compute bifurcation diagram over [2,4]
Nr= 200    # numb of steps over range + 1 [100]
x0= 0.1    # starting x value [0.1]
Nsettl= 50 # numb of runs allowed for 'settling' [50]
Nit= 100   # numb of iterations to plot for a given value of r [200]
# -----------

# ---
rmin = range[0]
rmax = range[1] 
h = (rmax-rmin)/Nr   # step-size
r = np.arange(rmin,rmax,h) 
r_indx = np.arange(0,Nr,1)
m_indx = np.arange(0,Nsettl+Nit,1)
fig1=plt.plot([])  # set up fig
# ---
# loop through each r value
for n in r_indx:
    rT = r[n]   # grab r-val
    x = x0      # reset to IC
    xS = []    # reset array to hold iterates (for given r)
    xS = np.insert(xS,0,x)
    indx = 2
    # loop through the iterations of the map
    for m in m_indx:
        x = rT*x*(1-x)    # deal w/ mapping
        xS = np.insert(xS,m+1,x)  # store vals
        indx = indx+1    # update indexer

    # plot points for a given iteration *past* the settling time
    xSb =  xS[Nsettl:Nsettl+Nit]
    rTa = rT*np.ones(len(xSb))
    fig1 = plt.plot(rTa,xSb,'k.')


# --- doll up the plot
fig1=plt.xlabel('r')  
fig1=plt.ylabel('x_n') 
fig1=plt.grid()
fig1=plt.title('Logistic Map Bifurcation Diagram [x_{n+1} = r*x_n*(1-x_n)]') 
#fig1=plt.ylim([0.494, 0.508])


