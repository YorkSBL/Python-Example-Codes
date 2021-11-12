#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EXestimatePI.py
Created on Wed Jun  9 09:13:55 2021
Purpose: Use a random # generator to estimate pi by considering the ratio of areas
 of a circle to a square

@author: CB
"""

import math
import numpy
import random
import matplotlib.pyplot as plt
# -----------
N=1000;      # num of points to use
# -----------
x=[]
y=[]
c=[]
# create array of random (x,y) pairs on square
for n in range(N):
    xt= 2*random.uniform(0,1)-1
    yt= 2*random.uniform(0,1)-1
    x.append(xt)
    y.append(yt)
    # determine if point lies within unit circle
    if math.sqrt(xt**2 + yt**2)<=1:
        c.append(1)

# --- determine pi estimate (and display)
pi_est= 4*len(c)/N
print(f'Estimate of pi = {str(pi_est)} ')
# --- visualize
plt.close("all")
#plt.figure(1)
plt.plot(x,y,'x')
plt.xlim([-1,1]) 
plt.xlabel('x')  
plt.ylabel('y') 
plt.grid()
# --- also draw unit circle
theta= numpy.linspace(0,2*math.pi,100); 
xC= numpy.cos(theta); 
yC= numpy.sin(theta);
plt.plot(xC,yC,'r-')
# --- to save figure
if 1==1:
    plt.savefig("test.pdf")
    
 