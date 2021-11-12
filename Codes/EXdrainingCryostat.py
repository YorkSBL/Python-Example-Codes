#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EXdrainingCryostat.py

Purpose: Code to find cubic roots as per the "draining cryostat" example (E6.10)
of Hill (2020). Here we have a spherical cryostat of radius R filled to
some height h0 that is draining at a constant rate F. Upon solving the
associated ODE, on gets a cubic of the form:
a*h^3 + b*h^2 +(0)*h + d =0
where d=d(t)


o basically a Python version of EXdrainingCryostat.m
o HOWEVER, am adding a bit to directly numerically integrate the ODE too

Created on Thu Jun 17 10:48:07 2021
@author: CB
"""
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import sys
Polynomial = np.polynomial.Polynomial

# ----------- [define RK4]
def rk4(r,t,ss):
    k1 = ss*f(r,t)
    k2 = ss*f(r+0.5*k1,t+0.5*h)
    k3 = ss*f(r+0.5*k2,t+0.5*h)
    k4 = ss*f(r+k3,t+h)
    r += (k1+2*k2+2*k3+k4)/6
    return r

# ----------- [define ODE system]
def f(r,t):
    hT = r[0]
    dhdt = F/(math.pi*hT*(hT-2*R))
    return np.array([dhdt],float)

# ---------------
R= 1.5;      # radius of cryostat [m]
h0= 0.1*R;   # initial height of fluid (must be between 0 <2*R [m]
F= 0.0002;      # const. flow rate out [m^3/s]
N= 200;      # # of time points
# ---------------
tMax= (math.pi*h0)/F * (R*h0-(1/3)*h0**2);  # total time to flow out [s]
t = np.arange(0,tMax,tMax/(N-1))  # create time array
t = np.append(t,[tMax])
# --- determine fixed coefficients (as per SOL to ODE)t
a= 1/3;  # cubic term
b= -R;   # quadratic term
c= 0;  # linear term
# ======= Method 1: roots of analytic sol 
h = []
h = np.insert(h,0,h0)  # initialize at initial val.
# --- loop to go through different time vals
for nn in np.arange(1,N,1):
    tV= t[nn]; # grab current time val
    #d= (4/3)*R^3 - (F/pi)*tV;  # const. coeffic. [** full tank case **]
    d= R*h0**2 -(1/3)*h0**3 - (F/math.pi)*tV;
    #coeff= [a, b, c, d]
    coeff= [d,c,b,a];
    p = Polynomial(coeff)
    roots = p.roots()
    ## NOTE: these "roots" are vals of h that satisfy the cubic polynomial.
    ## For us, we must have h be a val between 0 and h0
    new_h= roots[(0 <= roots) & (roots <= 2*R)]
    # when close to zero, a "double root" can occur (so just grab 1st)
    if len(new_h)>1:
        new_h = new_h[0]
        
    h = np.insert(h,nn,new_h)

#sys.exit()   # program kill line

# ======= Method 2: Numerically integrate the ODE
ss = (tMax)/N   # step-size
tpoints = np.arange(0,tMax,ss)   # time point array
hRK = []
r = [h0]
indx = 0
# --- integration loop
for tn in tpoints:
    hRK = np.insert(hRK,indx,r[0])
    r = rk4(r,tn,ss)
    indx = indx+1
    
# --- visualize
fig1=plt.plot([])  # set up fig
fig1 = plt.plot(t,h,'k.',label='Meth. 1: via roots')
fig1 = plt.plot(tpoints,hRK,'r-',label='Meth. 2: via RK4')
fig1=plt.xlabel('r')  
fig1=plt.ylabel('x_n') 
fig1=plt.grid()
fig1=plt.title('SOLs to draining cryostat problem')
fig1=plt.legend(loc=1)

