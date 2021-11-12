#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EXodeSIR.py

Purpose: Solve a modified SIR model using a hard-coded RK4 AND
the built-in Python solver scipy.integrate.solve_ivp
[this stems from 4030 F21 HW6 as per Hazkeel]

REF:
+ https://www.frontiersin.org/articles/10.3389/fpubh.2020.00230/full   

Created on Thu Oct 21 09:30:17 2021
@author: C. Bergevin
"""

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt


# ========= [User Vals]  ==================
# --- system params
N = 1000000  # population size (1 million)
R0 = 10   # Total number of people an infected person infects
d = 7   # [days]
gamma = 1.0/d   # [1/days]
beta= R0*gamma    # people infected per day by an infected person
delta = 1/7   # incubation period (1/days)
alpha= 0.15    # 15% death rate of the population
rho= 1/14      # 14 days from infection until death
# --- ICs
S0 = N-1  # initial # of suseptible
E0 = 0
I0 = 1
R0 = 0
D0 = 0
# --- integration params
tE = 100        # integration time (assumes tS=0) [days]
SR= 1*10**1     # "sample rate" (i.e., step-size = 1/SR)

# ===========================================

# ----------- [define ODE system]
def f(r,t):
    S = r[0]
    E = r[1]
    I = r[2]
    R = r[3]
    D = r[4]
    dSdt = -beta*I*S/N
    dEdt = beta*I*S/N - delta*E
    dIdt = delta*E - (1-alpha)*gamma*I - alpha*rho*I
    dRdt = (1-alpha)*gamma*I
    dDdt = alpha*rho*I
    return np.array([dSdt,dEdt,dIdt,dRdt,dDdt],float)

# ----------- [define ODE system]
def f2(t,y):
    S, E, I, R, D = y
    dSdt = -beta*I*S/N
    dEdt = beta*I*S/N - delta*E
    dIdt = delta*E - (1-alpha)*gamma*I - alpha*rho*I
    dRdt = (1-alpha)*gamma*I
    dDdt = alpha*rho*I
    out = np.array([dSdt, dEdt, dIdt, dRdt, dDdt])
    return out


# ----------- [define RK4]
def rk4(r,t,h):
    k1 = h*f(r,t)
    k2 = h*f(r+0.5*k1,t+0.5*h)
    k3 = h*f(r+0.5*k2,t+0.5*h)
    k4 = h*f(r+k3,t+h)
    r += (k1+2*k2+2*k3+k4)/6
    return r

# ===== bookkeeping
r0 = [S0,E0,I0,R0,D0]  # package ICs
tS = 0         # start time
Nt = tE*SR     # num of time points
h = (tE-tS)/Nt   # step-size
#tpoints = np.arange(tS,tE,h)   # time point array

tpoints = np.linspace(tS,tE,int(np.round(tE/h)))

S = []; E = []; I = []; R = []; D = [] 


# ==== integration loop for hard-coded RK4
indx = 0
r= r0
for t in tpoints:
    S = np.insert(S,indx,r[0])
    E = np.insert(E,indx,r[1])
    I = np.insert(I,indx,r[2])
    R = np.insert(R,indx,r[3])
    D = np.insert(D,indx,r[4])
    r = rk4(r,t,h)
    indx = indx+1

# ==== now use blackbox solver
#tpointsB = np.arange(tS,tE,1)
sol = sp.integrate.solve_ivp(f2,[0,tE],r0,method='RK45',\
                             t_eval=tpoints, rtol = 1e-6)  
rBB = np.array([sol.y[0],sol.y[1],sol.y[2],sol.y[3],sol.y[4]])
    
# ==== visualize RK4 SOL
plt.close("all")
fig = plt.figure()
plt.plot(tpoints,S,'r-',label='Susceptible')
plt.plot(tpoints,E,'b-',label='Exposed')
plt.plot(tpoints,I,'g-',label='Infected')
plt.plot(tpoints,R,'y-',label='Recovered')
plt.plot(tpoints,D,'k-',label='Dead')
plt.xlabel("t (days)")
plt.ylabel("Population")
plt.title("SEIRD: Hard-coded RK4")
plt.legend()
plt.grid()
plt.show()
# ==== visualize BB vers
fig = plt.figure()
plt.plot(tpoints,rBB[0],'r-',label='Susceptible')
plt.plot(tpoints,rBB[1],'b-',label='Exposed')
plt.plot(tpoints,rBB[2],'g-',label='Infected')
plt.plot(tpoints,rBB[3],'y-',label='Recovered')
plt.plot(tpoints,rBB[4],'k-',label='Dead')
plt.xlabel("t (days)")
plt.ylabel("Population")
plt.title("SEIRD: Blackbox solver (scipy.integrate.solve_ivp)")
plt.legend()
plt.grid()
plt.show()