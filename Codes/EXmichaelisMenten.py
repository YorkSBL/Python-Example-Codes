#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EXmichaelisMenten.py

Purpose: Solve the Michaelis-Menten system (MM) via
the built-in Python solver scipy.integrate.solve_ivp
[this stems from 4030 F21 Lect.13; see notes there]

Notes
o The typical MM curve is for a "steady-state" condition, which
translates to when dc/dt~0. One way to get this is let
kf >> kr > kc. That way, things quickly go to a const. c
and will remain as such until s peters out. If you then zoom
in for Fig.2 on the time region where c~const., you'll see
the typical saturation (e.g., let kf=10, kr=0.05, and kc=0.01)

REFS: 
+ https://chem.libretexts.org/Bookshelves/Biological_Chemistry/Supplemental_Modules_(Biological_Chemistry)/Enzymes/Enzymatic_Kinetics/Michaelis-Menten_Kinetics
+ https://depts.washington.edu/wmatkins/kinetics/index.html
+ http://fcostin.github.io/cmepy/example_mm.html

Created on Tue Oct 26 11:04:02 2021
@author: C. Bergevin
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as sp

# ========= [User Vals]  ==================
# --- rate consts
kf= 10      # forward rate const (i.e., E+S -> C)
kr = 0.05   # reverse rate const(i.e., C -> E+S)
kc= 0.01      # catalyzation rate const (i.e., C -> E+P)
#  --- ICs
S0 = 1         # initial substrate conc. {1}  
E0 = 0.3     # initial enzyme conc. {0.3?}  
C0 = 0       # initial complex conc. {0}  
P0 = 0       # initial product conc. {0}  
# --- 
tMax = 200    # max. time to integrate to {100?}
# ===========================================


# ==== define ODE system
def mm(t,y):
    S, E, C, P = y
    dSdt = -kf*S*E + kr*C 
    dEdt = -kf*S*E + kr*C + kc*C
    dCdt = kf*S*E - kr*C - kc*C
    dPdt = kc*C
    out = np.array([dSdt, dEdt, dCdt, dPdt])
    return out

# ==== bookeeping
t = np.linspace(0,tMax, 100*tMax) 
y0 = [S0,E0,C0,P0]

# ==== run solver & repackage
sol = sp.solve_ivp(mm, [t[0],t[-1]], y0, method='RK45', \
                   t_eval=t, rtol = 1e-6)
S= sol.y[0]
E= sol.y[1]
C= sol.y[2]
P= sol.y[3]
# ==== determine v (= dP/dt)
dt= t[1]-t[0]   # time-step size
v = np.gradient(P,dt)

# ==== Visualize timecourse
plt.close("all")
if 1==1:    
    fig = plt.figure()
    plt.plot(t,S,'r-',label='Substrate')
    plt.plot(t,E,'b-',label='Enzyme')
    plt.plot(t,C,'g-',label='Complex')
    plt.plot(t,P,'k-',label='Product')
    plt.xlabel("t [arb]")
    plt.ylabel("Concentration")
    plt.title("Michaelis-Menten kinetics")
    plt.legend()
    plt.grid()
    plt.show()
    
# ==== Visualize v(=dP/dt)
if 1==1:    
    fig = plt.figure()
    plt.plot(S,v,'r-',label='dP/dt')
    plt.xlabel("Substrate conc.")
    plt.ylabel("dP/dt")
    plt.title("Michaelis-Menten kinetics")
    plt.legend()
    plt.grid()
    plt.show()