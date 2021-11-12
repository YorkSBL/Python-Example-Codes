#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EXplotSOAE.py

[IN PROGRESS]
Purpose: Load in an SOAE txt file and plot

o Modeled after .../Cscript/Library/PlotSOAE5.m

Created on Wed Jun 23 10:07:44 2021
@author: CB
"""

import os
import numpy as np
import matplotlib.pyplot as plt

# -------------------------------- 
root= '../Files/'         # root path to file
fileN = 'ACsb18learSOAEsuppK3.txt'   # file name

# -------------------------------- 

# ==== bookeeping
fname = os.path.join(root,fileN)
data = np.loadtxt(fname)
# ==== 
freq = data[:,0]
mag = data[:,1]
# ====
fig1 = plt.plot(freq/1000,mag,'k-')
fig1=plt.xlabel('Frequency [kHz]')  
fig1=plt.ylabel('Magnitude [dB SPL]') 
fig1=plt.grid()
fig1=plt.title(os.path.join('SOAE spectrum: ',fileN)) 
fig1=plt.xlim([0, 5.5])
fig1=plt.ylim([-25, 10])
