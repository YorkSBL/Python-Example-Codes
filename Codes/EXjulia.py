#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EXjulia.py

Purpose: Creates a (color) fractal based upon the Julia set 
[http://mathworld.wolfram.com/JuliaSet.html]

o Uses the difference map (i.e., discrete ODE) 
   z_{n+1} = (z_n)^2 + c  (z,c are both complex)
where z = x+i*y (i.e., some point in the complex plane) and c is 
some complex constant
    
o Borrowed code bits from my EXfractal2.m 
o Had some trouble re translating Matlab indexing tricks to Python, so utilized:
  https://www.learnpythonwithrune.org/numpy-calculate-the-julia-set-with-vectorization/
   + Don't: divT = np.zeros((Npx,Npy))
   + Do:    divT = np.zeros(z.shape, dtype=int)
   + Don't: k = np.ones((Npx,Npy))
   + Do:    k = np.full(c.shape, True, dtype=bool)
o Another reference (not used here):
  https://scipython.com/book/chapter-7-matplotlib/problems/p72/the-julia-set/

o Useful way to check array size: >>> size = k.shape  # size of complex array
o A way to clear variables at the command line: >>> reset -f 

o Some c-vals to try:
--- (successively zooming in)
c= -0.297491-0.641051*i, xB=[-1.35 1.35], yB=[-1.05 1.05], scale=0.005 {BROAD VIEW}
c= -0.297491-0.641051*i, xB=[-0.7 -0.2], yB=[-0.05 0.4], scale=0.0005
c= -0.297491-0.641051*i, xB=[-0.6 -0.48], yB=[0.24 0.36], scale=0.00005
---
c= -0.726895347709+ 0.1888871290438*i, xB=[-1.5 1.5], yB=[-1 1], scale=0.005 {BROAD VIEW}
---
c= -0.123 + 0.745*i, xB=[-1.5 1.5], yB=[-1 1], scale=0.005 {Douady's Rabbit Fractal}

Created on Wed Jun 16 12:55:47 2021
@author: CB
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm


# -------------------------------------
# --- complex const. value (changing leads to different fractals)
#c= complex(-0.28,0.66)
c= complex(-0.8,0.156)
#c= complex(0.285,0.01)
# --- spatial scale of complex plane (x+i*y) to consider
xB= [-1.5, 1.5]; # x-range 
yB= [-1, 1]; # y-range 
scale = 0.005;  # spacing between adjacent points
# --- other params
Nmax = 200;      # # of iterations to run to check for divergence {30}
maxZ= 2;        # threshold |z| to determine divergence {2}
# -------------------------------------
# --- bookkeeping
i = complex(0,1)   # define "i" [i.e., sqrt(-1)]
Npx= round((xB[1]-xB[0])/scale)
Npy= round((yB[1]-yB[0])/scale)


# --- generate a grid of points in the specified complex 
# plane, and determine the associated c value (i.e., simply x+i*y)
x = np.linspace(xB[0],xB[1],Npx)
y = np.linspace(yB[0],yB[1],Npy)
xv, yv = np.meshgrid(x,y)
z = xv + i*yv  # Initialize matricies for iterative loop
#zI= z  
c = np.full(z.shape, c)  # create a z-shaped array of c vals
# --- other key arrays
k = np.full(c.shape, True, dtype=bool)  # boolean array re which points have diverged
divT = np.zeros(z.shape, dtype=int)   # initialize array of divergence vals (to plot!)


# --- loop to iterate the grid and check for divergence (flagging the array k)
# Any k locations for which abs(zI)>maxZ at this iteration and no
# previous iteration get assigned the value of nn (indicates rate of 
# divergence and therefore means to color code)
for nn in range(Nmax):
        z[k] = z[k]**2 + c[k]
        k[np.abs(z) > maxZ] = False
        divT[k] = nn
        
        
# --- visualize
fig1=plt.plot([])  # set up fig
plt.imshow(divT, interpolation='nearest', cmap=cm.jet)
plt.title('Julia set: c = '+ str(c[0,0].real) + ' + i*(' + str(c[0,0].imag) +')')
# save as pdf??
if 1==0:
    plt.savefig('fractal.pdf')

# Note: For other coloramaps, see:
# https://matplotlib.org/stable/tutorials/colors/colormaps.html

