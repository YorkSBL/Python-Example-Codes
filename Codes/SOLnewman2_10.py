#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  9 16:25:26 2021
@author: CB

o A complete version in Matlab can be found in SOLnewmanEX2_10.m

"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator

import seaborn as sns; sns.set_theme()

# ==================================
# [User Params]
A= 58  # mass #
Z= 28  # atomic #
# ==================================

# --- other consts
a1= 15.8
a2= 18.3
a3= 0.714
a4= 23.2
# ====== Parts A & B: determine B, approx. nuclear binding 
# energy [10^6 eV] for fixed Z and A
# --- determine a5
if (A%2 != 0): 
    a5 = 0
elif (A%2 == 0) and (Z%2 == 0):
    a5 = 12.0
else:
    a5 = -12.0
# --- determine associated B  
B = a1*A - a2*(A**(2/3)) - a3*(Z**2)/(A**(1/3)) - a4*((A-2*Z)**2)/A \
    + a5/(A**(1/2))
# --- print val to screen
print("=== Parts A & B ===")
print("For mass # A = ",A, " and atomic # Z = ",Z,":")
print("--> Approx. nuclear binding energy B= ",round(B,1)," [MeV]")
print("--> Binding energy per nucleon (B/A)= ",round(B/A,1)," [MeV]")
# ====== Part C: now determine B as a function of A
print("=== Part C ===")
print("For atomic # Z = ",Z,":")
# --- determine range of A vals to compute B over
Ar= np.linspace(Z,3*Z,num=2*Z+1);
# --- determine appropriate val. of a5 (since A varies here)
# NOTE: Python (Ar%2==0)*1 is equiv. to Matlab's ~isodd(Ar)
a5c= 12*np.ones(len(Ar))* ((Ar%2==0)*1)* ((Z%2==0)*2-1)
# --- now calculate corresponding B vals.
Br= a1*Ar - a2*(Ar**(2/3)) - a3*(Z**2)/(Ar**(1/3)) - \
    a4*((Ar-2*Z)**2)/Ar + a5c/(Ar**(1/2))
BrM= np.max(Br/Ar)  # determine max B val. PER NUCLEON
indx= np.where((Br/Ar) == np.amax(Br/Ar))  # find assoc. array index
Amax= Ar[indx]
# NOTE: Python indexing starts at 0 and Matlab at 1
print("--> Mass # re largest binding energy =  ",Amax[0])
print("--> Assoc. binding energy per nucleon = ",round(BrM,2)," [MeV]")
# ====== Part : now determine B as a function of Z
print("=== Part D ===")
Amax= 250  # # of possible atomic #s to scan through for each Z
Zmax= 100  # # of possible mass #s to scan through for each A
Zr= np.linspace(0,99,num=Zmax);
# --- initialize storage arrays
BrdM= [];  indxd= [];  AmaxD= []
Bplot= np.empty((Zmax,Amax)); 
Nnum= np.empty((Zmax,Amax)); 
Znum= np.empty((Zmax,Amax));  
# --- now loop through range of Z vals. and for each, determine 
# largest B across range of A
for n in Zr:
    # --- determine appropriate range of A vals for computation
    # note that Zr(n)=n+1 (so we could simplify below)
    Zt= Zr[int(n)]+1
    Ard= np.linspace(1,Amax,num=Amax);
   # --- determine appropriate val. of a5
    a5d= 12*np.ones(len(Ard))* ((Ard%2==0)*1)* ((Zt%2==0)*2-1)
   # --- calc. all assoc. binding energies
    Brd= a1*Ard - a2*(Ard**(2/3)) - a3*(Zt**2)/(Ard**(1/3)) - \
    a4*((Ard-2*Zt)**2)/Ard + a5d/(Ard**(1/2))
    BrdM.append(np.max(Brd/Ard))  # determine max B val. PER NUCLEON
    #BrdM[int(n)]= np.max(Brd/Ard)  # determine max B val. PER NUCLEON
    indxd.append(np.where((Brd/Ard) == np.amax(Brd/Ard))[0][0])  # find assoc. array index
    AmaxD.append(Ard[indxd[int(n)]])  # determine corresponding mass #
    # --- store away #s for later plotting (note that the row # keeps track
    # of Z as per the for loop, while the columns represent atomic # as per Ar)
    # NOTE: am thinking Matlab-style here, which makes this tricky...
    Bplot[int(n),:]= Brd/Ard  # bind. enrg. per nucleon for a given Z
    Nnum[int(n),:]= Ard-Zr[int(n)]  # keep track of # of neutrons re that Z
    Znum[int(n),:]= Zr[int(n)]*np.ones(len(Ard))
    


# --- find desired max val
maxB= np.max(BrdM)
indxd2= np.where((BrdM) == np.amax(BrdM))
print("--> Mass # re Max. stability = ",int(Zr[indxd2][0]+1))



# ==== plot bound. energy/nucleon vs N (i.e., # of neutrons)
# NOTE: Can't seem to get this last bit working (though I think 
# the #s are all correct)

# ------ CB solution (initial attempt)
"""
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
surf = ax.plot_surface(Nnum,Znum,Bplot, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

ax.set_xlim(0, 140)
ax.set_ylim(0, 100)
ax.set_zlim(0, 9)
c= plt.colorbar(surf, shrink=0.5, aspect=5)

ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')
#ax.view_init(elev=90., azim=0)
plt.show()
"""

# ------ plot bound. energy/nucleon vs N (i.e., # of neutrons)
# (thanks to Ben M. for the syntax!)
"""
fig, ax = plt.subplots()
im = ax.imshow(Bplot, origin='lower', cmap='viridis', vmin=0, vmax=9)
ax.set_xticks(range(0,140+1,20))
ax.set_xlim((0,140)) 
ax.set_xlabel('Number of Neutrons $N$')
ax.set_yticks(range(0,90+1,20))
ax.set_ylim((0,90))
ax.set_ylabel('Number of Protons $Z$')
cbar = plt.colorbar(im)
cbar.set_label('Binding Energy per Nucleon (MeV)', rotation=90) 
ax.grid(False)
plt.show()
"""

# ------ [alternative approach via Elham; ]

Zp = range(0, 100)
Np = range(0, 250)
plt.figure(figsize=(14, 9))
XX, YY = np.meshgrid(Np, Zp)
#be = binding_energy_per_nucleon(XX + YY, YY)
pcolor = plt.pcolormesh(XX,YY,Bplot,shading='auto',cmap='viridis')
plt.xlabel('Number of neutrons N')
plt.ylabel('Number of protons Z')
plt.clim(0, 10)
cbar = plt.colorbar(pcolor,label = 'Binding energy per nucleon (MeV)')
contour = plt.contour(XX,YY,Bplot, levels=[0, 0.8, 4.8, 6.8, 7.8, 8.2, 8.5, 8.7], vmin=0, cmap="inferno")
cbar.add_lines(contour)






