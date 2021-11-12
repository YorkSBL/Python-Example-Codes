#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EXconvert2binary.py
Created on Mon Sep 13 13:42:18 2021
@author: CB

Purpose: Take in a user-specified integer and determine it's digital 
form in binary as well as # of bits required to encode it

o REF
https://www.rapidtables.com/convert/number/decimal-to-binary.html
"""

import math
import numpy as np

# ==================================
# [User Params]
N= 24321   # integer to encode as binary
# ==================================
# ---
n=0  # initialize indexer re # of bits
Nn= N  # initialize val. for loop
bit=str()   # initialize string to keep track of binary rep
while (Nn//2)>0:
     bit+= str(Nn%2)  # append string
     Nn= Nn//2   # Python's "floor" division

# --- add a leading 1 (which would  correspond to the highest bit #, 
# or equiv. the 0 quotient (kludge as the loop should do this!)
bit= bit+str(1)
# --- need to flip it so "largest" bit is on the right (using "reverse
# slice" syntax of Python)
bit= bit[::-1]

# --- display out various quantities
print("------")
print("Integer = ",N)
print("binary version = ",bit)
print("# of required bits = ",len(bit))
print("binary vers. via bin = ",bin(N))
print("shortened vers. via bin = ",int(bin(N)[2:]))
print("decimal vers. converted back via int =  ",int(bit,2))

