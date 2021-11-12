#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EXreadSVG.py

[IN PROGRESS]
Purpose: Fiddle w/ svg files (as used for 3-D printing)

o REFs
+ https://pypi.org/project/svglib/
+ http://svg2stl.com/

NOTE: this code uses the package svglib

Created on Mon Oct 18 16:20:34 2021
@author: C. Bergevin
"""

from svglib.svglib import svg2rlg
from reportlab.graphics import renderPDF, renderPM


fileN= "./Images/rippled_surface.svg"

# ------------------

drawing = svg2rlg(fileN)
renderPDF.drawToFile(drawing, "file.pdf")
#renderPM.drawToFile(drawing, "file.png", fmt="PNG")