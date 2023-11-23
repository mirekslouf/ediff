'''
Package: EDIFF
--------------
Utilities for processing of electron diffraction patterns.

* Input:
    - 2D powder electron diffraction pattern
* Output:
    - 1D powder electron diffraction pattern
    - i.e. radial profile of the original pattern with background subtraction
    - the final profile is usually calibrated and compared with calculated PXRD

EDIFF modules:

* ediff.background = background correction (employs auxilliary package BGROUND)    
* ediff.center = find center of 2D-diffraction pattern
* ediff.io = input/output operations (read diffractogram, set plot params...)
* ediff.pxrd = calculate 1D-PXRD pattern for a known structure
* ediff.radial = calculate 1D-radial profile from 2D-diffraction pattern

EDIFF auxiliary package BGROUND and its modules:

* the package enables simple interactive background correction
* it is imported in this file so that it could be used as a part of EDIFF
'''

__version__ = "0.2.1"

# Import of modules so that we could use the package as follows:
# >>> import ediff as ed
# >>> ed.io.read_image ...
import ediff.center
import ediff.io
import ediff.pxrd
import ediff.radial

# This is a slightly special import:
# (import of a module that imports parts of external package bground
# (this "double import" enables us to use the ediff module as follows:
# >>> import ediff as ed
# >>> DATA  = ed.background.InputData ...
# >>> PPAR  = ed.background.PlotParams ...
# >>> IPLOT = ed.background.InteractivePlot ...
import ediff.background
