'''
Package: EDIFF
--------------
Utilities for processing of electron diffraction patterns.

* Input:
    - 2D powder electron diffraction pattern
* Output:
    - 1D powder electron diffraction pattern
    - the 1D pattern/profile is obtained by radial averaging of 2D pattern
    - the final 1D profile can be calibrated and compared with calculated PXRD

EDIFF modules:

* ediff.background = background correction (employs auxilliary package BGROUND)    
* ediff.center = find center of 2D-diffraction pattern
* ediff.io = input/output operations (read diffractogram, set plot params...)
* ediff.pxrd = calculate 1D-PXRD pattern for a known structure
* ediff.radial = calculate 1D-radial profile from 2D-diffraction pattern

Auxiliary package BGROUND:

* the package enables simple, semi-automated, interactive background correction
* it is imported during ediff initialization and accessible as ediff.background
'''

__version__ = "0.2.2"

# Import of modules so that we could use the package as follows:
# >>> import ediff as ed
# >>> ed.io.read_image ...
import ediff.center
import ediff.io
import ediff.pxrd
import ediff.radial

# This is a slightly special import:
# * ediff (1) imports ediff.background, which (2) imports bground package
# * see additional imports in ediff.background module to see what is done 
# * this "two-step import" enables us to use the ediff module as follows:
# >>> import ediff as ed
# >>> DATA  = ed.background.InputData ...
# >>> PPAR  = ed.background.PlotParams ...
# >>> IPLOT = ed.background.InteractivePlot ...
import ediff.background
