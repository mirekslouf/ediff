'''
Package: EDIFF
--------------
Utilities for processing of electron diffraction patterns.

* Input:
    - 2D powder electron diffraction pattern
* Output:
    - 1D powder electron diffraction pattern
* Notes:
    - The 1D pattern/profile is obtained by radial averaging of 2D pattern.
    - The 1D profile can be calibrated and compared with the theoretical PXRD.
    - The calculation of theoretical PXRD pattern is a part of this package. 

EDIFF modules:

* ediff.background = background correction (employs auxilliary package BGROUND)    
* ediff.calibration = calibration of SAED diffractograms (pixels -> q-vectors)
* ediff.center = find center of an arbitrary 2D-diffraction pattern
* ediff.io = input/output operations (read diffractogram, set plot params...)
* ediff.pxrd = calculate the 1D-PXRD pattern for a known structure
* ediff.radial = calculate the 1D-radial profile from 2D-diffraction pattern

Auxiliary package BGROUND:

* The is an external package, which enables a 1D background correction.
* It is imported during initialization dnd accesible as ediff.background.
'''

__version__ = "0.3"

# Import of modules so that we could use the package as follows:
# >>> import ediff as ed
# >>> ed.io.read_image ...
import ediff.calibration
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
