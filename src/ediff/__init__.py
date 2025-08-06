'''
Package: EDIFF
--------------
Processing of powder electron diffraction patterns.

* Input:  2D powder electron diffraction pattern (raw experimental data).
* Output: 1D powder electron diffraction pattern (final, compared with PXRD).
    - The 1D pattern/profile is obtained by radial averaging of 2D pattern.
    - The 1D profile is calibrated and compared with the theoretical PXRD.
    - The calculation of theoretical PXRD patterns is a part of this package. 

EDIFF modules:

* ediff.background = background correction (employs sub-package BGROUND)    
* ediff.calibration = calibration of SAED diffractograms (pixels -> q-vectors)
* ediff.center = find center of an arbitrary 2D-diffraction pattern
* ediff.io = input/output operations (read diffractogram, set plot params...)
* ediff.pxrd = calculate the 1D-PXRD pattern for a known structure
* ediff.radial = calculate the 1D-radial profile from 2D-diffraction pattern

Auxiliary package BGROUND:

* BGROUND is an external package, which enables a 1D background correction.
* It is imported during initialization to be accesible as ediff.background.
'''

__version__ = "0.6.2"


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


# Obligatory acknowledgement -- the development was co-funded by TACR.
#  TACR requires that the acknowledgement is printed when we run the program.
#  Nevertheless, Python packages run within other programs, not directly.
# The following code ensures that the acknowledgement is printed when:
#  (1) You run this file: __init__.py
#  (2) You run the package from command line: python -m ediff
# Technical notes:
#  To get item (2) above, we define __main__.py (next to __init__.py).
#  The usage of __main__.py is not very common, but still quite standard.

def acknowledgement():
    print('EDIFF package - process powder electron diffraction patterns.')
    print('------')
    print('The development of the package was co-funded by')
    print('the Technology agency of the Czech Republic,')
    print('program NCK, project TN02000020.')
    
if __name__ == '__main__':
    acknowledgement()
