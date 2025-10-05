'''
Package: EDIFF
--------------
Processing of electron diffraction patterns.

* Input:
    - Image file = an experimental 2D electron diffraction pattern.
    - CIF file = a text file describing the expected crystal structure.
* Processing:
    - Experimental diffractogram is the adjusted/calibrated image file.
    - Theoretical diffractogram is calculated from the CIF file.
* Output:
    - Comparison of the *experimental* and *theoretical* diffractogram.
    - If the two difractograms are equivalent, the sample is identified.

Technical notes:

* CIF files can be downloaded from open-acces databases:
    - Good source of CIF's: https://www.crystallography.net/cod
    - Alternativelly, CIF's can be created within EDIFF or found in www.
* EDIFF provides algorithms for calculation of theoretical diffractograms:
    - Polycrystal/powder diffractograms - finalized, fully working.
    - Monocrystal/spotty diffractograms - under development. 

EDIFF modules:

* ediff.background = background correction (employs sub-package BGROUND)    
* ediff.calibration = calibration of SAED diffractograms (pixels -> q-vectors)
* ediff.center = find the center of an arbitrary 2D-diffraction pattern
* ediff.io = input/output operations (read diffractogram, set plot params...)
* ediff.gcryst = selected useful functions from geometric crystallography
* ediff.mcryst = process monocrystal diffraction patterns
* ediff.pcryst = process polycrystal/powder diffraction patterns
* ediff.radial = calculate the 1D-radial profile from a 2D-diffraction pattern

Auxiliary package BGROUND:

* BGROUND is an external package, which enables a 1D background correction.
* It is imported during initialization to be accesible as ediff.background.
'''

__version__ = "0.8.0"


# Import of modules so that we could use the package as follows:
# >>> import ediff as ed
# >>> ed.io.read_image ...
import ediff.calibration
import ediff.center
import ediff.gcryst
import ediff.io
import ediff.pcryst
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
    print('EDIFF package - process electron diffraction patterns.')
    print('------')
    print('The development of the package was co-funded by')
    print('the Technology agency of the Czech Republic,')
    print('program NCK, project TN02000020.')
    
if __name__ == '__main__':
    acknowledgement()
