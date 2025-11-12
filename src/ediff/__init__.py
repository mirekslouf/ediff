'''
Package: EDIFF
--------------
Processing of electron diffraction patterns.

* Input:
    - Image file = an experimental 2D electron diffraction pattern.
    - CIF file = a description the expected/theoretical crystal structure.
* Output:
    - Comparison of the *experimental* and *theoretical* diffractogram.
    - If the two difractograms are equivalent, the sample has been identified.

EDIFF modules:

* ediff.bkg = background subtraction for 1D diffraction profiles
* ediff.bkg2d = background subtraction for 2D diffraction patterns  
* ediff.calibration = calibration of SAED diffractograms (pixels -> q-vectors)
* ediff.center = find the center of an arbitrary 2D-diffraction pattern
* ediff.io = input/output operations (read diffractogram, set plot params...)
* ediff.gcryst = functions from geometric/general crystallography
* ediff.mcryst = process monocrystal diffraction patterns
* ediff.pcryst = process polycrystal/powder diffraction patterns
* ediff.radial = calculate the 1D-radial profile from a 2D-diffraction pattern
'''

__version__ = "0.8.5"


# Import of modules so that we could use the package as follows:
# >>> import ediff as ed
# >>> ed.io.Diffractogram.read...
import ediff.calibration
import ediff.center
import ediff.gcryst
import ediff.io
import ediff.pcryst
import ediff.radial


# This is a slightly special import:
# * ediff (1) imports ediff.bkg, which (2) imports external bground package
# * see additional imports in ediff.bkg module to see how it is performed 
# * this "two-step import" enables us to use the ediff module as follows:
# >>> import ediff as ed
# >>> DATA  = ed.bkg.InputData ...
# >>> PPAR  = ed.bkg.PlotParams ...
# >>> IPLOT = ed.bkg.InteractivePlot ...
import ediff.bkg


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
