'''
Module: ediff.bkg
-----------------
Background subtraction in 1D diffraction profiles.    

* This module just imports key objects from external bground package.
* The source code and documentation are rather brief (basically just imports).
* The comments inside the code describe how it works (= how it can be used).
* Complete documentation of bground package: https://pypi.org/project/bground
'''

# Explanation of the following import commands
# 
# The 1st import command = all modules from bground.api to THIS module
#   - now ediff.background knows the same modules as bground.api
#   - but NOT yet the classes within bground.api - these are imported next
# The 2nd set of import commands = classes from bground.api to THIS module
#   - now ediff.bacground contains the key objects from bground.api
#   - THIS module now contains InputData, PlotParams, InteractivePlot ...
#
# Final conclusion => the users can do:
#
# >>> import ediff.background
# >>> DATA  = ediff.bkg.InputData ...
# >>> PPAR  = ediff.bkg.PlotParams ...
# >>> IPLOT = ediff.bkg.InteractivePlot ...

import bground.api
from bground.api import InputData, BkgParams 
from bground.api import InteractivePlot, RestoreFromPoints
from bground.api import SimpleFuncs, BaseLines, Wavelets
from bground.api import Run, Help
