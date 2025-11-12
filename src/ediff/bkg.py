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
#  - now ediff.background knows the same modules as bground.api
#   - but NOT yet the classes within bground.api - these are imported next
# The following import commands = classes from bground.api to THIS module
#   - now ediff.bacground contains the three objects from bground.api
#   - THIS module now contains InputData, PlotParams, InteractivePlot ...
#
# Final conclusion => the users can do:
#
# >>> import ediff.background
# >>> DATA  = ediff.background.InputData ...
# >>> PPAR  = ediff.background.PlotParams ...
# >>> IPLOT = ediff.background.InteractivePlot ...

import bground.api
from bground.api import InputData, PlotParams 
from bground.api import InteractivePlot, WaveletMethod
