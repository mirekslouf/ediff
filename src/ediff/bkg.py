'''
Module: ediff.bkg
-----------------
Background subtraction for XY-data.    

* This module just imports key objects from external bground package.
* Therefore, it is just a formal incorporation of bground package to ediff.

The source code is brief (= just imports),
but the comments describe how it works (= how it can be used).

* See the source code of ediff.bkg
  if you are interested in technical details concerning the import.
* See documentation of bground package at https://pypi.org/project/bground
  to find out how the background correction works.
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
