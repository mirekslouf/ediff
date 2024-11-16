'''
Module: ediff.background
------------------------
Semi-automated background correction.    

* This module just imports key objects from external bground package.
* Therefore, it is just a formal incorporation of bground package to ediff.

The source code is brief (just imports),
but it describes how it works (how it can be used).

* See the source code of ediff.background
  if you are interested in technical details concerning the import.
* See documentation of bground package at https://pypi.org/project/bground
  to find out how the background correction works.
'''

# Explanation of the following two import commands
# 
# The 1st import command = all modules from bground.ui to THIS module
#  - now ediff.background knows the same modules as bground.ui
#   - but NOT yet the classes within bground.ui - these are imported next
# The 2nd import command = three key classes from bground.ui to THIS module
#   - now ediff.bacground contains the three objects from bground.ui
#   - THIS module now contains InputData, PlotParams, InteractivePlot
#
# Final conclusion => the users can do:
#
# >>> import ediff.background
# >>> DATA  = ediff.background.InputData ...
# >>> PPAR  = ediff.background.PlotParams ...
# >>> IPLOT = ediff.background.InteractivePlot ...

import bground.ui
from bground.ui import InputData, PlotParams, InteractivePlot
