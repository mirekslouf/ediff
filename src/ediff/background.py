'''
Module ediff.background
-----------------------
Interactive background definition.    
'''

# The module just imports key objects from external bground module.
# This is a formal incorporation of bground module to ediff.background module.

# How it works? = comment to the next two import commands:
# The 1st import command = all modules from bground.ui to THIS module
#  => now ediff.background knows the same modules as bground.ui
#  => but NOT yet the classes within bground.ui - these are imported next
# The 2nd import command = three key classes from bground.ui to THIS module
#  => now ediff.bacground contains the three objects from bground.ui
#  => THIS module now contains InputData, PlotParams, InteractivePlot
# Now the users of this module can do:
#  >>> import ediff.background
#  >>> DATA  = ediff.background.InputData ...
#  >>> PPAR  = ediff.background.PlotParams ...
#  >>> IPLOT = ediff.background.InteractivePlot ...

import bground.ui
from bground.ui import InputData, PlotParams, InteractivePlot
