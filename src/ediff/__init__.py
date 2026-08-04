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
* Documentation, help, additional notes:
    - Quick start examples/demos - https://mirekslouf.github.io/ediff/docs
    - CIF files can be obtained from www - https://www.crystallography.net/cod

EDIFF modules:

* ediff.bkg = background subtraction for 1D diffraction profiles
* ediff.bkg2d = background subtraction for 2D diffraction patterns  
* ediff.calibration = calibration of SAED diffractograms (pixels -> q-vectors)
* ediff.center = find the center of an arbitrary 2D-diffraction pattern
* ediff.io = input/output operations (read diffractogram, set plot params...)
* ediff.gcryst = functions from geometric/general crystallography
* ediff.gui    = graphical user interface to ediff package
* ediff.mcryst = process monocrystal diffraction patterns
* ediff.pcryst = process polycrystal/powder diffraction patterns
* ediff.radial = calculate the 1D-radial profile from a 2D-diffraction pattern
'''

__version__ = "1.2.7"


# Import of modules so that we could use the package as follows:
# >>> import ediff as ed
# >>> ed.io.Diffractogram1D.show...
# from ediff import bkg
# from ediff import bkg2d
from ediff import center
# from ediff import gcryst
# from ediff import gui
from ediff import io
from ediff import mcryst
# from ediff import pcryst
from ediff import radial

# ediff.bkg, ediff.bkg2d, ediff.gcryst, ediff.pcryst, and ediff.gui are
# imported LAZILY -- only the first time they're actually accessed (e.g.
# `ed.bkg.something(...)` or `ed.pcryst.XRD_polycrystal(...)`) -- rather
# than unconditionally here. Profiling the GUI's interactive worker
# process (which only ever needs ediff.center/io for manual center
# detection, ellipse correction, and background point-picking) showed
# `import ediff` alone cost close to 8 seconds, and the overwhelming
# majority of that -- about 5.8 seconds, ~70% of the total -- came from
# ediff.bkg pulling in the external `bground` package (and its own heavy
# dependency `pybaselines`), completely unconditionally, even for
# operations that never touch background subtraction at all.
# ediff.gcryst and ediff.pcryst separately pull in pymatgen (itself a
# slow import), and ediff.gui pulls in the whole Qt/matplotlib GUI stack
# -- none of which detect_center or correct_ellipse need either.
# This is the standard PEP 562 module-level __getattr__ mechanism:
# `ed.bkg...`, `ed.bkg2d...`, `ed.gcryst...`, `ed.pcryst...`, and
# `ed.gui...` all keep working exactly as before -- only WHEN the
# underlying import happens changes. Operations that genuinely need one
# of these (e.g. subtract_background_interactive genuinely needs
# ediff.bkg, via pcryst.ELD_polycrystal.subtract_background) still pay
# that cost, but only when actually used, not on every cold start
# regardless of which operation is running.
import importlib

_LAZY_SUBMODULES = {"bkg", "bkg2d", "gcryst", "pcryst", "gui"}

def __getattr__(name):
    if name in _LAZY_SUBMODULES:
        module = importlib.import_module(f"ediff.{name}")
        globals()[name] = module   # cache it, so this only runs once
        return module
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

def __dir__():
    return sorted(list(globals().keys()) + list(_LAZY_SUBMODULES))

# List of all ediff modules, which:
# 1) enables >>> from ediff import *
# 2) defines all sub-modules clearly and helps tools like VSCode/Pylance
__all__ = [
    "bkg", "bkg2d", "calibration", "center", "gcryst",
    "gui", "io", "pcryst", "mcryst", "radial"]

# Technical notes - special imports
# ----------------------------------
# 1) ediff.bkg => imports external bground package
# 2) ediff.bkg2d => under development - call funcs in auxiliary IDIFF package
# 3) ediff.gui => graphical user interface - sub-package in its own subdir 

# Obligatory acknowledgement -- the development was co-funded by TACR
# -------------------------------------------------------------------
#  TACR requires that the acknowledgement is printed when we run the program.
#  Nevertheless, Python packages run within other programs, not directly.
# The following code ensures that the acknowledgement is printed when:
#  (1) You run this file: __init__.py
#  (2) You run the package from command line: python -m ediff
# Technical notes:
#  To get item (2) above, we define __main__.py (next to __init__.py).
#  The usage of __main__.py is not too common, but still quite standard.

def acknowledgement():
    print('EDIFF package - process electron diffraction patterns.')
    print('------')
    print('The development of the package was co-funded by')
    print('the Technology agency of the Czech Republic,')
    print('program NCK, project TN02000020.')

def acknowledgementGUI():
    print('EDIFF package GUI - process electron diffraction patterns.')
    print('------')
    print('The development of the package was co-funded by')
    print('the Technology agency of the Czech Republic,')
    print('program NCK, project TN02000020.')

if __name__ == '__main__':
    acknowledgement()
