EDIFF :: processing of powder electron diffraction patterns
-----------------------------------------------------------
* EDIFF is under development, but key modules do work:
    - io = input/output data treatment
	- pxrd = calculate PXRD pattern(s) for the known/expected crystal(s)
	- center = find center of 2D powder electron diffraction pattern
	- radial = calculate radial distribution (2D-pattern &rarr; 1D-pattern) 
	- background = semi-automated background subtraction (external module)
	- calibration = calibration of SAED diffractograms (pixels &rarr; q-vector)
* If you use EDIFF in your research, **please cite** the OpenAccess paper:
	- Materials 14 (2011) 7550.
	  [https://doi.org/10.3390/ma14247550](https://doi.org/10.3390/ma14247550)
	- The paper describes {stemdiff} package, {ediff} is a part of it.

Installation
------------
* Requirement: Python with sci-modules: numpy, matplotlib, scipy, pandas
* `pip install scikit-image` = 3rd party package for advanced image processing 
* `pip install pymatgen` = 3rd party package employed in PXRD calculation
* `pip install bground`= our package, interactive background subtraction
* `pip install ediff` = EDIFF package itself (uses all packages above)

Quick start
-----------
* See how it works:
	- Look at [worked example](https://www.dropbox.com/scl/fi/3hb78voxd17wb3fzh9n1p/01_ediff_au.nb.pdf?rlkey=qmbvwaw80o1gbe262hwgjvmgx&dl=0)
      in Jupyter.
* Try it yourself:
	- Download [complete examples with data](https://www.dropbox.com/scl/fo/td6rkdgp2usxosj1vqeku/h?rlkey=41carfdej5h2f8f4yscbuvagm&dl=0)
	  and scripts and basic instructions.
	- After downloading, unzip it and follow the instructions in *readme* file.

Documentation, help and examples
--------------------------------
* [PyPI](https://pypi.org/project/ediff) repository.
* [GitHub](https://github.com/mirekslouf/ediff) repository.
* [GitHub Pages](https://mirekslouf.github.io/ediff/)
  with [documentation](https://mirekslouf.github.io/ediff/docs).

Versions of EDIFF
-----------------

* Version 0.0.1 = the 1st draft
* Version 0.0.2 = pxrd module works
* Version 0.0.3 = pxrd module works including profiles
* Version 0.0.4 = bground module incorporated + docstrings improved
* Version 0.1   = 1st semi-complete version with basic documentation
* Version 0.2   = important improvements of center.py
* Version 0.3   = calibration module + improved docs + improved ediff template
* Version 0.4   = TODO: further improvements of center.py
