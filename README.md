EDIFF :: processing of powder electron diffraction patterns
-----------------------------------------------------------
* EDIFF is under development, but key modules do work:
    - io = input/output data treatment
	- background = background subtraction
	- center = find center of 2D powder diffraction pattern
	- radial = calculate radial distribution (2D-pattern &rArr; 1D-pattern) 
	- pxrd = calculation of theoretical powder X-ray diffraction patterns

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
	- Look at [worked example](https://www.dropbox.com/scl/fi/p8mi2jnhrn19qyzzsbue4/01_ediff.nb.html.pdf?rlkey=qny4agmi2osb0k8w6olpeh3wm&dl=0)
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

* Version 0.0.1 = just draft
* Version 0.0.2 = pxrd module works
* Version 0.0.3 = pxrd module works including profiles
* Version 0.0.4 = bground module incorporated + slightly improved docstrings
* Version 0.1.0 = 1st semi-complete version with basic documentation
* Version 0.1.1 = v.0.1.0 + improved/simplified outputs
* Version 0.1.2 = v.0.1.1 + small improvements of code and documentation
* Version 0.2   = important improvements of center.py
* Version 0.2.1 = consolidation, update of docs and examples on www
