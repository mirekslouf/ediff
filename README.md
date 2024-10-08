EDIFF :: processing of powder electron diffraction patterns
-----------------------------------------------------------
* EDIFF program:
	- converts 2D powder electron diffractogram to 1D diffraction profile (ELD)
	- calculates theoretical 1D powder X-ray diffraction profile (XRD)
	- compares the experimental ELD with theoretical XRD
* EDIFF modules:
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
* `pip install bground` = our package, interactive background subtraction
* `pip install ediff` = EDIFF package itself (uses all packages above)

Quick start
-----------
* Look at [worked example](https://www.dropbox.com/scl/fi/3hb78voxd17wb3fzh9n1p/01_ediff_au.nb.pdf?rlkey=qmbvwaw80o1gbe262hwgjvmgx&dl=0)
  to see how EDIFF works.
* Download [complete examples with data](https://www.dropbox.com/scl/fo/td6rkdgp2usxosj1vqeku/h?rlkey=41carfdej5h2f8f4yscbuvagm&dl=0)
  and try EDIFF yourself.

Documentation, help and examples
--------------------------------
* [PyPI](https://pypi.org/project/ediff) repository.
* [GitHub](https://github.com/mirekslouf/ediff) repository.
* [GitHub Pages](https://mirekslouf.github.io/ediff/)
  with [documentation](https://mirekslouf.github.io/ediff/docs).

Versions of EDIFF
-----------------

* Versions 0.0.x = the 1st drafts, testing of {pxrd} module
* Version 0.1 = the 1st semi-complete version with basic documentation
* Version 0.2 = important improvements of {center} module
* Version 0.3 = {calibration} module + improved docs + improved ediff template
* Version 0.4 = TODO: further improvements of {center} module
