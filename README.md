EDIFF :: processing powder electron diffraction patterns
--------------------------------------------------------

* EDIFF package:
	- converts a 2D powder electron diffractogram to a 1D diffractogram (ELD)
	- calculates a theoretical 1D powder X-ray diffractogram (XRD)
	- compares the experimental ELD with the theoretically calculated XRD
* If you use EDIFF in your research, **please cite** our recent paper:
	- *Microscopy and Microanalysis* 31, 2025, ozaf045. <br>
	  [https://doi.org/10.1093/mam/ozaf045](https://doi.org/10.1093/mam/ozaf045)


Principle
---------

<img src="https://mirekslouf.github.io/ediff/docs/assets/principle.pptx.png" alt="EDIFF principle" width="600"/>


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
  to see how EDIFF in action.
* Download [complete examples with data](https://www.dropbox.com/scl/fo/td6rkdgp2usxosj1vqeku/h?rlkey=41carfdej5h2f8f4yscbuvagm&dl=0)
  and try EDIFF yourself.


Documentation, help and examples
--------------------------------

* [PyPI](https://pypi.org/project/ediff) repository -
  the stable version to install.
* [GitHub](https://github.com/mirekslouf/ediff) repository - 
  the current version under development.
* [GitHub Pages](https://mirekslouf.github.io/ediff/)
  with [help](https://mirekslouf.github.io/ediff/docs)
  and [complete package documentation](https://mirekslouf.github.io/ediff/docs/pdoc.html/ediff.html).


Versions of EDIFF
-----------------

* Version 0.0 = the 1st drafts, testing of {pxrd} module
* Version 0.1 = the 1st functional version with basic documentation
* Version 0.2 = important improvements of {center} module
* Version 0.3 = {calibration} module + various updates + better ediff template
* Version 0.4 = {center} module: better structure + better funcs + saving coords
* Version 0.5 = consolidation of code and documentation - all stable and working
* Version 0.6 = finalization of {center} module
* Version 1.0 = TODO: complete release with finalized documentation


Acknowledgement
---------------

The development was co-funded by TACR, program NCK,
project [TN02000020](https://www.isibrno.cz/en/centre-advanced-electron-and-photonic-optics).
