EDIFF :: processing electron diffraction patterns
-------------------------------------------------

* EDIFF package:
	- processes an experimental 2D electron diffractogram (ELD)
	- calculates a theoretical X-ray diffractogram (XRD)
	- compares ELD and XRD to identify the structure
* Powder vs. monocrystal diffractograms:
	- powder/polycrystal diffractograms - finished, working (see examples below)
	- monocrystal diffractograms - under development (not fully implemented yet)
* If you use EDIFF in your research, **please cite** our recent paper:
	- *Microscopy and Microanalysis* 31, 2025, ozaf045. <br>
	  [https://doi.org/10.1093/mam/ozaf045](https://doi.org/10.1093/mam/ozaf045)


Principle
---------

<img src="https://mirekslouf.github.io/ediff/docs/assets/principle.pptx.png" alt="EDIFF principle" width="600"/>


Installation
------------

* Requirement: Python with sci-modules: numpy, matplotlib, scipy, pandas, pillow
* `pip install scikit-image` = 3rd party package for advanced image processing 
* `pip install pymatgen` = 3rd package employed in crystallographic calculations
* `pip install bground` = our package, interactive background subtraction
* `pip install ediff` = EDIFF package itself (uses all packages above)


Quick start
-----------

* [Worked example](https://drive.google.com/file/d/19C6ClDDepcpNNjl12OduE_d_WpFr0HMd/view?usp=sharing)
  shows the EDIFF package in action.
* [Help on GitHub](https://mirekslouf.github.io/ediff/docs/)
  with complete
  [online documentation](https://mirekslouf.github.io/ediff/docs/pdoc.html/ediff.html)
  and
  [additional examples](https://drive.google.com/drive/folders/1QXs6KyLu9JkMq7UQ-q2TDFPcPIICgoHE?usp=sharing).


Other resources
---------------

* [PyPI](https://pypi.org/project/ediff) repository -
  the stable version to install.
* [GitHub](https://github.com/mirekslouf/ediff) repository - 
  the current version under development.
* [GitHub Pages](https://mirekslouf.github.io/ediff/)
  with 
  the more user-friendly version of GitHub website.

  
Versions of EDIFF
-----------------

* Version 0.0 = the 1st draft, testing of {pxrd} module
* Version 0.1 = the 1st functional version with basic documentation
* Version 0.2 = important improvements of {center} module
* Version 0.3 = {calibration} module + various updates + better ediff template
* Version 0.4 = {center} module: better structure + better funcs + saving coords
* Version 0.5 = consolidation of code and documentation - all stable and working
* Version 0.6 = finalization of {center} module
* Version 0.7 = {gcryst} module: geometrical crystallography, zone axes ...
* Version 0.8 = rewrite of modules: {gcryst}, {pxrd &rArr; pcryst}, {mcryst}
* Version 1.0 = TODO: complete processing of powder diffractograms
* Version 2.0 = TODO: complete processing of powder + monocrystal diffractograms


Acknowledgement
---------------

The development was co-funded by TACR, program NCK,
project [TN02000020](https://www.isibrno.cz/en/centre-advanced-electron-and-photonic-optics).
