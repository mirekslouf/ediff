EDIFF :: processing of electron diffraction patterns
----------------------------------------------------

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

* Installation of `ediff` is simple:
	- `pip install ediff`
* Additional packages (auto-installed with ediff if not present):
	- Standard sci-libraries: `numpy`, `matplotlib`, `scipy`, `pandas`
	- Additional libs: `bground`, `scikit-image`, `python-opencv`, `pymatgen`


Quick start
-----------

* [Worked example](https://drive.google.com/file/d/1aI0h4YqQLmlQcQiOf4EXsGt-Bq8aANmy/view?usp=drive_link)
  shows the EDIFF package in action.
* [Help on GitHub](https://mirekslouf.github.io/ediff/docs/)
  with complete
  [package documentation](https://mirekslouf.github.io/ediff/docs/pdoc.html/ediff.html)
  and
  [additional examples](https://drive.google.com/drive/folders/1QXs6KyLu9JkMq7UQ-q2TDFPcPIICgoHE?usp=sharing).


Additional resources
--------------------

* [PyPI](https://pypi.org/project/ediff) repository -
  the stable version to install.
* [GitHub](https://github.com/mirekslouf/ediff) repository - 
  the current version under development.
* [GitHub Pages](https://mirekslouf.github.io/ediff/) -
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
* Version 0.8 = rewrite: {gcryst}, {pxrd &rArr; pcryst} + new: {mcryst}
* Version 0.9 = new: {bkg2d} + improved: {calibration} + better documentation
* Version 1.0 = complete processing of powder diffractograms
* Version 1.1 = improved and simplified interface for {pcryst}
* Version 1.2 = automated background subtraction in 1D {bkg}
* Version 1.3 = optional GUI interface
* TODO (upcomming versions):
	- automated background subtraction in 2D {bkg2d}
	- automated peak detection in 2D {peaks} 
	- automated peak identification in 2D, ELD-XRD {mcryst}
	- complete processing of monocrystal diffraction patterns
	

Acknowledgement
---------------

The development was co-funded by TACR, program NCK,
project [TN02000020](https://www.isibrno.cz/en/centre-advanced-electron-and-photonic-optics).
