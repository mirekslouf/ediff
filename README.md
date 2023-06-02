EDIFF :: processing of powder electron diffraction patterns
-----------------------------------------------------------

* EDIFF is under development, but key modules do work:
    - io = input/output data treatment
	- bkgr = background subtraction
	- center = find center of 2D powder diffraction pattern
	- radial = calculate radial distribution (2D-pattern to 1D-pattern) 
	- pxrd = calculation of theoretical powder X-ray diffraction patterns

Quick start
-----------

* See how it works:
	- Look at [worked example](https://mirekslouf.github.io/ediff/docs/examples/ex1_ediff.nb.html)
      in Jupyter.
* Try it yourself:
	- Download and unzip the [complete example with data](https://www.dropbox.com/scl/fo/nmsvdtef7xtmb7r2ku5aa/h?dl=0&rlkey=2evadkk009wp248rp2c3ij2nj).
	- Look at `00readme.txt` and run the example in Jupyter.

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
