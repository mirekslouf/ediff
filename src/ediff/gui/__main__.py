'''
Module: ediff.gui.__main__
---------------------------
Enables launching the GUI as:

    python -m ediff.gui
'''

from pathlib import Path

# Start the interactive worker's cold-start import (numpy/scipy/
# matplotlib/cv2/ediff.center, commonly several seconds) as early as
# absolutely possible -- BEFORE even importing ediff.gui.app, which
# itself has a substantial import chain of its own (matplotlib, and
# PyQt/PySide binding resolution via qtpy). Importing ediff.gui.pipeline
# alone is comparatively light (numpy, pathlib, and ediff's own already-
# lazy-loaded core modules), so this buys the worker extra head start
# equal to however long app.py's own imports take, on top of the
# pre-warm already happening at the very start of main().
from ediff.gui.pipeline import EdiffPipeline

_workdir = Path.home() / ".ediff_workdir"
_workdir.mkdir(exist_ok=True)
_prewarmed_pipeline = EdiffPipeline(_workdir)

from ediff.gui.app import main

if __name__ == "__main__":
    main(prewarmed_pipeline=_prewarmed_pipeline)