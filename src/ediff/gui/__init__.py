'''
Package: ediff.gui
-------------------
Graphical user interface for the EDIFF workflow (structure files + a
diffraction image, guided step by step through XRD calculation, ellipse
correction, center detection, radial averaging, background subtraction,
and calibration against theoretical XRD).

This is an *optional* part of the ediff package — install it with:

    pip install ediff[gui]

Launch the GUI in either of these equivalent ways once installed:

    ediff-gui                  # console-script entry point
    python -m ediff.gui        # module invocation

Or from Python:

    >>> import ediff.gui
    >>> ediff.gui.main()

Technical note for maintainers
-------------------------------
Three of this subpackage's steps (ellipse correction, center detection,
and interactive background subtraction) open their own interactive
matplotlib/Qt windows. Those specifically need to own the process's main
thread, so they run as separate subprocesses launched via
``python -m ediff.gui.interactive_center`` (etc.) rather than as plain
function calls inside the main GUI process — see pipeline.py's
_run_interactive() for the launching logic, and each interactive_*.py
module for the corresponding subprocess-side entry point.
'''

def __getattr__(name):
    # Lazy import: importing `ediff.gui` (or any of its submodules, like
    # `ediff.gui.pipeline`) should not require PyQt6 to be installed —
    # only actually using the GUI (i.e. touching `ediff.gui.main`) should.
    if name == "main":
        from ediff.gui.app import main
        return main
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

__all__ = ["main"]
