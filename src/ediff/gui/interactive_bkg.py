import sys, io, builtins as _bi
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace", line_buffering=True)
sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8", errors="replace", line_buffering=True)
_real_print = _bi.print
_amap = str.maketrans({"\u2190": "<-", "\u2192": "->", "\u2191": "^", "\u2193": "v"})
def _safe_print(*a, **k):
    _real_print(*[s.translate(_amap) if isinstance(s, str) else s for s in a], **k)
_bi.print = _safe_print

'''
Module: ediff.gui.interactive_bkg
------------------------------------
Interactive background subtraction, which needs its own matplotlib
window for the user to place background control points.

``subtract_background_interactive`` is the reusable entry point -- call
it directly from any Python code. It is also run as a standalone
subprocess by ``ediff.gui.pipeline`` (via
``python -m ediff.gui.interactive_bkg``), since an interactive
matplotlib window needs to own its process's main thread, which the
main GUI process's own Qt event loop is already using. The ``_cli()``
function below is only the thin command-line adapter for that specific
invocation; it carries no logic of its own.
'''

import matplotlib
matplotlib.use("QtAgg")


def subtract_background_interactive(eld_path, center_file, bkg_file):
    '''
    Interactively subtract background from a 1D diffractogram.

    Parameters
    ----------
    eld_path : str
        Path to the (2D) experimental diffraction pattern.
    center_file : str
        Path to a previously saved center-coordinates file.
    bkg_file : str
        Where to save the chosen background control points.

    Returns
    -------
    None
        The background points are written to bkg_file; nothing is
        returned.
    '''
    import ediff as ed

    eld = ed.pcryst.ELD_polycrystal(description="ELD", diffractogram2D=eld_path)

    eld.find_center(detection=None, refinement=None,
                    icut=200, verbose=0, final_print=False,
                    in_file=center_file)

    eld.calculate_1Ddiffractogram()
    eld.subtract_background("InteractivePlot", bkg_file=bkg_file)


def _cli():
    '''
    Command-line adapter for subprocess invocation.
    Usage: python -m ediff.gui.interactive_bkg <eld_path> <center_file>
                  <bkg_file>
    Prints: BKG_DONE on success, BKG_FAILED <reason> on error.
    '''
    import argparse
    import traceback

    parser = argparse.ArgumentParser()
    parser.add_argument("eld_path")
    parser.add_argument("center_file")
    parser.add_argument("bkg_file")
    args = parser.parse_args()

    try:
        subtract_background_interactive(args.eld_path, args.center_file, args.bkg_file)
        print("BKG_DONE", flush=True)
    except Exception as e:
        print(f"BKG_FAILED {e}", flush=True)
        traceback.print_exc(file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    _cli()
