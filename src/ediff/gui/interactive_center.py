import sys, io, builtins as _bi
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace", line_buffering=True)
sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8", errors="replace", line_buffering=True)
_real_print = _bi.print
_amap = str.maketrans({"\u2190": "<-", "\u2192": "->", "\u2191": "^", "\u2193": "v"})
def _safe_print(*a, **k):
    _real_print(*[s.translate(_amap) if isinstance(s, str) else s for s in a], **k)
_bi.print = _safe_print

'''
Module: ediff.gui.interactive_center
-------------------------------------
Center detection, including the (semi-)interactive detection/refinement
methods that need their own matplotlib/Qt window.

``detect_center`` is the reusable entry point -- call it directly from any
Python code. It is also run as a standalone subprocess by
``ediff.gui.pipeline`` (via ``python -m ediff.gui.interactive_center``),
since an interactive matplotlib window needs to own its process's main
thread, which the main GUI process's own Qt event loop is already using.
The ``_cli()`` function below is only the thin command-line adapter for
that specific invocation; it carries no logic of its own.
'''

import matplotlib
matplotlib.use("QtAgg")


def detect_center(eld_path, detection, refinement, icut,
                  in_file=None, out_file=None):
    '''
    Run center detection and return the detected coordinates.

    Parameters
    ----------
    eld_path : str
        Path to the (2D) experimental diffraction pattern.
    detection : str or None
        Detection method name (see ediff.center), or None to skip
        detection and load coordinates from in_file instead.
    refinement : str or None
        Refinement method name, or None to skip refinement. Passing
        'manual' opens an interactive matplotlib window for keyboard-based
        adjustment.
    icut : int
        Maximum intensity clip used when rendering the pattern.
    in_file : str, optional
        If given, load a previously saved center from this file instead
        of detecting one.
    out_file : str, optional
        If given, save the final center coordinates to this file.

    Returns
    -------
    (x, y) : tuple of float
        The detected/refined center coordinates.

    Raises
    ------
    RuntimeError
        If center detection does not produce a result.
    '''
    import ediff as ed

    eld = ed.pcryst.ELD_polycrystal(description="ELD", diffractogram2D=eld_path)

    kwargs = dict(detection=detection, refinement=refinement,
                  icut=icut, verbose=1, final_print=False)
    if in_file:
        kwargs["in_file"] = in_file
    if out_file:
        kwargs["out_file"] = out_file

    eld.find_center(**kwargs)

    if not (hasattr(eld, "center") and eld.center is not None):
        raise RuntimeError("no_center")
    try:
        return eld.center.x, eld.center.y
    except AttributeError:
        raise RuntimeError("no_coords")


def _cli():
    '''
    Command-line adapter for subprocess invocation.
    Usage: python -m ediff.gui.interactive_center <eld_path> <detection>
                  <refinement> <icut> [--in_file <path>] [--out_file <path>]
    Prints: CENTER x y on success, CENTER_FAILED <reason> on error.
    '''
    import argparse
    import traceback

    parser = argparse.ArgumentParser()
    parser.add_argument("eld_path")
    parser.add_argument("detection")
    parser.add_argument("refinement")
    parser.add_argument("icut", type=int)
    parser.add_argument("--in_file",  default=None)
    parser.add_argument("--out_file", default=None)
    args = parser.parse_args()

    detection  = None if args.detection  == "None" else args.detection
    refinement = None if args.refinement == "None" else args.refinement
    in_file  = args.in_file  if args.in_file  and args.in_file  != "None" else None
    out_file = args.out_file if args.out_file and args.out_file != "None" else None

    try:
        x, y = detect_center(args.eld_path, detection, refinement, args.icut,
                             in_file=in_file, out_file=out_file)
        print(f"CENTER {x} {y}", flush=True)
    except Exception as e:
        print(f"CENTER_FAILED {e}", flush=True)
        traceback.print_exc(file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    _cli()
