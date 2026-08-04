import sys, io, builtins as _bi
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace", line_buffering=True)
sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8", errors="replace", line_buffering=True)
_real_print = _bi.print
_amap = str.maketrans({"\u2190": "<-", "\u2192": "->", "\u2191": "^", "\u2193": "v"})
def _safe_print(*a, **k):
    _real_print(*[s.translate(_amap) if isinstance(s, str) else s for s in a], **k)
_bi.print = _safe_print

'''
Module: ediff.gui.interactive_ellipse
---------------------------------------
Elliptical distortion correction, which needs its own interactive
matplotlib window for the user to click 4 points on a diffraction ring.

``correct_ellipse`` is the reusable entry point -- call it directly from
any Python code. It is also run as a standalone subprocess by
``ediff.gui.pipeline`` (via ``python -m ediff.gui.interactive_ellipse``),
since an interactive matplotlib window needs to own its process's main
thread, which the main GUI process's own Qt event loop is already using.
The ``_cli()`` function below is only the thin command-line adapter for
that specific invocation; it carries no logic of its own.
'''

import numpy as np
import matplotlib
matplotlib.use("QtAgg")


def correct_ellipse(eld_path, icut, cmap, out_path):
    '''
    Interactively correct elliptical distortion and save the result.

    Parameters
    ----------
    eld_path : str
        Path to the (2D) experimental diffraction pattern to correct.
    icut : int
        Maximum intensity clip used when rendering the pattern for
        point-picking.
    cmap : str
        Matplotlib colormap name used to display the pattern.
    out_path : str
        Where to save the corrected image (16-bit PNG).

    Returns
    -------
    None
        The corrected image is written to out_path; nothing is returned.
    '''
    import cv2
    import ediff as ed
    import ediff.io
    import ediff.center

    img = ediff.io.Diffractogram2D(eld_path).astype(float)

    # Pass the FULL, unclipped image so the actual transform and saved
    # output preserve true intensity values. icut is only used as a
    # display-only color-scale clip (vmax) — passing it as before via
    # np.clip(img, 0, icut) baked that clip permanently into the data
    # itself, which then got stretched back to the full 16-bit range on
    # save, amplifying an already-crushed low-contrast signal into a
    # blown-out, mostly-binary-looking result once redisplayed.
    dummy = object.__new__(ed.center.CenterLocator)
    corrected = dummy.ellipse_distortion(
        img, show=False, cmap=cmap, method='manual', vmax=icut)

    # Save preserving the true dynamic range — no renormalization, just
    # clip to a valid 16-bit range (should already be a no-op for
    # legitimate data; this only guards against any out-of-range values
    # introduced by the warp's border fill or interpolation).
    cv2.imwrite(out_path, np.clip(corrected, 0, 65535).astype(np.uint16))


def _cli():
    '''
    Command-line adapter for subprocess invocation.
    Usage: python -m ediff.gui.interactive_ellipse <eld_path> <icut> <cmap>
                  <out_path>
    Prints: ELLIPSE_DONE on success, ELLIPSE_FAILED <reason> on error.
    '''
    import argparse
    import traceback

    parser = argparse.ArgumentParser()
    parser.add_argument("eld_path")
    parser.add_argument("icut",     type=int)
    parser.add_argument("cmap",     type=str)
    parser.add_argument("out_path", type=str)
    args = parser.parse_args()

    try:
        correct_ellipse(args.eld_path, args.icut, args.cmap, args.out_path)
        print("ELLIPSE_DONE", flush=True)
    except Exception as e:
        print(f"ELLIPSE_FAILED: {e}", flush=True)
        traceback.print_exc(file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    _cli()