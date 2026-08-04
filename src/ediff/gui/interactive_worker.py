'''
Module: ediff.gui.interactive_worker
--------------------------------------
Persistent worker process for ediff's manual/interactive operations
(center detection, ellipse correction, interactive background
point-picking).

Why this exists
----------------
Each of the three interactive operations needs its own PROCESS because
an interactive matplotlib/Qt window must own its process's main thread,
which the main GUI process's own Qt event loop is already using.

Previously, ``pipeline.py`` spawned a brand-new subprocess (a full,
fresh Python interpreter) for every single interactive click. Each of
those processes then had to ``import ediff`` from scratch -- pulling in
numpy, scipy, skimage, opencv, and matplotlib/Qt every time -- before
the window could even appear. That import cost (roughly one to a few
seconds) was being paid on *every* interactive action, which is exactly
why the app felt much slower than running ediff directly in a console
or Jupyter session: there, the heavy imports only happen once, at
kernel/interpreter startup, and every later call reuses what's already
loaded.

This module fixes that. It is started ONCE (lazily, on the first
interactive action -- see ``EdiffPipeline._ensure_worker`` in
pipeline.py) and stays alive for the rest of the session. It performs
the heavy imports exactly once, then sits in a loop reading one JSON
command per line from stdin, dispatches it to the matching function,
and writes one JSON response per line back on stdout. Every interactive
action after the very first one is then just as fast as calling these
functions directly in an already-warm Jupyter/console session.

Protocol
--------
Request  (stdin,  one JSON object per line):
    {"id": <int>, "op": "center" | "ellipse" | "bkg", "kwargs": {...}}

Response (stdout, one JSON object per line):
    {"id": <int>, "ok": true,  "result": <op-specific>}
    {"id": <int>, "ok": false, "error": "<message>"}

``kwargs`` are passed by keyword directly to the corresponding function:
    "center"  -> detect_center(**kwargs)                   -> [x, y]
    "ellipse" -> correct_ellipse(**kwargs)                 -> null
    "bkg"     -> subtract_background_interactive(**kwargs) -> null

The worker never exits on its own. The parent ends it either by closing
its stdin (clean shutdown, e.g. on app exit) or by killing the process
directly to abort whatever operation is currently in progress (see
``EdiffPipeline.terminate_subprocess``); in the latter case
``pipeline.py`` transparently spawns a fresh worker the next time an
interactive action is needed.
'''

import sys
import io
import json
import importlib
import builtins as _bi
import traceback

# ── Reserve the REAL stdio buffers before anything else touches them ────
# These are the actual pipes connecting this worker to the parent
# process (EdiffPipeline._call_worker). Everything below exists to
# protect them from being closed prematurely.
_real_stdout_buffer = sys.stdout.buffer
_real_stderr_buffer = sys.stderr.buffer


def _isolated_import(module_path, *attr_names):
    """
    Import `module_path` and return the requested attributes from it,
    giving its top-level code a throwaway, ISOLATED stdout/stderr that
    is never touched again afterward.

    Why this is needed
    -------------------
    interactive_center.py, interactive_ellipse.py, and interactive_bkg.py
    EACH wrap and reassign sys.stdout/sys.stderr at import time (for
    console-encoding safety when run standalone as their own process).
    That's harmless when only one of them is ever imported per process —
    the old one-shot-subprocess design. But this worker imports all
    three into ONE persistent process.

    A first attempt just gave the whole import block a single shared
    disposable buffer. That still broke: module A wraps buffer B (fine),
    but the moment module B's own top-level code does
    `sys.stdout = io.TextIOWrapper(sys.stdout.buffer, ...)`, it is
    re-wrapping module A's *same* buffer B a second time. Reassigning
    sys.stdout drops the only reference to module A's wrapper, which
    CPython deallocates immediately (pure refcounting) — closing B as a
    side effect, since TextIOWrapper.close() closes its underlying
    buffer by default. Module C then tries to wrap the now-closed B
    again, and either that construction or a later write raises
    "ValueError: I/O operation on closed file" — the exact crash seen
    in practice. Using one shared disposable buffer only moved the
    problem from the real pipe onto a fake one; it didn't remove the
    chain of repeated re-wrapping that causes it.

    The actual fix is to never let two modules' rewrap chains touch the
    same buffer at all: give EACH module's import its own fresh,
    independent BytesIO that no one else ever sees or wraps again. Its
    wrapper may still end up "poisoned" (wrapping an now-closed buffer)
    by the time this function returns, but that's fine — we only pull
    the plain function objects out of the module and never write
    through that wrapper again.
    """
    sys.stdout = io.TextIOWrapper(io.BytesIO(), encoding="utf-8", errors="replace")
    sys.stderr = io.TextIOWrapper(io.BytesIO(), encoding="utf-8", errors="replace")
    module = importlib.import_module(module_path)
    return tuple(getattr(module, name) for name in attr_names)


(detect_center,) = _isolated_import(
    "ediff.gui.interactive_center", "detect_center")
(correct_ellipse,) = _isolated_import(
    "ediff.gui.interactive_ellipse", "correct_ellipse")
(subtract_background_interactive,) = _isolated_import(
    "ediff.gui.interactive_bkg", "subtract_background_interactive")

# ── Now install our own real, single-owner stdio wrappers ───────────────
# These are built directly from the pristine buffers captured at the very
# top of this file, completely independent of whatever disposable-buffer
# wrapper chain any of the three imports above left behind in
# sys.stdout/sys.stderr — so this construction can never be affected by
# anything that happened during import.
sys.stdout = io.TextIOWrapper(_real_stdout_buffer, encoding="utf-8",
                               errors="replace", line_buffering=True)
sys.stderr = io.TextIOWrapper(_real_stderr_buffer, encoding="utf-8",
                               errors="replace", line_buffering=True)

# Reserve THIS stdout object exclusively as the JSON protocol channel.
# Nothing but our own explicit protocol writes below is ever allowed to
# reach it.
_PROTOCOL_OUT = sys.stdout

_amap = str.maketrans({"\u2190": "<-", "\u2192": "->", "\u2191": "^", "\u2193": "v"})
_real_print = _bi.print
def _safe_print(*a, **k):
    k.setdefault("file", sys.stderr)
    _real_print(*[s.translate(_amap) if isinstance(s, str) else s for s in a], **k)
_bi.print = _safe_print

# From here on, redirect stdout to stderr: any print() from ediff
# internals or from detect_center/correct_ellipse/
# subtract_background_interactive (e.g. verbose progress messages) must
# never land on the protocol channel.
sys.stdout = sys.stderr

_DISPATCH = {
    "center":  detect_center,
    "ellipse": correct_ellipse,
    "bkg":     subtract_background_interactive,
}


def _handle(req):
    '''Run one request dict and return the JSON-serializable response dict.'''
    req_id = req.get("id")
    op = req.get("op")
    kwargs = req.get("kwargs") or {}
    fn = _DISPATCH.get(op)
    if fn is None:
        return {"id": req_id, "ok": False, "error": f"Unknown operation: {op!r}"}
    try:
        result = fn(**kwargs)
        return {"id": req_id, "ok": True, "result": result}
    except Exception as e:
        # Full traceback goes to stderr for debugging; only a concise
        # message goes back over the protocol channel.
        traceback.print_exc(file=sys.stderr)
        return {"id": req_id, "ok": False, "error": str(e)}


def main():
    '''
    Main loop: read one JSON request per line from stdin, respond on the
    protocol channel, forever, until stdin is closed (EOF) -- at which
    point the worker exits cleanly.
    '''
    for line in sys.stdin:
        line = line.strip()
        if not line:
            continue
        try:
            req = json.loads(line)
        except json.JSONDecodeError as e:
            print(json.dumps({"id": None, "ok": False,
                               "error": f"Bad request: {e}"}),
                  file=_PROTOCOL_OUT, flush=True)
            continue
        resp = _handle(req)
        print(json.dumps(resp), file=_PROTOCOL_OUT, flush=True)


if __name__ == "__main__":
    main()
