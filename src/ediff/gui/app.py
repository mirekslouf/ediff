"""
EDIFF GUI - Qt desktop application (via qtpy, works with PyQt5, PyQt6, or
PySide6 - whichever is available in the current environment)
Run with: python -m ediff gui
"""

import sys, os, shutil, json, math, ast
from pathlib import Path

# No QT_API override here on purpose: qtpy auto-detects whichever Qt
# binding is actually installed (PyQt5, PyQt6, PySide2, PySide6). Forcing
# a specific one here would defeat the point of using qtpy at all - if the
# person's environment already has a working PyQt5 (e.g. alongside
# Spyder, which also uses qtpy), this lets the GUI use that instead of
# requiring PyQt6 specifically. If someone needs to force a particular
# binding, they can still set QT_API themselves before launching.

import matplotlib
matplotlib.use("QtAgg")
import matplotlib.pyplot as plt

from qtpy.QtWidgets import (
    QApplication, QMainWindow, QWidget, QSplitter,
    QVBoxLayout, QHBoxLayout, QScrollArea,
    QLabel, QPushButton, QComboBox, QLineEdit,
    QSpinBox, QDoubleSpinBox, QFileDialog,
    QSizePolicy, QFrame, QTabWidget, QMessageBox,
    QSplashScreen, QStackedWidget, QProgressBar, QDialog, QMenu,
    QListWidget,
)

from qtpy.QtCore import Qt, QThread, Signal, QObject, QTimer, QSize, QLocale
from qtpy.QtGui import QFont, QPixmap, QPainter, QColor, QLinearGradient, QShortcut, QKeySequence, QIcon

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavToolbar
from matplotlib.figure import Figure

from ediff.gui.pipeline import EdiffPipeline, default_workdir

# The app bunddles a handful of image files that need to be loaded as QPixmap
# regardless of how is the GUI/app/ediff running. A naive implementation breaks
# (e.g., open(hardoced_path))

def _resource_bytes(name: str) -> bytes:
    """
    Return the raw bytes of a bundled GUI resource file (e.g. 'ediff.png'),
    working in every context this module can run in:
      - pip-installed package (pip install ediff[gui])
      - running directly from a source checkout
      - frozen into a standalone executable via PyInstaller

    Bytes rather than a filesystem path are used deliberately: if the
    package ever ends up inside a zipped wheel, a path handed out by
    importlib.resources.as_file() can point to a temp file that's cleaned
    up the moment its context manager exits, which loading bytes avoids
    entirely — Qt's QPixmap.loadFromData() works identically either way.
    
    Parameters
    ----------
    name : str
        Filename of the resource to load, relative to the
        ``ediff.gui.resources`` package (e.g. ``'ediff.png'``). Must
        match a file that actually exists in that package/folder.
    
    Returns
    -------
    bytes
        The raw file contents, exactly as stored on disk (or inside
        the PyInstaller bundle) — no decoding or parsing is applied
        here; that's left to whatever consumes the bytes (e.g.
        QPixmap.loadFromData for image files).
    
    Raises
    ------
    FileNotFoundError
        If `name` doesn't exist in the resources package/folder, in
        any of the three contexts checked above.

    """
    # BRANCH 1: the frozen/PyInstaller case (bundled exe)
    if getattr(sys, 'frozen', False):
        # PyInstaller sets sys._MEIPASS to the path of a temp directory where 
        # it has already extracted every bundled data file onto a real 
        # filesystem => read the file normally from the known extraction folder
        with open(Path(sys._MEIPASS) / name, 'rb') as f:
            return f.read()
    
    # BRANCH 2: the normal, non-frozen case (pip-installed or from source)
    # correct way to read a data file that lives inside a package, as opposed 
    # to reading an arbitrary file from disk
    try:
        import importlib.resources as _res
        # this works correctly even if the package itself is sitting inside
        # a compressed wheel/zip
        return (_res.files('ediff.gui.resources') / name).read_bytes()
    
    # BRANCH 3: the except fallback
    # Path(__file__) is this very module's path, so .parent / 'resources' finds 
    # the sibling resources/ folder next to app.py itself. This only works 
    # if the module is sitting unpacked on a normal filesystem (not inside a zip)
    except Exception:
        with open(Path(__file__).parent / 'resources' / name, 'rb') as f:
            return f.read()


def _pixmap(name: str) -> QPixmap:
    """
    Load a bundled resource image as a QPixmap. 
    A convenience wrapper around `_resource_bytes`: fetches the named resource's 
    raw bytes (transparently handling whichever of the
               pip-installed / source-checkout / PyInstaller-frozen 
   contexts this module is running in — see `_resource_bytes` for the full
    explanation) and decodes them into a QPixmap via Qt's own ``loadFromData``, 
    entirely in memory. No file needs to exist on disk at the point this is 
    called; only the bytes matter.

    Parameters
    ----------
    name : str
        Filename of the image resource to load (e.g. ``'ediff.png'``,
        ``'tacr_logo.png'``) — same as the `name` argument to
        `_resource_bytes`.

    Returns
    -------
    QPixmap
        The decoded image, ready to assign directly to a QLabel,
        window icon, or any other Qt widget expecting a QPixmap. 
    """
    px = QPixmap()
    px.loadFromData(_resource_bytes(name))
    return px


# ── Config ────────────────────────────────────────────────────────────────────
# pointer to a hidden JSON file in the user's home directory (a fixed, predictable
# location that survives between separate runs of the app)
# a small quality-of-life memory that makes file dialogs reopen where you left
# off, rather than always starting from some default folder.
CONFIG_PATH = Path.home() / ".ediff_config.json"

def load_config():
    """ It tries to read and JSON-parse the file, and if anything at all goes 
    wrong — the file doesn't exist yet (first-ever run of the app), it's 
    corrupted, a permissions problem — it just returns an empty dict {} instead
    of raising. That means every other part of the app can call load_config() 
    without ever needing to check "does this file exist first" — it's 
    guaranteed to hand back something usable.
    """
    try: return json.loads(CONFIG_PATH.read_text())
    except: return {}

def save_config(cfg):
    """ It tries to write the dict out as JSON, and if that fails for any 
    reason (disk full, no write permission, whatever), silently give up.
    """
    try: CONFIG_PATH.write_text(json.dumps(cfg, indent=2))
    except: pass

# ── Plot style ────────────────────────────────────────────────────────────────
# Controls of the default appearance of every figure created afterward, 
# anywhere in the process, unless a specific plot explicitly overrides 
# something. Since this call sits at module level (not inside a function), it 
# runs exactly once, the moment this module is imported — before any of the 
# app's actual plots (XRD curves, 2D diffractograms, 1D profiles, 
# background/calibration comparisons) get generated. That's what lets every one
# of those very different plotting functions in pipeline.py share one 
# consistent, polished visual identity without repeating the same styling code 
# in each of them individually.

plt.rcParams.update({
    "figure.facecolor": "#ffffff", "axes.facecolor": "#ffffff",
    "axes.edgecolor": "#e5e5ea", "axes.linewidth": 0.8,
    "axes.grid": True, "grid.color": "#f0f0f5", "grid.linewidth": 0.6,
    "xtick.color": "#86868b", "ytick.color": "#86868b",
    "xtick.labelsize": 9, "ytick.labelsize": 9,
    "axes.labelcolor": "#1d1d1f", "axes.labelsize": 10,
    "axes.titlesize": 11, "axes.titlecolor": "#1d1d1f",
    "axes.titleweight": "bold", "lines.linewidth": 1.5,
    "font.family": "sans-serif", "font.size": 9,
    "legend.fontsize": 9, "legend.framealpha": 0.9,
    "savefig.dpi": 150, "savefig.bbox": "tight",
})

# ── Stylesheet ────────────────────────────────────────────────────────────────
# The app's global visual theme: a Qt Style Sheet applied once, app-wide, that
# gives every widget its consistent look.

# (a) Workaround for a bug (dropdown popups were rendering gray text on a black
#     background, unreadable) -- apply this exact same set of rules directly as
#     each individual combobox's own .setStyleSheet(COMBO_FIX) call. COMBO_FIX 
#     exists as a small, reusable string so that fix doesn't need to be retyped 
#     at every one of the app's ~11 comboboxes.
COMBO_FIX = (
    "QComboBox { background:#ffffff; color:#1d1d1f; }"
    "QComboBox QAbstractItemView { background:#ffffff; color:#1d1d1f; }"
    "QComboBox QAbstractItemView::item { color:#1d1d1f; background:#ffffff; }"
    "QComboBox QAbstractItemView::item:selected { background:#e8f0fe; color:#0071e3; }"
)

# (b) The app's global QSS theme
#     This long string is assigned once to the whole QApplication and it 
#     cascades down to every widget in the app unless something more specific
#     overrides it (e.g., COMBO_FIX).
#  >> QMainWindow, QWidget :: base window styling; sets the default light-gray
#     background, near-black text, and font for every widget in the app
#  >> QScrollArea, QScrollBars :: thin (6px), transparent track, rounded 
#     gray handle, no visible up/down arrow buttons (height: 0 on add-line/
#     sub-line) — the flat, minimal scrollbar style 
#  >> QLineEdit, QSpinBox, QDoubleSpinBox, QComboBox :: input fields; shared 
#     rounded-rectangle styling, white background, subtle border
#  >> QComboBox QAbstractItemView block :: the global version of the exact 
#     same fix COMBO_FIX provides individually. Both exist for the reason
#     explained above: belt-and-suspenders, since the global one alone 
#     wasn't reliable.
#  >> QSpinBox :: spinbox up/down buttons
#  >> QPushButton :: buttons; the default button is solid blue with white 
#     text (the primary action style buttons in the left panel), the 
#     secondary buttons (white background, blue border and text) are used 
#     for less-emphasized actions, and stop buttons (red, shown while a step
#     is running)
#  >> QTabWidget/QTabBar :: flat, underline-style tabs; transparent background, 
#     gray inactive text, blue text with a blue underline when selected, 
#     lighter gray when disabled, and a subtle darkening on hover for tabs that
#     are neither selected nor disabled.
#  >> QStatusBar :: status bar; white background, thin top border separating it
#     from the content above, small muted-gray text — deliberately unobtrusive 
#     since it's just carrying transient status messages.
#  >> QToolTip/QMessageBox :: these are explicit overrides that exist because 
#     relying purely on inherited/palette-based styling wasn't reliable enough 
#     in practice

STYLE = """
QMainWindow, QWidget {
    background-color: #f5f5f7;
    color: #1d1d1f;
    font-family: "Segoe UI", Arial, sans-serif;
    font-size: 13px;
}
QScrollArea { border: none; background: transparent; }
QScrollBar:vertical { width: 6px; background: transparent; }
QScrollBar::handle:vertical { background: #c7c7cc; border-radius: 3px; min-height: 30px; }
QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical { height: 0; }

QScrollArea#sidebar {
    background: #ffffff;
    border-right: 1px solid #e5e5ea;
    border-top: none; border-bottom: none; border-left: none;
}

QLineEdit, QSpinBox, QDoubleSpinBox, QComboBox {
    background: #ffffff; border: 1px solid #d1d1d6;
    border-radius: 8px; padding: 5px 10px;
    color: #1d1d1f; font-size: 12px; min-height: 28px;
}
QLineEdit:focus, QSpinBox:focus, QDoubleSpinBox:focus, QComboBox:focus {
    border-color: #0071e3;
}
QLineEdit:disabled, QSpinBox:disabled, QDoubleSpinBox:disabled, QComboBox:disabled {
    background: #f5f5f7; color: #c7c7cc; border-color: #e5e5ea;
}
QComboBox::drop-down { border: none; width: 24px; }
QComboBox::down-arrow { image: none; width: 0; }
QComboBox QAbstractItemView {
    background: #ffffff;
    color: #1d1d1f;
    alternate-background-color: #f9f9fb;
    selection-background-color: #e8f0fe;
    selection-color: #0071e3;
    border: 1px solid #d1d1d6;
    border-radius: 8px;
    padding: 4px;
    outline: none;
}
QComboBox QAbstractItemView::item {
    padding: 5px 10px;
    min-height: 24px;
    color: #1d1d1f;
    background: #ffffff;
}
QComboBox QAbstractItemView::item:selected {
    background: #e8f0fe;
    color: #0071e3;
}
QComboBox QAbstractItemView::item:hover {
    background: #f0f6ff;
    color: #0071e3;
}
QSpinBox::up-button, QSpinBox::down-button,
QDoubleSpinBox::up-button, QDoubleSpinBox::down-button {
    width: 16px; border: none; background: transparent;
}
QSpinBox, QDoubleSpinBox {
    padding-right: 20px;
}
QPushButton {
    background-color: #0071e3; color: #ffffff; border: none;
    border-radius: 8px; padding: 8px 16px;
    font-size: 13px; font-weight: 600; min-height: 34px;
}
QPushButton:hover   { background-color: #0077ed; }
QPushButton:pressed { background-color: #005bbf; }
QPushButton:disabled { background-color: #e5e5ea; color: #aeaeb2; }
QPushButton#secondary {
    background-color: #ffffff; color: #0071e3;
    border: 1.5px solid #0071e3;
}
QPushButton#secondary:hover    { background-color: #e8f0fe; }
QPushButton#secondary:disabled { color: #aeaeb2; border-color: #e5e5ea; }
QPushButton#stop {
    background-color: #ff3b30; color: #ffffff; border: none;
}
QPushButton#stop:hover { background-color: #ff453a; }

QTabWidget::pane  { border: none; background: #f5f5f7; }
QTabBar           { background: #ffffff; border-bottom: 1px solid #e5e5ea; }
QTabBar::tab {
    background: transparent; color: #86868b;
    padding: 9px 20px; font-size: 12px;
    border-bottom: 2px solid transparent; margin-bottom: -1px;
}
QTabBar::tab:selected  { color: #0071e3; border-bottom: 2px solid #0071e3; font-weight: 600; }
QTabBar::tab:disabled  { color: #c7c7cc; }
QTabBar::tab:hover:!selected:!disabled { color: #1d1d1f; }
QStatusBar {
    background: #ffffff; border-top: 1px solid #e5e5ea;
    color: #86868b; font-size: 11px; padding: 2px 12px;
}
QToolTip {
    background-color: #ffffff;
    color: #1d1d1f;
    border: 1px solid #d1d1d6;
    border-radius: 6px;
    padding: 5px 9px;
    font-size: 11px;
    opacity: 255;
}
QMessageBox {
    background-color: #ffffff;
}
QMessageBox QLabel {
    background: transparent;
    color: #1d1d1f;
    font-size: 12px;
}
"""

# ── Option maps ───────────────────────────────────────────────────────────────
# This whole section is the shared "dictionary and toolkit" that everything else 
# in the app draws from — pick a detection method by its friendly name, get a 
# tooltip by its key, build a properly-styled field or button in one line, 
# rather than repeating any of that construction logic at each of the dozens 
# of places it's needed.

# (a) Option maps — translating human-readable labels to API values
#       - the dictionary where the keys are the long, descriptive strings shown
#         in the Center Detection method dropdown, and the values are the short
#         strings ediff.center.CenterLocator actually expects for its 
#         detection= parameter.
#       - refinement options - same idea for the refinement dropdown
DETECTION_OPTIONS = {
    "phase  - auto (phase correlation)":    "phase",
    "ccorr  - auto (cross-correlation)":    "ccorr",
    "hough  - auto (Hough transform)":      "hough",
    "intensity - auto (intensity center)":  "intensity",
    "dirbeam - auto (direct beam)":         "dirbeam",
    "curvefit - auto (curve fit)":          "curvefit",
    "manual - interactive (click 3 pts)":   "manual",
    "from file - read saved center":        None,
}

REFINEMENT_OPTIONS = {
    "sum - maximize intensity sum":         "sum",
    "var - minimize intensity variance":    "var",
    "manual - interactive adjustment":      "manual",
    "none - skip refinement":               None,
}

# BKG_METHODS and CAL_METHODS are simpler — just plain lists of the valid 
# method names for their dropdowns (Background subtraction and Calibration 
# respectively), since those don't need the "friendly label → short code" 
# translation the detection/refinement ones do.
BKG_METHODS = [
    "InteractivePlot", "RestoreFromPoints", "SimpleFuncs", "Baselines"
    ]
CAL_METHODS = [
    "MaxPeaks", "MaxPeaksInRange", "Microscope", "MicroscopeParams"
    ]

# Algorithms available inside the two automatic background methods, and the
# dedicated parameter(s) each one exposes in the GUI. Each parameter is
# (name, kind, default, tooltip): `name` is the exact keyword passed to
# bground; `kind` is "int" / "float" / "bool". Defaults are the values
# recommended in bground's own documentation (verified to run correctly),
# NOT the library defaults - e.g. pybaselines' own `lam` default differs
# per algorithm (10 / 1e6 / None), which is why the default is stored
# per algorithm here rather than once per parameter.
#
# TopHat deliberately has no dedicated parameter: it exists only in the
# development version of bground.sfunc and its keyword arguments are not
# documented, so none are invented here. Anything it (or any algorithm)
# needs beyond these fields can be given in the "Extra parameters" line.
BKG_ALGORITHMS = {
    "SimpleFuncs": {
        "RollingBall": [
            ("radius", "int", 70,
             "Radius of the rolling ball, in pixels (default 70).\n"
             "Smaller values follow the background more closely; larger "
             "values give a smoother background. Judge it on the "
             "'With background' plot."),
        ],
        "TopHat": [],
    },
    "Baselines": {
        "peak_filling": [
            ("half_window", "int", 1,
             "Half-size of the moving-average window "
             "(window = 2 x half_window + 1). Recommended: 1."),
        ],
        "snip": [
            ("decreasing", "bool", True,
             "True iterates the window from large to small, which gives a "
             "smoother baseline (recommended). False goes small to large."),
        ],
        "pspline_iasls": [
            ("lam", "float", 5.0,
             "Smoothing parameter - larger values give a smoother "
             "baseline. Recommended: 5."),
        ],
        "irsqr": [
            ("lam", "float", 1000.0,
             "Smoothing parameter - larger values give a smoother "
             "baseline. Recommended: 1000."),
        ],
        "rubberband": [
            ("lam", "float", 2.0,
             "Smoothing used when interpolating between the baseline "
             "points. Set to 0 for plain linear interpolation. "
             "Recommended: 2."),
        ],
    },
}

# (b) Tooltips — the app's tooltip text library, indexed by key
#        - a single dictionary holding every tooltip string used across the 
#          entire sidebar, keyed by a short identifier ("wavelength", "icut",
#          "finetune", etc.).
TOOLTIPS = {
    "wavelength": "X-ray wavelength in Angstrom (e.g. 0.71 for Mo Ka)",
    "ttmin":      "Minimum 2-theta angle for the theoretical diffractogram",
    "ttmax":      "Maximum 2-theta angle for the theoretical diffractogram",
    "sigma":      "Width of the Gaussian peak profile (smaller = sharper peaks)",
    "detection":  "Algorithm used to find the center of the 2D diffraction pattern",
    "refinement": "Method to refine the initially detected center position",
    "icut":       "Intensity cutoff for display - pixels above this value are clipped",
    "cmap":       "Colormap for displaying the 2D diffraction image",
    "cfile":      "File to read/write the detected center coordinates",
    "csq":        "Size of the zoomed detail view around the detected center (pixels)",
    "px_min":     "Start of the pixel range shown in the 1D profile",
    "px_max":     "End of the pixel range shown in the 1D profile",
    "bkg_method": "Background subtraction algorithm",
    "bkg_file":   "File storing the background control points",
    "bkg_algo":   "Algorithm used to calculate the background",
    "bkg_range":  "Pixel range in which the background is calculated. 'auto' "
                  "uses the start / end of the profile. NOTE: the corrected "
                  "intensity is set to zero outside this range.",
    "bkg_extra":  "Optional extra algorithm parameters as name=value pairs, "
                  "e.g.  size=50, max_half_window=40\n"
                  "Values are numbers, True/False, None or quoted text. They "
                  "override the fields above if both set the same name.",
    "cal_method": "Calibration method to convert pixels to q values",
    "cal_ref":    "CIF phase used as the reference for calibration",
    "finetune":   "Small multiplicative correction applied to the q axis after calibration",
    "qmin":       "Minimum q value shown in the final comparison plot [1/A]",
    "qmax":       "Maximum q value shown in the final comparison plot [1/A]",
    "ellipse_method": (
        "Method for detecting the elliptical distortion.\n"
        "interactive: click points on a ring to fit the ellipse.\n"
        "auto: automatic ring detection and ellipse fitting."
    ),
    "ellipse_icut":  "Intensity cutoff for the ellipse fitting display.",
    "ellipse_npts":  "Number of points to click when using interactive ellipse fitting.",
    "ellipse_file":  "File to save/restore the ellipse correction parameters.",
}

# ── Helpers ───────────────────────────────────────────────────────────────────
def parse_extra_params(text):
    """Parse a "name=value, name=value" line into a dict of keyword arguments.

    Used for the free-form "Extra parameters" field of the background step,
    so any algorithm keyword that has no dedicated field (e.g. an option of
    TopHat, or of a pybaselines algorithm) can still be passed through.

    Values are parsed with ast.literal_eval, so only Python literals are
    accepted: numbers, True/False, None, quoted strings, and tuples/lists/
    dicts of those. Nothing is ever executed, so a stray or malicious entry
    cannot run code. Commas inside brackets or quotes are handled correctly
    (e.g. ``bounds=(1, 5), size=50``).

    Parameters
    ----------
    text : str
        The raw text. Blank / whitespace-only returns an empty dict.

    Returns
    -------
    dict
        Parameter name -> parsed value.

    Raises
    ------
    ValueError
        If the text is not a valid list of ``name=value`` pairs, or a value
        is not a plain literal.
    """
    text = (text or "").strip()
    if not text:
        return {}
    msg = ("Extra parameters must be written as  name=value, name=value  "
           "(for example:  size=50, decreasing=True).")
    try:
        call = ast.parse(f"dict({text})", mode="eval").body
        # Reject positional arguments and **unpacking - only name=value
        if call.args or any(k.arg is None for k in call.keywords):
            raise ValueError(msg)
        return {k.arg: ast.literal_eval(k.value) for k in call.keywords}
    except (SyntaxError, ValueError):
        raise ValueError(msg) from None


def divider():
    """ 
    Create a thin horizontal separator line.

    QFrame with HLine shape is Qt's built-in way to draw a horizontal
    rule. This wraps that setup in one function so any part of the UI
    that wants a subtle visual break between sections can just call
    `divider()`, rather than repeating the frame construction and
    styling at every call site.

    Returns
    -------
    QFrame
        A styled horizontal-line frame — a barely-visible, pale gray
        1px rule (matching the muted styling used elsewhere in the
        app, e.g. the plot gridline color) — ready to add directly to
        a layout.
    """
    f = QFrame()
    f.setFrameShape(QFrame.Shape.HLine)
    f.setStyleSheet(
        "color:#f0f0f5; background:#f0f0f5; border:none; max-height:1px;"
        )
    return f

def _glabel2(name, unit, tip=None):
    """
    Build a compact two-line label: bold name above a smaller unit.

    Stacks two QLabels vertically inside one widget — a bold, dark,
    larger name on top (e.g. "Wavelength") and a smaller, muted-gray
    unit underneath (e.g. "[Å]"). Used for the grid-based, side-by-side
    parameter rows (e.g. "2θ min | 2θ max" sharing one row), where each
    field needs its own compact name+unit label rather than the
    full-width, single-line label `field()` uses below.

    Parameters
    ----------
    name : str
        The parameter's name, shown bold on the top line.
        
    unit : str
        The parameter's unit (e.g. "[Å]", "[deg]"), shown smaller and
        gray on the line below.
        
    tip : str or None, optional, default=None
        Tooltip text. If given, it's applied to both the name and unit
        labels, so hovering over either part of the two-line widget
        shows the same explanation. If None, no tooltip is set.

    Returns
    -------
    QWidget
        A widget containing the stacked name/unit labels, with a
        transparent background, ready to add to a grid layout.
    """
    w = QWidget(); w.setStyleSheet("background:transparent;")
    v = QVBoxLayout(w); v.setContentsMargins(0, 0, 0, 0); v.setSpacing(0)
    n = QLabel(name); n.setStyleSheet(
        "font-size:11px; font-weight:700; color:#1d1d1f;"
        )
    u = QLabel(unit); u.setStyleSheet("font-size:10px; color:#86868b;")
    v.addWidget(n); v.addWidget(u)
    if tip:
        n.setToolTip(tip); u.setToolTip(tip)
    return w

def field(label_text, widget, tip_key=None):
    """
    Wrap a label and an input widget into one labeled sidebar field.

    Given a label string and an input widget, stacks them vertically
    (label on top, input below) with consistent spacing and styling,
    and returns the whole pair as one widget ready to drop into a
    sidebar layout — the standard, single-field building block used
    throughout the app's step panels.

    Sizing differs by widget type: QSpinBox/QDoubleSpinBox widgets are
    capped at 140px width and left-aligned, since a spinbox only ever
    shows a short number — letting it stretch to fill a wide sidebar
    just leaves empty box for no visual reason. Dropdowns and text
    fields are left free to expand instead, since a combobox's
    selected option, or a file path, can genuinely need the extra
    width.

    Parameters
    ----------
    label_text : str
        The label displayed above the input widget.
        
    widget : QWidget
        The input widget itself (QLineEdit, QSpinBox, QDoubleSpinBox,
        QComboBox, etc.) — its size policy and, for spinboxes, width
        and focus policy are set by this function, so it should not
        already be placed in another layout.
        
    tip_key : str or None, optional, default=None
        Key into the module-level `TOOLTIPS` dict. If given and
        present in `TOOLTIPS`, the corresponding text is set as the
        widget's tooltip (not the label's — tooltips belong on the
        thing being interacted with, not duplicated onto the label
        too). If None, or not found in `TOOLTIPS`, no tooltip is set.

    Returns
    -------
    QWidget
        A widget containing the label stacked above the input widget,
        ready to add to a sidebar layout.
    
    """
    w = QWidget()
    v = QVBoxLayout(w)
    v.setContentsMargins(0, 0, 0, 0)
    v.setSpacing(3)
    lbl = QLabel(label_text)
    lbl.setStyleSheet("font-size:11px; color:#1d1d1f; font-weight:600;")
    v.addWidget(lbl)
    
    # QSpinBox/QDoubleSpinBox widgets get capped at 140px and left-aligned 
    # (since a spinbox only ever shows a short number — letting it stretch to 
    # fill a wide sidebar just leaves empty box for no reason), while dropdowns 
    # and text fields are left free to expand (since a combobox's selected 
    # option, or a file path, can genuinely need the extra width). The tip_key 
    # parameter, if given and present in TOOLTIPS, sets that tooltip on the 
    # widget itself.
    
    if isinstance(widget, (QSpinBox, QDoubleSpinBox)):
        widget.setSizePolicy(
            QSizePolicy.Policy.Maximum, QSizePolicy.Policy.Fixed
            )
        widget.setMaximumWidth(140)
        widget.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        v.addWidget(widget, alignment=Qt.AlignmentFlag.AlignLeft)
    else:
        widget.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed
            )
        v.addWidget(widget)
        
    if tip_key and tip_key in TOOLTIPS:
        widget.setToolTip(TOOLTIPS[tip_key])
        
    return w

def make_btn(text, primary=True, tip=None, obj_name=None):
    """ Create a QPushButton with the app's standard visual variants.

    Small factory wrapping QPushButton construction so every button in
    the app picks up its look through one consistent mechanism: Qt
    object names, which is exactly what the STYLE stylesheet's
    `#secondary` and `#stop` selectors key off of. Centralizing this
    here means a button's visual variant is chosen with one keyword
    argument at the call site, rather than every caller needing to
    know and set the right object name directly.
    
    Parameters
    ----------
    text : str
        Label text displayed on the button.
        
    primary : bool, optional, default=True
        If True, the button is left unnamed and picks up the default
        QPushButton rule from STYLE — solid blue, used for the main
        action on a given step (e.g. "Calculate XRD", "Find center").
        If False, sets the object name to "secondary", giving the
        white-background/blue-border look used for less prominent
        actions (e.g. "Browse structure files...").
        
    tip : str or None, optional, default=None
        Tooltip text to show on hover. If None, no tooltip is set.
        
    obj_name : str or None, optional, default=None
        Explicit object name, for variants other than the
        primary/secondary pair — currently only "stop" (the red Stop
        buttons shown while a step is running) uses this. If given,
        this takes precedence over whatever `primary` would have set,
        since it's applied after the `primary` check.

    Returns
    -------
    QPushButton
        The constructed button, styled and tooltipped as requested,
        ready to add to a layout.
    """
    b = QPushButton(text)
    if not primary:
        b.setObjectName("secondary")
    if obj_name:
        b.setObjectName(obj_name)
    if tip:
        b.setToolTip(tip)
    return b

# ── Worker ────────────────────────────────────────────────────────────────────
# Qt apps have exactly one thread that's allowed to touch the GUI — update labels, 
# redraw widgets, respond to clicks. If you run a slow computation directly
# on that thread, the entire application window freezes solid until it's done: 
# no clicking, no repainting, the OS might even flag it as "not responding."

# Worker by itself doesn't create a thread — it's a plain QObject holding 
# a function and its arguments, waiting to be told to run. The actual threading 
# happens in run_in_thread (the function that constructs a Worker, then does 
# worker.moveToThread(thread) and thread.start()). 
# The pattern is: create a QThread, create a Worker, move the worker onto that 
# thread, then hook thread.started up to worker.run — so when the thread starts, 
# it's run() that actually executes there, not just some arbitrary code.

class Worker(QObject):
    """
    Run a function on a background thread, with cooperative cancellation.

    A QObject wrapper around an arbitrary callable, designed to be
    moved to a QThread (see `run_in_thread`) so a long-running
    computation doesn't block the GUI's main/event thread. Results,
    errors, and cancellation are reported back via Qt signals rather
    than a return value, since the caller lives on a different thread
    and can't safely receive a plain return value directly.

    Cancellation here is cooperative, not forceful, and that
    distinction matters: calling `cancel()` does NOT interrupt
    `fn(*args, **kwargs)` if it's already running — Python threads
    cannot be safely killed mid-execution, so a computation that has
    already started will run to completion regardless of when
    `cancel()` is called. What `cancel()` actually controls is only
    which signal gets emitted once the call resolves: `cancelled()`
    instead of `finished()`/`error()`. A caller that needs to reflect
    this honestly in its UI — e.g. not re-enabling a "run" button as
    if the step were immediately ready to retry — must account for the
    fact that the underlying work may still be quietly running in the
    background even after `cancel()` returns.

    Parameters
    ----------
    fn : callable
        The function to run on the background thread.
        
    *args, **kwargs
        Positional and keyword arguments passed through to `fn` when
        `run()` executes.

    """

    finished  = Signal(object)
    error     = Signal(str)
    cancelled = Signal()

    def __init__(self, fn, *args, **kwargs):
        """Store the callable and its arguments; does not start it.

        Parameters
        ----------
        fn : callable
            The function to run when `run()` is called.
            
        *args, **kwargs
            Arguments to pass to `fn`.
        """
        super().__init__()
        self._fn, self._args, self._kwargs = fn, args, kwargs
        self._cancel = False

    def cancel(self):
        """Request cancellation.

        Sets an internal flag checked by `run()` before starting and
        again after `fn` returns or raises. Does not interrupt `fn`
        if it is already executing — see the class docstring for why
        that's not possible here, and what this means in practice.
        """
        self._cancel = True

    def run(self):
        """Execute `fn` and emit the appropriate result signal.

        Intended to be connected to a QThread's `started` signal (see
        `run_in_thread`), not called directly from the main thread.

        Three possible outcomes (the cancellation flag gets checked at three 
        separate points, not once):

        (1) Before starting at all — if cancel() was already called before 
            the thread even got around to running this (e.g. the user clicked 
            Stop in the split second between clicking "Run" and the thread 
            actually starting), fn is never called. 
            Cheapest possible cancellation.
        
        (2) fn raises an exception — normally this emits error, but if 
            cancellation was requested, it emits cancelled instead — the caller 
            doesn't need to see an "error" for something the user themselves 
            stopped.
        
        (3) fn succeeds — normally emits finished with the result, but again, 
            if cancel() was called at any point, the result gets thrown away 
            and cancelled is emitted instead.
        """
        if self._cancel:
            self.cancelled.emit(); return
        try:
            result = self._fn(*self._args, **self._kwargs)
        except Exception as e:
            if self._cancel:
                self.cancelled.emit()
            else:
                self.error.emit(str(e))
            return
        # Only emit finished if not cancelled — suppresses the done callback
        if self._cancel:
            self.cancelled.emit()
        else:
            self.finished.emit(result)
            

class _ResultDispatcher(QObject):
    """Receives a worker's result on whatever thread it finishes on, and
    safely hands it off to the plain Python on_done/on_error callbacks.

    This class is a one-purpose adapter — it exists solely to convert "a plain 
    function Qt can't identify a thread for" into "a bound method of an object 
    Qt can definitely place on the main thread," closing a real gap where 
    an explicit QueuedConnection request could otherwise be silently ignored, 
    and GUI-touching code could end up running on the wrong thread without any 
    obvious error to point at why.
    
    Parameters
    ----------
    on_done : callable or None
        Plain function to call with the result once the worker
        finishes successfully. May be None if the caller doesn't need
        a completion callback.
        
    on_error : callable or None
        Plain function to call with the error message (str) if the
        worker's function raises. May be None if the caller doesn't
        need an error callback.
        
    parent : QObject or None, optional, default=None
        Standard Qt parent, for object lifetime management.
        
    """
    def __init__(self, on_done, on_error, parent=None):
        super().__init__(parent)
        self._on_done  = on_done
        self._on_error = on_error

    def handle_done(self, result):
        """Forward a successful result to the stored on_done callback.

       Intended to be connected to a Worker's `finished` signal via a
       QueuedConnection, so this runs safely on the main thread even
       though the signal was emitted from a worker thread.

       Parameters
       ----------
       result : object
           The value emitted by Worker.finished — whatever the
           wrapped function returned.
       """
        if self._on_done:
            self._on_done(result)

    def handle_error(self, err):
        """Forward an error message to the stored on_error callback.

        Intended to be connected to a Worker's `error` signal via a
        QueuedConnection, so this runs safely on the main thread even
        though the signal was emitted from a worker thread.

        Parameters
        ----------
        err : str
            The error message emitted by Worker.error (the string
            representation of the exception the wrapped function
            raised).
        """
        if self._on_error:
            self._on_error(err)


def run_in_thread(parent, fn, *args, on_done=None, on_error=None, busy_key=None, **kwargs):
    """ Function that every "run this without freezing the GUI" call in the app 
    goes through — it assembles the thread, the worker, and the dispatcher into 
    one working pipeline, and its two most consequential design choices 
    (routing through _ResultDispatcher, and tracking busy state by actual 
    thread completion rather than by cancellation request) are both direct 
    responses to real bugs this project already hit and fixed.
    
    Parameters
    ----------
    parent : QObject
        Owning object for the new QThread and the _ResultDispatcher,
        and the object whose `_busy_keys`/`_threads` attributes are
        used for busy-tracking and reference-keeping (created on
        `parent` on first use if they don't already exist). Typically
        the main window instance.
        
    fn : callable
        The function to execute on the background thread.
        
    *args
        Positional arguments passed through to `fn`.
        
    on_done : callable or None, optional, default=None
        Called with `fn`'s return value on the main thread if `fn`
        completes successfully and no cancellation was requested.
        
    on_error : callable or None, optional, default=None
        Called with the exception's string representation on the main
        thread if `fn` raises and no cancellation was requested.
        
    busy_key : hashable or None, optional, default=None
        If given, added to `parent._busy_keys` for the duration of the
        computation (created as an empty set on `parent` if it doesn't
        exist yet), and removed only once the thread has genuinely
        finished — see point 2 above. If None, no busy-tracking is
        done for this call.
        
    **kwargs
        Keyword arguments passed through to `fn`.

    Returns
    -------
    Worker
        The Worker instance driving the computation. Callers typically
        keep a reference to this (e.g. `self._active_worker_2 = w`) so
        `.cancel()` can be called on it later, or so the corresponding
        QThread can be located via `parent._threads` when needed (e.g.
        to wait for its `finished` signal before re-enabling a button
        after an abort).
        
    """
    # Create a genuine OS-level thread
    thread = QThread(parent)
    
    # Wrap the target function fn and its arguments in a Worker
    worker = Worker(fn, *args, **kwargs)
    
    # Move the worker onto the new thread
    worker.moveToThread(thread)
    
    # Run the thread's event loop
    thread.started.connect(worker.run)
    
    # Route through a dispatcher QObject (see class above) rather than
    # connecting worker signals directly to on_done/on_error.
    dispatcher = _ResultDispatcher(on_done, on_error, parent=parent)
    if on_done:  worker.finished.connect(
            dispatcher.handle_done, 
            Qt.ConnectionType.QueuedConnection
            )
    if on_error: worker.error.connect(
            dispatcher.handle_error, 
            Qt.ConnectionType.QueuedConnection
            )
    
    # Cleanup wiring: whichever of these three outcomes actually happens,
    # the thread's event loop is told to stop, and once the thread has genuinly
    # finished, the worker object itself is scheduled for deletion
    worker.finished.connect(thread.quit)
    worker.error.connect(thread.quit)
    worker.cancelled.connect(thread.quit)
    thread.finished.connect(worker.deleteLater)
    
    # The key gets added the moment the thread launches, and is only ever 
    # removed when thread.finished fires => the underlying computation has
    # genuinely completed.
    if busy_key is not None:
        # a pure-Python computation can't be interrupted mid-execution
        # prevents launching a second computation for the same step while 
        # the first is still silently running in the background, which could 
        # otherwise corrupt shared pipeline state.
        if not hasattr(parent, "_busy_keys"): parent._busy_keys = set()
        parent._busy_keys.add(busy_key)
        thread.finished.connect(lambda: parent._busy_keys.discard(busy_key))
    
    # Keeping references alive 
    # thread and worker are local variables inside run_in_threa — without 
    # storing them somewhere that outlives this function call, Python's garbage 
    # collector could reclaim them the moment run_in_thread returns (potentially 
    # while the thread is still running) causing the whole thing  to crash or 
    # silently stop. Appending (thread, worker) to a list attached to parent 
    # (the EdiffWindow instance, which lives for the whole session) keeps 
    # a strong reference around for as long as needed.
    if not hasattr(parent, "_threads"): parent._threads = []
    parent._threads.append((thread, worker))
    
    # Start the thread
    thread.start()
    
    return worker

# ── Splash ───────────────────────────────────────────────────────────────────
# This builds the app's startup splash screen — the branded loading image shown 
# briefly while the main window is being constructed. 
# Rather than loading a pre-made splash image file, it paints the whole thing 
# programmatically, pixel by pixel, using QPainter onto a blank QPixmap canvas.

def make_splash():
    """Build the app's startup splash screen.

    Paints a branded loading screen entirely with QPainter, rather
    than loading a pre-made splash image file: a soft diagonal
    gradient background, the bundled 'ediff.png' logo scaled to fit
    and vertically centered on the left, and the app name/tagline/
    "Loading..." text laid out to its right.
    
    Degrades gracefully if the logo resource fails to load (e.g. an
    unsupported/corrupt image — see `_pixmap`'s docstring for why that
    doesn't raise) — the text is simply positioned closer to the left
    edge instead of failing outright.
    
    Returns
    -------
    QSplashScreen
        A ready-to-show splash screen, already flagged to stay on top
        of other windows so it's visible during startup rather than
        potentially hidden behind them. Caller is responsible for
        calling `.show()` and eventually `.finish(main_window)` or
        `.close()` once startup completes.
    """
    import os
    
    # a blank white 520×300 canvas
    W, H = 520, 300
    px = QPixmap(W, H)
    px.fill(QColor("#ffffff"))
    p = QPainter(px)

    # Subtle gradient background
    grad = QLinearGradient(0, 0, W, H)
    grad.setColorAt(0, QColor("#f0f6ff"))
    grad.setColorAt(1, QColor("#ffffff"))
    p.fillRect(0, 0, W, H, grad)

    # Logo image — left side
    logo_px = _pixmap("ediff.png")
    
    # If the logo somehow fails to load (logo_px.isNull()), it degrades 
    # gracefully — no crash, just less left margin for the text.
    if not logo_px.isNull():
        logo_scaled = logo_px.scaled(
            200, 200,
            Qt.AspectRatioMode.KeepAspectRatio,
            Qt.TransformationMode.SmoothTransformation)
        p.drawPixmap(24, (H - logo_scaled.height()) // 2, logo_scaled)
        text_x = 24 + logo_scaled.width() + 20
    else:
        text_x = 40

    # Text — right of logo
    f = QFont("Segoe UI", 32, QFont.Weight.Bold)
    p.setFont(f); p.setPen(QColor("#1d2d44"))
    p.drawText(text_x, 110, "EDIFF")

    f2 = QFont("Segoe UI", 11); p.setFont(f2)
    p.setPen(QColor("#86868b"))
    p.drawText(text_x, 136, "Crystallographic Analysis")
    p.drawText(text_x, 154, "from 4D STEM")

    p.fillRect(text_x, 168, 64, 3, QColor("#0071e3"))

    f3 = QFont("Segoe UI", 10); p.setFont(f3)
    p.setPen(QColor("#aeaeb2"))
    p.drawText(text_x, 256, "Loading...")
    p.end()
    
    # Wraps the finished painted image in a QSplashScreen and sets it to always 
    # stay on top of other windows, so it's actually visible during startup 
    # rather than potentially hidden behind something.
    splash = QSplashScreen(px)
    splash.setWindowFlag(Qt.WindowType.WindowStaysOnTopHint)
    return splash

# ── Spinner / empty state ─────────────────────────────────────────────────────
# These are the two "nothing to show yet" placeholder widgets used inside every 
# PlotCard's content area, covering the two moments when there's no real plot 
# to display:
#    >> SpinnerWidget — shown while a computation is running: an indeterminate 
#       progress bar plus a message label above it, both centered vertically 
#       and horizontally in whatever space they're given.
#    >> EmptyState — shown before anything has been run at all: a simple icon, 
#       a bold title ("No results yet"), and a muted subtitle nudging the user 
#       toward the sidebar steps — the very first thing you'd see in a tab 
#       before its corresponding step has ever executed.    

class SpinnerWidget(QWidget):
    """Centered "computing" placeholder shown while a step is running.

    Displays a message label above an indeterminate progress bar (a
    QProgressBar with range (0, 0), which puts it into marquee/
    indeterminate mode rather than showing a percentage, since there's
    no way to know in advance how long a given computation will take).
    Both are vertically and horizontally centered within whatever
    space this widget is given — intended to be shown inside a
    PlotCard's content area while its step's background computation is
    in progress.
    """
    def __init__(self, parent=None):
        """
        Build the spinner layout with a default "Computing..." message.

        Parameters
        ----------
        parent : QWidget or None, optional, default=None
            Standard Qt parent widget.
        """
        # initialize the base QWidget
        super().__init__(parent)
        
        # vertical box layout; installing it as this widget's layout
        lay = QVBoxLayout(self)
        
        # center its group of contents within whatever extra space is available
        lay.setAlignment(Qt.AlignmentFlag.AlignCenter)
        
        # message label
        self._lbl = QLabel("Computing...")
        self._lbl.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self._lbl.setStyleSheet("color:#86868b; font-size:13px;")
        
        # progress bar
        bar = QProgressBar()
        
        # put QProhressBar into intermediatee mode (it shows a continuously 
        # animated marquee/pulsing bar)
        bar.setRange(0, 0); bar.setFixedWidth(220); bar.setFixedHeight(4)
        
        # hide the percentage label QProgressBar would otherwise draw inside 
        # itself, which wouldn't mean anything for an indeterminate bar anyway.
        bar.setTextVisible(False)
        
        # QSS styling
        bar.setStyleSheet(
            "QProgressBar{background:#e5e5ea;border-radius:2px;border:none;}"
            "QProgressBar::chunk{background:#0071e3;border-radius:2px;}")
        
        # insert an invisible, flexible spacer that expands to consume any
        # available extra room
        lay.addStretch()
        lay.addWidget(self._lbl, alignment=Qt.AlignmentFlag.AlignCenter)
        lay.addSpacing(12)
        lay.addWidget(bar, alignment=Qt.AlignmentFlag.AlignCenter)
        lay.addStretch()

    def set_message(self, msg): self._lbl.setText(msg)

class EmptyState(QWidget):
    """Centered placeholder shown before a step has produced any result.

    A simple, static microscope icon, a bold "No results yet" title,
    and a muted subtitle prompting the user to use the sidebar steps —
    all vertically and horizontally centered. Intended as the initial
    content of a PlotCard's tab before its corresponding step has ever
    been run.
    """
    def __init__(self, parent=None):
        """
        Build the empty-state layout.

        Parameters
        ----------
        parent : QWidget or None, optional, default=None
            Standard Qt parent widget.
        """
        # initialize the base QWidget
        super().__init__(parent)
        
        # create a vertical layout and install it on self in one step
        lay = QVBoxLayout(self)
        
        # center itself within available space, rather than the layout 
        # stretching its contents to fill the widget edge-to-edge
        lay.setAlignment(Qt.AlignmentFlag.AlignCenter)
        
        # The emoji itself is just placed as ordinary text inside a QLabel,
        # Qt renders emoji characters using the system font, so no separate 
        # image asset is needed here
        icon = QLabel("🔬"); icon.setAlignment(Qt.AlignmentFlag.AlignCenter)
        icon.setStyleSheet("font-size:48px;")
        
        # the heading
        title = QLabel("No results yet")
        title.setAlignment(Qt.AlignmentFlag.AlignCenter)
        title.setStyleSheet("font-size:15px; font-weight:600; color:#1d1d1f;")
        
        # the subtitle
        sub = QLabel("Follow the steps on the left to analyse your data.")
        sub.setAlignment(Qt.AlignmentFlag.AlignCenter)
        sub.setStyleSheet("font-size:12px; color:#86868b;")
        
        # stretch-to-stretch vertical centering pattern; insert an invisible, 
        # flexible spacer that expands to consume any available extra room
        lay.addStretch()
        lay.addWidget(icon)
        lay.addSpacing(8)
        lay.addWidget(title)
        lay.addSpacing(6)
        lay.addWidget(sub)
        lay.addStretch()

# ── Plot card ─────────────────────────────────────────────────────────────────
# RunPanel == one tab's worth of content inside a PlotCard (a single attempt or
# run of a step). Every card you have seen with multiple tabs is stacking 
# multiple RunPanel instances, one per attempt, letting you to compare them
# side by side and pick which one to accept

class RunPanel(QWidget):
    """
    Single run result: 
        matplotlib canvas + toolbar + 'Close others' + 'Use this run'.
    """
    # two signals, both emitting self (a child widget needs to tell its 
    # parent/container something happened, without needing a direct reference 
    # back to that parent)
    use_requested          = Signal(object)   # emits self when Use-this-run 
    close_others_requested = Signal(object)   # emits self when Close-others 

    def __init__(self, label, parent=None, show_use_button=True):
        # building the tab's layout
        super().__init__(parent)
        
        # A vertical stack with zero margins/spacing
        lay = QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0); lay.setSpacing(0)

        # Top bar: toolbar + "Close others" + "Use this run"
        # a horizontal strip across the top of the panel
        top = QWidget()
        top.setStyleSheet("background:#ffffff;")
        top_lay = QHBoxLayout(top); top_lay.setContentsMargins(8, 4, 8, 4); top_lay.setSpacing(8)
        
        # embedded matplotlib figure 
        self.fig    = Figure(facecolor="#ffffff", tight_layout=True)
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        self.toolbar = NavToolbar(self.canvas, self)
        self._style_toolbar()
        top_lay.addWidget(self.toolbar, stretch=1)

        # "Close others" — same size/shape as "Use this run", neutral color
        # so it doesn't compete with it. Hidden until a sibling tab exists 
        # (PlotCard turns it on for every panel once there's something 
        # to close).
        self._close_others_btn = QPushButton("Close others")
        self._close_others_btn.setToolTip(
            "Close every other run tab, keeping only this one")
        self._close_others_btn.setStyleSheet(
            "QPushButton {"
            "  background:#f5f5f7; color:#6e6e73; border:1px solid #d1d1d6;"
            "  border-radius:8px; padding:5px 14px;"
            "  font-size:12px; font-weight:600;"
            "}"
            "QPushButton:hover { background:#e8e8ed; }"
            "QPushButton:pressed { background:#d8d8dd; }"
        )
        self._close_others_btn.clicked.connect(
            lambda: self.close_others_requested.emit(self)
            )
        
        self._close_others_btn.setVisible(False)
        top_lay.addWidget(self._close_others_btn)
        
        # the Accept button
        self._use_btn = QPushButton("✓  Use this run")
        self._use_btn.setToolTip("Apply this result to subsequent steps")
        self._use_btn.setStyleSheet(
            "QPushButton {"
            "  background:#0071e3; color:#ffffff; border:none;"
            "  border-radius:8px; padding:5px 14px;"
            "  font-size:12px; font-weight:600;"
            "}"
            "QPushButton:hover { background:#0077ed; }"
            "QPushButton:pressed { background:#005bb5; }"
            "QPushButton[active='true'] {"
            "  background:#34c759; color:#ffffff;"
            "}"
        )
        self._use_btn.clicked.connect(lambda: self.use_requested.emit(self))
        # The run card on the right carries the "Use this run" button now; this
        # one stays only so the old code that reads it keeps working.
        self._show_use = show_use_button
        self._use_btn.setVisible(False)
        top_lay.addWidget(self._use_btn)
        
        # Stacks the top bar above the canvas, and initializes the "am I the 
        # currently selected run" flag to false.
        lay.addWidget(top)
        lay.addWidget(self.canvas)
        self._active = False

    def _style_toolbar(self):
        """Make matplotlib's navigation toolbar small, flat and light.

        Smaller icons in a soft grey (instead of heavy black ones) with a
        rounded highlight on hover / when pressed or checked, and without the
        two buttons that mean little here ("Subplots" and "Customize" - the card's
        "Axis limits" control covers the useful part). Separators left next to
        each other by the hidden buttons are hidden too. If tinting the icons
        fails for any reason the original icons are kept.
        """
        tb = self.toolbar
        tb.setStyleSheet(
            "QToolBar { background:transparent; border:none; padding:0px; spacing:2px; }"
            "QToolButton { background:transparent; border:none; border-radius:7px;"
            "  padding:5px; }"
            "QToolButton:hover { background:#eef1f6; }"
            "QToolButton:checked { background:#dfeafb; }"
            "QToolButton:pressed { background:#d3e1f7; }")
        tb.setIconSize(QSize(18, 18))
        for act in tb.actions():
            if act.text() in ("Subplots", "Customize"):
                act.setVisible(False)
        was_separator = True          # also hides a separator that would come first
        for act in tb.actions():
            if not act.isVisible():
                continue
            if act.isSeparator():
                if was_separator:
                    act.setVisible(False)
                was_separator = True
            else:
                was_separator = False
        try:
            for act in tb.actions():
                icon = act.icon()
                if icon.isNull():
                    continue
                pm = icon.pixmap(QSize(36, 36))
                painter = QPainter(pm)
                painter.setCompositionMode(QPainter.CompositionMode.CompositionMode_SourceIn)
                painter.fillRect(pm.rect(), QColor("#5a6270"))
                painter.end()
                act.setIcon(QIcon(pm))
        except Exception:
            pass

    def set_close_others_visible(self, visible: bool):
        """
        PlotCard calls this on every panel whenever the number of tabs changes, 
        to keep "Close others" shown/hidden consistently across all of them.
        """
        # The sidebar's "Close others" link replaces this button (it stays hidden).
        pass

    def set_active(self, active: bool):
        """Mark this run as the currently active/selected one.
        
        When a panel becomes the accepted/selected run, its button's label 
        changes ("Use this run" → "In use") and its color changes (blue → green, 
        via that active property feeding the [active='true'] QSS selector 
        from __init__).
        """
        self._active = active
        if active:
            self._use_btn.setText("✓  In use")
            self._use_btn.setProperty("active", "true")
        else:
            self._use_btn.setText("✓  Use this run")
            self._use_btn.setProperty("active", "false")
            
        # Force stylesheet re-evaluation
        # Qt style sheets are evaluated once and cached; changing a dynamic 
        # property after the fact doesn't automatically trigger Qt to re-check
        # which QSS rules now apply. 
        #  => unpolish() followed by polish() is the standard, documented 
        #     workaround — it forces the widget's style engine to discard its 
        #     cached appearance and fully re-evaluate the stylesheet from 
        #     scratch, picking up the new property value this time. 
        self._use_btn.style().unpolish(self._use_btn)
        self._use_btn.style().polish(self._use_btn)

    def show_buffers(self, bufs_and_titles):
        """ Method for displaying plots. """
        # import libraries
        from PIL import Image
        import numpy as np
        
        # Figures out a grid layout for however many images are being shown 
        # in this one panel:
        #     <3 :: a single row, they sit side by side rather than stacking
        #     >3 :: a square grid (ceil(sqrt(n)) columns), so a larger number 
        #           of images doesn't end up as one absurdly long row.
        n = len(bufs_and_titles)
        cols = n if n <= 3 else math.ceil(math.sqrt(n))
        rows = 1 if n <= 3 else math.ceil(n / cols)
        
        # wipe out whatever was there before
        self.fig.clear(); self.fig.set_facecolor("#ffffff")
        
        # for each label, PNG-buffer pair passed in it decodes the PNG back into
        # a raw pixel array via Pillow, and displays that array as an image 
        # inside a fresh subplot
        for i, (lbl, buf) in enumerate(bufs_and_titles):
            buf.seek(0)
            img = np.array(Image.open(buf))
            ax = self.fig.add_subplot(rows, cols, i + 1)
            ax.imshow(img); ax.axis("off")
            
            # Titles are only added when there's more than one image and a label
            # was actually given
            if lbl and n > 1:
                ax.set_title(lbl, fontsize=10, color="#1d1d1f", pad=8)
                
        # if the grid has more cells than images, the leftover empty cells are
        # explicitly created and then hidden (otherwise matplotlib would leave 
        # behind empty axis boxes with visible borders/ticks in the unused grid 
        # positions.)                
        for j in range(n, rows * cols):
            self.fig.add_subplot(rows, cols, j + 1).set_visible(False)
            
        self.fig.tight_layout(pad=1.5)
        self.canvas.draw()


# ── Modern plot card: run cards + runs sidebar ────────────────────────────────
# (names are imported here, next to the classes that use them)
from qtpy.QtCore import QPoint, QSize
from qtpy.QtGui import QColor, QIcon, QPainter
from qtpy.QtWidgets import (QFrame, QGridLayout, QGraphicsDropShadowEffect,
                            QScrollArea, QToolButton)


class RunCard(QFrame):
    """One finished run, drawn as a card in the runs list beside the plot.

    It replaces a tab. The card shows the run's name and the time it finished,
    a close button, and either a "Use this run" button or - for the run that
    the next steps will use - a green "In use" chip. Clicking the card itself
    makes the plot show that run.

    The card holds no state of its own. ``PlotCard`` rebuilds all cards from its
    tab widget whenever anything changes, so the cards can never disagree with
    the tabs.

    Signals
    -------
    clicked, use_clicked, close_clicked : Signal(object)
        Each carries the ``RunPanel`` this card stands for.
    """
    clicked       = Signal(object)
    use_clicked   = Signal(object)
    close_clicked = Signal(object)

    STYLE = (
        "QFrame#runCard { background:#ffffff; border:1px solid #e5e5ea;"
        "  border-radius:10px; }"
        "QFrame#runCard:hover { border-color:#c7c7cc; }"
        "QFrame#runCard[current='true'] { background:#f4f8ff; border:1px solid #0071e3; }"
        "QFrame#runCard[inuse='true'] { border-left:4px solid #34c759; }"
        "QFrame#runCard QLabel { background:transparent; border:none; }"
    )

    def __init__(self, panel, title, stamp, in_use, current, can_use, parent=None):
        """
        Parameters
        ----------
        panel : RunPanel
            The run this card stands for (passed back in every signal).
        title : str
            Run name, e.g. ``"phase / sum  icut=200"``.
        stamp : str
            Time the run finished (may be empty).
        in_use : bool
            This is the run the next steps use.
        current : bool
            This is the run the plot is showing right now.
        can_use : bool
            Show a "Use this run" button (some steps have nothing to choose).
        """
        super().__init__(parent)
        self._panel = panel
        self.setObjectName("runCard")
        self.setStyleSheet(self.STYLE)
        self.setProperty("current", "true" if current else "false")
        self.setProperty("inuse", "true" if in_use else "false")
        self.setCursor(Qt.CursorShape.PointingHandCursor)
        self.setToolTip(f"{title}\nfinished {stamp}" if stamp else title)
        lay = QVBoxLayout(self)
        lay.setContentsMargins(12, 10, 8, 10); lay.setSpacing(6)

        top = QHBoxLayout(); top.setSpacing(4)
        name = QLabel(title); name.setWordWrap(True)
        name.setStyleSheet("font-size:12px; font-weight:600; color:#1d1d1f;")
        close = QToolButton(); close.setText("\u2715")
        close.setToolTip("Close this run")
        close.setCursor(Qt.CursorShape.PointingHandCursor)
        close.setStyleSheet(
            "QToolButton { border:none; background:transparent; color:#aeaeb2;"
            "  font-size:11px; padding:2px 5px; border-radius:5px; }"
            "QToolButton:hover { background:#f0f0f3; color:#ff3b30; }")
        # `_=False` absorbs the bool that clicked() passes
        close.clicked.connect(lambda _=False: self.close_clicked.emit(self._panel))
        top.addWidget(name, 1); top.addWidget(close, 0, Qt.AlignmentFlag.AlignTop)
        lay.addLayout(top)

        meta = QHBoxLayout(); meta.setSpacing(6)
        when = QLabel(stamp); when.setStyleSheet("font-size:11px; color:#86868b;")
        meta.addWidget(when); meta.addStretch()
        if in_use:
            chip = QLabel("\u2713 In use")
            chip.setStyleSheet(
                "font-size:11px; font-weight:600; color:#1e8e3e; background:#e6f6ea;"
                "border-radius:9px; padding:2px 9px;")
            meta.addWidget(chip)
        elif can_use:
            use = QPushButton("Use this run")
            use.setCursor(Qt.CursorShape.PointingHandCursor)
            use.setStyleSheet(
                "QPushButton { background:#0071e3; color:#ffffff; border:none; min-height:24px; max-height:24px;"
                "  border-radius:12px; padding:0px 12px; font-size:11px; font-weight:600; }"
                "QPushButton:hover { background:#0077ed; }"
                "QPushButton:pressed { background:#005bb5; }")
            use.setFixedHeight(24)
            use.clicked.connect(lambda _=False: self.use_clicked.emit(self._panel))
            meta.addWidget(use)
        lay.addLayout(meta)

    def mousePressEvent(self, ev):
        """A click anywhere on the card (not on its buttons) shows this run."""
        if ev.button() == Qt.MouseButton.LeftButton:
            self.clicked.emit(self._panel)
        super().mousePressEvent(ev)


class RunsSidebar(QFrame):
    """The column of run cards on the right of a plot card.

    A header ("Runs", a count, and a "Close others" link) above a scrolling list
    of ``RunCard`` widgets. It only displays: ``PlotCard`` fills it (``clear_cards``
    then ``add_card`` per run) and decides when it is visible.

    Signals
    -------
    close_others_clicked : Signal()
        The "Close others" link was clicked.
    """
    close_others_clicked = Signal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setObjectName("runsSidebar")
        self.setFixedWidth(250)
        self.setStyleSheet(
            "QFrame#runsSidebar { background:#f7f8fa; border-left:1px solid #e5e5ea; }")
        lay = QVBoxLayout(self); lay.setContentsMargins(0, 0, 0, 0); lay.setSpacing(0)

        head = QWidget(); hl = QHBoxLayout(head)
        hl.setContentsMargins(14, 12, 10, 6)
        self._title = QLabel("RUNS")
        self._title.setStyleSheet("font-size:11px; font-weight:700; color:#86868b;")
        self._others = QPushButton("Close others")
        self._others.setCursor(Qt.CursorShape.PointingHandCursor)
        self._others.setToolTip("Close every run except the one being shown")
        self._others.setStyleSheet(
            "QPushButton { border:none; background:transparent; color:#0071e3; min-height:0px;"
            "  font-size:11px; padding:2px 4px; }"
            "QPushButton:hover { text-decoration:underline; }")
        self._others.clicked.connect(lambda _=False: self.close_others_clicked.emit())
        self._others.setVisible(False)
        hl.addWidget(self._title); hl.addStretch(); hl.addWidget(self._others)
        lay.addWidget(head)

        scroll = QScrollArea(); scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.Shape.NoFrame)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        scroll.setStyleSheet(
            "QScrollArea { background:transparent; border:none; }"
            "QScrollArea > QWidget > QWidget { background:transparent; }")
        inner = QWidget(); inner.setObjectName("runsInner")
        self._list = QVBoxLayout(inner)
        self._list.setContentsMargins(10, 2, 10, 12); self._list.setSpacing(8)
        self._list.addStretch()
        scroll.setWidget(inner)
        lay.addWidget(scroll, 1)

    def clear_cards(self):
        """Remove every card (the trailing stretch stays)."""
        while self._list.count() > 1:
            item = self._list.takeAt(0)
            w = item.widget()
            if w is not None:
                w.setParent(None); w.deleteLater()

    def add_card(self, card):
        """Append a card above the trailing stretch."""
        self._list.insertWidget(self._list.count() - 1, card)

    def set_count(self, n):
        """Show how many runs there are in the header."""
        self._title.setText(f"RUNS  \u00b7  {n}")

    def set_close_others_visible(self, visible):
        """Show the "Close others" link (only useful with more than one run)."""
        self._others.setVisible(visible)


class PlotCard(QWidget):
    """Container for one step's plot area: empty/loading/results states.

    Manages the three states a step's plot display can be in — nothing
    run yet (EmptyState), a computation currently in progress
    (SpinnerWidget), or one or more completed runs shown as tabs (a
    QTabWidget of RunPanel instances) — switching between them via an
    internal QStackedWidget. Every step in the app that shows a plot
    (XRD, Ellipse, Center, 1D, Background, Comparison) owns one
    PlotCard, and drives its state transitions from that step's
    `_on_XXX` handler (show_spinner while running, show_buffers once a
    result is ready).

    Also owns the "Axis limits" row — a small X min/max, Y min/max
    input row with Apply/Reset buttons, hidden until a refresh callback
    is registered via `set_refresh_callback`. Applying new limits
    re-renders the currently active plot with those limits by calling
    back into whatever step-specific logic knows how to regenerate
    that particular plot (e.g. re-calling `pipeline.plot_xrd(xlim=...,
    ylim=...)`), rather than relying on the embedded matplotlib
    toolbar's own zoom — which, since each RunPanel displays a
    rasterized PNG rather than a live data axis (see RunPanel's
    `show_buffers`), can only zoom into a picture of the plot, not
    the underlying data at a different resolution.

    Signals
    -------
    run_selected(object)
        Emitted when the user clicks "Use this run" on one of this
        card's panels. Payload is the RunPanel instance that was
        selected — connected to by each step's `_on_XXX_run_selected`
        handler to know which specific attempt to accept.

    Attributes
    ----------
    _panels : list of RunPanel
        Every run panel added to this card, in tab order.
    _active_idx : int
        Index into `_panels` of the currently "in use" run, or -1 if
        none has been selected yet.
    _refresh_callback : callable or None
        The function(xlim, ylim) -> list[(label, buf)] registered via
        `set_refresh_callback`, used by the Axis limits row's Apply
        button to regenerate the current plot(s). None (and the axis
        row hidden) until a step registers one.
    """
    # Emitted when the user clicks "Use this run" on a panel.
    # Payload: the RunPanel that was selected.
    run_selected = Signal(object)

    def __init__(self, subtitle="", parent=None):
        """Build the empty/spinner/tabs stack and the Axis limits row.

        Parameters
        ----------
        subtitle : str, optional, default=""
            Short descriptive text shown in a thin bar above the plot
            area (e.g. "Radially averaged intensity profile"). If
            empty, no subtitle bar is created at all.
            
        parent : QWidget or None, optional, default=None
            Standard Qt parent widget.
        """
        super().__init__(parent)
        
        # Setting up state and the outer layout
        # Initializes the bookkeeping this card needs:
        self._run_count  = 0    # counter
        self._panels     = []   # list of RunPanel in tab order
        self._active_idx = -1   # index of the currently "in use" panel
        self._refresh_callback = None  # set via set_refresh_callback()
        
        # a vertical layout with zero margins/spacing, installed onto self
        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0); outer.setSpacing(0)
        
        # ── Header: the subtitle on the left, the axis-limits control on the right ──
        self._header = QFrame()
        self._header.setObjectName("plotHeader")
        self._header.setStyleSheet(
            "QFrame#plotHeader { background:#ffffff; border:none;"
            "  border-bottom:1px solid #f0f0f3; }")
        self._header.setFixedHeight(46)
        hl = QHBoxLayout(self._header)
        hl.setContentsMargins(18, 0, 14, 0); hl.setSpacing(8)
        self._subtitle_lbl = QLabel(subtitle)
        self._subtitle_lbl.setStyleSheet("font-size:12px; color:#86868b; background:transparent;")
        hl.addWidget(self._subtitle_lbl); hl.addStretch()
        self._view_btn = QPushButton("\u2699  Axis limits")
        self._view_btn.setCursor(Qt.CursorShape.PointingHandCursor)
        self._view_btn.setToolTip(
            "Set the X / Y range of the plot (highlighted when custom limits are set)")
        self._view_btn.setStyleSheet(
            "QPushButton { background:#f5f5f7; color:#3a3a3c; border:none; min-height:30px; max-height:30px;"
            "  border-radius:15px; padding:0px 14px; font-size:12px; font-weight:600; }"
            "QPushButton:hover { background:#e8e8ed; }"
            "QPushButton[custom='true'] { background:#e8f0fe; color:#0071e3; }")
        self._view_btn.setFixedHeight(30)
        self._view_btn.clicked.connect(self._toggle_axis_popover)
        self._view_btn.setVisible(False)       # shown once a refresh callback is set
        hl.addWidget(self._view_btn)
        self._header.setVisible(bool(subtitle))
        outer.addWidget(self._header)

        # ── Axis-limits popover: floats under the button, closes on a click outside ──
        # (it keeps the old name _axis_row; the same fields and buttons as before)
        self._axis_row = QFrame(self, Qt.WindowType.Popup | Qt.WindowType.FramelessWindowHint)
        self._axis_row.setAttribute(Qt.WidgetAttribute.WA_TranslucentBackground, True)
        # so a click on the button that closes the popover does not instantly reopen it
        self._axis_row.setAttribute(Qt.WidgetAttribute.WA_NoMouseReplay, True)
        pop_lay = QVBoxLayout(self._axis_row)
        pop_lay.setContentsMargins(12, 6, 12, 18)      # room for the shadow
        pop = QFrame(); pop.setObjectName("axisCard")
        pop.setStyleSheet(
            "QFrame#axisCard { background:#ffffff; border:1px solid #e0e0e5;"
            "  border-radius:14px; }"
            "QFrame#axisCard QLabel { background:transparent; border:none; }")
        shadow = QGraphicsDropShadowEffect(pop)
        shadow.setBlurRadius(28); shadow.setOffset(0, 6); shadow.setColor(QColor(0, 0, 0, 55))
        pop.setGraphicsEffect(shadow)
        pop_lay.addWidget(pop)
        grid = QGridLayout(pop)
        grid.setContentsMargins(18, 16, 18, 16)
        grid.setHorizontalSpacing(8); grid.setVerticalSpacing(10)

        def _mini_field(placeholder):
            """A small, rounded numeric box for one axis limit."""
            w = QLineEdit()
            w.setPlaceholderText(placeholder)
            w.setFixedWidth(92)
            w.setStyleSheet(
                "QLineEdit { font-size:12px; padding:6px 9px; border:1px solid #d1d1d6;"
                "  border-radius:8px; background:#ffffff; }"
                "QLineEdit:focus { border:1px solid #0071e3; }")
            return w

        axis_title = QLabel("Axis limits")
        axis_title.setStyleSheet("font-size:13px; font-weight:700; color:#1d1d1f;")
        axis_label_style = "font-size:12px; color:#6e6e73; font-weight:600;"
        xl = QLabel("X"); xl.setStyleSheet(axis_label_style)
        self._axis_xmin = _mini_field("min")
        self._axis_xmax = _mini_field("max")
        yl = QLabel("Y"); yl.setStyleSheet(axis_label_style)
        self._axis_ymin = _mini_field("min")
        self._axis_ymax = _mini_field("max")
        hint = QLabel("Fill both min and max to limit an axis;\nan axis left empty auto-scales.")
        hint.setStyleSheet("font-size:11px; color:#86868b;")
        self._axis_apply_btn = QPushButton("Apply")
        self._axis_apply_btn.setToolTip(
            "Re-render the current plot with these axis limits. An axis is "
            "limited only when both its min and max are filled; otherwise it auto-scales.")
        self._axis_apply_btn.setStyleSheet(
            "QPushButton { background:#0071e3; color:#ffffff; border:none; min-height:32px; max-height:32px;"
            "  border-radius:10px; padding:0px 20px; font-size:12px; font-weight:600; }"
            "QPushButton:hover { background:#0077ed; }"
            "QPushButton:pressed { background:#005bb5; }")
        self._axis_reset_btn = QPushButton("Reset")
        self._axis_reset_btn.setToolTip("Clear all four limits (auto-scale) and redraw")
        self._axis_reset_btn.setStyleSheet(
            "QPushButton { background:#f5f5f7; color:#3a3a3c; border:none; min-height:32px; max-height:32px;"
            "  border-radius:10px; padding:0px 18px; font-size:12px; font-weight:600; }"
            "QPushButton:hover { background:#e8e8ed; }")
        self._axis_apply_btn.setFixedHeight(32); self._axis_reset_btn.setFixedHeight(32)
        # the original handlers do the work; the extra ones tidy the popover afterwards
        self._axis_apply_btn.clicked.connect(self._on_apply_axis_limits)
        self._axis_reset_btn.clicked.connect(self._on_reset_axis_limits)
        self._axis_apply_btn.clicked.connect(self._after_axis_action)
        self._axis_reset_btn.clicked.connect(self._after_axis_action)
        for _f in (self._axis_xmin, self._axis_xmax, self._axis_ymin, self._axis_ymax):
            _f.textChanged.connect(lambda _t: self._refresh_view_badge())
        grid.addWidget(axis_title, 0, 0, 1, 3)
        grid.addWidget(xl, 1, 0); grid.addWidget(self._axis_xmin, 1, 1); grid.addWidget(self._axis_xmax, 1, 2)
        grid.addWidget(yl, 2, 0); grid.addWidget(self._axis_ymin, 2, 1); grid.addWidget(self._axis_ymax, 2, 2)
        grid.addWidget(hint, 3, 0, 1, 3)
        btn_row = QHBoxLayout(); btn_row.setSpacing(8)
        btn_row.addStretch(); btn_row.addWidget(self._axis_reset_btn); btn_row.addWidget(self._axis_apply_btn)
        grid.addLayout(btn_row, 4, 0, 1, 3)

        # ── Body: the plot area (empty / spinner / runs) with the runs sidebar on its right ──
        body = QWidget()
        body_lay = QHBoxLayout(body)
        body_lay.setContentsMargins(0, 0, 0, 0); body_lay.setSpacing(0)
        self._stack = QStackedWidget()
        self._empty   = EmptyState()
        self._spinner = SpinnerWidget()
        # The tab widget stays the single source of truth for the runs (pages, order,
        # labels, current run); only its tab bar is hidden, because the cards on the
        # right replace it.
        self._run_tabs = QTabWidget()
        self._run_tabs.setTabsClosable(True)
        self._run_tabs.setDocumentMode(True)
        self._run_tabs.tabBar().hide()
        self._run_tabs.tabCloseRequested.connect(self._close_run)
        self._run_tabs.currentChanged.connect(lambda _i: self._refresh_cards())
        self._stack.addWidget(self._empty)    # 0
        self._stack.addWidget(self._spinner)  # 1
        self._stack.addWidget(self._run_tabs) # 2
        self._stack.setCurrentIndex(0)
        self._runs_bar = RunsSidebar()
        self._runs_bar.setVisible(False)      # shown once there is at least one run
        self._runs_bar.close_others_clicked.connect(lambda: self._close_other_runs())
        body_lay.addWidget(self._stack, 1)
        body_lay.addWidget(self._runs_bar)
        outer.addWidget(body, 1)

        # The compatibility placeholders
        self.fig    = None
        self.canvas = None
        self.toolbar = None

    def _close_run(self, idx):
        """ Remove one run tab and keep all bookkeeping consistent.

        Removes both the internal RunPanel record and the corresponding
        visible tab at `idx`, then corrects `_active_idx` for the shift
        this causes (removing an item from a list moves every later item
        down by one position, so a stored index pointing past the removed
        item would otherwise silently point at the wrong panel), and falls
        back to the empty-state placeholder if no tabs remain.
    
        Parameters
        ----------
        idx : int
            Index of the run tab to close, matching its position in both
            `self._panels` and the underlying QTabWidget.
    
        Returns
        -------
        None
    
        Notes
        -----
        Three cases are handled for `_active_idx` after removal:
    
        - If the closed tab *was* the active one, there's no active
          selection anymore (reset to -1).
        - If the active one sat *after* the closed index, its position is
          decremented by one so it still points at the same panel it did
          before removal.
        - If the active one sat *before* the closed index, it's
          unaffected and left unchanged.
    
        If this was the last remaining tab, the stack falls back to the
        empty-state placeholder (index 0) and all bookkeeping is reset.
        Either way, `_refresh_close_others_visibility` is called afterward
        to keep the "Close others" button's visibility correct on whatever
        panels remain.
        """
        # Remove both the internal record and the actual visible tab at that
        # index. Panel is captured but never used afterward
        panel = self._panels.pop(idx) if idx < len(self._panels) else None
        self._run_tabs.removeTab(idx)
        
        # If we closed the active run, clear the active selection
        # _active_idx == plain integer position into _panels => if you remove 
        # an item from a list, every item after it shifts down by one position,
        # so any stored index pointing past the removed item is now silently 
        # wrong unless corrected. 
        # Two cases:
        #   >> if the closed tab was the active one, there's no active selection 
        #      anymore (-1)
        #   >> if the active one was somewhere after the closed one, its 
        #      position needs to shift down by one to still point at the same 
        #      panel it did before. 
        if self._active_idx == idx:
            self._active_idx = -1
        elif self._active_idx > idx:
            self._active_idx -= 1
            
        # If that was the last remaining tab, fall back to the empty-state 
        # placeholder (index 0), and reset everything cleanly. Either way, 
        # re-check whether "Close others" should still be shown on the remaining
        # panels (see next method).
        if self._run_tabs.count() == 0:
            self._stack.setCurrentIndex(0)
            self._panels.clear()
            self._active_idx = -1
        self._refresh_close_others_visibility()
        self._refresh_cards()

    def _refresh_close_others_visibility(self):
        """ Show or hide "Close others" on every panel based on tab count.

        Shows the "Close others" button on every panel once there is more
        than one run tab open; hides it again when only one remains, since
        there is nothing else left to close.
        """
        show = len(self._panels) > 1
        for p in self._panels:
            p.set_close_others_visible(show)

    def _on_close_others_requested(self, panel):
        """ Handle a RunPanel's close_others_requested signal.

        Translates the RunPanel object emitted by the signal into its
        current index within `self._panels` (since that's what
        `_close_other_runs` actually operates on), then delegates to it.
    
        Parameters
        ----------
        panel : RunPanel
            The panel whose "Close others" button was clicked — i.e. the
            one panel that should be kept.
        """
        idx = self._panels.index(panel) if panel in self._panels else -1
        if idx >= 0:
            self._close_other_runs(idx)

    def _close_other_runs(self, keep_idx=None):
        """ Close every run tab except one.

        Parameters
        ----------
        keep_idx : int or None, optional, default=None
            Index of the tab to keep open; every other tab is closed. If
            None, keeps whichever tab is currently being viewed
            (`self._run_tabs.currentIndex()`). If that resolves to a
            negative index (nothing currently selected), this is a no-op.
        """
        if keep_idx is None:
            keep_idx = self._run_tabs.currentIndex()
            
        if keep_idx < 0:
            return
        
        # Remove in descending order of original index, skipping keep_idx.
        # This keeps every index we're about to act on valid: removing a tab
        # only shifts the positions of tabs AFTER it, and we always work our
        # way down from the highest index, so anything at-or-below the one
        # we're currently closing hasn't moved yet.
        
        for idx in reversed(range(self._run_tabs.count())):
            if idx != keep_idx:
                self._close_run(idx)

    def _refresh_cards(self):
        """Rebuild the run cards on the right from the tab widget.

        The tab widget is the single source of truth (which runs exist, their
        order and labels, which one is shown); the "in use" run is
        ``_active_idx``. Rebuilding everything is cheap (a handful of cards) and
        means the cards can never be out of step with the tabs. Safe to call at
        any time, including while widgets are being torn down.
        """
        try:
            bar = self._runs_bar
            bar.clear_cards()
            shown = self._run_tabs.currentIndex()
            for i, panel in enumerate(self._panels):
                card = RunCard(
                    panel,
                    self._run_tabs.tabText(i).lstrip("\u2713 "),
                    self._run_tabs.tabToolTip(i),
                    in_use=(i == self._active_idx),
                    current=(i == shown),
                    can_use=getattr(panel, "_show_use", True))
                card.clicked.connect(self._on_card_clicked)
                card.use_clicked.connect(self._on_use_requested)
                card.close_clicked.connect(self._on_card_close)
                bar.add_card(card)
            bar.set_count(len(self._panels))
            bar.set_close_others_visible(len(self._panels) > 1)
            bar.setVisible(len(self._panels) > 0)
        except RuntimeError:
            pass      # a widget was already deleted (window closing)

    def _on_card_clicked(self, panel):
        """Show the run whose card was clicked."""
        if panel in self._panels:
            self._run_tabs.setCurrentIndex(self._panels.index(panel))

    def _on_card_close(self, panel):
        """Close the run whose card's x was clicked."""
        if panel in self._panels:
            self._close_run(self._panels.index(panel))

    def _toggle_axis_popover(self):
        """Open the axis-limits popover under its button, or close it."""
        pop = self._axis_row
        if pop.isVisible():
            pop.hide()
            return
        pop.adjustSize()
        btn = self._view_btn
        pop.move(btn.mapToGlobal(QPoint(btn.width() - pop.width() + 12, btn.height() - 2)))
        pop.show()
        self._axis_xmin.setFocus()

    def _after_axis_action(self, _checked=False):
        """After Apply or Reset: close the popover and update the button's highlight."""
        self._axis_row.hide()
        self._refresh_view_badge()

    def _refresh_view_badge(self):
        """Highlight the "Axis limits" button while any limit is set."""
        custom = any(f.text().strip() for f in (
            self._axis_xmin, self._axis_xmax, self._axis_ymin, self._axis_ymax))
        self._view_btn.setProperty("custom", "true" if custom else "false")
        self._view_btn.style().unpolish(self._view_btn)
        self._view_btn.style().polish(self._view_btn)

    def show_spinner(self, msg="Computing..."):
        """ Switch this card to the loading/spinner state.

        Parameters
        ----------
        msg : str, optional, default="Computing..."
            Message displayed above the indeterminate progress bar (e.g.
            "Fitting elliptical distortion — follow the window...",
            "Radially averaging...").
        """
        self._spinner.set_message(msg)
        self._stack.setCurrentIndex(1)

    def reset_idle(self):
        """Clear the spinner without discarding any existing successful runs.

        `show_spinner()` switches the stack to the loading state, but
        nothing automatically switches it back. When a computation
        succeeds, `show_buffers()` naturally switches to the results tab
        as a side effect of adding one, so that case is already covered.
        But if a computation fails, nothing calls `show_buffers()`, and
        the spinner would otherwise sit there indefinitely — this is the
        explicit undo needed for that case, called from a step's error
        handler.
    
        Distinguishes two different "back to idle" scenarios: if any
        successful runs already exist from before this failed attempt,
        shows those (index 2) rather than reverting to a blank
        empty-state as if nothing had ever worked; if this was the very
        first attempt and it failed, there is genuinely nothing to show,
        so falls back to the empty-state placeholder (index 0).
        """
        self._stack.setCurrentIndex(2 if self._panels else 0)

    def show_buffers(self, bufs_and_titles, run_label=None, auto_select=False,
                     show_use_button=True):
        """ Add a new run tab showing the given plot image(s).

        Parameters
        ----------
        bufs_and_titles : list of (str, BytesIO)
            One or more (label, PNG-buffer) pairs to display in this run's
            panel — passed straight through to RunPanel.show_buffers.
            
        run_label : str or None, optional, default=None
            Label shown on the tab itself. If None, falls back to a
            generic "Run N" using this card's internal run counter.
            
        auto_select : bool, optional, default=False
            If True, immediately marks this run as "in use" and emits
            `run_selected` — appropriate when there is only one possible
            result (e.g. a step just completed with no alternatives to
            compare against). If False, the "Use this run" button is shown
            instead and the user picks explicitly.
            
        show_use_button : bool, optional, default=True
            If False, hides the "Use this run" affordance entirely.
            Appropriate for cards where selecting an older run has no real
            effect (no snapshot/restore behavior behind it), so the button
            doesn't invite a click that wouldn't actually do anything.
    
        Returns
        -------
        RunPanel
            The newly created panel, so callers can attach extra data to
            it afterward (e.g. `panel._center_snapshot = ...`) before it
            might be selected.
    
        Notes
        -----
        Also updates `self.fig`/`self.canvas`/`self.toolbar` to point at
        the newly added panel's figure/canvas/toolbar — kept only for
        external compatibility (see the `__init__` comment); note these
        always reflect the most *recently added* run, not necessarily
        whichever tab is currently being viewed if the user has since
        switched tabs.
        """
        from datetime import datetime
        
        # Increments the run counter (used only as a fallback label if the caller 
        # didn't supply a more descriptive run_label), and captures a timestamp 
        # string — used as a tooltip on the tab a few lines down, presumably so 
        # hovering over a tab shows exactly when that attempt was run.
        self._run_count += 1
        ts  = datetime.now().strftime("%H:%M:%S")
        tab_label = run_label if run_label else f"Run {self._run_count}"
        
        # Build a fresh RunPanel, render the given image buffers, wire its two
        # signals to this card's own handlers, closing the loop
        panel = RunPanel(tab_label, show_use_button=show_use_button)
        panel.show_buffers(bufs_and_titles)
        panel.use_requested.connect(self._on_use_requested)
        panel.close_others_requested.connect(self._on_close_others_requested)
        
        # Empty placeholders from __init__ get populated here:
        self.fig    = panel.fig
        self.canvas = panel.canvas
        self.toolbar = panel.toolbar
        
        # Register the panel
        self._panels.append(panel)
        
        # add it as an actual tab (based on the timestamp)
        idx = self._run_tabs.addTab(panel, tab_label)
        
        # switch focus to the new tab
        self._run_tabs.setTabToolTip(idx, ts)
        self._run_tabs.setCurrentIndex(idx)
        self._stack.setCurrentIndex(2)
        
        # refresh visibility across all panels
        self._refresh_close_others_visibility()
        self._refresh_cards()
        
        # immediately marks this run as the accepted one via _select_panel
        if auto_select:
            self._select_panel(idx)
            
        return panel

    def set_refresh_callback(self, callback):
        """ Register the function used to re-render this card's plot(s).

        Parameters
        ----------
        callback : callable or None
            A function(xlim, ylim) -> list[(label, buf)] that re-renders
            the current plot(s) with new axis limits, using whatever
            underlying data/parameters the owning step's handler already
            has in scope. Passing a callback here also reveals the "Axis
            limits" row; pass None to hide it again (e.g. for a card with
            no meaningful numeric axes to set).
        """
        # Stores the callback and shows/hides the row based on whether one was 
        # actually given — exactly matching what was described when this 
        # feature was designed.
        self._refresh_callback = callback
        # The axis limits live behind the header's button; it only makes sense
        # once there is something to redraw.
        self._view_btn.setVisible(callback is not None)
        self._header.setVisible(bool(self._subtitle_lbl.text()) or callback is not None)
        if callback is None:
            self._axis_row.hide()

    def set_axis_fields(self, xlim, ylim):
        """Pre-fill the axis-limit input boxes without triggering a refresh.

        Purely cosmetic — used to show the user which carried-over limits
        (from an earlier step with the same kind of axes) are about to be,
        or already have been, applied to a freshly-generated plot.
    
        Parameters
        ----------
        xlim : tuple of (float, float) or None
            (min, max) to show in the X fields, or None/falsy to leave
            them blank.
            
        ylim : tuple of (float, float) or None
            (min, max) to show in the Y fields, or None/falsy to leave
            them blank.
        """
        # Pre-fills the four boxes to reflect carried-over limits, without
        # triggering an actual refresh — purely cosmetic, so the user can see 
        # what's about to be (or already was) applied.
        self._axis_xmin.setText("" if not xlim else str(xlim[0]))
        self._axis_xmax.setText("" if not xlim else str(xlim[1]))
        self._axis_ymin.setText("" if not ylim else str(ylim[0]))
        self._axis_ymax.setText("" if not ylim else str(ylim[1]))

    def _parse_axis_fields(self):
        """ Read the four axis-limit text boxes into (xlim, ylim) tuples.

        Returns
        -------
        xlim : tuple of (float, float) or None
            (min, max) parsed from the X fields, or None if either field
            is blank or fails to parse as a float — a lone filled-in min
            with a blank max produces None (auto-scale), not a broken
            one-sided range.
            
        ylim : tuple of (float, float) or None
            Same logic as `xlim`, independently, for the Y fields.
        """
        
        def _val(w):
            """
            Parse one field's text as a float, or None if blank/invalid.

            Treats an empty/blank field and text that fails to parse as a
            float identically — both become None, rather than raising or
            silently defaulting to zero.
    
            Parameters
            ----------
            w : QLineEdit
                The field to read.
    
            Returns
            -------
            float or None
                The parsed value, or None if the field was blank or not a
                valid float.
            """
            t = w.text().strip()
            if not t:
                return None
            try:
                return float(t)
            except ValueError:
                return None
        
        # xlim is only constructed as an actual (min, max) tuple if both 
        # the min and max fields parsed successfully — a lone filled-in min with 
        # a blank max produces xlim = None (auto-scale), not a broken one-sided 
        # range. Same logic for ylim, independently.
        xmin, xmax = _val(self._axis_xmin), _val(self._axis_xmax)
        ymin, ymax = _val(self._axis_ymin), _val(self._axis_ymax)
        
        xlim = (xmin, xmax) if xmin is not None and xmax is not None else None
        ylim = (ymin, ymax) if ymin is not None and ymax is not None else None
        
        return xlim, ylim

    def _on_apply_axis_limits(self):
        """ Handle the Apply button: re-render the current plot with new limits.

        Parses the axis-limit fields and runs the registered refresh
        callback on a background thread (via run_in_thread, the same
        mechanism every other step in this app already uses), then
        applies the result once it's ready. Regenerating a plot means a
        full matplotlib render + PNG encode -- for the comparison view,
        one such render per XRD phase -- which previously ran directly
        on the main thread and froze the whole GUI for its duration on
        every single "Apply" click, on every step that has a plot. This
        is one of the most frequently-clicked controls in the app (axis
        tweaking is naturally iterative), so that synchronous call was
        a much bigger source of perceived GUI lag than any one-time
        startup cost.

        Does nothing if no callback is registered. If the callback
        raises (e.g. the underlying recalculation itself fails), shows
        a warning dialog rather than letting the exception propagate
        and crash the app. On success, replaces the content of
        whichever tab was being viewed at the moment "Apply" was
        clicked (`self._run_tabs.currentIndex()`, captured before
        dispatching to the background thread -- not whatever tab
        happens to be current once the computation finishes, since the
        user could switch tabs while it's still running) in place via
        `RunPanel.show_buffers` on that same panel, rather than
        creating a new tab.
        """
        # guard against no callback being registered
        if self._refresh_callback is None:
            return

        # guard against overlapping refreshes: the previous, always-
        # synchronous version couldn't have two of these running at
        # once (each blocked the whole GUI until it finished); moving
        # this to a background thread makes that newly possible if the
        # user clicks "Apply" again before the first render completes.
        if getattr(self, '_axis_apply_busy', False):
            return
        self._axis_apply_busy = True
        if hasattr(self, '_axis_apply_btn'):
            self._axis_apply_btn.setEnabled(False)

        # parse the fields (fast; Qt widget access, main thread only)
        xlim, ylim = self._parse_axis_fields()

        # Capture the callback and target tab now, on the main thread,
        # before dispatching to the background -- the same "capture
        # GUI state now, use it later on the worker thread" pattern
        # used by every other long-running step in this app. This also
        # means the result always applies to the tab that was active
        # when the user clicked Apply, regardless of whatever tab is
        # current by the time the background render finishes.
        refresh_callback = self._refresh_callback
        idx = self._run_tabs.currentIndex()

        def do():
            return refresh_callback(xlim, ylim)

        def done(bufs_and_titles):
            self._axis_apply_busy = False
            if hasattr(self, '_axis_apply_btn'):
                self._axis_apply_btn.setEnabled(True)
            if not bufs_and_titles:
                return
            if 0 <= idx < len(self._panels):
                self._panels[idx].show_buffers(bufs_and_titles)

        def on_err(err):
            self._axis_apply_busy = False
            if hasattr(self, '_axis_apply_btn'):
                self._axis_apply_btn.setEnabled(True)
            QMessageBox.warning(
                self,
                "EDIFF",
                f"Could not update the plot:\n{err}"
                )

        run_in_thread(self, do, on_done=done, on_error=on_err)

    def _on_reset_axis_limits(self):
        """ Clear all four axis-limit fields and re-render with auto-scaled 
        axes.

        Blanks the X/Y min/max fields, then reuses `_on_apply_axis_limits`
        directly — since blank fields parse to None limits, this naturally
        triggers a re-render with fully auto-scaled axes, without needing
        any separate reset-specific rendering logic.
        """
        self._axis_xmin.clear(); self._axis_xmax.clear()
        self._axis_ymin.clear(); self._axis_ymax.clear()
        self._on_apply_axis_limits()

    def select_run(self, panel):
        """Publicly select a specific panel as the accepted run.

        Distinct from clicking the "Use this run" button directly, but
        triggers the exact same selection logic. Used to auto-accept a run
        programmatically after attaching extra data to it (e.g. a
        snapshot).
    
        Parameters
        ----------
        panel : RunPanel
            The panel to select. Must already be attached to this card
            (i.e. present in `self._panels`).
    
        Notes
        -----
        Any data a caller needs available to `run_selected`'s handler
        (e.g. `panel._center_snapshot`) must be attached to `panel` BEFORE
        calling this — selection emits `run_selected` immediately and
        synchronously, and the connected handler reads that data at the
        moment the signal fires, not afterward.
        """
        # a public entry point, distinct from clicking the button
        self._on_use_requested(panel)

    def _on_use_requested(self, panel):
        """Handle a RunPanel's use_requested signal.

        Translates the RunPanel object emitted by the signal into its
        current index within `self._panels`, then delegates the actual
        selection work to `_select_panel`.
    
        Parameters
        ----------
        panel : RunPanel
            The panel whose "Use this run" button was clicked.
        """
        idx = self._panels.index(panel) if panel in self._panels else -1
        if idx >= 0:
            self._select_panel(idx)

    def _select_panel(self, idx):
        """Mark the panel at `idx` as the accepted/active run.

        If a different panel was previously active, deselects it (reverts
        its button's styling via `set_active(False)`) and strips the "✓ "
        checkmark prefix from its tab label. Marks the new panel active,
        adds the "✓ " prefix to its own tab label (only if not already
        present, so re-selecting an already-checked tab doesn't stack up
        multiple checkmarks), switches the visible tab to it, and finally
        emits `run_selected` — the signal every step's
        `_on_XXX_run_selected` handler listens for, to actually restore
        that specific run's underlying data/snapshot.
    
        Parameters
        ----------
        idx : int
            Index, into `self._panels`, of the panel to select.
        """
        # Deselect previous
        if 0 <= self._active_idx < len(self._panels):
            self._panels[self._active_idx].set_active(False)
            self._run_tabs.setTabText(
                self._active_idx,
                self._run_tabs.tabText(self._active_idx).lstrip("✓ "))
        # Select new
        self._active_idx = idx
        panel = self._panels[idx]
        panel.set_active(True)
        tab_text = self._run_tabs.tabText(idx)
        if not tab_text.startswith("✓"):
            self._run_tabs.setTabText(idx, f"✓ {tab_text}")
        self._run_tabs.setCurrentIndex(idx)
        self._refresh_cards()
        self.run_selected.emit(panel)
        
    # Read only status properties; accessed like plain attributes
    @property
    def has_selection(self) -> bool:
        return self._active_idx >= 0

    @property
    def active_panel(self):
        if 0 <= self._active_idx < len(self._panels):
            return self._panels[self._active_idx]
        return None
    
    # full reset, wipes everything
    def clear_runs(self):
        self._run_tabs.clear()
        self._panels.clear()
        self._run_count  = 0
        self._active_idx = -1
        self._stack.setCurrentIndex(0)
        self._refresh_cards()


class SaveSummaryCard(QWidget):
    """
    Container for Step 9's (Save) result area: empty/loading/summary states.

    A much simpler sibling of PlotCard: Step 9 doesn't produce a plot to
    compare across multiple runs, it produces a one-shot listing of what
    was written to disk. Reuses the same EmptyState/SpinnerWidget
    building blocks and the same three-state QStackedWidget pattern
    (nothing run yet / computing / result shown) as PlotCard, but the
    "result" state is a simple output-folder label plus a file list
    rather than a QTabWidget of RunPanel plot instances, and there's no
    axis-limits row (nothing to zoom into for a file listing).
    """

    def __init__(self, subtitle="", parent=None):
        """Build the empty/spinner/summary stack.

        Parameters
        ----------
        subtitle : str, optional, default=""
            Short descriptive text shown in a thin bar above the
            summary area, matching PlotCard's own subtitle bar. If
            empty, no subtitle bar is created at all.

        parent : QWidget or None, optional, default=None
            Standard Qt parent widget.
        """
        super().__init__(parent)

        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0); outer.setSpacing(0)

        if subtitle:
            sb = QWidget()
            sb.setStyleSheet(
                "background:#ffffff; border-bottom:1px solid #f0f0f3;"
            )
            sb.setFixedHeight(46)
            sbl = QHBoxLayout(sb)
            sbl.setContentsMargins(18, 0, 14, 0)
            lbl = QLabel(subtitle)
            lbl.setStyleSheet("font-size:12px; color:#86868b;")
            sbl.addWidget(lbl)
            outer.addWidget(sb)

        self._stack = QStackedWidget()
        self._empty = EmptyState()
        self._spinner = SpinnerWidget()

        # Result state: a simple listing of the output folder + every
        # file written, rather than a QTabWidget of RunPanel/plot
        # instances -- Save has nothing to compare across multiple
        # attempts, just one listing of the most recent save.
        self._summary = QWidget()
        self._summary.setStyleSheet("background:#ffffff;")
        sum_lay = QVBoxLayout(self._summary)
        sum_lay.setContentsMargins(24, 24, 24, 24)
        sum_lay.setSpacing(12)
        sum_lay.setAlignment(Qt.AlignmentFlag.AlignTop)

        header = QLabel("\u2713 Saved successfully")
        header.setStyleSheet(
            "font-size:15px; font-weight:600; color:#1d1d1f;"
        )
        sum_lay.addWidget(header)

        self._folder_lbl = QLabel("")
        self._folder_lbl.setStyleSheet("font-size:12px; color:#86868b;")
        self._folder_lbl.setWordWrap(True)
        sum_lay.addWidget(self._folder_lbl)

        self._file_list = QListWidget()
        self._file_list.setStyleSheet(
            "QListWidget { border:1px solid #e5e5ea; border-radius:6px; "
            "font-size:12px; } "
            "QListWidget::item { padding:4px 8px; }")
        sum_lay.addWidget(self._file_list, stretch=1)

        self._stack.addWidget(self._empty)    # 0
        self._stack.addWidget(self._spinner)  # 1
        self._stack.addWidget(self._summary)  # 2

        self._stack.setCurrentIndex(0)
        outer.addWidget(self._stack)

    def show_spinner(self, msg="Saving..."):
        """
        Switch to the loading/spinner state.

        Parameters
        ----------
        msg : str, optional, default="Saving..."
            Message shown above the indeterminate progress bar.
        """
        self._spinner.set_message(msg)
        self._stack.setCurrentIndex(1)

    def reset_idle(self):
        """
        Clear the spinner without discarding a previous successful save.

        Mirrors PlotCard.reset_idle: if a save already succeeded once
        this session, an error on a later attempt shouldn't wipe out
        that previous summary -- falls back to it (index 2) rather
        than the blank empty-state (index 0) if one exists.
        """
        self._stack.setCurrentIndex(
            2 if self._file_list.count() else 0
        )

    def show_summary(self, folder, filenames):
        """
        Show the result state: output folder + list of saved files.

        Parameters
        ----------
        folder : str or PathLike
            The output folder the files were saved into.

        filenames : list of str
            Names of the files that were written, in the order to
            display them.

        Returns
        -------
        None
            Updates the summary state and switches the stack to show it.
        """
        self._folder_lbl.setText(f"Output folder: {folder}")
        self._file_list.clear()
        self._file_list.addItems(filenames)
        self._stack.setCurrentIndex(2)


class StepBlock(QWidget):
    """ One collapsible, numbered step section in the sidebar.

    Represents a single workflow step (e.g. "1  Input files",
    "2  XRD parameters") as a circular numbered badge plus title in a
    header, with that step's actual input fields and buttons stacked
    below in a body — the same building block used to construct every
    numbered section in the sidebar. Tracks and visually reflects one
    of three states (LOCKED, ACTIVE, DONE), enabling/disabling its
    body's contents and restyling the header/icon/title accordingly.

    Attributes
    ----------
    LOCKED, ACTIVE, DONE : int
        The three possible states, as simple integer constants (0, 1,
        2 respectively) used as an enum via `set_state`.
    """
    # The three states: simple integer constants used as an enum
    LOCKED = 0; ACTIVE = 1; DONE = 2

    def __init__(self, number, title, parent=None):
        """Build the header (numbered badge + title) and empty body.

        The block starts in the LOCKED state, with the same styling
        `set_state(LOCKED)` would apply, since the constructor needs
        some initial appearance before `set_state` is ever explicitly
        called.

        Parameters
        ----------
        number : int
            The step's number, shown in the circular badge and
            prefixed to the title text (e.g. `1` for "1  Input
            files").
            
        title : str
            The step's title text, shown next to the numbered badge.
            
        parent : QWidget or None, optional, default=None
            Standard Qt parent widget.
        """
        # building the header and body
        super().__init__(parent)
        
        # store the step's number/title 
        self._number = number
        self._title  = title
        
        # Tracks every widget added via add(), so set_state can loop
        # over them to enable/disable each one individually.
        self._children = []

        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0); outer.setSpacing(0)

        # Header — styled via setStyleSheet only on this widget (no children)
        self._header = QWidget(); self._header.setFixedHeight(48)
        self._header.setStyleSheet("background:transparent;")
        hl = QHBoxLayout(self._header); hl.setContentsMargins(16, 0, 16, 0)
        
        # The circular numbered badge - a QLabel showing just the step number
        # as text, turining a square label into a perfect circle. 
        self._icon = QLabel(str(number))
        self._icon.setFixedSize(22, 22)
        self._icon.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self._icon.setStyleSheet(
            "font-size:10px; font-weight:700; color:#c7c7cc;"
            "background:#f5f5f7; border-radius:11px; border:1.5px solid #e5e5ea;")
        
        # the step's title text, starting in the same muted locked-state gray
        self._title_lbl = QLabel(f"{number}  {title}")
        self._title_lbl.setStyleSheet("font-size:13px; color:#aeaeb2;")
        
        # lays out icon, small gap, title, then a stretch to push everything 
        # left and leave the remaining header space empty
        hl.addWidget(self._icon); hl.addSpacing(6)
        hl.addWidget(self._title_lbl); hl.addStretch()
        outer.addWidget(self._header)

        # Body
        # NEVER call setEnabled or setStyleSheet on this container :: every 
        # child gets setEnabled() called on it individually, and self._body 
        # itself is never touched that way => setting setEnabled on each child 
        # directly is more verbose, but fully predictable — each widget's own 
        # enabled/disabled QSS rules apply exactly as designed, with no 
        # container-level side effects layered on top.
        self._body = QWidget()
        self._body_lay = QVBoxLayout(self._body)
        self._body_lay.setContentsMargins(44, 6, 16, 20)
        self._body_lay.setSpacing(8)
        outer.addWidget(self._body)
        
        # a thin separator line below each block, and the initial state set to
        # locked - matching the icon/title styling already applied above at 
        # construction time
        outer.addWidget(divider())
        self._state = self.LOCKED

    def add(self, widget):
        """Add a child widget to this step's body, starting disabled.

        Adds `widget` to the body's layout, tracks it in
        `self._children` so `set_state` can find it later, and
        disables it immediately — every block starts locked, so the
        widget stays disabled until a subsequent `set_state(ACTIVE)`
        (or `DONE`) call re-enables it.

        Parameters
        ----------
        widget : QWidget
            The widget to add to this step's body (e.g. a field built
            via `field()`, or a button built via `make_btn()`).
        """
        # Three things happen for every widget added to a step's body:
        # (1) Add it into the body's layout
        self._body_lay.addWidget(widget)
        
        # (2) Track it in _children (so set_state can find it later)
        self._children.append(widget)
        
        # (3) Disable it immediately 
        widget.setEnabled(False)  # start locked; set_state will enable

    def set_state(self, state):
        """ Transition this block to a new state and update its appearance.

        Enables every tracked child widget if the new state is
        anything other than LOCKED (both ACTIVE and DONE mean
        "usable" — a completed step's fields stay editable in case the
        user wants to go back and adjust something), enabling/
        disabling each child individually rather than the body
        container itself, consistent with the warning in `__init__`.
        Then restyles the numbered badge, title, and header background
        to match the new state.

        Parameters
        ----------
        state : int
            One of `StepBlock.LOCKED`, `StepBlock.ACTIVE`, or
            `StepBlock.DONE`.

        Notes
        -----
        Visual appearance per state:

        - LOCKED: badge shows the plain step number, styled light
          gray-on-light-gray with a thin border; title in muted gray,
          normal weight; header background transparent.
        - ACTIVE: badge shows the step number, solid blue background
          with white text (the app's accent color); title also blue
          and semi-bold; header gets a light blue tint — the highlight
          marking "you are here" in the sidebar.
        - DONE: badge's text is replaced with a checkmark, solid green
          background (the same green used elsewhere in the app for
          accepted/in-use runs); title reverts to near-black, normal
          weight; header background back to transparent.
        """
        # Store the new state
        self._state = state
        
        # Flag :: if the new state is anything other than locked
        enabled = state != self.LOCKED
        
        # Enable/disable each child individually — never the body container
        for w in self._children:
            w.setEnabled(enabled)
                        
        # Locked appearance :: style only labels, never containers
        if state == self.LOCKED:
            self._icon.setText(str(self._number))
            self._icon.setStyleSheet(
                "font-size:10px; font-weight:700; color:#c7c7cc;"
                "background:#f5f5f7; border-radius:10px; border:1.5px solid #e5e5ea;"
                )
            self._title_lbl.setStyleSheet(
                "font-size:13px; color:#aeaeb2; font-weight:normal;"
                )
            self._header.setStyleSheet("background:transparent;")
        
        # Active appearance :: this is the subtle highlight that visually marks
        # "you are here" in the sidebar.
        elif state == self.ACTIVE:
            self._icon.setText(str(self._number))
            self._icon.setStyleSheet(
                "font-size:10px; font-weight:700; color:#ffffff;"
                "background:#0071e3; border-radius:10px; border:none;")
            self._title_lbl.setStyleSheet(
                "font-size:13px; color:#0071e3; font-weight:600;"
                )
            self._header.setStyleSheet("background:#f0f6ff;")
            
        # Done appearance
        elif state == self.DONE:
            self._icon.setText("✓")
            self._icon.setStyleSheet(
                "font-size:10px; font-weight:700; color:#ffffff;"
                "background:#34c759; border-radius:10px; border:none;")
            self._title_lbl.setStyleSheet(
                "font-size:13px; color:#1d1d1f; font-weight:normal;")
            self._header.setStyleSheet("background:transparent;")

# ── Tutorial texts ────────────────────────────────────────────────────────────
# Extracted from ediff.center verbose=1 printouts

# A dictionary of dictionaries, one entry per interactive method that has its 
# own keyboard shortcuts to explain — manual center detection, manual refinement, 
# ellipse correction, interactive background subtraction. Each entry has a title 
# and a steps list of (badge_text, description) pairs. 

# This content originates from ediff's own console/terminal help text, 
# repackaged here as structured data so it can be rendered as a proper GUI 
# dialog instead of scrolling past in a console log. 

# Keeping this as pure data, separate from TutorialDialog itself, means adding 
# a new tutorial for some future interactive method is just adding a new 
# dictionary entry — no widget code to write.

TUTORIAL_TEXTS = {
    "manual_detection": {
        "title": "Manual center detection — how to use",
        "steps": [
            ("Goal", "Click 3 points on a diffraction ring to define a circle. The center of that circle becomes the detected center."),
            ("Key  1", "Press '1' to mark the point at the current cursor position."),
            ("Key  2", "Press '2' to delete the last point."),
            ("Key  3", "Press '3' to delete the point closest to the cursor."),
            ("Key  d", "Press 'd' when done (3 points selected) to confirm and close."),
            ("Close", "Closing the figure without pressing 'd' cancels detection."),
        ],
    },
    "manual_refinement": {
        "title": "Manual center refinement — how to use",
        "steps": [
            ("Goal", "Nudge the circle to align it perfectly with a bright diffraction ring."),
            ("← → ↑ ↓", "Arrow keys move the center left / right / up / down."),
            ("Key  +", "Press '+' to increase the circle radius."),
            ("Key  -", "Press '-' to decrease the circle radius."),
            ("Key  b", "Press 'b' to increase the step size (coarse movement)."),
            ("Key  l", "Press 'l' to decrease the step size (fine movement)."),
            ("Key  d", "Press 'd' when satisfied to confirm and close."),
            ("Note", "Arrow key navigation in the matplotlib toolbar is disabled while this dialog is open — use the keyboard shortcuts above."),
        ],
    },
    "ellipse_correction": {
        "title": "Ellipse correction — how to use",
        "steps": [
            ("Goal", "Click 4 points on an elliptical diffraction ring. The software fits a perspective transform that maps those 4 points to a perfect circle."),
            ("Choosing points", "Pick points evenly spaced around the ring (roughly at 12, 3, 6 and 9 o'clock positions) for the best fit."),
            ("Accuracy", "After clicking each point it is automatically snapped to the nearest local intensity maximum."),
            ("Result", "The corrected image replaces the raw pattern for all subsequent steps. You can still skip and use the original."),
        ],
    },
    "interactive_bkg": {
        "title": "Interactive background subtraction — keyboard shortcuts",
        "steps": [
            ("Goal", "Place background control points along the baseline of the 1D diffractogram. The program fits a spline through them and subtracts it."),
            ("Key  1", "Add a background point at the current mouse cursor position."),
            ("Key  2", "Delete the background point closest to the mouse cursor."),
            ("Key  3", "Show the plot with all current background points marked."),
            ("Key  4", "Preview background as a linear spline through your points."),
            ("Key  5", "Preview background as a quadratic spline through your points."),
            ("Key  6", "Preview background as a cubic spline (recommended)."),
            ("Key  a", "Load previously saved background points from file."),
            ("Key  b", "Save current background points to the .bp file."),
            ("Key  t", "Subtract current background and save result to TXT file."),
            ("Key  u", "Subtract current background and update the profile in memory."),
            ("Close", "Close the window when done. Results are auto-saved."),
        ],
    },
}

class TutorialDialog(QWidget):
    """Floating, always-on-top tutorial panel for an interactive step.

    Renders one entry from the module-level TUTORIAL_TEXTS dict as a
    small reference window: a title, an accent bar, a grid of
    (badge, description) rows describing each available keyboard
    shortcut, and a Close button. Sized dynamically based on how many
    steps the requested tutorial has, since different interactive
    methods (manual detection, manual refinement, ellipse correction,
    interactive background subtraction) have very different numbers of
    shortcuts to explain.

    Deliberately a real top-level window with WindowStaysOnTopHint
    (not embedded in the main app window), since its purpose is to
    stay visible as a live reference while the user simultaneously
    interacts with a separate interactive matplotlib window (e.g.
    clicking points on a diffraction ring) — it would be useless as a
    reference if it could be buried behind that window.

    Attributes
    ----------
    _KEY_PREFIXES : tuple of str
        Prefixes used to classify each step: a step whose badge text
        starts with any of these ("Key", or an arrow symbol) is styled
        as a keyboard-keycap badge; anything else (e.g. "Goal", "Note")
        is styled as a plain blue label badge instead.
    """

    # Classify step keys: keyboard keys get a 'key' badge, others get a 'label'
    # badge
    _KEY_PREFIXES = ("Key", "←", "→", "↑", "↓")

    def __init__(self, key, parent=None):
        """Build the tutorial dialog for one entry in TUTORIAL_TEXTS.

        Parameters
        ----------
        key : str
            Key into the module-level TUTORIAL_TEXTS dict (e.g.
            "manual_detection", "ellipse_correction"). If not found,
            the dialog is built with an empty title and no steps
            rather than raising.
        parent : QWidget or None, optional, default=None
            Standard Qt parent widget.
        """
        # two combined window flags:
        #    - Window = genuine top-level window, not embedded inside the main
        #               app window
        #    - WindowStaysOnTopHint = keeps it floating above other windows
        super().__init__(
            parent, Qt.WindowType.Window | Qt.WindowType.WindowStaysOnTopHint
            )
        
        # Looks up the requested tutorial by its key, defaulting to an empty
        # dict/list if the key is not found (passing an unrecognized key 
        # produces an empty, mostly-blank dialog rather than a KeyError crash.)
        data = TUTORIAL_TEXTS.get(key, {})
        steps = data.get('steps', [])
        self.setWindowTitle(f"Tutorial — {data.get('title','')}")

        # Width fixed; height: header ~70px + per-row ~58px + close btn ~44px
        # height is computed, since different tutorials have different numbers
        # of steps
        fixed_w = 560
        fixed_h = 80 + len(steps) * 62 + 44
        self.setFixedSize(fixed_w, fixed_h)
        self.setStyleSheet(
            "QWidget { background:#ffffff; font-family:'Segoe UI',Arial,sans-serif; }"
        )

        lay = QVBoxLayout(self)
        lay.setContentsMargins(20, 14, 20, 14)
        lay.setSpacing(8)

        # Heading
        heading = QLabel(data.get('title', ''))
        heading.setStyleSheet("font-size:14px; font-weight:700; color:#1d2d44;")
        heading.setWordWrap(True)
        lay.addWidget(heading)

        accent = QFrame()
        accent.setFixedHeight(2)
        accent.setStyleSheet("background:#0071e3; border:none;")
        lay.addWidget(accent)

        # Use a grid so the badge column auto-sizes to the widest badge
        # there are genuinely two columns with different sizing needs: a narrow 
        # badge column on the left, a wide description column on the right.
        from qtpy.QtWidgets import QGridLayout
        grid_w = QWidget(); grid_w.setStyleSheet("background:transparent;")
        grid = QGridLayout(grid_w)
        grid.setContentsMargins(0, 4, 0, 4)
        grid.setHorizontalSpacing(12)
        grid.setVerticalSpacing(10)
        grid.setColumnStretch(1, 1)   # description column stretches
        grid.setColumnMinimumWidth(0, 110)  # badge column: minimum but can grow
        
        # For each step, classify it using the constant from earlier 
        # is_key becomes True for things like "Key  1" or "← → ↑ ↓", 
        #                False for things like "Goal" or "Note"
        for row_idx, (step_key, desc) in enumerate(steps):
            is_key = any(step_key.startswith(p) for p in self._KEY_PREFIXES)
            
            # Strips the literal "Key  " prefix text (note the double space, 
            # matching exactly how it's written in TUTORIAL_TEXTS) before 
            # displaying it — so the data dictionary can keep "Key  1" as 
            # a clear, self-documenting string, while the badge itself just 
            # shows the short "1".
            kl = QLabel(step_key.replace("Key  ", ""))
            kl.setAlignment(Qt.AlignmentFlag.AlignCenter)
            kl.setSizePolicy(
                QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Fixed
                )
            
            if is_key:
                kl.setStyleSheet(
                    "font-size:12px; font-weight:700; color:#1d1d1f;"
                    "background:#f0f0f0;"
                    "border:1px solid #aaaaaa;"
                    "border-bottom:3px solid #888888;"
                    "border-radius:5px;"
                    "padding:3px 10px;")
            else:
                kl.setStyleSheet(
                    "font-size:11px; font-weight:700; color:#0071e3;"
                    "background:#e8f0fe;"
                    "border-radius:5px;"
                    "padding:3px 10px;")
            
            # The description label 
            dl = QLabel(desc)
            dl.setStyleSheet(
                "font-size:12px; color:#1d1d1f; background:transparent;"
                )
            dl.setWordWrap(True)
            dl.setAlignment(
                Qt.AlignmentFlag.AlignVCenter | Qt.AlignmentFlag.AlignLeft
                )

            grid.addWidget(kl, row_idx, 0, Qt.AlignmentFlag.AlignVCenter)
            grid.addWidget(dl, row_idx, 1, Qt.AlignmentFlag.AlignVCenter)

        lay.addWidget(grid_w)
        lay.addStretch()
        close_btn = QPushButton("Close")
        close_btn.setFixedWidth(80)
        close_btn.clicked.connect(self.close)
        lay.addWidget(close_btn, alignment=Qt.AlignmentFlag.AlignRight)


def make_tutorial_btn(key):
    """Build a "? Help" button that opens (or refocuses) a TutorialDialog.

    Returns a small, styled button that, when clicked, shows the
    TutorialDialog for `key` — creating a new one only if none exists
    yet or the previous one has been closed, and simply bringing an
    already-open one to the front and giving it focus otherwise. This
    guarantees at most one tutorial window per button is ever open at
    a time, no matter how many times it's clicked.

    Parameters
    ----------
    key : str
        Key into the module-level TUTORIAL_TEXTS dict, passed through
        to TutorialDialog when it's created.

    Returns
    -------
    QPushButton
        The "? Help" button, ready to add to a layout. Its click
        handler manages the associated TutorialDialog's lifecycle
        internally via a closure — no external tracking of the dialog
        is needed by the caller.

    Notes
    -----
    The dialog reference is held in a single-element list (`_dialog = [None]`) 
    rather than a plain variable, specifically so the nested `_show()` function 
    can reassign it across multiple calls without needing the `nonlocal` 
    keyword — mutating an already-referenced list's contents doesn't require it, 
    while reassigning an enclosing-scope variable directly would.
    """
    btn = QPushButton("?  Help")
    btn.setObjectName("secondary")
    btn.setToolTip(
        "Show step-by-step instructions for this interactive method"
        )
    btn.setFixedHeight(28)
    btn.setStyleSheet(
        "QPushButton { background:#f0f6ff; color:#0071e3; border:1px solid #b3d0f7;"
        "  border-radius:6px; font-size:11px; font-weight:600; padding:0 10px; }"
        "QPushButton:hover { background:#e0ecff; }"
    )
    
    # Python's closures can read a variable from an enclosing scope freely, but 
    # reassigning a plain variable from within a nested function requires 
    # the nonlocal keyword (nonlocal _dialog) — without it, 
    # _dialog = TutorialDialog(key) inside _show() would create a brand-new 
    # local variable inside _show, shadowing the outer one instead of actually 
    # updating it, and the "remembered" dialog reference would silently never 
    # persist between clicks. Wrapping it in a single-element list sidesteps 
    # that entirely: _dialog itself (the list object) is never reassigned — only 
    # its contents are mutated (_dialog[0] = ...), and mutating an object you 
    # already have a reference to doesn't require nonlocal at all.
    # => each call to make_tutorial_btn(key) creates its own independent 
    #    _dialog slot, remembered across every click of that specific button 
    #    (since _show is a closure capturing its own _dialog list)
    _dialog = [None]
    
    def _show():
        # if a dialog was created dbefore but has since been closed, forget it
        if _dialog[0] and not _dialog[0].isVisible():
            _dialog[0] = None
            
        # if there is no remembered dialog at all (either never created, or 
        # just reset), build a fresh one and show it
        if _dialog[0] is None:
            _dialog[0] = TutorialDialog(key)
            _dialog[0].show()
            
        # otherwise (dialog exists and is still open) do not create a second
        # one, instead bring the existing one to the front and give it 
        # keyboard focus
        else:
            _dialog[0].raise_()
            _dialog[0].activateWindow()

    # Net effect: clicking the help button repeatedly never spawns multiple 
    # copies of the same tutorial window pilling up on screen. It always shows
    # exactly one, creating it fresh only when needed, and simply refocusing 
    # the existing one otherwise  
          
    btn.clicked.connect(_show)
    return btn



# =============================================================================
# About dialog

# QDdialog is QT's purpose-built class for modal, short-lived informational
# or confirmation windows
# =============================================================================
class AboutDialog(QDialog):
    """Simple modal About dialog: app description, PyPI link, funding.

    Shows the bundled funding-body logo, the app name and tagline, a
    clickable "View on PyPI" hyperlink (opened in the system's default
    browser), a thin divider, and the project's funding acknowledgement
    text, followed by a Close button.

    Deliberately a QDialog rather than a plain QWidget (unlike
    TutorialDialog): this is a short-lived, modal informational popup
    with a single dismissal action, which is exactly what QDialog's
    built-in accept()/reject() semantics are designed for — the Close
    button below is wired directly to `self.accept()`.
    """

    def __init__(self, parent=None):
        """
        Build the About dialog's fixed-size layout.

        Parameters
        ----------
        parent : QWidget or None, optional, default=None
            Standard Qt parent widget.
        """
        super().__init__(parent)
        
        # Window settings
        self.setWindowTitle("About EDIFF")
        self.setFixedSize(380, 400)
        self.setStyleSheet(
            "QDialog { background:#ffffff; font-family:'Segoe UI',Arial,sans-serif; }"
        )

        lay = QVBoxLayout(self)
        lay.setContentsMargins(32, 30, 32, 24)
        lay.setSpacing(4)
        lay.setAlignment(Qt.AlignmentFlag.AlignHCenter)
        
        # Load the funding-body logo via _pixmap
        logo_lbl = QLabel()
        logo_lbl.setStyleSheet("background:transparent;")
        logo_px = _pixmap("tacr_logo.png")
        if not logo_px.isNull():
            logo_px = logo_px.scaled(
                72, 72,
                Qt.AspectRatioMode.KeepAspectRatio,
                Qt.TransformationMode.SmoothTransformation)
            logo_lbl.setPixmap(logo_px)
        logo_lbl.setFixedSize(72, 72)
        logo_lbl.setAlignment(Qt.AlignmentFlag.AlignCenter)
        lay.addWidget(logo_lbl, alignment=Qt.AlignmentFlag.AlignHCenter)
        lay.addSpacing(14)
        
        # The app name, tagline
        title_lbl = QLabel("EDIFF")
        title_lbl.setStyleSheet(
            "background:transparent; font-size:19px; font-weight:700; color:#1d1d1f;")
        title_lbl.setAlignment(Qt.AlignmentFlag.AlignHCenter)
        lay.addWidget(title_lbl)

        subtitle_lbl = QLabel("Process electron diffraction patterns")
        subtitle_lbl.setStyleSheet(
            "background:transparent; font-size:12px; color:#86868b;")
        subtitle_lbl.setAlignment(Qt.AlignmentFlag.AlignHCenter)
        lay.addWidget(subtitle_lbl)
        lay.addSpacing(8)
        
        # Clickable hyperlin
        link_lbl = QLabel(
            '<a href="https://pypi.org/project/ediff/" '
            'style="color:#0071e3; text-decoration:none;">View on PyPI</a>')
        link_lbl.setStyleSheet("background:transparent; font-size:12px;")
        link_lbl.setAlignment(Qt.AlignmentFlag.AlignHCenter)
        link_lbl.setOpenExternalLinks(True)
        lay.addWidget(link_lbl)
        lay.addSpacing(20)
        
        # Functionally identical to the module-level divider() helper, but it's
        # built inline here with its own separate code
        accent = QFrame()
        accent.setFixedHeight(1)
        accent.setStyleSheet("background:#e5e5ea; border:none;")
        lay.addWidget(accent)
        lay.addSpacing(16)
        
        # The funding acknowledgement text itself
        funding_lbl = QLabel(
            "The development of the package was co-funded by "
            "the Technology Agency of the Czech Republic, "
            "program NCK, project TN02000020."
        )
        funding_lbl.setStyleSheet(
            "background:transparent; font-size:12px; color:#3a3a3c;")
        funding_lbl.setAlignment(Qt.AlignmentFlag.AlignHCenter)
        funding_lbl.setWordWrap(True)
        lay.addWidget(funding_lbl)
        
        # A stretch to push the Close button down to the bottom of the fixed
        # size dialog
        lay.addStretch()

        close_btn = make_btn("Close", primary=True)
        close_btn.clicked.connect(self.accept)
        lay.addWidget(close_btn)



# =============================================================================
# Startup dialog

# This is the very first screen a user sees — the mode-selection dialog shown 
# before the main window even exists. 
# =============================================================================
class StartupDialog(QWidget):
    """First screen shown at app launch: pick polycrystalline or 
    monocrystalline mode.

    A branded splash-like dialog (logo, app name/tagline/accent bar,
    matching the visual language of `make_splash()`'s painted splash
    screen) followed by two large, image-based mode-selection buttons.
    Only "Polycrystalline" is currently enabled; "Monocrystalline" is
    shown grayed out with a "Coming soon" tooltip and an explanatory
    note beneath the buttons, since that mode is not yet implemented.

    Emits `mode_selected` with the chosen mode string and closes itself
    once the user picks an enabled option — the caller (see `main()`)
    is responsible for constructing the actual main application window
    in response.

    Signals
    -------
    mode_selected(str)
        Emitted with `"polycrystal"` (currently the only reachable
        value, since the other option is disabled) when the user
        clicks an enabled mode button.
    """
    # signal carrying a string — this is what the outer code uses to find out 
    # which mode the user picked
    mode_selected = Signal(str)

    def __init__(self):
        """Build the logo header, mode prompt, and the two mode buttons.

        Takes no arguments beyond `self` — this dialog is always constructed 
        the same way, with no parent (it's a genuine top-level window shown 
        before any other part of the application exists).
        """
        super().__init__()
        self.setWindowTitle("EDIFF")
        self.setFixedSize(520, 440)
        self.setStyleSheet(
            "QWidget { background:#ffffff; font-family:'Segoe UI',Arial,sans-serif; }"
            )
        self.setWindowFlag(Qt.WindowType.WindowContextHelpButtonHint, False)
        
        # Startup window settings
        lay = QVBoxLayout(self)
        lay.setContentsMargins(40, 24, 40, 28)
        lay.setSpacing(0)

        # Logo header
        header_row = QWidget()
        header_row.setStyleSheet("background:transparent;")
        hr_lay = QHBoxLayout(header_row)
        hr_lay.setContentsMargins(0, 0, 0, 0); hr_lay.setSpacing(16)
        hr_lay.setAlignment(Qt.AlignmentFlag.AlignVCenter)
        
        # load-scale-fallback
        logo_lbl = QLabel()
        logo_px = _pixmap("ediff.png")
        if not logo_px.isNull():
            logo_px = logo_px.scaled(72, 72,
                Qt.AspectRatioMode.KeepAspectRatio,
                Qt.TransformationMode.SmoothTransformation)
            logo_lbl.setPixmap(logo_px)
        logo_lbl.setFixedSize(72, 72)
        hr_lay.addWidget(logo_lbl)
        
        # A vertical column stacking the large app name, a smaller gray tagline, 
        # then a small 40×3px blue accent bar underneath
        text_col = QWidget(); text_col.setStyleSheet("background:transparent;")
        tc_lay = QVBoxLayout(text_col)
        tc_lay.setContentsMargins(0, 0, 0, 0); tc_lay.setSpacing(2)
        title = QLabel("EDIFF")
        title.setStyleSheet(
            "font-size:28px; font-weight:700; color:#1d2d44; letter-spacing:1px;"
            )
        sub = QLabel("Crystallographic Analysis from 4D STEM")
        sub.setStyleSheet(
            "font-size:11px; color:#86868b; letter-spacing:0.3px;"
            )
        accent = QFrame()
        accent.setFixedSize(40, 3)
        accent.setStyleSheet("background:#0071e3; border:none;")
        tc_lay.addWidget(title)
        tc_lay.addWidget(sub)
        tc_lay.addSpacing(6)
        tc_lay.addWidget(accent)
        
        # Logo, then text column, then a stretch absorbing any remaining 
        # horizontal space
        hr_lay.addWidget(text_col)
        hr_lay.addStretch()
        lay.addWidget(header_row)
        lay.addSpacing(20)
        
        # The mode-selection prompt and buttons
        lbl = QLabel("Select analysis mode:")
        lbl.setStyleSheet("font-size:13px; font-weight:600; color:#1d1d1f;")
        lay.addWidget(lbl)
        lay.addSpacing(16)

        # Image buttons row
        btn_row = QWidget(); btn_row.setStyleSheet("background:transparent;")
        brl = QHBoxLayout(btn_row)
        brl.setContentsMargins(0, 0, 0, 0); brl.setSpacing(16)

        def img_btn(img_name, label, enabled, mode):
            """ Build one large, image-based mode-selection button.

            A small factory, defined locally since it's only needed here,
            building a 160×160 square button whose entire visible surface is
            a scaled icon (rather than a QLabel-plus-pixmap sitting inside a
            button, so the whole square is clickable as one unit, image
            included), with a text caption beneath it. Styling, tooltip, and
            click behavior all differ depending on whether this option is
            currently selectable.
        
            Parameters
            ----------
            img_name : str
                Filename of the bundled resource image to load and scale onto
                the button (e.g. "pcryst.png"), passed to `_pixmap`. If the
                image fails to load, the button is shown with no icon rather
                than raising.
                
            label : str
                Caption text shown centered beneath the button (e.g.
                "Polycrystalline").
                
            enabled : bool
                If True, the button is clickable, styled with the normal
                light card look plus hover/pressed states, given a
                "Start {label} analysis" tooltip, and wired to call
                `self._pick(mode)` when clicked. If False, the button is
                disabled, styled with a flatter, lower-contrast, semi-
                transparent look, given a "Coming soon" tooltip instead, and
                no click handler is connected at all.
                
            mode : str
                The mode string reported to `self._pick` if this button is
                clicked (irrelevant when `enabled` is False, since no click
                handler is attached in that case).
        
            Returns
            -------
            QWidget
                A widget containing the image button stacked above its
                caption label, both horizontally centered, ready to add to
                the mode-selection row.

            """
            w = QWidget()
            w.setStyleSheet("background:transparent;")
            wl = QVBoxLayout(w)
            wl.setContentsMargins(0,0,0,0)
            wl.setSpacing(8)
            wl.setAlignment(Qt.AlignmentFlag.AlignHCenter)
            
            # a large, square button whose entire visible content is an icon
            # the whole square is clickable as one unit, image included
            btn = QPushButton()
            btn.setFixedSize(160, 160)
            btn.setEnabled(enabled)
            
            # Load and scale image
            pix = _pixmap(img_name)
            if not pix.isNull(): # guard pattern
                pix = pix.scaled(140, 140,
                    Qt.AspectRatioMode.KeepAspectRatio,
                    Qt.TransformationMode.SmoothTransformation)
                btn.setIcon(QIcon(pix))
                btn.setIconSize(QSize(140, 140))
                
            # clear hover and pressed states, wired to actually call a mode
            # when clicked, with a tooltip describing what it does
            if enabled:
                btn.setStyleSheet(
                    "QPushButton{background:#f5f5f7; border:2px solid #e5e5ea;"
                    "border-radius:16px; padding:4px;}"
                    "QPushButton:hover{background:#e8f0fe; border-color:#0071e3;}"
                    "QPushButton:pressed{background:#dce4ed; border-color:#005bbf;}")
                btn.clicked.connect(lambda checked, m=mode: self._pick(m))
                btn.setToolTip(f"Start {label.lower()} analysis")
                
            # a flatter, grayer, lower-contrast look visually communicating 
            # unavailable beyond just the border/background colors. There is no
            # live handler waiting behind a button the user cannot actually 
            # press
            else:
                btn.setStyleSheet(
                    "QPushButton{background:#fafafa; border:2px solid #f0f0f0;"
                    "border-radius:16px; padding:4px; opacity:0.5;}")
                btn.setToolTip("Coming soon")
            
            # The text label under the image button
            lbl_w = QLabel(label)
            lbl_w.setAlignment(Qt.AlignmentFlag.AlignCenter)
            if enabled:
                lbl_w.setStyleSheet(
                    "font-size:12px; font-weight:600; color:#1d1d1f;"
                    )
            else:
                lbl_w.setStyleSheet(
                    "font-size:12px; font-weight:600; color:#aeaeb2;"
                    )

            wl.addWidget(btn, alignment=Qt.AlignmentFlag.AlignHCenter)
            wl.addWidget(lbl_w)
            return w
        
        # Two calls to the factory: Polycrystalline is enabled and reports 
        # "polycrystal"; Monocrystalline is disabled and reports "monocrystal"
        # — a mode string that, if it were ever clickable, would presumably be 
        # caught by the _reopen_with_mode logic that already explicitly handles
        # it with an "not yet implemented" message.
        poly_img = "pcryst.png"
        mono_img = "mcryst.png"

        brl.addStretch()
        brl.addWidget(
            img_btn(
                poly_img, 
                "Polycrystalline", 
                True,  
                "polycrystal"
                )
            )
        brl.addWidget(
            img_btn(
                mono_img, 
                "Monocrystalline", 
                False, 
                "monocrystal"
                )
            )
        brl.addStretch()
        lay.addWidget(btn_row)
        lay.addSpacing(14)
        
        # A small, muted, centered clarifying note beneath the buttons — 
        # explicitly telling the user why one of the two options looks grayed 
        # out
        note = QLabel("Monocrystalline mode coming soon.")
        note.setStyleSheet("font-size:10px; color:#aeaeb2;")
        note.setAlignment(Qt.AlignmentFlag.AlignCenter)
        lay.addWidget(note)

    def _pick(self, mode):
        """Report the chosen mode and close this dialog.

        Parameters
        ----------
        mode : str
            The mode string to emit via `mode_selected` (e.g.
            `"polycrystal"`).
        """
        self.mode_selected.emit(mode)
        self.close()


class EdiffWindow(QMainWindow):
    """Main application window: sidebar steps plus tabbed result plots.

    The top-level window for the entire session — constructed once at
    startup (after StartupDialog reports a mode) and torn down/rebuilt
    only via New Session or Change Mode. Owns the EdiffPipeline
    instance driving all actual ediff-library work, the sidebar of
    numbered StepBlocks, and the tabbed PlotCards showing each step's
    results.

    Attributes
    ----------
    TAB_XRD, TAB_ELLIPSE, TAB_ELD2D, TAB_1D, TAB_BKG, TAB_COMP : int
        Indices into `self._tabs` matching the order tabs are added in
        `_build_ui`, used throughout instead of bare integers (e.g.
        `self._tabs.setCurrentIndex(self.TAB_XRD)`).
        
    N_STEPS : int
        Total number of sidebar steps (7), used by step-progression
        logic (e.g. `_complete_step`) so it doesn't need editing if a
        step is ever added or removed.
    """

    # The main application window — the one that exists for the entire session, 
    # holding the sidebar of StepBlocks and the result PlotCards. 
    # QMainWindow is Qt's dedicated class for a top-level window with a central 
    # content area plus optional toolbars/status bar/menu bar
    
    # integer indices matching tab order in self._tabs (built late)
    TAB_PREVIEW=0
    TAB_XRD=1
    TAB_ELLIPSE=2
    TAB_ELD2D=3
    TAB_1D=4
    TAB_BKG=5
    TAB_COMP=6
    TAB_SAVE=7

    # the total step count => using this constant rather than a hard-coded 7
    # means the step-progression logic doesn't need editing if a step is ever 
    # added or removed.
    N_STEPS = 7

    # Interface sizes: UI_SCALES is what Ctrl+= / Ctrl+- step through, and
    # UI_SCALE_PRESETS is what the "Size" menu offers directly (see
    # _choose_ui_scale). 1.0 = the original size.
    UI_SCALES = (0.75, 0.9, 1.0, 1.15, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0)
    UI_SCALE_PRESETS = (1.0, 1.25, 1.5, 1.75, 2.0, 2.5)

    def __init__(self, pipeline=None):
        """Construct the window: state, pipeline, UI, then shortcuts.

        Initializes thse working directory (defaulting to the user's home folder,
        later changeable via Step 1's workdir picker) and the EdiffPipeline
        instance, all session-tracking state (loaded files, axis-limit memory,
        busy-step tracking, etc.), then builds every widget via `_build_ui`,
        activates Step 1, and finally wires up global keyboard shortcuts — in
        that order, since the shortcut handlers depend on widgets `_build_ui`
        creates.

        Parameters
        ----------
        pipeline : EdiffPipeline or None, optional
            A pipeline instance created earlier (e.g. in `main()`, right
            after QApplication is constructed, so its interactive worker
            subprocess gets the earliest possible head start on its
            cold-start imports). If None, a fresh EdiffPipeline is created
            here instead, using the default workdir.
        """
        super().__init__()
        self.setWindowTitle("EDIFF")
        self.resize(1440, 900)
        
        # Window icon (taskbar + title bar)
        _icon_px = _pixmap("ediff.png")
        if not _icon_px.isNull():
            self.setWindowIcon(QIcon(_icon_px))

        # Reuse a pipeline created earlier in main() (with its worker
        # already pre-warming in the background) if one was given,
        # instead of constructing a brand new one here -- that would
        # throw away the head start entirely.
        if pipeline is not None:
            self.pipeline = pipeline
            self.workdir = pipeline.workdir
        else:
            self.workdir = default_workdir()
            self.pipeline = EdiffPipeline(self.workdir)

        # The three pieces of file-loading state
        # (1) the list of loaded structure files
        self.cif_paths = []
        
        # (2) the currently selected diffraction image path
        self.eld_path  = None
        
        # (3) which path the pipeline's own loaded data currently reflects, as 
        #     distinct from eld_path => these can genuinely differ, and that 
        #     gap is deliberate: it's exactly the mechanism that lets the app 
        #     avoid reloading/reprocessing an image that's already loaded, 
        #     while still correctly detecting when the user has picked 
        #     a genuinely different one since the last load.
        
        self._loaded_eld_path = None   
        
        # Remembers the last axis limits explicitly applied, grouped by
        # what the axes actually mean — so setting xlim/ylim on one
        # step's plot carries over automatically to the next step that
        # shows a plot with the same kind of axes (e.g. 1D profile ->
        # Background, both "Pixels" vs "Intensity"), without carrying
        # over to unrelated axis types (e.g. XRD/Comparison's "q [1/Å]").
        # Keys: 'pixels_intensity', 'q_intensity', 'image_2d'.
        self._axis_memory = {}
        
        # Load the small persisted-preferences dict from disk (guaranteed to 
        # return something usable (an empty dict at worst) rather than ever 
        # raising)
        self._cfg      = load_config()

        # Every start begins in the user's own directory. Folders remembered
        # from earlier sessions are not used: the "last folder" settings the
        # load/save dialogs read are set to the starting directory instead.
        # Within this session they update as usual (each dialog opens where
        # the previous one ended).
        _start_dir = str(default_workdir())
        for _key in ("last_cif_dir", "last_eld_dir", "last_center_dir",
                     "last_save_dir"):
            self._cfg[_key] = _start_dir
        # When this session began. "Save all results" and "Download results"
        # copy only files created or changed since then, because the working
        # directory is now the user's own folder, full of unrelated files.
        import time
        self._session_started = time.time()
        
        # Session state:
        # (1) whether the ellipse-correction preview has been shown yet
        # self._ellipse_preview_loaded = False

        # (2) which step (if any) currently has a background computation running
        self._current_step_n = None
        
        # (3) whether Tutorial mode is toggled on
        self._tutorial_mode  = False
        
        # (4) the busy-tracking set that run_in_thread/_step_busy rely on to 
        #     prevent a second computation launching for a step while the first 
        #     is still silently finishing
        self._busy_keys      = set()
        
        # Ordered construction steps:
        # (1) Build every widget first
        self._build_ui()
        
        # (2) Activate step 1
        self._go_to_step(1)
        
        # (3) Wire up the global keyboard shortcuts
        self._add_shortcuts()

    def _add_shortcuts(self):
        """Wire up the app's global keyboard shortcuts.

        Binds Ctrl+O/E/S/D/N (browse structure files, browse diffraction image, 
        save, download, new session respectively) to their handlers via 
        QShortcut, each attached to the whole main window so they fire 
        regardless of which widget currently has focus. Must be called after 
        `_build_ui`, since the connected handlers depend on widgets it creates.
        """
        QShortcut(QKeySequence("Ctrl+O"), self).activated.connect(
            self._browse_cif
            )
        QShortcut(QKeySequence("Ctrl+E"), self).activated.connect(
            self._browse_eld
            )
        QShortcut(QKeySequence("Ctrl+S"), self).activated.connect(
            self._on_save
            )
        QShortcut(QKeySequence("Ctrl+D"), self).activated.connect(
            self._on_download
            )
        QShortcut(QKeySequence("Ctrl+N"), self).activated.connect(
            self._on_new_session
            )

        # Interface size: larger / smaller / reset. "Ctrl+=" is the same
        # key as "Ctrl++" without Shift, which is how most keyboards type it.
        # These only save the size; use the Size menu to restart and apply.
        for seq in (QKeySequence.StandardKey.ZoomIn, "Ctrl+="):
            QShortcut(QKeySequence(seq), self).activated.connect(
                lambda: self._change_ui_scale(+1))
        for seq in (QKeySequence.StandardKey.ZoomOut, "Ctrl+-"):
            QShortcut(QKeySequence(seq), self).activated.connect(
                lambda: self._change_ui_scale(-1))
        QShortcut(QKeySequence("Ctrl+0"), self).activated.connect(
            lambda: self._set_ui_scale(DEFAULT_UI_SCALE))

    # ── Build UI ──────────────────────────────────────────────────────────────
    def _build_ui(self):
        """Construct the entire main window: top bar, sidebar, result tabs.

        Assembles the window's central widget as a vertical stack of the top
        bar (see `_build_topbar`) above a horizontal QSplitter, whose left pane
        is a scrollable sidebar containing all seven StepBlocks (built via 
        `_build_step1` through `_build_step7`, stored in `self._blocks` keyed 
        by step number) and whose right pane is a QTabWidget of six PlotCards
        — one per result view (XRD, Ellipse, ELD 2D, ELD 1D, Background, 
        Comparison) — in the exact order matching the TAB_* class constants.
    
        Every step block starts locked (per StepBlock's own default state)
        and every result tab starts disabled; nothing becomes usable until
        `_go_to_step(1)` (called immediately after this, in `__init__`)
        and subsequent step completions unlock them.
    
        Must run before `_add_shortcuts`, since that method's handlers
        depend on widgets created here — and the Ellipse/Center colormap
        sync (`self._ell_cmap`/`self._ctr_cmap`) can only be wired up
        partway through this method, after the step-building loop that
        creates those specific widgets has already run.
        """
        # Create central widget (a plain container holding a vertical layout)
        central = QWidget()
        self.setCentralWidget(central)
        root = QVBoxLayout(central)
        root.setContentsMargins(0, 0, 0, 0); root.setSpacing(0)
        root.addWidget(self._build_topbar())
        
        # Widget for two or more panes with a user-draggable divider between 
        # them
        splitter = QSplitter(Qt.Orientation.Horizontal)
        splitter.setHandleWidth(1)
        splitter.setStyleSheet("QSplitter::handle { background: #e5e5ea; }")

        # Sidebar scroll area 
        # necessary because seven step blocks, each with its own fields, will 
        # almost certainly not all fit vertically at once; this gives it 
        # a vertical scrollbar for free.
        sidebar = QScrollArea()
        sidebar.setObjectName("sidebar")
        
        # resize its contained widget to fill the available width
        sidebar.setWidgetResizable(True) 
        sidebar.setMinimumWidth(360); sidebar.setMaximumWidth(560)
        sidebar.setHorizontalScrollBarPolicy(
            Qt.ScrollBarPolicy.ScrollBarAlwaysOff
            )
        
        # The actual scrollable content: 
        # a plain widget holding a vertical layout
        steps_w = QWidget()
        s_lay = QVBoxLayout(steps_w)
        
        # small top margin, larger bottom margin, zero spacing
        s_lay.setContentsMargins(0, 8, 0, 24); s_lay.setSpacing(0)
        sidebar.setWidget(steps_w)
        
        # Loop that builds all seven StepBlocks
        self._blocks = {}
        for n, builder in [
            # iterating over a list of (step number, builder method) pairs
            (1, self._build_step1), (2, self._build_step2),
            (3, self._build_step3), (4, self._build_step4),
            (5, self._build_step5), (6, self._build_step6), 
            (7, self._build_step7),
        ]:
            # calling each builder method
            b = builder()
            
            # storing the results keyed by step number
            self._blocks[n] = b
            
            # adding it to the sidebar's layout
            s_lay.addWidget(b)
        s_lay.addStretch()

        # Result tabs
        # Six PlotCard instances, each given its own subtitle string
        self._tabs = QTabWidget(); self._tabs.setDocumentMode(True)

        self._card_preview = PlotCard(
            "Raw diffraction pattern — inspect before proceeding"
        )
        self._card_xrd      = PlotCard(
            "Calculated from crystal structure files"
            )
        self._card_ellipse  = PlotCard(
            "Elliptical distortion correction (optional)"
            )
        self._card_eld2d    = PlotCard(
            "Experimental electron diffraction pattern"
            )
        self._card_1d       = PlotCard(
            "Radially averaged intensity profile"
            )
        self._card_bkg      = PlotCard(
            "Original / with background / corrected"
            )
        self._card_comp     = PlotCard(
            "Calibrated ELD vs theoretical XRD"
            )
        self._card_save = SaveSummaryCard(
            "Output files written to disk"
        )
        
        # Adds all six cards as tabs in the exact order that matches TAB_XRD
        # through TAB_COMP's numeric values from the class-level constants.
        for lbl, card in [
            ("Preview", self._card_preview),
            ("XRD", self._card_xrd),
            ("Ellipse", self._card_ellipse),
            ("ELD 2D", self._card_eld2d),
            ("ELD 1D", self._card_1d),
            ("Background", self._card_bkg),
            ("Comparison", self._card_comp),
            ("Save", self._card_save),
        ]:
            self._tabs.addTab(card, lbl)
            
        
        # Every tab starts disabled, matching the sidebar's own "everything
        # starts locked" philosophy from StepBlock.
        for i in range(self._tabs.count()):
            self._tabs.setTabEnabled(i, False)
        
        # Assembling the splitter and finishing up
        # splitter on the left, tabs on the right
        splitter.addWidget(sidebar); splitter.addWidget(self._tabs)
        
        # control of what happens when the whole window is resized
        splitter.setStretchFactor(0, 0); splitter.setStretchFactor(1, 1)
        
        # initial pixel widths of each pane at startup
        splitter.setSizes([520, 920])
        splitter.setMinimumWidth(420)
        
        # add the fully-assembled splitter to the main layout, and sets 
        # the initial status bar message
        root.addWidget(splitter)
        self.statusBar().showMessage(
            "Load files to begin.   Ctrl+O = CIF  Ctrl+E = ELD  Ctrl+S = Save  Ctrl+N = New session")

    def _build_topbar(self):
        """Build the top bar: logo, step progress indicator, and action buttons
    
        Assembles the horizontal strip shown across the top of the main
        window: the app logo and name, a live step-progress breadcrumb
        ("Files — XRD — Ellipse — Center — 1D — Background — Calibrate",
        each step's badge/label pair stored in `self._step_indicators` so
        other code can recolor them as steps become active/complete), and
        a row of action buttons — Download results, Abort, Tutorial mode
        toggle, a Session dropdown menu (New/Save/Load/Run all steps),
        Change mode, and About.
    
        Returns
        -------
        QWidget
            The fully assembled top bar, ready to add to the main window's
            layout.
        """
        # Basic setup
        import os
        bar = QWidget(); bar.setFixedHeight(52)
        bar.setStyleSheet(
            "background:#ffffff; border-bottom:1px solid #e5e5ea;"
            )
        lay = QHBoxLayout(bar)
        lay.setContentsMargins(16, 0, 20, 0)
        lay.setSpacing(10)

        # Logo image in topbar
        logo_lbl = QLabel()
        logo_px = _pixmap("ediff.png")
        if not logo_px.isNull():
            logo_px = logo_px.scaled(32, 32,
                Qt.AspectRatioMode.KeepAspectRatio,
                Qt.TransformationMode.SmoothTransformation)
            logo_lbl.setPixmap(logo_px)
        logo_lbl.setFixedSize(32, 32)
        logo_lbl.setStyleSheet("background:transparent;")
        lay.addWidget(logo_lbl)
    
        app_name = QLabel("EDIFF")
        app_name.setStyleSheet(
            "font-size:15px; font-weight:700; color:#1d2d44; letter-spacing:0.5px;"
            )
        lay.addWidget(app_name)
        
        # A vertical line
        sep = QFrame(); sep.setFrameShape(QFrame.Shape.VLine)
        sep.setFixedHeight(20)
        sep.setStyleSheet("background:#e5e5ea; border:none; max-width:1px;")
        lay.addWidget(sep); lay.addSpacing(4)
        
        # The step process indicator (a live-updating breadcrumb trail)
        # the row of small numbered labels visible across the very top 
        # of the project; building it from a plain list of names (rather than 
        # one widget per name hard-coded separately) keeps this compact and 
        # easy to maintain.
        self._step_indicators = []
        
        for i, name in enumerate(
                [
                    "Files",
                     "XRD",
                     "Ellipse",
                     "Center",
                     "1D",
                     "Background",
                     "Calibrate",
                     "Comparison",
                     "Save"
                 ]
                ):
            # it's False only for i == 0 (the very first item, "Files"), 
            # so a "—" separator gets inserted before every step except 
            # the first — exactly the right number of separators (6) for 7 items.
            if i:
                arr = QLabel("—")
                arr.setStyleSheet("color:#e5e5ea; font-size:11px; margin: 0 2px;")
                lay.addWidget(arr)
                
            cell = QWidget()
            cell.setStyleSheet("background:transparent;")
            cl = QHBoxLayout(cell)
            cl.setContentsMargins(0,0,0,0)
            cl.setSpacing(6)
            
            # For each step name, builds a small circular numbered badge   
            badge = QLabel(str(i+1))
            badge.setFixedSize(20, 20)
            badge.setAlignment(Qt.AlignmentFlag.AlignCenter)
            badge.setStyleSheet(
                "background:#e5e5ea; color:#aeaeb2; border-radius:10px;"
                "font-size:10px; font-weight:700;")
            lbl = QLabel(name)
            lbl.setStyleSheet("font-size:12px; color:#aeaeb2; font-weight:400;")
            cl.addWidget(badge); cl.addWidget(lbl)
            lay.addWidget(cell)
            
            # badge and lbl are stored as a tuple, which lets some other method 
            # elsewhere, reach back into this exact list later and recolor a 
            # specific step's badge/label once that step becomes active or 
            # completes, mirroring the same LOCKED/ACTIVE/DONE state concept
            # from StepBlock, just rendered as a compact horizontal strip 
            # instead of the full sidebar block.
            self._step_indicators.append((badge, lbl))
        
        # Push everything built so far to the left, leaving the remaining space 
        # on the right for the action buttons that follow.
        lay.addStretch()
        
        # Download button
        # built via the shared make_btn() factory, wired to _on_download
        self._btn_dl = make_btn(
            "Download results", 
            tip="Copy result files to a folder  (Ctrl+D)"
            )
        self._btn_dl.clicked.connect(self._on_download)
        
        # add to widget
        lay.addWidget(self._btn_dl)
        lay.addSpacing(8)
        
        # Abort button
        # becomes clickable once a step is actually running
        abort_style = (
            "QPushButton {"
            "  background: #fff0ee; color: #c0392b; border: 1px solid #f5c6c0;"
            "  border-radius: 14px; padding: 0px 14px; font-size: 12px; font-weight: 600;"
            "  min-height: 28px; max-height: 28px;"
            "}"
            "QPushButton:hover { background: #ffe0db; border-color: #e8a09a; }"
            "QPushButton:pressed { background: #ffd0ca; }"
            "QPushButton:disabled { background: #f5f5f7; color: #c7c7cc; border-color: #e5e5ea; }"
        )

        self._btn_abort = QPushButton("Abort")
        self._btn_abort.setToolTip(
            "Stop the current running step without clearing results.\n"
            "Use to interrupt a stuck or unwanted calculation.")
        self._btn_abort.setStyleSheet(abort_style)
        self._btn_abort.setEnabled(False)
        self._btn_abort.clicked.connect(self._on_abort)
        
        # add to widget
        lay.addWidget(self._btn_abort)
        lay.addSpacing(4)
        
        # Tutorial mode toggle -- a manually-managed two-state button
        tutorial_off_style = (
            "QPushButton {"
            "  background: #e8edf2; color: #3a3a3c; border: 1px solid #d1d1d6;"
            "  border-radius: 14px; padding: 0px 14px; font-size: 12px; font-weight: 500;"
            "  min-height: 28px; max-height: 28px;"
            "}"
            "QPushButton:hover { background: #dce4ed; }"
        )
        tutorial_on_style = (
            "QPushButton {"
            "  background: #0071e3; color: #ffffff; border: 1px solid #005bbf;"
            "  border-radius: 14px; padding: 0px 14px; font-size: 12px; font-weight: 600;"
            "  min-height: 28px; max-height: 28px;"
            "}"
            "QPushButton:hover { background: #0077ed; }"
        )
        self._tutorial_on_style  = tutorial_on_style
        self._tutorial_off_style = tutorial_off_style

        self._btn_tutorial_mode = QPushButton("Tutorial  ○")
        self._btn_tutorial_mode.setToolTip(
            "Toggle Tutorial mode.\n"
            "ON: step-by-step instructions open automatically\n"
            "    when an interactive method is selected.\n"
            "OFF: help buttons are still available on demand.")
        self._btn_tutorial_mode.setStyleSheet(tutorial_off_style)
        
        # We manage state manually -- Qt buttons have a built-in toggle 
        # mechanism
        self._btn_tutorial_mode.setCheckable(False)   
        self._btn_tutorial_mode.clicked.connect(self._on_toggle_tutorial)
        self._tutorial_mode = False
        lay.addWidget(self._btn_tutorial_mode)
        lay.addSpacing(6)

        pill_style = (
            "QPushButton {"
            "  background: #e8edf2; color: #3a3a3c; border: 1px solid #d1d1d6;"
            "  border-radius: 14px; padding: 0px 14px; font-size: 12px; font-weight: 500;"
            "  min-height: 28px; max-height: 28px;"
            "}"
            "QPushButton:hover { background: #dce4ed; border-color: #b0b8c4; }"
            "QPushButton:pressed { background: #cdd6e0; }"
            # The menu buttons (Session, Size) already show a "\u25be" in their
            # label; hide the second arrow Qt draws for any button with a menu.
            "QPushButton::menu-indicator { image: none; width: 0px; }"
        )
        
        # Session menu -- a button with a dropdown, not a click handler
        # the small downward-triangle character (▾) appended to the label
        # (a standard visual cue that this button opens a menu rather than 
        # performing an action directly)
        self._btn_session = QPushButton("Session  \u25be")
        
        # tooltip
        self._btn_session.setToolTip("New, save, load, or run a full session")
        self._btn_session.setStyleSheet(pill_style)
        
        # show the dropdown whenever this button is clicked
        session_menu = QMenu(self._btn_session)
        session_menu.setStyleSheet(
            "QMenu { background:#ffffff; color:#1d1d1f; "
            "border:1px solid #d1d1d6; border-radius:6px; padding:4px; }"
            "QMenu::item { background:#ffffff; color:#1d1d1f; "
            "padding:6px 20px; border-radius:4px; }"
            "QMenu::item:selected { background:#e8f0fe; color:#0071e3; }"
            "QMenu::separator { height:1px; background:#e5e5ea; margin:4px 8px; }"
        )
        
        # create a menu item and connect its own triggered signal to the given 
        # handler in one line — a convenience overload distinct from the more
        # verbose two-step "create a QAction, then connect its triggered signal" 
        # pattern
        session_menu.addAction("New session", self._on_new_session)
        session_menu.addSeparator()
        session_menu.addAction("Save session...", self._on_save_session)
        session_menu.addAction("Load session...", self._on_load_session)
        session_menu.addSeparator()
        session_menu.addAction("Run all steps", self._on_run_all)
        self._btn_session.setMenu(session_menu)
        
        lay.addWidget(self._btn_session)
        lay.addSpacing(6)
        
        # Change mode and About buttons
        # (1) Change mode :: switch between polycrystal and monocrystal modes
        #                    (opens the initial menu with mode selection)
       
        self._btn_mode = QPushButton("Change mode")
        self._btn_mode.setToolTip(
            "Go back to mode selection (polycrystal / monocrystal)"
            )
        self._btn_mode.setStyleSheet(pill_style)
        self._btn_mode.clicked.connect(self._on_change_mode)
        lay.addWidget(self._btn_mode)
        lay.addSpacing(6)

        # Interface size menu (also Ctrl+= / Ctrl+- / Ctrl+0). Qt can only
        # change its scale factor at startup, so a new size is applied by
        # restarting EDIFF - the menu offers to do that (see _choose_ui_scale).
        self._btn_zoom = QPushButton("Size  \u25be")
        self._btn_zoom.setToolTip(
            "Make the whole interface larger or smaller.\n"
            "EDIFF restarts to apply the new size.")
        self._btn_zoom.setStyleSheet(pill_style)
        zoom_menu = QMenu(self._btn_zoom)
        zoom_menu.setStyleSheet(
            "QMenu { background:#ffffff; color:#1d1d1f; "
            "border:1px solid #d1d1d6; border-radius:6px; padding:4px; }"
            "QMenu::item { background:#ffffff; color:#1d1d1f; "
            "padding:6px 20px; border-radius:4px; }"
            "QMenu::item:selected { background:#e8f0fe; color:#0071e3; }"
            "QMenu::item:disabled { color:#86868b; }"
            "QMenu::separator { height:1px; background:#e5e5ea; margin:4px 8px; }"
        )
        # First line: what is running now / what is saved (filled in when
        # the menu opens, see _refresh_size_menu).
        self._size_info = zoom_menu.addAction("")
        self._size_info.setEnabled(False)
        zoom_menu.addSeparator()
        # Direct choices. `checked=False` absorbs the bool that QAction's
        # triggered signal passes, so it cannot overwrite the bound `s`.
        self._size_actions = []
        for s in self.UI_SCALE_PRESETS:
            act = zoom_menu.addAction(
                f"{round(s * 100)}%" + ("  (default)" if s == DEFAULT_UI_SCALE else ""),
                lambda checked=False, s=s: self._choose_ui_scale(s))
            act.setCheckable(True)
            self._size_actions.append((s, act))
        zoom_menu.addSeparator()
        zoom_menu.addAction("Larger\tCtrl++", lambda: self._change_ui_scale(+1))
        zoom_menu.addAction("Smaller\tCtrl+-", lambda: self._change_ui_scale(-1))
        zoom_menu.addSeparator()
        self._size_restart = zoom_menu.addAction(
            "Restart EDIFF to apply...", self._restart_from_menu)
        zoom_menu.aboutToShow.connect(self._refresh_size_menu)
        self._btn_zoom.setMenu(zoom_menu)
        lay.addWidget(self._btn_zoom)
        lay.addSpacing(6)
        
        # (2) About :: displays the acknowledgment
        about_style = (
            "QPushButton {"
            "  background: #e8f0fe; color: #0071e3; border: 1px solid #cfe1fb;"
            "  border-radius: 14px; padding: 0px 14px; font-size: 12px; font-weight: 500;"
            "  min-height: 28px; max-height: 28px;"
            "}"
            "QPushButton:hover { background: #dbe9fd; border-color: #b7d3f9; }"
            "QPushButton:pressed { background: #cce0fb; }"
        )
        self._btn_about = QPushButton("About")
        self._btn_about.setToolTip("About EDIFF and its funding")
        self._btn_about.setStyleSheet(about_style)
        self._btn_about.clicked.connect(self._on_about)
        lay.addWidget(self._btn_about)
        
        # Returns the fully assembled, top bar widget, to be added into root
        # by _build_ui().
        return bar

    # ── Step builders ─────────────────────────────────────────────────────────
    def _build_step1(self):
        """
        Build Step 1: working directory, structure files, diffraction image.

        Constructs the sidebar's entry-point step block: a working-
        directory display with a "Change..." button, an empty (populated
        via `_refresh_cif_list`) list of loaded structure files with its
        own "Browse structure files..." button, and an empty (populated
        via `_refresh_eld_list`) single-slot diffraction-image display
        with its own browse button.
    
        Unlike every other step, this block's fields are explicitly
        force-enabled immediately after construction (overriding
        StepBlock.add's default "start locked" behavior), since Step 1 is
        the app's entry point and must be usable before any step-
        activation logic has run.
    
        Returns
        -------
        StepBlock
            The fully constructed Step 1 block, ready to add to the
            sidebar.
        """
    
        # Create the StepBlock container for step 1 - everything added from here
        # becomes part of its body
        b = StepBlock(1, "Input files")

        # Working directory — where all output files (corrected images,
        # center/background files, 1D profiles, etc.) get saved. Defaults
        # to a sensible location, changeable at any point.
        wd_lbl_row = QLabel("Working directory:")
        wd_lbl_row.setStyleSheet("font-size:11px; font-weight:600; color:#1d1d1f;")
        b.add(wd_lbl_row)
        
        # A read-only-looking path display, styled like a disabled/inactive
        # input field (it is meant to be changed via the button next to it)
        wd_row = QWidget(); wd_row.setStyleSheet("background:transparent;")
        wd_lay = QHBoxLayout(wd_row)
        wd_lay.setContentsMargins(0, 4, 0, 10); wd_lay.setSpacing(8)

        self._wd_path_lbl = QLabel()
        self._wd_path_lbl.setStyleSheet(
            "font-size:11px; color:#3a3a3c; background:#f5f5f7;"
            "border:1px solid #e5e5ea; border-radius:6px; padding:5px 8px;")
        wd_lay.addWidget(self._wd_path_lbl, 1) 
        # 1 is a stretch factor, which gives this label first claim on available
        # horizontal space in the row, so it grows to show as much of a long
        # path as possible, while the button beside it stays its natural size
        
        # the button that opens a folder-picker
        wd_browse_btn = QPushButton("Change...")
        wd_browse_btn.setStyleSheet(
            "QPushButton { background:#f5f5f7; color:#0071e3; border:1px solid #d1d1d6;"
            "border-radius:6px; padding:5px 12px; font-size:11px; }"
            "QPushButton:hover { background:#e8f0fe; }"
        )
        wd_browse_btn.setToolTip(
            "Choose where output files (corrected images, center/background\n"
            "files, 1D profiles, etc.) get saved.")
        wd_browse_btn.clicked.connect(self._browse_workdir)
        wd_lay.addWidget(wd_browse_btn)

        b.add(wd_row)
        
        # populate self._wd_path_lbl with the current working directory path
        # right at construction time
        self._update_workdir_label()

        # Structure files — one row per file, each with its own red remove
        # button; re-browsing adds to this list rather than replacing it.
        # an empty container widget with its own layout
        self._cif_list_w = QWidget()
        self._cif_list_w.setStyleSheet("background:transparent;")
        self._cif_list_lay = QVBoxLayout(self._cif_list_w)
        self._cif_list_lay.setContentsMargins(0, 0, 0, 0)
        self._cif_list_lay.setSpacing(4)
        b.add(self._cif_list_w)
        
        # populate the list (container) -- the actual contents of the list 
        # change dynamically throughout the app's life as files re added and
        # removed, while the container widget itself is only ever created once,
        # at construction time
        self._refresh_cif_list()
        
        # button triggering the file dialog
        btn = make_btn("Browse structure files...", primary=False,
            tip="Select structure files (CIF, POSCAR, VASP, etc.) — adds to "
                "the list below, up to 3 total  (Ctrl+O)")
        btn.clicked.connect(self._browse_cif); b.add(btn)

        # Diffraction image — only one at a time, but still removable
        # an empty container populated via its own refresh method, plus its own 
        # browse button. 
        # structurally identical to the CIF section above, but conceptually 
        # different: as the tooltip implies and the helper methods below 
        # confirm, only one ELD image can ever be selected at a time, unlike 
        # CIF files (up to 3)
        self._eld_list_w = QWidget()
        self._eld_list_w.setStyleSheet("background:transparent;")
        self._eld_list_lay = QVBoxLayout(self._eld_list_w)
        self._eld_list_lay.setContentsMargins(0, 0, 0, 0)
        self._eld_list_lay.setSpacing(4)
        b.add(self._eld_list_w)
        self._refresh_eld_list()

        btn2 = make_btn("Browse diffraction image...", primary=False,
            tip="Select experimental diffraction image (PNG, TIFF)  (Ctrl+E)")
        btn2.clicked.connect(self._browse_eld); b.add(btn2)

        # Global display settings — icut (intensity clip) and colormap,
        # shared across every step that displays the 2D diffraction
        # pattern (Preview, Ellipse correction, Center detection), so
        # there's a single source of truth instead of several
        # independent, easily-out-of-sync copies.
        from qtpy.QtWidgets import QGridLayout

        self._icut = QSpinBox()
        self._icut.setRange(1, 9999);
        self._icut.setValue(200)
        self._icut.setToolTip(TOOLTIPS["icut"])
        self._cmap = QComboBox()
        self._cmap.addItems(
            [
                "gray",
                "viridis",
                "inferno",
                "plasma",
                "hot"
            ]
        )
        self._cmap.setToolTip(TOOLTIPS["cmap"])
        self._cmap.setStyleSheet(COMBO_FIX)

        row_disp = QWidget();
        row_disp.setStyleSheet("background:transparent;")
        gl_disp = QGridLayout(row_disp)
        gl_disp.setContentsMargins(0, 8, 0, 0);
        gl_disp.setSpacing(6)
        gl_disp.addWidget(
            _glabel2(
                "icut",
                "max clip",
                TOOLTIPS["icut"]
            ),
            0,
            0
        )

        gl_disp.addWidget(self._icut, 0, 1)
        gl_disp.addWidget(
            _glabel2(
                "Colormap",
                "display",
                TOOLTIPS["cmap"]
            ),
            0,
            2
        )

        gl_disp.addWidget(self._cmap, 0, 3)
        self._icut.setMaximumWidth(100)
        self._icut.setFocusPolicy(
            Qt.FocusPolicy.StrongFocus
        )
        gl_disp.setColumnStretch(1, 1);
        gl_disp.setColumnStretch(3, 1)
        b.add(row_disp)

        # Preview button — the only way to trigger the raw-pattern
        # preview now. Re-click it after changing icut/colormap to
        # update the preview; there is no automatic preview to keep
        # in sync anymore.

        btn_preview = make_btn(
            "Preview",
            primary=False,
            tip="Show the raw diffraction pattern using the icut/colormap "
                "settings above — click again after changing them to "
                "update the preview")
        btn_preview.clicked.connect(self._on_preview_clicked)
        b.add(btn_preview)

        # Step 1 starts enabled (it's the entry point)
        # this directly overrides stepblock's own default behaviour
        for w in b._children: w.setEnabled(True)
        
        return b

    def _on_preview_clicked(self):
        """
        Handle the Preview button in Step 1.

        Shows the raw diffraction pattern using the current global
        icut/colormap settings. Re-clicking after changing those
        settings re-renders the preview, since there is no automatic
        preview to keep in sync anymore.

        Returns
        -------
        None
            Shows the preview, or a status-bar message if no
            diffraction image has been loaded yet.
        """
        if self.pipeline.eld is None:
            self.statusBar().showMessage(
                "Load a diffraction image first, then click Preview.")
            return
        self._show_raw_eld_preview(switch_tab=True)


    def _clear_layout(self, layout):
        """Remove and fully delete every widget currently in a layout.

        Merely removing a widget from a layout (`layout.takeAt`) does not
        destroy it — the widget object remains alive, still parented, unless
        explicitly reparented and scheduled for deletion. This helper does 
        both: detaches each widget from its parent (`setParent(None)`) and 
        schedules its actual destruction (`deleteLater()`, which defers 
        deletion to a safe point in Qt's event loop rather than deleting 
        synchronously). 
        
        Used by  `_refresh_cif_list`/`_refresh_eld_list` to fully wipe a file 
        list before rebuilding it from scratch, avoiding a widget leak on every
        refresh.
    
        Parameters
        ----------
        layout : QLayout
            The layout to empty.
        """
        # simply removing widgets from a layout only detaches it from that
        # layout, leaving the widget object itself still alive, still parented
        # to whatever it was before, and still potentially visible or at least
        # present in memory
        
        while layout.count():
            item = layout.takeAt(0)
            w = item.widget()
            
            if w is not None:
                # sever the parent-child relationship entirely
                w.setParent(None)
                
                # schedule the actial destruction later, because Qt objects
                # should not be deleted synchronously while their own event 
                # handlers might still be executing
                w.deleteLater()

    def _make_file_row(self, filename, full_path, on_remove):
        """Build one file-list row: a filename label plus a remove button.

        Parameters
        ----------
        filename : str
            Short filename shown as the row's main text, styled green (the
            app's "accepted/loaded" color).
            
        full_path : str
            Full path shown as the label's tooltip on hover.
            
        on_remove : callable
            Handler connected to the row's small red "✕" button's clicked
            signal. This method has no built-in knowledge of what should
            actually be removed — that logic belongs entirely to the
            caller-supplied callback.
    
        Returns
        -------
        QWidget
            The assembled row, ready to add to a file-list layout.
        """
        # display which files were successfully loaded during the first step
        row = QWidget()
        row.setStyleSheet("background:transparent;")
        rl = QHBoxLayout(row)
        rl.setContentsMargins(0, 0, 0, 0)
        rl.setSpacing(6)
        lbl = QLabel(filename)
        lbl.setStyleSheet("font-size:11px; color:#34c759;")
        lbl.setWordWrap(True)
        lbl.setToolTip(full_path)
        rl.addWidget(lbl, stretch=1)
        
        # a small solid red "x" button to remove selected files
        remove_btn = QPushButton("✕")
        remove_btn.setFixedSize(18, 18)
        remove_btn.setToolTip(f"Remove {filename}")
        remove_btn.setStyleSheet(
            "QPushButton {"
            "  background:#ff3b30; color:#ffffff; border:none;"
            "  border-radius:9px; font-size:10px; font-weight:700; padding:0px;"
            "}"
            "QPushButton:hover { background:#ff453a; }"
            "QPushButton:pressed { background:#d70015; }"
        )
        
        # callback passed in by whichever caller builds this row
        remove_btn.clicked.connect(on_remove)
        rl.addWidget(remove_btn)
        return row


    def _refresh_cif_list(self):
        """Rebuild the structure-files list to match self.cif_paths.

        Fully clears the list container (via `_clear_layout`) and rebuilds
        it from scratch: a muted placeholder if `self.cif_paths` is empty,
        or one row per path otherwise, each wired to remove that specific
        file via `_remove_cif`.
        """
        # start by fully clearing the container
        self._clear_layout(self._cif_list_lay)
        
        # if there are no CIF paths, shows a muted gray placeholder instead 
        # of an empty, blank-looking space
        if not self.cif_paths:
            placeholder = QLabel("No structure files selected")
            placeholder.setStyleSheet("font-size:11px; color:#aeaeb2;")
            self._cif_list_lay.addWidget(placeholder)
            return
        
        # build one row per file
        for p in self.cif_paths:
            row = self._make_file_row(
                Path(p).name, p,
                on_remove=lambda checked=False, path=p: self._remove_cif(path))
            self._cif_list_lay.addWidget(row)


    def _refresh_eld_list(self):
        """Rebuild the diffraction-image display to match self.eld_path.

        Fully clears the container and rebuilds it: a muted placeholder if
        `self.eld_path` is None, or a single row for the loaded image
        otherwise, wired to remove it via `_remove_eld`.
        """
        # start by fully clearing the container
        self._clear_layout(self._eld_list_lay)
        
        # if there are no images, shows a muted gray placeholder instead 
        # of an empty, blank-looking space
        if not self.eld_path:
            placeholder = QLabel("No diffraction image selected")
            placeholder.setStyleSheet("font-size:11px; color:#aeaeb2;")
            self._eld_list_lay.addWidget(placeholder)
            return
        
        # build the row
        row = self._make_file_row(
            Path(self.eld_path).name, self.eld_path,
            on_remove=lambda checked=False: self._remove_eld())
        self._eld_list_lay.addWidget(row)


    def _remove_cif(self, path):
        """Remove one structure file from the loaded list.

        Parameters
        ----------
        path : str
            Path of the CIF file to remove. If not currently in
            `self.cif_paths`, this is a no-op.
    
        Notes
        -----
        Distinguishes two different consequences of removal: if this was
        the last remaining structure file (count goes from >=1 to 0), the
        pipeline's overall mode changes (e.g. full pipeline -> ELD-only),
        so `_check_files_ready()` re-evaluates the whole step-unlocking
        state machine from scratch. If other structure files remain, the
        mode is unaffected — only a status message is shown, since
        whatever step the user is currently on shouldn't be disturbed by
        trimming the candidate list.
        """
        # if the path is not within the list, nothing happens
        if path not in self.cif_paths:
            return
        
        # check whether there is at least 1 cif file
        was_has_cif = bool(self.cif_paths)
        
        # remove a selected cif file from the input path
        self.cif_paths.remove(path)
        
        # populate the containers
        self._refresh_cif_list()
        
        # check whether there is still at least 1 cif file
        now_has_cif = bool(self.cif_paths)
        
        if was_has_cif and not now_has_cif:
            # Dropped from "has structure files" to "none" — the pipeline
            # mode itself changes (full pipeline -> ELD-only, or -> nothing),
            # so let the state machine re-evaluate from scratch.
            self._check_files_ready()
        else:
            # Still have at least one — just trimmed the candidate list.
            # Don't disturb whatever step the user is currently on.
            self.statusBar().showMessage(
                f"Removed {Path(path).name}. Recalculate XRD to update the phase list.")

    def _remove_eld(self):
        """Remove the currently loaded diffraction image, if any.

        Resets `self.eld_path` to None, restores the window title to its
        plain default, refreshes the display, and always re-runs
        `_check_files_ready()` — since only one ELD image can ever be
        loaded at a time, removing it is always a transition from "have
        one" to "have none", unlike the CIF case where some files can
        remain after a removal.
        """
        if not self.eld_path:
            return
        
        # since there's only ever at most one ELD image, removing it is always 
        # a "went from having one to having none" transition — there's no 
        # analogous "still have some left, just trimmed the list" middle case 
        # the way there is for CIF files.
        self.eld_path = None
        self.setWindowTitle("EDIFF")
        
        # populate the containers
        self._refresh_eld_list()
        
        # Only one ELD image can ever be active, so removing it always
        # changes the pipeline mode — always re-evaluate from scratch.
        self._check_files_ready()

    def _build_step2(self):
        """Build Step 2: XRD calculation parameters.

        Constructs five spinboxes (wavelength, 2θ min/max, peak sigma,
        temperature factor) laid out two-per-row in a grid — matching the
        same side-by-side pattern used elsewhere for related short-value
        fields — plus a full-width "Peak profile" dropdown (pseudo_voigt /
        gaussian / lorentzian, all three of which `ediff.pcryst` actually
        supports), and the Calculate/Stop button pair wired to `_on_xrd`
        and `_stop_step` respectively.
    
        Returns
        -------
        StepBlock
            The fully constructed Step 2 block, ready to add to the
            sidebar.
    
        Notes
        -----
        The peak profile combobox's popup styling is written out inline
        here, duplicating the module-level `COMBO_FIX` constant's exact
        rules rather than referencing it directly — likely a holdover from
        before `COMBO_FIX` was extracted as a shared constant, since this
        combobox was the one originally patched individually. Functionally
        equivalent either way, but worth being aware of if this styling
        ever needs to change, since it currently exists in two separate
        places.
        """
        b = StepBlock(2, "XRD parameters")

        # Spinboxes
        # (1) Wavelength
        self._xrd_wl = QDoubleSpinBox(); self._xrd_wl.setDecimals(4)
        self._xrd_wl.setRange(0.001,10); self._xrd_wl.setValue(0.71)
        self._xrd_wl.setToolTip(TOOLTIPS["wavelength"])
        
        # (2) Angle
        self._xrd_ttmin = QSpinBox()
        self._xrd_ttmin.setRange(1,90)
        self._xrd_ttmin.setValue(2)
        self._xrd_ttmin.setToolTip(TOOLTIPS["ttmin"])
        
        self._xrd_ttmax = QSpinBox()
        self._xrd_ttmax.setRange(10,180)
        self._xrd_ttmax.setValue(80)
        self._xrd_ttmax.setToolTip(TOOLTIPS["ttmax"])
        
        # (3) Peak sigma (the Gaussian peak width parameter)
        self._xrd_sigma = QDoubleSpinBox(); self._xrd_sigma.setDecimals(4)
        self._xrd_sigma.setRange(0.001,5); self._xrd_sigma.setValue(0.05)
        self._xrd_sigma.setToolTip(TOOLTIPS["sigma"])
        
        # (4) Temperature factor
        self._xrd_temp  = QDoubleSpinBox(); self._xrd_temp.setDecimals(2)
        self._xrd_temp.setRange(0.01,5.0); self._xrd_temp.setValue(0.8)
        self._xrd_temp.setToolTip(
            "Debye-Waller temperature factor (default 0.8).\n"
            "Higher values broaden and reduce peak intensities.")
        
        # The peak profile dropdown
        # TODO: the same rule set as COMBO_FIX, but note it's written out again 
        # inline here rather than actually using COMBO_FIX itself. Functionally 
        # it's identical, but it's a small maintenance inconsistency: if this 
        # exact styling ever needs to change in the future, someone would have 
        # to remember to update it in two places rather than one. 
        self._xrd_profile = QComboBox()
        self._xrd_profile.setStyleSheet(
            "QComboBox { background:#ffffff; color:#1d1d1f; }"
            "QComboBox QAbstractItemView { background:#ffffff; color:#1d1d1f; }"
            "QComboBox QAbstractItemView::item { color:#1d1d1f; background:#ffffff; }"
            "QComboBox QAbstractItemView::item:selected { background:#e8f0fe; color:#0071e3; }"
        )
        
        # three peak profile shapes within a dropdown menu
        self._xrd_profile.addItems(["pseudo_voigt", "gaussian", "lorentzian"])
        self._xrd_profile.setToolTip(
            "Shape of calculated diffraction peaks.\n"
            "pseudo_voigt is the standard choice.")

        # Two-column grid layout
        grid_w = QWidget(); grid_w.setStyleSheet("background:transparent;")
        
        # QGridLayout imported locally right here (rather than at the top of 
        # the file, alongside the other qtpy imports) — this is the only place 
        # it's needed, so it's imported right where it's used rather than 
        # cluttering the module-level import list for a widget type used in only 
        # a handful of methods.
        from qtpy.QtWidgets import QGridLayout
        
        gl = QGridLayout(grid_w)
        gl.setContentsMargins(0,0,0,0)
        gl.setSpacing(6)

        # four-row-by-four-column grid
        # (0) wavelengths
        gl.addWidget(_glabel2("Wavelength","[Å]",TOOLTIPS["wavelength"]), 0, 0)
        gl.addWidget(self._xrd_wl, 0, 1)
        
        # (1) 2θ min/max (sharing 1 row)
        gl.addWidget(_glabel2("2θ min","[deg]",TOOLTIPS["ttmin"]), 1, 0)
        gl.addWidget(self._xrd_ttmin,  1, 1)
        gl.addWidget(_glabel2("2θ max","[deg]",TOOLTIPS["ttmax"]),1, 2)
        gl.addWidget(self._xrd_ttmax,1, 3)
        
        # (2) Peak sigma/Temp. factor (sharing 1 row)
        gl.addWidget(_glabel2("Peak sigma","[deg]",TOOLTIPS["sigma"]),2, 0)
        gl.addWidget(self._xrd_sigma, 2, 1)
        gl.addWidget(_glabel2("Temp. factor","[B]"), 2, 2)
        gl.addWidget(self._xrd_temp, 2, 3)
        
        # bug fix to allow typing into spinboxes
        for _sb in (self._xrd_wl, self._xrd_ttmin, self._xrd_ttmax,
                    self._xrd_sigma, self._xrd_temp):
            _sb.setMaximumWidth(130)
            _sb.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        
        # (3) Peak profile dropdown menu
        gl.addWidget(_glabel2("Peak profile",""),  3, 0)
        gl.addWidget(self._xrd_profile,  3, 1, 1, 3)
        
        # horizontal stretch
        gl.setColumnStretch(1, 1); gl.setColumnStretch(3, 1)
        
        # the assembled grid widget is added as a single unit to the StepBlock's 
        # body
        b.add(grid_w)
        
        # Buttons
        # calculate XRD
        self._btn_xrd = make_btn("Calculate XRD",
            tip="Calculate theoretical diffractograms from structure files")
        self._btn_xrd.clicked.connect(self._on_xrd); b.add(self._btn_xrd)
        
        # stop the calculation 
        self._btn_stop2 = make_btn("Stop", obj_name="stop", 
            tip="Stop calculation and edit parameters"
            )
        
        # buttons shown only if step 2 is running
        self._btn_stop2.setVisible(False) 
        self._btn_stop2.clicked.connect(
            lambda: self._stop_step(
                2, self._btn_xrd, self._btn_stop2, self._card_xrd
                )
            )
        b.add(self._btn_stop2)
        
        return b

    def _build_step4(self):
        """Build Step 4: center detection method and its many sub-options.

        Constructs the detection method dropdown (driving everything else
        in this step via `_update_center_fields`), the refinement method,
        icut/colormap, a center-file field that switches between a text
        input (saving a new file) and a Browse button (reading an
        existing one), a "Method-specific options" section covering ring
        radius method, histogram equalization, Sobel filter, beam-stop
        masking, central square size, and minimum intensity threshold,
        plus a detail-view size field, the Find center/Stop buttons, and
        two tutorial buttons (one for manual detection, one for manual
        refinement — independently visible based on which is selected).
    
        Wires both the detection and refinement dropdowns to
        `_update_center_fields`, and calls it once immediately so the
        correct fields are shown from the start rather than everything
        being visible until the user first changes a dropdown.
    
        Returns
        -------
        StepBlock
            The fully constructed Step 4 block, ready to add to the
            sidebar.
        """
    
        b = StepBlock(4, "Center detection")
        from qtpy.QtWidgets import QGridLayout

        def glabel(txt, tip=None):
            l = QLabel(txt)
            l.setStyleSheet("font-size:11px; color:#1d1d1f; font-weight:600;")
            if tip: l.setToolTip(tip)
            return l

        # Detection method — full width
        # populated directly from DETECTION_OPTIONS's keys 
        self._ctr_det = QComboBox()
        self._ctr_det.addItems(DETECTION_OPTIONS.keys())
        self._ctr_det.setToolTip(TOOLTIPS["detection"])
        self._ctr_det.setStyleSheet(COMBO_FIX)
        b.add(field("Detection method", self._ctr_det, "detection"))

        # Refinement — full width (hidden for manual/from file)
        # stored as an instance attribute specifically to be hidden later
        self._ctr_ref = QComboBox()
        self._ctr_ref.addItems(REFINEMENT_OPTIONS.keys())
        self._ctr_ref.setToolTip(TOOLTIPS["refinement"])
        self._ctr_ref.setStyleSheet(COMBO_FIX)
        self._f_ref = field("Refinement", self._ctr_ref)
        b.add(self._f_ref)

        # Center file — full width, contents change depending on mode:
        # a Browse button when READING a previously-saved file (so the
        # user picks an actual existing file rather than typing a
        # filename and hoping it matches something in the workdir), or
        # the plain text field when SAVING a new one (where there's
        # nothing to browse to yet — the file doesn't exist until this
        # step creates it).
        self._ctr_read_file_path = None  # absolute path once browsed

        self._ctr_file = QLineEdit("center.txt")
        self._ctr_file.setToolTip(TOOLTIPS["cfile"])
        
        self._f_file_lbl = QLabel("Save center file")
        self._f_file_lbl.setStyleSheet(
            "font-size:11px; color:#1d1d1f; font-weight:600;"
            )
        
        self._f_file_w = QWidget()
        self._f_file_w.setStyleSheet("background:transparent;")
        
        _fv = QVBoxLayout(self._f_file_w)
        _fv.setContentsMargins(0,0,0,0); _fv.setSpacing(3)
        _fv.addWidget(self._f_file_lbl); _fv.addWidget(self._ctr_file)
        b.add(self._f_file_w)

        self._ctr_read_file_lbl = QLabel("No file selected")
        self._ctr_read_file_lbl.setStyleSheet("color:#86868b; font-size:12px;")
        self._btn_browse_center_file = make_btn(
            "Browse center file...", 
            primary=False
            )
        self._btn_browse_center_file.clicked.connect(self._browse_center_file)
        self._f_read_file_w = QWidget()
        self._f_read_file_w.setStyleSheet("background:transparent;")
        _frv = QVBoxLayout(self._f_read_file_w)
        _frv.setContentsMargins(0,0,0,0); _frv.setSpacing(6)
        _frv.addWidget(self._ctr_read_file_lbl)
        _frv.addWidget(self._btn_browse_center_file)
        self._f_read_file_w.setVisible(False)
        b.add(self._f_read_file_w)

        # The method-specific section separator
        self._ctr_extra_lbl = QLabel("Method-specific options:")
        self._ctr_extra_lbl.setStyleSheet(
            "font-size:10px; color:#86868b; font-weight:600; letter-spacing:0.05em;")
        b.add(self._ctr_extra_lbl)

        # Ring radius method — full width
        self._ctr_rtype = QComboBox()
        self._ctr_rtype.addItems([
            "1 - radial distribution (recommended)",
            "2 - masked radial distribution (for beam stop)",
            "0 - peak matching"])
        self._ctr_rtype.setToolTip(
            "Method to determine radius of diffraction ring for refinement.\n"
            "'Masked radial distribution' excludes the beam-stop shadow "
            "from the calculation — prefer it over the plain radial "
            "distribution method whenever a beam stop is visible in the "
            "pattern, since the beam stop can otherwise make the wrong "
            "ring appear brightest.")
        self._ctr_rtype.setStyleSheet(COMBO_FIX)
        self._f_rtype = field("Ring radius method", self._ctr_rtype)
        b.add(self._f_rtype)

        # Histogram eq + Sobel — side by side (bool dropdowns)
        self._ctr_heq = QComboBox()
        self._ctr_heq.addItems(["False","True"])
        self._ctr_heq.setToolTip("Apply histogram equalization before detection.")
        self._ctr_heq.setStyleSheet(COMBO_FIX)
        self._ctr_sobel = QComboBox()
        self._ctr_sobel.addItems(["False","True"])
        self._ctr_sobel.setToolTip("Apply Sobel edge filter before detection.")
        self._ctr_sobel.setStyleSheet(COMBO_FIX)

        row_hs = QWidget(); row_hs.setStyleSheet("background:transparent;")
        gl_hs = QGridLayout(row_hs)
        gl_hs.setContentsMargins(0,0,0,0)
        gl_hs.setSpacing(6)
        
        gl_hs.addWidget(glabel("Histogram eq."),  0, 0)
        gl_hs.addWidget(self._ctr_heq,            0, 1)
        gl_hs.addWidget(glabel("Sobel filter"),   0, 2)
        gl_hs.addWidget(self._ctr_sobel,          0, 3)
        gl_hs.setColumnStretch(1, 1); gl_hs.setColumnStretch(3, 1)
        self._f_heq_row = row_hs
        b.add(row_hs)

        # Beam-stop masking — full width, defaults to True (matching
        # CenterLocator's own default) since it's usually helpful
        self._ctr_masking = QComboBox()
        self._ctr_masking.addItems(["True", "False"])
        self._ctr_masking.setToolTip(
            "Automatically mask out the beam-stop shadow before detection "
            "so it doesn't interfere with finding the ring center. "
            "Disable only if it's masking out real diffraction data.")
        self._ctr_masking.setStyleSheet(COMBO_FIX)
        self._f_masking = field("Beam-stop masking", self._ctr_masking)
        b.add(self._f_masking)

        # Normalize background — relevant for phase/ccorr methods only.
        # Divides by a Gaussian-smoothed version of the image before
        # correlation. Usually helps, but can backfire on images with
        # an unusually large/asymmetric beam-stop shadow, since the
        # smoothing then picks up the shadow's edge as a strong,
        # spurious feature that can dominate the correlation. Disable
        # it if phase/ccorr detection gives an implausible result on
        # such an image.
        self._ctr_normalize_bg = QComboBox()
        self._ctr_normalize_bg.addItems(["True", "False"])
        self._ctr_normalize_bg.setToolTip(
            "Normalize background before phase correlation. Usually "
            "helps, but can mislead detection on images with a large, "
            "asymmetric beam-stop shadow -- try disabling if the "
            "detected center looks clearly wrong.")
        self._ctr_normalize_bg.setStyleSheet(COMBO_FIX)
        self._f_normalize_bg = field(
            "Normalize background", self._ctr_normalize_bg)
        b.add(self._f_normalize_bg)

        # Live plot — relevant for Hough detection only. Shows a
        # real-time visualization window while the Hough transform runs.
        self._ctr_live_plot = QComboBox()
        self._ctr_live_plot.addItems(["False", "True"])
        self._ctr_live_plot.setToolTip(
            "Show a real-time visualization window while Hough "
            "detection runs.")
        self._ctr_live_plot.setStyleSheet(COMBO_FIX)
        self._f_live_plot = field("Live plot", self._ctr_live_plot)
        b.add(self._f_live_plot)

        # Central square + min intensity — side by side
        self._ctr_csq = QSpinBox()
        self._ctr_csq.setRange(10, 2000)
        self._ctr_csq.setValue(50)
        self._ctr_csq.setToolTip(TOOLTIPS["csq"])
        self._ctr_cintensity = QDoubleSpinBox()
        self._ctr_cintensity.setDecimals(2)
        self._ctr_cintensity.setRange(0.0, 1.0)
        self._ctr_cintensity.setValue(0.8)
        self._ctr_cintensity.setSingleStep(0.05)
        self._ctr_cintensity.setToolTip(
            "Intensity threshold (0-1): pixels below this fraction of max intensity are ignored.")

        row_ci = QWidget(); row_ci.setStyleSheet("background:transparent;")
        gl_ci = QGridLayout(row_ci)
        gl_ci.setContentsMargins(0,0,0,0)
        gl_ci.setSpacing(6)
        gl_ci.addWidget(glabel("Central sq. (px)", TOOLTIPS["csq"]), 0, 0)
        gl_ci.addWidget(self._ctr_csq, 0, 1)
        gl_ci.addWidget(glabel("Min. intensity"), 0, 2)
        gl_ci.addWidget(self._ctr_cintensity, 0, 3)
        
        self._ctr_csq.setMaximumWidth(100)
        self._ctr_csq.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self._ctr_cintensity.setMaximumWidth(100)
        self._ctr_cintensity.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        gl_ci.setColumnStretch(1, 1); gl_ci.setColumnStretch(3, 1)
        self._f_csq_row = row_ci
        b.add(row_ci)

        # Detail view — full width
        # parameter feeding pipeline.plot_center(csquare=csq, ...), referenced
        # when the axis-limit carry-over feature was built for the center step
        self._ctr_detail_csq = QSpinBox()
        self._ctr_detail_csq.setRange(50, 2000)
        self._ctr_detail_csq.setValue(300)
        self._ctr_detail_csq.setToolTip(
            "Size of zoomed center detail view in pixels."
            )
        self._f_detail_csq = field("Detail view (px)", self._ctr_detail_csq)
        b.add(self._f_detail_csq)
        
        # Buttons - main action, stop, and two tutorial buttons
        # (a) run center detection
        self._btn_center = make_btn("Find center",
            tip="Detect center of 2D diffraction pattern")
        self._btn_center.clicked.connect(self._on_center)
        b.add(self._btn_center)
        
        # (b) stop the center detction
        self._btn_stop4 = make_btn("Stop", obj_name="stop")
        self._btn_stop4.setVisible(False)
        self._btn_stop4.clicked.connect(
            lambda: self._stop_step(
                4, self._btn_center, self._btn_stop4, self._card_eld2d
                )
            )
        b.add(self._btn_stop4)
        
        # (c) tutorial buttons - separate popup window for detection/refinement
        self._btn_tutorial_center_det = make_tutorial_btn("manual_detection")
        self._btn_tutorial_center_ref = make_tutorial_btn("manual_refinement")
        b.add(self._btn_tutorial_center_det)
        b.add(self._btn_tutorial_center_ref)
        
        # Wiring the reactive update
        # Both the detection and refinement dropdowns trigger the same update 
        # function -- the manual-refinement tutorial button's visibility 
        # depends on the refinement choice specifically, not just detection. 
        # The final explicit call, right before return b, runs the visibility 
        # logic once immediately at construction time — without this, the block 
        # would show every field simultaneously until the user first touched 
        # either dropdown, since nothing would have triggered
        # _update_center_fields yet.
        self._ctr_det.currentTextChanged.connect(self._update_center_fields)
        self._ctr_ref.currentTextChanged.connect(self._update_center_fields)
        self._update_center_fields(self._ctr_det.currentText())
        return b


    def _browse_center_file(self):
        """
        Handle "Browse center file...": pick an existing center-coordinates file.

        Opens a file dialog (defaulting to whichever folder a center file
        was last browsed from, remembered via the persisted config) for
        reading a previously-saved center-coordinates text file. On
        selection, stores the full path, updates the read-file display
        label (styled green, matching the app's "successfully loaded"
        convention), and remembers the containing folder for next time.
        Does nothing if the dialog is cancelled.
        """
    
        last = self._cfg.get("last_center_dir", "")
        
        # obtain selected data path 
        path, _ = QFileDialog.getOpenFileName(
            self, "Select saved center file", last,
            "Text files (*.txt);;All files (*.*)")
        
        # handles the user cancelling the dialog (an empty string is returned 
        # in that case)
        if not path:
            return
        
        # on success: stores the full path, updates the display label's text
        self._ctr_read_file_path = path
        self._ctr_read_file_lbl.setText(Path(path).name)
        self._ctr_read_file_lbl.setStyleSheet(
            "color:#1d7a3c; font-size:12px; font-weight:600;")
        
        # remember the containing folder as "last_center_dir" in the persisted 
        # config, so the next time this dialog opens, it starts in the same 
        # folder rather than some default location.
        self._cfg["last_center_dir"] = str(Path(path).parent)
        save_config(self._cfg)

    def _update_center_fields(self, _=None):
        """Show or hide Step 4's fields based on the current method selection.

        Reads the current Detection method and Refinement dropdown values, 
        classifies them into five booleans (is_file, is_manual, is_auto, 
        is_intensity, is_transform, plus ref_is_manual for the refinement 
        dropdown independently), and toggles the visibility of every conditional 
        field in the step accordingly — refinement,icut/colormap, 
        the save-file/read-file field pair, the "Method-specific options" 
        section header, ring radius method, beam-stop masking, histogram 
        equalization/Sobel filter, central square/minimum intensity, and the 
        two manual-method tutorial buttons.
    
        Connected to both the Detection method and Refinement dropdowns'
        currentTextChanged signals in `_build_step4`, and also called once
        directly at construction time so the correct fields are shown
        immediately rather than everything being visible until the user
        first touches a dropdown.
    
        Parameters
        ----------
        _ : str or None, optional, default=None
            Unused. Accepts and discards the new text value Qt's
            currentTextChanged signal passes, since this method always
            re-reads both dropdowns' current state directly rather than
            relying on whichever single value triggered the call.
    
        Notes
        -----
        `is_auto` is computed as "not is_file and not is_manual" rather
        than listing every automatic detection method explicitly, so it
        stays correct without modification if a new automatic method is
        ever added to DETECTION_OPTIONS.
    
        "Detail view (px)" (`_f_detail_csq`) is the one field always shown
        regardless of method — unlike every other line in this method, its
        visibility is set unconditionally to True, not tied to any of the
        computed booleans, since the zoomed detail view is useful
        regardless of how the center was found.
        """
        
        # Reads the current friendly-label text from both dropdowns and 
        # translates each through its options dictionary back to the short 
        # internal code
        det_text     = self._ctr_det.currentText()
        ref_text     = self._ctr_ref.currentText()
        detection    = DETECTION_OPTIONS.get(det_text)
        refinement   = REFINEMENT_OPTIONS.get(ref_text)
        
        # Five boolean classifications computed once, up front, and reused 
        # throughout the rest of the method rather than being re-derived or 
        # checked inline repeatedly.
        is_file      = detection is None
        is_manual    = detection == "manual"
        is_auto      = not is_file and not is_manual
        is_intensity = detection in ("intensity", "curvefit")
        is_transform = detection in ("phase", "ccorr", "hough")
        is_phase_or_ccorr = detection in ("phase", "ccorr")
        is_hough = detection == "hough"
        ref_is_manual = refinement == "manual"

        # Common fields
        # shown or hidden based only on the coarse file/auto/manual split, not 
        # on which specific automatic method is chosen. Refinement only makes 
        # sense for automatic detection (a manually-clicked point doesn't need 
        # a separate refinement algorithm, and reading from a file has no 
        # "detection" step to refine at all).
        self._f_ref.setVisible(is_auto)
        # self._f_icut_row.setVisible(not is_file)
        self._f_file_w.setVisible(not is_file)
        self._f_read_file_w.setVisible(is_file)

        # Method-specific section header
        # shown whenever either automatic or manual detection is selected 
        # — i.e., whenever at least one of the fields below it might actually 
        # be visible — but correctly hidden when reading from a file, since 
        # none of ring radius method, masking, histogram eq., Sobel, or central 
        # square/intensity apply to that mode at all; nothing would be below 
        # the header in that case.
        has_extras = is_auto or is_manual
        self._ctr_extra_lbl.setVisible(has_extras)

        # Method-specific fields
        # each keyed to exactly the classification that determines its 
        # relevance: 
        #    >> ring radius and beam-stop masking apply to any automatic method 
        #       broadly
        #    >> histogram eq./Sobel apply only to the three transform-based 
        #       methods specifically
        #    >> central square/minimum intensity apply to the two 
        #       intensity-centroid-style automatic methods and to manual 
        #       detection (since a manually-clicked point still benefits from 
        #       a central-square-constrained refinement search
        self._f_rtype.setVisible(is_auto)
        self._f_masking.setVisible(is_auto)
        self._f_heq_row.setVisible(is_transform)
        self._f_csq_row.setVisible(is_intensity or is_manual)
        self._f_detail_csq.setVisible(True)
        self._f_normalize_bg.setVisible(is_phase_or_ccorr)
        self._f_live_plot.setVisible(is_hough)

        # Tutorial buttons — show only when the relevant interactive method 
        # is selected
        self._btn_tutorial_center_det.setVisible(is_manual)
        self._btn_tutorial_center_ref.setVisible(ref_is_manual)


    def _build_step3(self):
        """Build Step 3: optional ellipse (elliptical distortion) correction.

        Constructs the icut/colormap row, the "Correct ellipse" button
        (which opens an interactive 4-point-click window), its tutorial
        button, a Stop button, and three outcome buttons that only appear
        once a correction has actually been attempted: "Use corrected
        image" (accept), "Show difference" (before/after comparison), and
        "Skip" (always visible, since skipping is valid at any time).
    
        Unlike most other steps, this one has no single linear "run then
        accept" flow — see `_on_use_ellipse` and `_on_skip_ellipse` for the
        branching logic handling both outcomes, including the case where a
        center was already detected in an earlier pass and step 4 can be
        skipped entirely.
    
        Returns
        -------
        StepBlock
            The fully constructed Step 3 block, ready to add to the
            sidebar.
        """
        b = StepBlock(3, "Ellipse correction  (optional)")

        from qtpy.QtWidgets import QGridLayout

        def glabel(txt, tip=None):
            l = QLabel(txt)
            l.setStyleSheet("font-size:11px; color:#1d1d1f; font-weight:600;")
            if tip: l.setToolTip(tip)
            return l

        # Buttons
        # the main action button
        self._btn_ellipse = make_btn(
            "Correct ellipse",
            tip="Opens a window — click 4 points on an elliptical ring to correct distortion")
        self._btn_ellipse.clicked.connect(self._on_ellipse)
        b.add(self._btn_ellipse)
        
        # tutorial button
        self._btn_tutorial_ellipse = make_tutorial_btn("ellipse_correction")
        b.add(self._btn_tutorial_ellipse)
        
        # stop button
        self._btn_stop3 = make_btn("Stop", obj_name="stop")
        self._btn_stop3.setVisible(False)
        self._btn_stop3.clicked.connect(lambda: self._stop_step_3b())
        b.add(self._btn_stop3)
        
        # accept correction and use the correcteed image
        self._btn_use_ellipse = make_btn(
            "✓  Use corrected image", primary=True,
            tip="Proceed to center detection using the ellipse-corrected pattern")
        self._btn_use_ellipse.setVisible(False)
        self._btn_use_ellipse.clicked.connect(self._on_use_ellipse)
        b.add(self._btn_use_ellipse)
        
        # show difference (before/after comparison)
        self._btn_show_diff = make_btn(
            "Show difference", primary=False,
            tip="Show |before - after| to see exactly how much the correction changed")
        self._btn_show_diff.setVisible(False)
        self._btn_show_diff.clicked.connect(self._on_show_ellipse_diff)
        b.add(self._btn_show_diff)
        
        # skip (always visible, since it is always a valid choice regardless 
        # of whether a correction has been attempted yet)
        self._btn_skip_ellipse = make_btn(
            "Skip — proceed to step 4", primary=False,
            tip="Skip ellipse correction and continue with the original pattern")
        self._btn_skip_ellipse.clicked.connect(self._on_skip_ellipse)
        b.add(self._btn_skip_ellipse)

        return b

    def _update_ellipse_fields(self, method_text=None):
        """No-op kept for compatibility.

        Retained so any existing signal connection or external reference
        expecting this method to exist doesn't break, without needing to
        track down and update every caller. Currently does nothing.
    
        Parameters
        ----------
        method_text : str or None, optional, default=None
            Unused.
        """
        pass

    def _stop_step_3b(self):
        """Stop the ellipse correction step.

        Delegates to the shared `_stop_step` (closing any open subprocess,
        figures, and tutorial dialogs) with this step's specific widgets —
        kept as its own named method, rather than an inline lambda like
        other steps' Stop buttons use, as a stable hook for any
        ellipse-specific stop behavior that might be added later.
        """
        self._stop_step(3, self._btn_ellipse, self._btn_stop3, self._card_ellipse)

    def _on_use_ellipse(self):
        """Accept the corrected image and advance past ellipse correction.

        Hides the "Use corrected image"/"Show difference" buttons and
        marks step 3 complete. If a center was already detected in an
        earlier pass (`_center_already_found`), skips step 4 entirely and
        jumps straight to step 5 — informing the auto-run queue (via
        `_advance_auto_run(4)`) that both steps 3 and 4 are now resolved,
        not just 3. Otherwise, proceeds normally to step 4.
        """
        self._btn_use_ellipse.setVisible(False)
        self._btn_show_diff.setVisible(False)
        self._complete_step(3)
        # If center already exists on the pipeline, jump straight to step 5
        if self._center_already_found():
            self.statusBar().showMessage(
                "Using corrected image. Center already detected — proceed to 1D conversion.")
            self._go_to_step(5)
            self._advance_auto_run(4)  # step 4 was skipped too, not just 3
        else:
            self.statusBar().showMessage(
                "Using corrected image — detect center to continue."
                )

    def _on_skip_ellipse(self):
        """Skip (or discard an already-applied) ellipse correction.

        Reloads the original, uncorrected diffraction image — first
        capturing any center already attached to the pipeline's current
        `eld` object, then reattaching it to the freshly-reloaded one
        afterward, so an existing detection result isn't lost purely
        because the image itself was swapped back to its uncorrected
        form.
    
        Sets step 3's state directly (rather than via `_complete_step`,
        since the branching that follows needs finer control than that
        generic path provides) and then, like `_on_use_ellipse`, checks
        whether a center was already found: if so, skips step 4 and jumps
        to step 5, otherwise proceeds to step 4 normally. Handles
        "ELD-only" mode (no structure file loaded) as a distinct case
        throughout, since `_go_to_step`'s generic progression logic
        doesn't correctly apply when steps 2 and 7 are permanently
        unavailable — block states and top-bar indicators are set directly
        via `_update_topbar_indicators_partial` in that case instead.
        """
        # Save the existing center (if any) before reloading the image
        # before reloading the original, uncorrected image (which would normally
        # wipe out any center coordinates associated with the pipeline's current,
        # possibly-corrected eld object), the existing center is captured first, 
        # then manually reattached to the freshly-reloaded object afterward.
        saved_center = None
        if hasattr(self.pipeline, 'eld') and self.pipeline.eld is not None:
            saved_center = getattr(self.pipeline.eld, 'center', None)

        # Reload original image if correction was applied
        if hasattr(self.pipeline, '_original_eld_path') and self.pipeline._original_eld_path:
            try:
                self.pipeline.load_eld(self.pipeline._original_eld_path)
            except Exception:
                pass

        # Restore center onto the freshly loaded ELD object
        if saved_center is not None and self.pipeline.eld is not None:
            self.pipeline.eld.center = saved_center
        
        
        self._btn_use_ellipse.setVisible(False)
        self._btn_show_diff.setVisible(False)
        
        # directly called - this method manages step 3's own completion state 
        # manually rather than delegating it, because the branching logic that
        # follows (deciding whether to advance to step 4 or jump to step 5, and 
        # doing so differently depending on eld_only mode) is specific enough 
        # to this particular path that it couldn't cleanly reuse whatever 
        # _complete_step does generically for other steps
        self._blocks[3].set_state(StepBlock.DONE)
        has_cif  = bool(self.cif_paths)
        has_eld  = bool(self.eld_path)
        eld_only = has_eld and not has_cif

        # If center already found, skip straight to step 5
        if self._center_already_found():
        
            if eld_only:
                # No structure file loaded, so the generic step-progression
                # logic (_go_to_step) doesn't correctly apply here — it
                # assumes the full 7-step sequence, but ELD-only mode
                # permanently excludes steps 2 (XRD calculation, needs a
                # structure file) and 7 (final comparison, needs XRD
                # results). Step 6 (Background) is included in `locked` too,
                # but only because it hasn't been reached yet in the normal
                # sequence — not because it's excluded by eld_only mode.
                # Block states and top-bar indicators are set directly here
                # instead, via the lower-level
                # _update_topbar_indicators_partial rather than _go_to_step.
                self._blocks[5].set_state(StepBlock.ACTIVE)
                self._blocks[4].set_state(StepBlock.DONE)
                self._update_topbar_indicators_partial(
                    done=[1, 3, 4], active=5, locked=[2, 6, 7])
     
            else:
                # A structure file is loaded, so the generic step-progression
                # logic applies normally — go to step 5 and proceed through
                # the rest of the workflow as usual.
                self._go_to_step(5)
        
            self.statusBar().showMessage(
                "Using original image. Center already detected — proceed to 1D conversion.")
            self._advance_auto_run(4)  # both 3 and 4 are effectively done
        
        # Center not yet found
        else:
            if eld_only:
                # Same reasoning as above: step 2 (XRD calculation) and step 7
                # (final comparison) are permanently excluded since there's no
                # structure file, while steps 5/6/7 in `locked` are simply not
                # reached yet, not structurally excluded.
                self._blocks[4].set_state(StepBlock.ACTIVE)
                self._update_topbar_indicators_partial(
                    done=[1, 3], active=4, locked=[2, 5, 6, 7])
            else:
                # A structure file is loaded — go to step 4 and proceed
                # through the rest of the workflow as usual.
                self._go_to_step(4)

            self.statusBar().showMessage("Using original image — detect center to continue.")
            self._advance_auto_run(3)

    def _center_already_found(self):
        """Check whether the pipeline currently holds a valid center result.

        Returns
        -------
        bool
            True if `self.pipeline.eld.center` exists and carries a
            non-None `x` coordinate (the signal that a real detection
            result is present, not just an empty placeholder object).
            False if any part of that chain is missing, absent, or raises
            (e.g. `self.pipeline.eld` itself not yet existing).
        """
        # checking not just that some center object exists, but that it 
        # specifically carries a non-None x coordinate, which is the actual 
        # signal that a real detection result (not just an empty placeholder 
        # object) is present.
        try:
            c = self.pipeline.eld.center
            return c is not None and hasattr(c, 'x') and c.x is not None
        except Exception:
            return False

    def _build_step5(self):
        """Build Step 5: 2D-to-1D radial averaging parameters.

        Constructs the Pixel min/max spinboxes (the primary x-range
        control for the 1D profile — see `_on_1d`, where these values take
        precedence over anything carried over via the axis-limit memory
        system), an optional "Save 1D profile" output filename field, and
        the Convert/Stop button pair wired to `_on_1d`.
    
        Returns
        -------
        StepBlock
            The fully constructed Step 5 block, ready to add to the
            sidebar.
        """
        from qtpy.QtWidgets import QGridLayout
        b = StepBlock(5, "Convert to 1D")
        
        # Spinboxes
        self._px_min = QSpinBox()
        self._px_min.setRange(0,9999); self._px_min.setValue(0)
        self._px_max = QSpinBox()
        self._px_max.setRange(0,9999); self._px_max.setValue(500)
        self._px_min.setMaximumWidth(100)
        self._px_min.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self._px_max.setMaximumWidth(100)
        self._px_max.setFocusPolicy(Qt.FocusPolicy.StrongFocus)

        def glabel(txt, tip=None):
            l = QLabel(txt)
            l.setStyleSheet("font-size:11px; color:#1d1d1f; font-weight:600;")
            if tip: l.setToolTip(tip)
            return l
        
        # Create grid
        row_px = QWidget(); row_px.setStyleSheet("background:transparent;")
        gl_px = QGridLayout(row_px)
        gl_px.setContentsMargins(0,0,0,0)
        gl_px.setSpacing(6)
        gl_px.addWidget(glabel("Pixel min", TOOLTIPS.get("px_min")), 0, 0)
        gl_px.addWidget(self._px_min,                                 0, 1)
        gl_px.addWidget(glabel("Pixel max", TOOLTIPS.get("px_max")), 0, 2)
        gl_px.addWidget(self._px_max,                                 0, 3)
        gl_px.setColumnStretch(1, 1); gl_px.setColumnStretch(3, 1)
        b.add(row_px)
        
        # Optional entry of the filename
        self._1d_outfile = QLineEdit("")
        self._1d_outfile.setPlaceholderText("optional — e.g. profile_1d.txt")
        self._1d_outfile.setToolTip(
            "If given, also saves the raw 1D profile (pixels, intensity) "
            "to a text file in the working directory. Leave blank to skip.")
        b.add(field("Save 1D profile", self._1d_outfile))
        
        # Calculate/STOP button pair
        self._btn_1d = make_btn("Convert 2D -> 1D", tip="Radially average 2D pattern to 1D profile")
        self._btn_1d.clicked.connect(self._on_1d); b.add(self._btn_1d)
        self._btn_stop5 = make_btn("Stop", obj_name="stop")
        self._btn_stop5.setVisible(False)
        self._btn_stop5.clicked.connect(lambda: self._stop_step(5, self._btn_1d, self._btn_stop5, self._card_1d))
        b.add(self._btn_stop5)
        
        return b

    def _build_step6(self):
        """Build Step 6: background subtraction method and its options.

        Constructs the Method dropdown (one of BKG_METHODS) and, below it,
        controls that appear depending on that choice:

        - InteractivePlot / RestoreFromPoints: the background-points file
          field (the only methods that read or write one), plus - for
          InteractivePlot - the tutorial button.
        - SimpleFuncs / Baselines: an Algorithm dropdown (see
          BKG_ALGORITHMS), the dedicated parameter field(s) of the chosen
          algorithm (e.g. ball radius, smoothing lam), an optional pixel
          range, and a free-form "Extra parameters" line for anything
          without a dedicated field. These methods ignore the points file,
          so its field is hidden for them rather than left looking
          meaningful.

        One widget is created per distinct parameter name (radius,
        half_window, decreasing, lam) and shared by every algorithm that
        uses it; `_update_bkg_algo_fields` shows the relevant one(s) and
        resets them to the chosen algorithm's recommended default.

        Wires the Method dropdown to `_update_bkg_fields` and the Algorithm
        dropdown to `_update_bkg_algo_fields`, and calls the former once at
        construction so the correct controls show from the start.

        Returns
        -------
        StepBlock
            The fully constructed Step 6 block, ready to add to the
            sidebar.
        """
        from qtpy.QtWidgets import QGridLayout
        b = StepBlock(6, "Background subtraction")
        self._bkg_method = QComboBox(); self._bkg_method.addItems(BKG_METHODS)
        self._bkg_method.setStyleSheet(COMBO_FIX)
        b.add(field("Method", self._bkg_method, "bkg_method"))

        # ── Algorithm (SimpleFuncs / Baselines only) ──────────────────────
        # Items are filled in by _update_bkg_fields, since they depend on
        # which of the two methods is selected.
        self._bkg_algo = QComboBox()
        self._bkg_algo.setStyleSheet(COMBO_FIX)
        self._f_bkg_algo = field("Algorithm", self._bkg_algo, "bkg_algo")
        b.add(self._f_bkg_algo)

        # ── Dedicated parameter widgets: one per distinct name ───────────
        # Shared across algorithms (e.g. `lam` serves three of them);
        # which is visible, its tooltip and its default are set per
        # algorithm in _update_bkg_algo_fields.
        labels = {"radius": "Ball radius [px]", "half_window": "Half-window",
                  "decreasing": "Decreasing window", "lam": "Smoothing (lam)"}
        kinds = {}
        for algos in BKG_ALGORITHMS.values():
            for plist in algos.values():
                for name, kind, _default, _tip in plist:
                    kinds[name] = kind
        self._bkg_param_widgets = {}
        self._bkg_param_fields  = {}
        for name, kind in kinds.items():
            if kind == "int":
                w = QSpinBox(); w.setRange(1, 100000)
            elif kind == "float":
                # min 0: rubberband treats lam=0 as plain linear interpolation
                w = QDoubleSpinBox(); w.setDecimals(2)
                w.setRange(0.0, 1e12)
            else:  # bool - a dropdown, like the other True/False options
                w = QComboBox(); w.addItems(["True", "False"])
                w.setStyleSheet(COMBO_FIX)
            self._bkg_param_widgets[name] = w
            f = field(labels.get(name, name), w)
            self._bkg_param_fields[name] = f
            b.add(f)

        # ── Pixel range (SimpleFuncs / Baselines only) ────────────────────
        # -1 is the minimum, shown as "auto": the start / end of the profile.
        self._bkg_xmin = QSpinBox(); self._bkg_xmax = QSpinBox()
        for sb in (self._bkg_xmin, self._bkg_xmax):
            sb.setRange(-1, 99999); sb.setValue(-1)
            sb.setSpecialValueText("auto")
            sb.setMaximumWidth(100)
            sb.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self._f_bkg_range = QWidget()
        self._f_bkg_range.setStyleSheet("background:transparent;")
        gl_rg = QGridLayout(self._f_bkg_range)
        gl_rg.setContentsMargins(0, 0, 0, 0); gl_rg.setSpacing(6)
        gl_rg.addWidget(_glabel2("Range start", "[px]", TOOLTIPS["bkg_range"]), 0, 0)
        gl_rg.addWidget(self._bkg_xmin,                                          0, 1)
        gl_rg.addWidget(_glabel2("Range end",   "[px]", TOOLTIPS["bkg_range"]), 0, 2)
        gl_rg.addWidget(self._bkg_xmax,                                          0, 3)
        gl_rg.setColumnStretch(1, 1); gl_rg.setColumnStretch(3, 1)
        b.add(self._f_bkg_range)

        # ── Extra, free-form algorithm parameters ─────────────────────────
        self._bkg_extra = QLineEdit()
        self._bkg_extra.setPlaceholderText("optional - e.g. size=50")
        self._f_bkg_extra = field("Extra parameters", self._bkg_extra, "bkg_extra")
        b.add(self._f_bkg_extra)

        # ── Points file (InteractivePlot / RestoreFromPoints only) ────────
        self._bkg_file = QLineEdit("background.txt.bp")
        self._f_bkg_file = field("Background file", self._bkg_file, "bkg_file")
        b.add(self._f_bkg_file)

        self._btn_bkg = make_btn("Subtract background", tip="Remove background from 1D diffractogram")
        self._btn_bkg.clicked.connect(self._on_bkg); b.add(self._btn_bkg)
        self._btn_stop6 = make_btn("Stop", obj_name="stop")
        self._btn_stop6.setVisible(False)
        self._btn_stop6.clicked.connect(lambda: self._stop_step(6, self._btn_bkg, self._btn_stop6, self._card_bkg))
        b.add(self._btn_stop6)
        self._btn_tutorial_bkg = make_tutorial_btn("interactive_bkg")
        b.add(self._btn_tutorial_bkg)
        self._bkg_method.currentTextChanged.connect(self._update_bkg_fields)
        self._bkg_algo.currentTextChanged.connect(self._update_bkg_algo_fields)
        self._update_bkg_fields(self._bkg_method.currentText())
        return b

    def _update_bkg_fields(self, method):
        """Show the Step 6 controls that apply to the selected method.

        Connected to the Method dropdown's currentTextChanged signal, and
        also called once at construction.

        For SimpleFuncs and Baselines (the methods listed in
        BKG_ALGORITHMS) this repopulates the Algorithm dropdown with that
        method's algorithms and shows the algorithm, pixel-range and
        extra-parameter controls. For InteractivePlot and RestoreFromPoints
        it hides those and shows the background-points file field instead -
        the only methods that actually use a points file. The interactive
        tutorial button is shown for InteractivePlot only.

        Parameters
        ----------
        method : str
            The currently selected background-subtraction method (one of
            BKG_METHODS).
        """
        algos = list(BKG_ALGORITHMS.get(method, {}))
        automatic = bool(algos)
        self._f_bkg_algo.setVisible(automatic)
        self._f_bkg_range.setVisible(automatic)
        self._f_bkg_extra.setVisible(automatic)
        self._f_bkg_file.setVisible(not automatic)
        self._btn_tutorial_bkg.setVisible(method == "InteractivePlot")
        # Repopulate silently, then update the parameter fields once
        # explicitly - otherwise they'd update for a half-built list.
        self._bkg_algo.blockSignals(True)
        self._bkg_algo.clear()
        self._bkg_algo.addItems(algos)
        self._bkg_algo.blockSignals(False)
        self._update_bkg_algo_fields(self._bkg_algo.currentText())

    def _update_bkg_algo_fields(self, algo):
        """Show the parameter field(s) of the chosen algorithm, at their defaults.

        Connected to the Algorithm dropdown's currentTextChanged signal,
        and called by `_update_bkg_fields` whenever the algorithm list is
        rebuilt. Shows only the parameter widgets that `algo` uses (per
        BKG_ALGORITHMS), gives each the algorithm-specific tooltip, and
        resets it to that algorithm's recommended default.

        The reset is intentional: the same widget (e.g. `lam`) serves
        several algorithms whose sensible values differ by orders of
        magnitude (5 / 1000 / 2), so carrying a value across a switch would
        be wrong. Callers restoring saved values (session loading) must
        therefore set the algorithm FIRST and the parameter values after.

        Parameters
        ----------
        algo : str
            The selected algorithm name, or "" when the current method has
            no algorithms (all parameter fields are then hidden).
        """
        method = self._bkg_method.currentText()
        spec = BKG_ALGORITHMS.get(method, {}).get(algo, [])
        active = {name for name, _k, _d, _t in spec}
        for name, f in self._bkg_param_fields.items():
            f.setVisible(name in active)
        for name, kind, default, tip in spec:
            w = self._bkg_param_widgets[name]
            w.setToolTip(tip)
            if kind == "bool":
                w.setCurrentText("True" if default else "False")
            else:
                w.setValue(default)

    def _read_bkg_params(self):
        """Read the dedicated parameter fields of the current algorithm.

        Returns
        -------
        dict
            Parameter name -> value, for the parameters the selected
            algorithm uses (e.g. {"radius": 70}, {"decreasing": True},
            {"lam": 1000.0}). Empty for an algorithm with none (TopHat) or
            a method with no algorithms. Values are typed correctly for the
            library: int for spin boxes, float for the lam box, bool for
            True/False dropdowns.
        """
        method = self._bkg_method.currentText()
        algo = self._bkg_algo.currentText()
        out = {}
        for name, kind, _default, _tip in BKG_ALGORITHMS.get(method, {}).get(algo, []):
            w = self._bkg_param_widgets[name]
            out[name] = (w.currentText() == "True") if kind == "bool" else w.value()
        return out

    def _bkg_algorithm_args(self):
        """Collect the SimpleFuncs/Baselines options, validated, for a run.

        Must be called on the main (GUI) thread, before any background
        thread starts, because it reads live widgets - and so that bad
        input is reported before the step is marked as running.

        Returns
        -------
        (algorithm, xrange, params) : (str, tuple or None, dict)
            For SimpleFuncs/Baselines: the algorithm name; the pixel range
            as (start, end) with None for an "auto" side, or None if both
            are auto; and all algorithm keyword arguments - the dedicated
            fields overlaid with the "Extra parameters" line (extras win on
            a name clash). For the point-based methods: (None, None, None),
            since they use none of this.

        Raises
        ------
        ValueError
            If the Extra parameters text is malformed (see
            `parse_extra_params`), or the range start is not smaller than
            its end.
        """
        if self._bkg_method.currentText() not in BKG_ALGORITHMS:
            return None, None, None
        params = self._read_bkg_params()
        params.update(parse_extra_params(self._bkg_extra.text()))
        lo, hi = self._bkg_xmin.value(), self._bkg_xmax.value()
        lo = None if lo < 0 else lo
        hi = None if hi < 0 else hi
        if lo is not None and hi is not None and lo >= hi:
            raise ValueError(
                f"Background range start ({lo}) must be smaller than its "
                f"end ({hi}).")
        xrange = None if (lo is None and hi is None) else (lo, hi)
        return self._bkg_algo.currentText(), xrange, params

    def _build_step7(self):
        """
        Build Step 7 of the GUI: calibration and comparison.
     
        This step lets the user define how diffraction-pattern calibration 
        should be performed, select the reference phase, optionally enter 
        microscope geometry parameters, calibrate the pixel-to-q conversion, 
        update comparison plots, and save the resulting data/figures.
     
        The calibration controls are split into several logical groups:
     
        1. Calibration method
           Selects the calibration approach. Some methods may require 
           microscope parameters such as accelerating voltage, camera length, 
           pixel size, and binning.
     
        2. Reference phase
           Selects which phase should be treated as the reference for 
           calibration.
     
        3. Microscope parameters
           Additional input fields used by calibration methods that rely 
           on nominal microscope/camera geometry.
     
        4. Fine-tuning and q-range
           Allows the user to manually adjust the calibration scale factor 
           and define the displayed q-range for comparison plots.
     
        5. Action buttons
           Provides controls for calibration, stopping a running operation, 
           updating plots, and saving results.
     
        Returns
        -------
        StepBlock
            Fully constructed Step 7 widget block containing all calibration 
            and comparison controls.
        """
        # Create the container block for this workflow step.
        b = StepBlock(7, "Calibrate & compare")
        from qtpy.QtWidgets import QGridLayout

        def glabel(txt, tip=None):
            """
            Create a small bold label used inside compact grid layouts.
    
            Parameters
            ----------
            txt : str
                Text shown in the label.
    
            tip : str, optional
                Tooltip shown when hovering over the label.
    
            Returns
            -------
            QLabel
                Styled QLabel used as a row/field label.
            """
            l = QLabel(txt)
            l.setStyleSheet("font-size:11px; color:#1d1d1f; font-weight:600;")
            
            # Tooltips are optional; only attach one if it is explicitly 
            # provided.
            if tip: l.setToolTip(tip)
            
            return l
        
        # Calibration method selection
        # Combo box listing all available calibration methods.
        self._cal_method = QComboBox()
        self._cal_method.addItems(CAL_METHODS)
        self._cal_method.setStyleSheet(COMBO_FIX)
        
        # Wrap the widget in your app-specific `field(...)` helper
        b.add(field("Calibration method", self._cal_method, "cal_method"))
        
        # Reference phase selection
        # Combo box for selecting the phase used as the calibration reference.
        self._cal_ref = QComboBox(); self._cal_ref.addItem("Phase 1")
        self._cal_ref.setStyleSheet(COMBO_FIX)
        self._cal_ref.setToolTip(TOOLTIPS["cal_ref"])
        b.add(field("Reference phase", self._cal_ref, "cal_ref"))

        # ── Microscope / MicroscopeParams extra fields ─────────────────────
        # These fields are shown only for calibration methods that need 
        # microscope geometry. 
        # Visibility is controlled later by `_update_cal_fields(...)`.
        
        # Container widget for all microscope-parameter rows.
        self._cal_mic_params = QWidget()
        self._cal_mic_params.setStyleSheet("background:transparent;")
        
        # Vertical layout because we want several compact horizontal/grid rows
        # stacked under each other.
        mic_lay = QVBoxLayout(self._cal_mic_params)
        mic_lay.setContentsMargins(0, 0, 0, 0); mic_lay.setSpacing(6)

        # Row 1: accelerating voltage 1 camera length
        # accelerating voltage in kV
        self._cal_voltage  = QDoubleSpinBox()
        self._cal_voltage.setDecimals(0); self._cal_voltage.setRange(1, 1000)
        self._cal_voltage.setValue(120)
        self._cal_voltage.setToolTip(
            "Accelerating voltage of the electron microscope [kV]"
            )
        
        # nominal camera length in mm
        self._cal_cl = QDoubleSpinBox()
        self._cal_cl.setDecimals(1); self._cal_cl.setRange(0.1, 5000)
        self._cal_cl.setValue(170)
        self._cal_cl.setToolTip(
            "Camera length — nominal value read from the microscope software [mm]"
            )
        # transparent row widget used only as a layout container.
        row_vc = QWidget(); row_vc.setStyleSheet("background:transparent;")
        
        # Grid layout allows labels and spin boxes to be aligned in columns:
        # label | value | label | value
        gl_vc = QGridLayout(row_vc)
        gl_vc.setContentsMargins(0,0,0,0)
        gl_vc.setSpacing(6)
        
        # add labels with units
        gl_vc.addWidget(_glabel2("Voltage","[kV]"), 0, 0)
        gl_vc.addWidget(self._cal_voltage, 0, 1)
        gl_vc.addWidget(_glabel2("Cam. length", "[mm]"), 0, 2)
        gl_vc.addWidget(self._cal_cl, 0, 3)
        
        # keep numeric fields compact so they do not dominate the row
        self._cal_voltage.setMaximumWidth(100)
        self._cal_cl.setMaximumWidth(100)
        
        # StrongFocus allows both tab focus and mouse focus.
        # This is useful for spin boxes because users often type values directly.
        self._cal_voltage.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self._cal_cl.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        
        # let the value columns absorb extra horizontal space
        gl_vc.setColumnStretch(1, 1); gl_vc.setColumnStretch(3, 1)
        mic_lay.addWidget(row_vc)

        # Row 2: physical pixel size + camera binning
        # Physical detector pixel size before binning, in micrometres.
        self._cal_pixsize = QDoubleSpinBox()
        self._cal_pixsize.setDecimals(2); self._cal_pixsize.setRange(0.1, 1000)
        self._cal_pixsize.setValue(13.2)
        self._cal_pixsize.setToolTip("Physical size of one camera pixel [µm]")
        
        # Camera binning factor.
        # Effective pixel size is usually pixel_size * binning.
        self._cal_binning = QSpinBox()
        self._cal_binning.setRange(1, 16); self._cal_binning.setValue(2)
        self._cal_binning.setToolTip("Camera binning factor (1 = no binning)")

        row_pb = QWidget()
        row_pb.setStyleSheet("background:transparent;")
        
        gl_pb = QGridLayout(row_pb)
        gl_pb.setContentsMargins(0,0,0,0)
        gl_pb.setSpacing(6)
        
        gl_pb.addWidget(_glabel2("Pixel size",  "[µm]"),       0, 0)
        gl_pb.addWidget(self._cal_pixsize,                     0, 1)
        gl_pb.addWidget(_glabel2("Binning",     "×"),             0, 2)
        gl_pb.addWidget(self._cal_binning,                     0, 3)
        
        self._cal_pixsize.setMaximumWidth(100)
        self._cal_binning.setMaximumWidth(100)
        
        self._cal_pixsize.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self._cal_binning.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        
        gl_pb.setColumnStretch(1, 1)
        gl_pb.setColumnStretch(3, 1)
        mic_lay.addWidget(row_pb)
        
        # Add the complete microscope-parameter section to Step 7.
        # Its visibility is managed by `_update_cal_fields(...)`.
        b.add(self._cal_mic_params)
        
        # Fine-tuning factor and q-range controls
        # Multiplicative correction factor applied after initial calibration.
        # This is useful when the nominal calibration is slightly off.
        self._cal_finetune = QDoubleSpinBox()
        self._cal_finetune.setDecimals(5)
        self._cal_finetune.setValue(1.003)
        self._cal_finetune.setRange(0.5,2.0)
        self._cal_finetune.setSingleStep(0.001)
        self._cal_finetune.setToolTip(TOOLTIPS.get("finetune", ""))
        
        # Minimum q value shown/used in comparison plots.
        self._cal_qmin = QDoubleSpinBox(); self._cal_qmin.setDecimals(2)
        self._cal_qmin.setValue(1.0); self._cal_qmin.setRange(0,20)
        self._cal_qmin.setToolTip(TOOLTIPS.get("qmin", ""))
        
        # Maximum q value shown/used in comparison plots.
        self._cal_qmax = QDoubleSpinBox(); self._cal_qmax.setDecimals(2)
        self._cal_qmax.setValue(8.0); self._cal_qmax.setRange(0,20)
        self._cal_qmax.setToolTip(TOOLTIPS.get("qmax", ""))
        
        # Compact grid row for fine-tune and q-range controls.
        row_q = QWidget()
        row_q.setStyleSheet("background:transparent;")
        
        gl_q = QGridLayout(row_q)
        gl_q.setContentsMargins(0,0,0,0)
        gl_q.setSpacing(6)
        
        # Fine-tune occupies the first row.
        gl_q.addWidget(glabel("Fine-tune"),    0, 0)
        gl_q.addWidget(self._cal_finetune,     0, 1)
        
        # q-range occupies the second row.
        gl_q.addWidget(glabel("q min [1/A]"), 1, 0)
        gl_q.addWidget(self._cal_qmin,        1, 1)
        gl_q.addWidget(glabel("q max [1/A]"), 1, 2)
        gl_q.addWidget(self._cal_qmax,        1, 3)
        
        self._cal_finetune.setMaximumWidth(100)
        self._cal_finetune.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self._cal_qmin.setMaximumWidth(100)
        self._cal_qmin.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self._cal_qmax.setMaximumWidth(100)
        self._cal_qmax.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        gl_q.setColumnStretch(1, 1); gl_q.setColumnStretch(3, 1)
        b.add(row_q)
        
        # Action buttons
        # (1) Main calibration button. Starts calibration using the selected method
        # and current parameter values.
        self._btn_cal = make_btn(
            "Calibrate", 
            tip="Calibrate pixel to q conversion"
            )
        
        self._btn_cal.clicked.connect(self._on_calibrate); b.add(self._btn_cal)
        
        # (2) Stop button for interrupting Step 7 processing.
        # Hidden by default and shown only while calibration is running.
        self._btn_stop7 = make_btn("Stop", obj_name="stop")
        self._btn_stop7.setVisible(False)
        
        # The lambda passes Step-7-specific controls to a generic stop handler.
        # `_card_comp` is presumably the result/comparison card affected by 
        # this step.
        self._btn_stop7.clicked.connect(
            lambda: self._stop_step(
                7, 
                self._btn_cal, 
                self._btn_stop7, 
                self._card_comp
                )
            )
        b.add(self._btn_stop7)
        
        # (3) Button for replotting comparison results without fully 
        # recalibrating. Initially disabled because there is nothing to update 
        # before calibration.
        self._btn_update = make_btn(
            "Update plots", 
            primary=False, 
            tip="Replot with current fine-tune / q range"
            )
        
        self._btn_update.clicked.connect(self._on_update_compare)
        self._btn_update.setEnabled(False); b.add(self._btn_update)
        
        # (4) Save button for exporting all final outputs.
        # Initially disabled until valid calibration/comparison results exist. 
        self._btn_save = make_btn(
            "Save all results  (Ctrl+S)", 
            tip="Save plots and data to a folder"
            )
        self._btn_save.clicked.connect(self._on_save)
        self._btn_save.setEnabled(False); b.add(self._btn_save)

        # Dynamic visibility of calibration fields
        # Whenever the calibration method changes, show/hide method-specific 
        # fields. For example, microscope parameters should probably be visible 
        # only for calibration methods based on microscope geometry.
        self._cal_method.currentTextChanged.connect(self._update_cal_fields)
        
        # Apply the correct visibility immediately, so the initial GUI state 
        # matches the default calibration method.
        self._update_cal_fields(self._cal_method.currentText())
        return b

    def _update_cal_fields(self, method):
        """
        Update visibility of calibration-method-specific input fields.
    
        Some calibration methods require additional microscope/camera 
        parameters,such as accelerating voltage, camera length, detector pixel 
        size, and binning. These controls are grouped 
        in ``self._cal_mic_params``.
    
        This function shows the microscope-parameter section only when the 
        selected calibration method needs those values. For all other 
        calibration methods, the extra fields are hidden to keep the Step 7 
        interface compact and avoid implying that unused parameters affect 
        the calibration.
    
        Parameters
        ----------
        method : str
            Name of the currently selected calibration method, usually taken 
            from ``self._cal_method.currentText()``.
    
        Returns
        -------
        None
            The function updates the GUI state in place.
        """
        # Only these calibration methods depend on the microscope/camera 
        # geometry fields contained in `self._cal_mic_params`.
        show = method in ("Microscope", "MicroscopeParams")
        
        # Show or hide the whole microscope-parameter widget group at once.
        # This is simpler than toggling each individual spin box separately.
        self._cal_mic_params.setVisible(show)

    # ── Step state machine ────────────────────────────────────────────────────
    # Defines the logical order of workflow steps in the EDIFF GUI.
    # This is used instead of simple numeric comparison so that the order can 
    # be changed later without rewriting all step-comparison logic.

    _STEP_ORDER = [1, 2, 3, 4, 5, 6, 7]

    def _step_lt(self, a, b):
        """
        Return whether workflow step ``a`` comes before workflow step ``b``.
    
        The comparison is based on ``self._STEP_ORDER`` rather than directly 
        using numeric comparison. This makes the state-machine logic robust to 
        future changes in the pipeline order.
    
        Parameters
        ----------
        a : int
            First step number.
    
        b : int
            Second step number.
    
        Returns
        -------
        bool
            True if step ``a`` appears before step ``b`` in ``self._STEP_ORDER``.
            False if either step is not present in the defined workflow order.
        """
        try:
            # Compare the positions of the two steps in the explicit pipeline 
            # order.
            return self._STEP_ORDER.index(a) < self._STEP_ORDER.index(b)
        
        except ValueError:
            # If either step is not in _STEP_ORDER, treat it as not ordered.
            # Returning False is safer than raising inside GUI-state updates.
            return False


    def _go_to_step(self, step):
        """
        Move the GUI workflow to a selected active step.
    
        All steps before ``step`` are marked as done, the requested step is 
        marked active, and all later steps are locked. The visual top-bar 
        indicators are updated to match the same state.
    
        Step 1 is treated specially: its input widgets remain enabled even 
        after later steps become active. This lets the user still access 
        initial loading controls without the first block becoming unusable.
    
        Parameters
        ----------
        step : int
            Step number that should become the active workflow step.
    
        Returns
        -------
        None
            Updates widget states and top-bar indicator styles in place.
        """
        # Update the state of each StepBlock according to its position relative
        # to the requested active step.
        for n, block in self._blocks.items():
            if self._step_lt(n, step):   
                block.set_state(StepBlock.DONE)
            elif n == step:              
                block.set_state(StepBlock.ACTIVE)
            else:                        
                block.set_state(StepBlock.LOCKED)
            
            
        # Step 1 always stays enabled, even when marked as done.
        # This is useful because Step 1 usually contains file/session loading
        # controls that the user may need again later.
        for w in self._blocks[1]._children:
            w.setEnabled(True)
            
        # Update the top progress indicators using a simple three-state visual
        # scheme: done, active, locked.
        for i, (badge, lbl) in enumerate(self._step_indicators):
            n = i + 1
            
            if n < step:
                # Done state:
                #   > green filled badge with a checkmark and green label.
                badge.setText("✓")
                badge.setStyleSheet(
                    "background:#34c759; color:#ffffff; border-radius:10px;"
                    "font-size:10px; font-weight:700;"
                    )
                lbl.setStyleSheet(
                    "font-size:12px; color:#34c759; font-weight:600;"
                    )
                
            elif n == step:
                # Active state:
                #   > blue filled badge showing the current step number.
                badge.setText(str(n))
                badge.setStyleSheet(
                    "background:#0071e3; color:#ffffff; border-radius:10px;"
                    "font-size:10px; font-weight:700;"
                    )
                lbl.setStyleSheet(
                    "font-size:12px; color:#0071e3; font-weight:600;"
                    )
                
            else:
                # Locked state:
                #   > pale grey outlined badge and muted label.
                badge.setText(str(n))
                badge.setStyleSheet(
                    "background:#f5f5f7; color:#c7c7cc; border-radius:10px;"
                    "border: 1.5px solid #e5e5ea;"
                    "font-size:10px; font-weight:700;"
                    )
                lbl.setStyleSheet(
                    "font-size:12px; color:#c7c7cc; font-weight:400;"
                    )

    def _complete_step(self, step):
        """
        Mark a workflow step as complete and advance the GUI state.
    
        This is the central state-machine function called after a step finishes
        successfully. It enables and switches to the corresponding result tab,
        updates StepBlock states, updates the top-bar indicators, and advances
        the automatic "Run all" queue if one is active.
    
        The function handles three different data-availability modes:
    
        1. ELD-only mode
           A diffraction image is loaded, but no structure files are loaded.
           Steps 3-6 can run, but XRD and calibration remain locked.
    
        2. CIF-only mode
           Structure files are loaded, but no diffraction image is loaded.
           Only XRD can run, because later image-processing steps require ELD 
           data.
    
        3. Full mode
           Both structure files and diffraction image are available.
           The complete workflow can progress through calibration.
    
        Parameters
        ----------
        step : int
            Step number that has just completed.
    
        Returns
        -------
        None
            Updates GUI state in place and may trigger the next automatic step.
        """
        # Map workflow steps to their corresponding result tabs.
        # Step 1 is not included because it is the loading/setup step.
        TAB = {
            2: self.TAB_XRD,
            3: self.TAB_ELLIPSE,
            4: self.TAB_ELD2D,
            5: self.TAB_1D,
            6: self.TAB_BKG,
            7: self.TAB_COMP,
        }
        
        # If this step has a result tab, enable it and switch to it immediately.
        if step in TAB:
            self._tabs.setTabEnabled(TAB[step], True)
            self._tabs.setCurrentIndex(TAB[step])
        
        # Determine which input data are currently available.
        has_cif  = bool(self.cif_paths)
        has_eld  = bool(self.eld_path)
        
        # These modes decide which subset of the pipeline is valid.
        eld_only = has_eld and not has_cif
        cif_only = has_cif and not has_eld

        if eld_only:
            # In ELD-only mode, image-processing steps 3-6 are available.
            # XRD step 2 and calibration step 7 stay locked because they need
            # structure files.
            self._blocks[step].set_state(StepBlock.DONE)
            
            next_step = step + 1
            
            # Step 2 is excluded because it cannot be completed without CIF data.
            done_steps = [s for s in range(1, step + 1) if s != 2]
            
            
            if next_step in (3, 4, 5, 6):
                # Continue through the available image-processing chain.
                self._blocks[next_step].set_state(StepBlock.ACTIVE)
                
                # Keep step 2 locked, and lock anything after the next active 
                # step.
                locked_steps = [2] + [
                    s for s in range(next_step + 1, 8) if s != 2
                    ]
                
                self._update_topbar_indicators_partial(
                    done=done_steps,
                    active=next_step,
                    locked=locked_steps
                    )
                
                self.statusBar().showMessage(
                    f"Step {step} done. Load structure files to enable XRD and calibration.")
            else:
                # We reached the final ELD-only step, usually background 
                # correction. Without structure files, calibration still cannot 
                # run.
                self._update_topbar_indicators_partial(
                    done=done_steps, 
                    active=0, 
                    locked=[2, 7]
                    )
                
                self.statusBar().showMessage(
                    f"Step {step} done. Load structure files to enable calibration."
                    )
                
        elif cif_only:
            # In CIF-only mode, only XRD can run.
            # Later steps require a diffraction image, so completing XRD does 
            # not
            # unlock the image-processing or calibration chain.
            self._blocks[step].set_state(StepBlock.DONE)
            self._update_topbar_indicators_partial(
                done=[1, step], 
                active=0,
                locked=[n for n in range(2, 8) if n != step]
                )
            
            self.statusBar().showMessage(
                "XRD calculated. Load a diffraction image to unlock the remaining steps.")
        
        elif step == 2:
            # XRD is independent of the image-processing chain.
            # Important case:
            # The user may first process an ELD-only dataset through steps 3-6,
            # then load CIF files and run XRD. In that case, completing step 2
            # should not reset or re-lock the already completed ELD steps.
            # Instead, mark XRD done and merge the current state into the full
            # pipeline, unlocking calibration if both sides are now available.
            self._blocks[2].set_state(StepBlock.DONE)
            self._merge_into_full_pipeline()
            
        elif step < self.N_STEPS:
            # Normal sequential progression:
            # mark current step done and activate the next one.
            self._go_to_step(step + 1)
            
        else:
            # Final step completed.
            # There is no next step to activate.
            self._blocks[step].set_state(StepBlock.DONE)

        # If a "Run all" sequence is in progress, automatically move on to
        # the next queued step now that this one has finished.
        self._advance_auto_run(step)


    def _on_run_all(self):
        """
        Start automatic execution of all currently applicable workflow steps.
    
        The queue is built from the data that are currently loaded:
    
        - CIF loaded: include XRD calculation, step 2. 
        - ELD loaded: include image-processing steps 3-6.
        - CIF and ELD both loaded: include calibration/comparison, step 7.
    
        This is especially useful after loading a saved session where all 
        necessary parameters are already filled, but it can also run from a 
        fresh setup.
    
        Returns
        -------
        None
            Initializes ``self._auto_run_queue`` and starts the first queued 
            step.
        """
        # Determine which branches of the workflow are available.
        has_cif = bool(self.cif_paths)
        has_eld = bool(self.eld_path)
        
        # At least one input source is needed before automatic execution makes 
        # sense.
        if not has_cif and not has_eld:
            QMessageBox.information(
                self, 
                "EDIFF",
                "Load a structure file and/or a diffraction image first."
                )
            return
        
        # Do not start a new automatic sequence while another computation is 
        # active. This prevents two workers from modifying the same GUI state 
        # or data model.
        if getattr(self, "_busy_keys", set()):
            QMessageBox.information(
                self, 
                "EDIFF",
                "A step is still running — please wait for it to finish "
                "before running all steps.")
            return
        
        # Build the automatic execution queue.
        queue = []
        
        # Structure branch: XRD calculation.
        if has_cif:
            queue.append(2)
            
        # Image branch: ellipse, center, 1D integration, background.            
        if has_eld:
            queue += [3, 4, 5, 6]
            
            # Calibration requires both calculated/reference structure 
            # information and processed diffraction data.
            if has_cif:
                queue.append(7)
                
        # Store queue on self so it can be consumed step by step after each
        # asynchronous operation completes.        
        self._auto_run_queue = queue
        self.statusBar().showMessage("Running all steps automatically...")
        
        # Start the first queued step. None means no previous step has 
        # completed yet.
        self._advance_auto_run(None)


    def _advance_auto_run(self, completed_step):
        """
        Advance the automatic "Run all" queue.
    
        This function is called once at the beginning of "Run all" and again 
        after each completed step. It removes completed steps from the queue 
        and schedules the next handler using ``QTimer.singleShot``.
    
        Parameters
        ----------
        completed_step : int or None
            Step that has just completed. Use None for the initial call, before 
            any queued step has run.
    
        Returns
        -------
        None
            Starts the next queued step or clears the queue when finished.
        """
        # Fetch the active queue, if one exists.
        queue = getattr(self, '_auto_run_queue', None)
        
        # No queue means no automatic sequence is currently running.
        if not queue:
            return
        
        
        if completed_step is not None:
            # Remove all queued steps up to and including the completed step.
            # This deliberately handles cases where one action skips over 
            # multiple logical steps. For example, skipping ellipse correction 
            # may effectively move the workflow past step 3 without a normal 
            # one-step transition.
            while queue and queue[0] <= completed_step:
                queue.pop(0)
                
        if not queue:
            # All queued steps have finished.
            self._auto_run_queue = None
            self.statusBar().showMessage(
                "Run all: finished — every step completed."
                )
            return
        
        # The next queued step determines which handler should be triggered.
        next_step = queue[0]
        
        if next_step == 3:
            # Ellipse correction is interactive.
            # If a saved session had an ellipse correction applied, reopen the
            # manual correction tool so the correction can be reproduced/
            # confirmed. Otherwise skip ellipse correction automatically.
            # In both cases, auto-run waits until _complete_step(3) is 
            # eventually called by the interactive workflow.
            
            handler = (
                self._on_ellipse
                if getattr(self, '_session_ellipse_was_corrected', False)
                else self._on_skip_ellipse
                )
            
        else:
            # Map non-interactive or internally managed steps to their handlers
            handler = {
                2: self._on_xrd,
                4: self._on_center,
                5: self._on_1d,
                6: self._on_bkg,
                7: self._on_calibrate,
            }.get(next_step)
            
        if handler:
            # Delay the next step slightly so the GUI has time to finish 
            # repainting and processing signals from the previous step.
            QTimer.singleShot(300, handler)

    def _close_external_windows(self):
        """
        Close external/intermediate UI elements created during interactive 
        steps.
    
        This includes:
    
        - interactive subprocess windows, for example ellipse or center 
          picking;
        - matplotlib figures opened in the current process;
        - tutorial dialogs opened by Tutorial mode.
    
        The method is used by Stop and Abort actions to clean up visible 
        windows before the workflow state is reset.
    
        Returns
        -------
        None
            Attempts cleanup in place. Individual cleanup failures are ignored.
        """
        # Terminate any running interactive subprocess first.
        # This is important for tools whose Qt/matplotlib windows live outside 
        # the main GUI process.
        try:
            self.pipeline.terminate_subprocess()
        except Exception:
            # Cleanup should be best-effort. Failure here should not crash 
            # Stop/Abort.
            pass
        
        # Close every open matplotlib figure in this process
        try:
            import matplotlib.pyplot as plt
            plt.close('all')
            
        except Exception:
            # Matplotlib may be unavailable or already in a partially closed 
            # state.
            pass
        
        # Close tutorial dialogs that were opened by Tutorial mode.
        if hasattr(self, '_tutorial_dialogs'):
            for d in list(self._tutorial_dialogs):
                try:
                    if d.isVisible():
                        d.close()
                except RuntimeError:
                    # The Python wrapper may still exist after the underlying C++
                    # Qt object has already been deleted.
                    pass
                
            # Reset the list after pruning/closing all known dialogs.    
            self._tutorial_dialogs = []
            

    def _stop_step(self, step_n, run_btn, stop_btn, card):
        """
        Stop a running workflow step and reset its visible GUI controls.
    
        This is used by step-specific Stop buttons. It terminates external windows,
        cancels the active worker if possible, restores the result card to its empty
        state, and makes the step active again so the user can rerun it.
    
        Parameters
        ----------
        step_n : int
            Number of the workflow step being stopped.
    
        run_btn : QPushButton
            The normal "Run" button for the step.
    
        stop_btn : QPushButton
            The corresponding "Stop" button for the step.
    
        card : QWidget-like
            Result card associated with the step. Expected to contain ``_stack``.
    
        Returns
        -------
        None
            Updates GUI state in place.
        """
        # Kill external windows before canceling the worker. This matters for
        # interactive subprocesses because their windows should disappear 
        # immediately.
        self._close_external_windows()
        
        # Fetch and cancel the worker associated with this step, if one exists.
        # The worker may not stop instantly if its underlying computation is 
        # not interruptible, but its result should be suppressed by the worker 
        # logic.
        worker = getattr(self, f"_active_worker_{step_n}", None)
        
        if worker: worker.cancel()
        
        # Return the result card to its initial/empty state 
        card._stack.setCurrentIndex(0)
        
        # Restore the normal button state:
        # hide Stop, show Run, and allow the user to run the step again.
        stop_btn.setVisible(False)
        run_btn.setVisible(True)
        run_btn.setEnabled(True)
        
        # Keep this step active so the user can adjust parameters and rerun it.
        self._blocks[step_n].set_state(StepBlock.ACTIVE)
        
        # Step 1 remains editable throughout the workflow.
        for w in self._blocks[1]._children: w.setEnabled(True)
        
        # Step-specific Stop has handled the interruption, so the global Abort
        # button no longer needs to be active.
        self._btn_abort.setEnabled(False)
        
        
        if step_n == 3:
            # Ellipse correction is special because the user can either rerun 
            # it or explicitly skip it to continue the workflow.
            self.statusBar().showMessage(
                "Ellipse correction stopped. Rerun or skip to continue.")
        else:
            self.statusBar().showMessage(
                f"Step {step_n} stopped. Adjust parameters and rerun.")


    def _step_busy(self, step_n):
        """
        Return whether a step still has a computation finishing 
        in the background.
    
        This is used after Stop/Abort actions. Some computations cannot be 
        killed immediately once started, so the GUI tracks them as busy until 
        their worker or thread actually finishes.
    
        Parameters
        ----------
        step_n : int
            Workflow step number.
    
        Returns
        -------
        bool
            True if ``step_n`` is present in ``self._busy_keys``.
        """
        # _busy_keys may not exist early in startup, so getattr defaults to an
        # empty set.
        return step_n in getattr(self, "_busy_keys", set())

    def _start_step(self, step_n, run_btn, stop_btn):
        """
        Put a workflow step into its running UI state.
    
        This hides the normal run button, shows the stop button, enables 
        the global Abort button, and records which step is currently active.
    
        Parameters
        ----------
        step_n : int
            Step number that is starting.
    
        run_btn : QPushButton
            Normal run button for the step.
    
        stop_btn : QPushButton
            Stop button for the step.
    
        Returns
        -------
        None
            Updates button visibility and active-step tracking in place.
        """
        # While the step is running, the user should not be able to click Run 
        # again.
        run_btn.setVisible(False)
        
        # The step-specific Stop button becomes available during execution.
        stop_btn.setVisible(True)
        stop_btn.setEnabled(True)
        
        # The global Abort button is also enabled while any step is active.
        self._btn_abort.setEnabled(True)
        
        # Record the active step so _on_abort can find the correct 
        # worker/buttons.
        self._current_step_n = step_n

    def _finish_step(self, step_n, run_btn, stop_btn):
        """
        Restore normal button state after a workflow step finishes.
    
        Parameters
        ----------
        step_n : int
            Step number that has finished. Currently unused directly, but kept in
            the signature for consistency with ``_start_step`` and call sites.
    
        run_btn : QPushButton
            Normal run button for the step.
    
        stop_btn : QPushButton
            Stop button for the step.
    
        Returns
        -------
        None
            Updates button visibility and active-step tracking in place.
        """
        # Hide Stop and show Run again.
        stop_btn.setVisible(False)
        run_btn.setVisible(True)
        run_btn.setEnabled(True)
        
        # No step is actively running, so the global Abort button can be 
        # disabled.
        self._btn_abort.setEnabled(False)
        
        # Clear active-step tracking.
        self._current_step_n = None

    def _on_abort(self):
        """
        Abort the currently running workflow step without clearing previous 
        results.
    
        Abort cancels the active worker and closes external windows, but it 
        does not delete already available results. If the underlying computation 
        cannot be interrupted immediately, the step remains internally busy 
        until the worker actually finishes. In that case, the run button stays 
        disabled and is re-enabled only when the background thread emits 
        ``finished``.
    
        Returns
        -------
        None
            Cancels the current step if one is active and updates GUI state.
        """
        # Identify the currently running step.
        step_n = getattr(self, '_current_step_n', None)
        
        if step_n is None:
            # Nothing is currently running
            # make sure Abort is disabled and exit.
            self._btn_abort.setEnabled(False)
            return
        
        # Close external/intermediate windows before canceling the worker.
        self._close_external_windows()
        
        # Cancel the worker associated with the current step.
        worker = getattr(self, f"_active_worker_{step_n}", None)
        if worker:
            worker.cancel()
            
        # Clear global running state.            
        self._btn_abort.setEnabled(False)
        self._current_step_n = None

        # Map step numbers to their corresponding buttons and result cards.
        # This allows one generic abort handler to restore the correct UI 
        # controls for whichever step was interrupted.
        btn_map = {
            2: self._btn_xrd,    3: self._btn_ellipse,
            4: self._btn_center, 5: self._btn_1d,
            6: self._btn_bkg,    7: self._btn_cal,
        }
        stop_map = {
            2: self._btn_stop2, 3: self._btn_stop3,
            4: self._btn_stop4, 5: self._btn_stop5,
            6: self._btn_stop6, 7: self._btn_stop7,
        }
        card_map = {
            2: self._card_xrd,     3: self._card_ellipse,
            4: self._card_eld2d,   5: self._card_1d,
            6: self._card_bkg,     7: self._card_comp,
        }
        
        if step_n in btn_map:
            # Hide the step-specific Stop button.
            stop_map[step_n].setVisible(False)
            
            # Return the result card to either the last completed run or 
            # the empty state. This avoids leaving the spinner/progress state 
            # visible after aborting.
            card = card_map[step_n]
            
            if card._run_tabs.count() > 0:
                card._stack.setCurrentIndex(2)   # show last run
            else:
                card._stack.setCurrentIndex(0)   # show empty state

            # If the computation is still busy internally, keep the run button
            # disabled until the worker thread genuinely finishes.            
            still_running = self._step_busy(step_n)
            
            btn_map[step_n].setVisible(True)
            btn_map[step_n].setEnabled(not still_running)
            
            if still_running:
                # Find the thread behind this worker and re-enable the button
                # for real once it actually finishes in the background, instead
                # of guessing at a delay.
                thread = next(
                    (t for t, w in getattr(self,"_threads",[]) if w is worker),
                    None,
                    )
                
                if thread is not None:
                    btn = btn_map[step_n]
                    thread.finished.connect(lambda: btn.setEnabled(True))
        
        # Keep the aborted step active so the user can adjust parameters 
        # and rerun.
        self._blocks[step_n].set_state(StepBlock.ACTIVE)
        
        if self._step_busy(step_n):
            # Honest UI message: the user aborted the result flow, but 
            # the low-level calculation may still be returning in the background.
            self.statusBar().showMessage(
                f"Step {step_n} aborted — the current calculation can't be "
                "interrupted mid-run, so it's still finishing in the "
                "background. The button will re-enable automatically once "
                "it's done; results so far are preserved.")
        else:
            self.statusBar().showMessage(
                f"Step {step_n} aborted. Results so far are preserved — adjust and rerun.")


    def _on_toggle_tutorial(self):
        """
        Toggle Tutorial mode and update the toolbar button appearance.
    
        When Tutorial mode is enabled, interactive methods can automatically 
        open explanatory help panels. Turning the mode off does not remove 
        normal Help buttons; it only disables automatic tutorial popups.
    
        Returns
        -------
        None
            Updates tutorial-mode state and button styling in place.
        """
        # Flip the current tutorial-mode boolean.
        self._tutorial_mode = not self._tutorial_mode
                
        if self._tutorial_mode:
            # Active appearance: filled/blue tutorial button.
            self._btn_tutorial_mode.setText("Tutorial  ●")
            self._btn_tutorial_mode.setStyleSheet(self._tutorial_on_style)
            
            # Show a short intro dialog slightly after the button updates, 
            # so the GUI has time to repaint first..
            QTimer.singleShot(50, self._show_tutorial_intro)
            
        else:
            # Inactive appearance: outline/off tutorial button.
            self._btn_tutorial_mode.setText("Tutorial  ○")
            self._btn_tutorial_mode.setStyleSheet(self._tutorial_off_style)
            self.statusBar().showMessage(
                "Tutorial mode off — help buttons still available on demand.")


    def _show_tutorial_intro(self):
        """
        Show a small introduction dialog explaining Tutorial mode.
    
        The dialog briefly tells the user what automatic tutorial panels do and
        visually distinguishes the normal Help button from the Tutorial-mode 
        toggle. It is displayed when Tutorial mode is switched on.
    
        Returns
        -------
        None
            Creates and shows a transient Qt dialog.
        """
        # Create a small top-level child window.
        # WA_DeleteOnClose ensures Qt deletes the underlying C++ object when 
        # closed.
        dlg = QWidget(
            self, 
            Qt.WindowType.Window | Qt.WindowType.WindowStaysOnTopHint
            )
        
        dlg.setWindowTitle("Tutorial mode is ON")
        dlg.setFixedSize(480, 280)
        dlg.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose)
        
        # Apply a local style so the dialog looks consistent regardless of the
        # parent widget styling.
        dlg.setStyleSheet(
            "QWidget { background:#ffffff; font-family:'Segoe UI',Arial,sans-serif; }"
            )
        
        # Main vertical layout for the dialog contents.
        lay = QVBoxLayout(dlg)
        lay.setContentsMargins(20, 16, 20, 16); lay.setSpacing(8)

        # Heading
        h = QLabel("Tutorial mode is now ON")
        h.setStyleSheet("font-size:13px; font-weight:700; color:#0071e3;")
        lay.addWidget(h)
        
        # Thin blue accent line under the heading.
        accent = QFrame()
        accent.setFixedHeight(2)
        accent.setStyleSheet("background:#0071e3; border:none;")
        lay.addWidget(accent)

        # Main explanatory paragraph.
        body = QLabel(
            "When you run an interactive method (manual center detection, "
            "manual refinement, or ellipse correction), a help panel opens "
            "automatically alongside the interactive window.")
        body.setStyleSheet("font-size:12px; color:#1d1d1f; background:transparent;")
        body.setWordWrap(True)
        lay.addWidget(body)

        lay.addSpacing(6)

        # Two-column explanation table:
        # left column shows a visual replica of a button,
        # right column describes what that button means.
        from qtpy.QtWidgets import QGridLayout
        
        table_w = QWidget(); table_w.setStyleSheet("background:transparent;")
        grid = QGridLayout(table_w)
        grid.setContentsMargins(0, 0, 0, 0)
        grid.setHorizontalSpacing(12)
        grid.setVerticalSpacing(8)
        grid.setColumnStretch(1, 1)

        def btn_replica(label, style):
            """
            Create a non-clickable QLabel that visually imitates a toolbar 
            button.
    
            This is used only inside the tutorial intro dialog so the user can
            recognize the real toolbar buttons without interacting with them.
    
            Parameters
            ----------
            label : str
                Text shown in the replica.
    
            style : str
                Qt stylesheet applied to the replica.
    
            Returns
            -------
            QLabel
                Styled label that looks like a small button.
            """
            b = QLabel(label)
            b.setFixedHeight(26)
            b.setAlignment(Qt.AlignmentFlag.AlignCenter)
            b.setStyleSheet(style)
            
            # Fixed size policy prevents the replica from stretching like 
            # a normal text label.
            b.setSizePolicy(
                QSizePolicy.Policy.Fixed, 
                QSizePolicy.Policy.Fixed
                )
            
            # Compute a reasonable width from the text and add horizontal 
            # padding.
            b.adjustSize()
            b.setMinimumWidth(b.sizeHint().width() + 16)
            return b

        def desc(text):
            """
            Create a description label for the tutorial-intro table.

            Parameters
            ----------
            text : str
                Description text.

            Returns
            -------
            QLabel
                Left-aligned label vertically centered in the grid row.
            """
            l = QLabel(text)
            l.setStyleSheet(
                "font-size:12px; color:#1d1d1f; background:transparent;"
            )
            l.setAlignment(
                Qt.AlignmentFlag.AlignVCenter | Qt.AlignmentFlag.AlignLeft
            )

            return l

        # Replica of a normal Help button.
        help_btn = btn_replica("?  Help",
            "font-size:11px; font-weight:600; color:#0071e3;"
            "background:#f0f6ff; border:1px solid #b3d0f7;"
            "border-radius:6px; padding:0 8px;")

        # Replica of the active Tutorial-mode button.
        tutorial_btn = btn_replica("Tutorial  ●",
            "font-size:11px; font-weight:600; color:#ffffff;"
            "background:#0071e3; border-radius:12px; padding:0 10px;")

        # Add first row: normal Help button.
        grid.addWidget(help_btn, 0, 0, Qt.AlignmentFlag.AlignVCenter)
        grid.addWidget(
            desc(
                "Available in each step — click any time for help."
            ),
            0,
            1,
            Qt.AlignmentFlag.AlignVCenter
        )

        # Add second row: Tutorial-mode toggle
        grid.addWidget(tutorial_btn, 1, 0, Qt.AlignmentFlag.AlignVCenter)
        grid.addWidget(
            desc("Click again in the toolbar to turn Tutorial mode off."),
            1,
            1,
            Qt.AlignmentFlag.AlignVCenter
        )

        lay.addWidget(table_w)

        # Push the OK button to the bottom.
        lay.addStretch()

        # Close button.
        ok = QPushButton("Got it")
        ok.setFixedWidth(80)
        ok.clicked.connect(dlg.close)
        lay.addWidget(ok, alignment=Qt.AlignmentFlag.AlignRight)

        # Center the dialog over the main application window.
        geo = self.geometry()
        dlg.move(geo.center().x() - dlg.width() // 2,
                 geo.center().y() - dlg.height() // 2
                 )

        dlg.show()

        self.statusBar().showMessage(
            "Tutorial mode ON — help panels open automatically for interactive methods."
        )


    def _maybe_show_tutorial(self, key):
        """
        Open a tutorial dialog for a given help topic if Tutorial mode
        is enabled.

        This helper is called by interactive or semi-interactive workflow
        methods. It silently does nothing when Tutorial mode is off.

        Parameters
        ----------
        key : str
            Tutorial/help-topic key passed to ``TutorialDialog``.
            The dialog uses this key to decide which help content
            to display.

        Returns
        -------
        None
            Shows a tutorial dialog if Tutorial mode is active.
        """
        # Do nothing unless Tutorial mode is currently enabled.
        if not self._tutorial_mode:
            return

        # Lazily create the dialog-tracking list.
        if not hasattr(self, '_tutorial_dialogs'):
            self._tutorial_dialogs = []

        # Remove stale dialog references before opening a new one.
        # Important Qt detail:
        # Calling isVisible() can raise RuntimeError if the Python object
        # still exists but the underlying C++ QWidget has already been
        # deleted.

        alive = []
        for d in self._tutorial_dialogs:
            try:
                if d.isVisible():
                    alive.append(d)
            except RuntimeError:
                # Underlying Qt object is already gone
                # discard this reference.
                pass

        self._tutorial_dialogs = alive

        # Create and show the tutorial dialog for the requested topic.
        dlg = TutorialDialog(key)

        # Keep a Python reference so the dialog is not garbage-collected
        # while open.
        self._tutorial_dialogs.append(dlg)
        dlg.show()

    # ── Step 1 ────────────────────────────────────────────────────────────────
    def _update_workdir_label(self):
        """
        Update the displayed working-directory label.

        The full working-directory path can be too long to fit
        nicely into the GUI. This helper shortens the displayed
        text using middle elision while keeping the complete path
        available as a tooltip.

        Returns
        -------
        None
            Updates ``self._wd_path_lbl`` in place.
        """
        # QFontMetrics knows how wide text will be with the label's
        # current font. This lets Qt produce an elided string that fits
        # into the chosen width.
        from qtpy.QtGui import QFontMetrics

        fm = QFontMetrics(self._wd_path_lbl.font())

        # Shorten the path in the middle, e.g.
        # C:/Users/.../experiment/session
        # Middle elision is preferable for paths because the start and
        # the filename or final folder are usually both useful.

        elided = fm.elidedText(
            str(self.workdir),
            Qt.TextElideMode.ElideMiddle,
            280
        )

        # Show the shortened path in the GUI.
        self._wd_path_lbl.setText(elided)

        # Keep the full path accessible on hover.
        self._wd_path_lbl.setToolTip(str(self.workdir))

    # ── Session save / load ──────────────────────────────────────────────────
    def _collect_session_dict(self):
        """
        Collect the current GUI/session state into a JSON-serializable
        dictionary.

        This function gathers all file paths and user-editable parameters
        currently stored in the GUI. The returned dictionary is intended
        to be written to a session JSON file and later restored using
        ``_apply_session_dict``.

        The saved state includes:
            - loaded CIF file paths;
            - loaded diffraction-image path;
            - current working directory;
            - XRD simulation parameters;
            - ellipse-correction settings;
            - center-detection settings;
            - 1D integration settings;
            - background-subtraction settings;
            - calibration/comparison settings.

        Important implementation detail
        -------------------------------
        For the ELD path, the function prefers ``self.pipeline._eld_path``
        if it exists. This means that if ellipse correction has already
        produced a corrected image, the session saves the currently active
        corrected image rather than only the originally browsed raw image.

        Returns
        -------
        dict
            Plain JSON-compatible session dictionary.
        """
        return {
            # Session format version. This gives you a hook for future
            # migration if the saved-session structure changes later.
            "version": 1,

            # Store all selected structure files.
            # Convert to list explicitly in case self.cif_paths is
            # another iterable.
            "cif_paths": list(self.cif_paths),

            # Save the image path that the pipeline is currently
            # using
            # Priority:
            #   1. pipeline._eld_path, if present;
            #   2. self.eld_path, if present;
            #   3. empty string.
            # This matters after ellipse correction: the corrected image
            # may be the active pipeline input, even though the original
            # Step-1 file path was araw ELD image.
            "eld_path": str(
                getattr(
                    self.pipeline,
                    '_eld_path',
                    None
                ) or self.eld_path or ""
            ),

            # Save current working directory as a string for JSON
            # compatibility
            "workdir": str(self.workdir),

            # All step-specific GUI parameters are grouped under
            # "params".
            "params": {
                # Global display settings shared across every step that
                # shows the 2D pattern (Preview, Ellipse, Center).
                "display": {
                    "icut": self._icut.value(),
                    "cmap": self._cmap.currentText(),
                },

                # Step 2: XRD calculation parameters
                "xrd": {
                    # Radiation wavelength used for simulated diffraction.
                    "wavelength": self._xrd_wl.value(),

                    # Two-theta range for the calculated XRD pattern.
                    "tt_min": self._xrd_ttmin.value(),
                    "tt_max": self._xrd_ttmax.value(),

                    # Peak broadening / smoothing parameter.
                    "sigma": self._xrd_sigma.value(),

                    # Temperature factor used in the XRD model.
                    "temp_factor": self._xrd_temp.value(),

                    # Selected peak/profile model.
                    "profile": self._xrd_profile.currentText(),
                },

                # Step 3: ellipse-correction parameters
                "ellipse": {
                    # Intensity cutoff used for display or processing.
                    # "icut": self._ell_icut.value(),

                    # Colormap used when showing the ELD image.
                    # "cmap": self._ell_cmap.currentText(),



                    # Remember whether the currently loaded ELD path differs
                    # from the raw ELD path. If yes, ellipse correction has
                    # probably already been applied and baked into the active
                    # image. This flag is used during session reload to inform
                    # the user that the correction was part of the saved
                    # workflow.
                    "was_corrected": bool(
                        getattr(self.pipeline, '_eld_path', None)
                        and getattr(self.pipeline, '_raw_eld_path', None)
                        and str(self.pipeline._eld_path)
                        != str(self.pipeline._raw_eld_path)
                    ),
                },

                # Step 4: center-detection parameters
                "center": {
                    # Main center-detection method.
                    "detection": self._ctr_det.currentText(),

                    # Optional refinement method applied after initial
                    # detection.
                    "refinement": self._ctr_ref.currentText(),

                    # Display/processing cutoff and visualization colormap.
                    # "icut": self._ctr_icut.value(),
                    # "cmap": self._ctr_cmap.currentText(),

                    # Method-specific parameters for center finding.
                    "csquare": self._ctr_csq.value(),
                    "cintensity": self._ctr_cintensity.value(),
                    "rtype": self._ctr_rtype.currentText(),
                    "heq": self._ctr_heq.currentText(),
                    "sobel": self._ctr_sobel.currentText(),
                    "masking": self._ctr_masking.currentText(),
                    "normalize_bg": self._ctr_normalize_bg.currentText(),
                    "live_plot": self._ctr_live_plot.currentText(),

                    # Method-specific parameters for center finding.
                    "save_file": self._ctr_file.text(),

                    # If detection reads from a file, save the actual path
                    # used, so "Run all" can point straight at it again.
                    "read_file_path": self._ctr_read_file_path,
                },

                # Step 5: 1D integration parameters
                "1d": {
                    # Pixel-radius range used for 1D integration.
                    "px_min": self._px_min.value(),
                    "px_max": self._px_max.value(),

                    # Output filename for the 1D curve.
                    "outfile": self._1d_outfile.text(),
                },

                # Step 6: background-subtraction parameters
                "background": {
                    # Selected background-removal method.
                    "method": self._bkg_method.currentText(),

                    # Path/name of background file, if the chosen method
                    # uses one.
                    "bkg_file": self._bkg_file.text(),

                    # SimpleFuncs/Baselines options: the algorithm, the
                    # active algorithm's dedicated parameters, the pixel
                    # range (-1 = "auto") and the free-form extras.
                    "algorithm": self._bkg_algo.currentText(),
                    "params": self._read_bkg_params(),
                    "range_start": self._bkg_xmin.value(),
                    "range_end": self._bkg_xmax.value(),
                    "extra": self._bkg_extra.text(),
                },

                # Step 7: calibration/comparison parameters
                "calibration": {
                    # Calibration method and selected reference phase.
                    "method": self._cal_method.currentText(),
                    "ref_idx": self._cal_ref.currentIndex(),

                    # Manual calibration scale correction and q-display
                    # range.
                    "finetune": self._cal_finetune.value(),
                    "qmin": self._cal_qmin.value(),
                    "qmax": self._cal_qmax.value(),

                    # Microscope/camera parameters used by Microscope or
                    # MicroscopeParams calibration methods.
                    "voltage": self._cal_voltage.value(),
                    "cl": self._cal_cl.value(),
                    "pixsize": self._cal_pixsize.value(),
                    "binning": self._cal_binning.value(),
                },
            },
        }

    def _apply_session_dict(self, data):
        """
        Restore GUI fields from a saved session dictionary.

        This function updates the GUI controls from a dictionary
        previously created by ``_collect_session_dict``. It only
        populates fields; it does not run any processing steps.

        The function is intentionally tolerant of missing keys.
        Each value is applied only if present, which makes older
        or partially written session files easier to load.

        Some saved methods are interactive and cannot be fully
        replayed without user input. For those cases, this function
        adds human-readable notes to the returned list so the caller
        can inform the user.

        Parameters
        ----------
        data : dict
            Session dictionary loaded from JSON.

        Returns
        -------
        list of str
            Notes describing anything that requires user attention
            during replay, such as manual point selection or
            interactive background picking.
        """
        # Notes are returned to the caller and can be shown in a message
        # box or status message after loading the session.
        notes = []

        # NOTE: this must come BEFORE the display block below, which reads
        # `p`; it used to follow it, so loading any session raised
        # UnboundLocalError.
        # Step parameters are grouped under "params".
        # Missing "params" is allowed, in which case defaults remain
        # unchanged.
        p = data.get("params", {})

        # Global display settings
        disp = p.get("display", {})
        if "icut" in disp:
            self._icut.setValue(disp["icut"])
        if "cmap" in disp:
            self._cmap.setCurrentText(disp["cmap"])


        # Step 2: restore XRD parameters
        xrd = p.get("xrd", {})

        # Each key is checked independently so older session files with
        # fewer parameters can still be loaded.
        if "wavelength" in xrd:
            self._xrd_wl.setValue(xrd["wavelength"])
        if "tt_min" in xrd:
            self._xrd_ttmin.setValue(xrd["tt_min"])
        if "tt_max" in xrd:
            self._xrd_ttmax.setValue(xrd["tt_max"])
        if "sigma" in xrd:
            self._xrd_sigma.setValue(xrd["sigma"])
        if "temp_factor" in xrd:
            self._xrd_temp.setValue(xrd["temp_factor"])
        if "profile" in xrd:
            self._xrd_profile.setCurrentText(xrd["profile"])

        # Step 3: restore ellipse-correction parameters
        ell = p.get("ellipse", {})
        # if "icut" in ell:
        #    self._ell_icut.setValue(ell["icut"])

        # if "cmap" in ell:
        #     self._ell_cmap.setCurrentText(ell["cmap"])

        # Remember whether the saved session had ellipse correction
        # applied. This flag is later used by Run all / ellipse workflow
        # logic.
        self._session_ellipse_was_corrected = bool(
            ell.get("was_corrected", False)
        )

        if self._session_ellipse_was_corrected:
            # The actual four manually selected points are not saved here.
            # Therefore, rerunning the step still requires user interaction.
            notes.append(
                "Ellipse correction was used — you'll need to click the "
                "4 points again when this step runs.")

        # Step 4: restore center-detection parameters
        ctr = p.get("center", {})

        if "csquare" in ctr:
            self._ctr_csq.setValue(ctr["csquare"])
        if "cintensity" in ctr:
            self._ctr_cintensity.setValue(ctr["cintensity"])
        if "rtype" in ctr:
            self._ctr_rtype.setCurrentText(ctr["rtype"])
        if "heq" in ctr:
            self._ctr_heq.setCurrentText(ctr["heq"])
        if "sobel" in ctr:
            self._ctr_sobel.setCurrentText(ctr["sobel"])
        if "masking" in ctr:
            self._ctr_masking.setCurrentText(ctr["masking"])
        if "normalize_bg" in ctr:
            self._ctr_normalize_bg.setCurrentText(ctr["normalize_bg"])
        if "live_plot" in ctr:
            self._ctr_live_plot.setCurrentText(ctr["live_plot"])
        if "save_file" in ctr:
            self._ctr_file.setText(ctr["save_file"])

        # Detection method and optional read-from-file path need
        # slightly more logic because manual detection requires
        # a note, while read-file detection may need to restore
        # a file label.
        detection = ctr.get("detection", "")
        read_path = ctr.get("read_file_path")

        if detection:
            self._ctr_det.setCurrentText(detection)
        if detection == "manual":
            # Manual center detection cannot be replayed automatically
            # unless the clicked center was saved somewhere else.
            notes.append(
                "Center detection was 'manual' — you'll need to click the "
                "center again when this step runs."
            )

        if read_path:
            # Preserve the saved read-file path for later use by the
            # center detection step.
            self._ctr_read_file_path = read_path

            # If the file still exists, update the GUI label to show
            # the basename and mark it visually as valid.
            if Path(read_path).exists():
                self._ctr_read_file_lbl.setText(Path(read_path).name)
                self._ctr_read_file_lbl.setStyleSheet(
                    "color:#1d7a3c; font-size:12px; font-weight:600;"
                )

        if "refinement" in ctr:
            self._ctr_ref.setCurrentText(ctr["refinement"])

        # Step 5: restore 1D integration parameters
        d1 = p.get("1d", {})
        if "px_min" in d1:
            self._px_min.setValue(d1["px_min"])
        if "px_max" in d1:
            self._px_max.setValue(d1["px_max"])
        if "outfile" in d1:
            self._1d_outfile.setText(d1["outfile"])

        # Step 6: restore background-subtraction parameters
        bkg = p.get("background", {})
        method = bkg.get("method", "")

        if "bkg_file" in bkg:
            self._bkg_file.setText(bkg["bkg_file"])

        if method:
            self._bkg_method.setCurrentText(method)

        # Order matters: selecting the algorithm resets the parameter
        # fields to that algorithm's defaults, so it must come BEFORE the
        # saved parameter values are applied (else they'd be overwritten).
        # Older session files lack all of these keys - each is optional.
        if bkg.get("algorithm"):
            self._bkg_algo.setCurrentText(bkg["algorithm"])
        for pname, pval in (bkg.get("params") or {}).items():
            w = self._bkg_param_widgets.get(pname)
            if w is None:
                continue
            try:
                if isinstance(w, QComboBox):
                    w.setCurrentText("True" if pval in (True, "True") else "False")
                else:
                    w.setValue(pval)
            except (TypeError, ValueError):
                pass  # hand-edited / corrupt value: keep the default
        if "range_start" in bkg:
            self._bkg_xmin.setValue(bkg["range_start"])
        if "range_end" in bkg:
            self._bkg_xmax.setValue(bkg["range_end"])
        if "extra" in bkg:
            self._bkg_extra.setText(bkg["extra"])

        if method == "InteractivePlot":
            # Interactive background selection cannot be fully restored
            # from just the method name, because the clicked background
            # points are not stored.
            notes.append(
                "Background method was 'InteractivePlot' — you'll need "
                "to click the background points again when this step runs."
            )

        # Step 7: restore calibration/comparison parameters
        cal = p.get("calibration", {})
        if "method" in cal:
            self._cal_method.setCurrentText(cal["method"])
        if "ref_idx" in cal:
            self._cal_ref.setCurrentIndex(cal["ref_idx"])
        if "finetune" in cal:
            self._cal_finetune.setValue(cal["finetune"])
        if "qmin" in cal:
            self._cal_qmin.setValue(cal["qmin"])
        if "qmax" in cal:
            self._cal_qmax.setValue(cal["qmax"])
        if "voltage" in cal:
            self._cal_voltage.setValue(cal["voltage"])
        if "cl" in cal:
            self._cal_cl.setValue(cal["cl"])
        if "pixsize" in cal:
            self._cal_pixsize.setValue(cal["pixsize"])
        if "binning" in cal:
            self._cal_binning.setValue(cal["binning"])

        # Return all user-facing notes to the caller.
        return notes

    def _on_save_session(self):
        """
        Save the current EDIFF GUI state to a session JSON file.

        The session file stores file paths and all relevant GUI
        parameters collected by ``_collect_session_dict``. It does
        not save computed arrays or figures directly; instead, it
        saves enough information to restore the setup and rerun
        the pipeline later.

        The default filename is ``session.ediff.json`` inside the
        current working directory.

        Returns
        -------
        None
            Writes a session file if the user chooses a path
            successfully.
        """
        # Ask the user where to save the session file.
        # The default location is the current EDIFF working directory.
        path, _ = QFileDialog.getSaveFileName(
            self,
            "Save session",
            str(self.workdir / "session.ediff.json"),
            "EDIFF session (*.ediff.json)"
        )

        # User cancelled the dialog.
        if not path:
            return

        # Ensure the saved file has a JSON-compatible extension.
        # Note:
        # This currently appends ".ediff.json" to anything that
        # does not end with ".json". For example, "my_session"
        # becomes "my_session.ediff.json".
        if not path.endswith(".json"):
            path += ".ediff.json"
        try:
            # Collect the current session state into a plain
            # dictionary.
            data = self._collect_session_dict()

            # Write the session using UTF-8 and indentation so
            # the file remains human-readable and easy to debug.
            with open(path, "w", encoding="utf-8") as f:
                json.dump(data, f, indent=2)

        except Exception as e:
            # Saving can fail because of permissions, invalid paths,
            # disk issues, serialization problems, etc. Show the exception
            # to help debugging.
            QMessageBox.warning(
                self,
                "EDIFF",
                f"Could not save session:\n{e}"
            )
            return

        # Confirm successful save in the status bar.
        self.statusBar().showMessage(f"Session saved to: {path}")

    def _on_load_session(self):
        """
        Load an EDIFF session JSON file and optionally rerun
        the pipeline.

        This function restores a previously saved session
        in three stages:

        1. Ask the user to choose a session JSON file.
        2. Load referenced CIF and ELD files, if they still exist.
        3. Apply all saved GUI parameters using ``_apply_session_dict``.

        Missing referenced files are reported to the user, but they
        do not prevent the remaining session settings from being
        loaded. After loading, the user is asked whether the whole
        applicable pipeline should be run immediately with
        the restored settings.

        Returns
        -------
        None
            Updates the GUI state from the selected session file.
        """
        # Ask the user to select a saved EDIFF session file.
        # The dialog starts in the current working directory.
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Load session",
            str(self.workdir),
            "EDIFF session (*.ediff.json *.json)"
        )

        # User cancelled the dialog.
        if not path:
            return

        try:
            # Read and parse the JSON session file.
            with open(path, "r", encoding="utf-8") as f:
                data = json.load(f)

        except Exception as e:
            # Loading can fail if the file is missing, unreadable,
            # malformed JSON, or not actually an EDIFF session file.
            QMessageBox.warning(
                self,
                "EDIFF",
                f"Could not read session file:\n{e}"
            )
            return

        # Restore referenced input files
        # Load referenced files similarly to Step 1 browsing.
        # Missing files are collected and reported later, but they
        # do not block restoring the rest of the GUI parameters.
        missing = []

        # Restore CIF/structure files.
        cif_paths = data.get("cif_paths") or []

        # Keep only paths that still exist on disk.
        valid_cifs = [p for p in cif_paths if Path(p).exists()]

        # Track paths from the session that could not be found.
        missing += [p for p in cif_paths if p not in valid_cifs]

        if valid_cifs:
            # Update internal CIF path list and refresh the Step-1
            # GUI display.
            self.cif_paths = valid_cifs
            self._refresh_cif_list()

        # Restore ELD/diffraction-image path.
        eld_path = data.get("eld_path")

        if eld_path and Path(eld_path).exists():
            # If the image still exists, load it into the GUI state.
            self.eld_path = eld_path
            self._refresh_eld_list()

            # Show the loaded image name in the main window title.
            self.setWindowTitle(f"EDIFF - {Path(eld_path).name}")

        elif eld_path:
            # The session referenced an image file, but it no
            # longer exists.
            missing.append(eld_path)

        # If at least one valid input file was restored, update the
        # pipeline state. This should unlock the appropriate workflow
        # steps depending on whether CIF, ELD, or both are available.
        if valid_cifs or (eld_path and Path(eld_path).exists()):
            self._check_files_ready()

        # Restore all saved GUI fields and collect notes about
        # interactive steps that cannot be replayed fully without
        # user input.
        notes = self._apply_session_dict(data)

        # Build user-facing summary message
        msg_lines = ["Session loaded."]
        if missing:
            # Report missing files in one message. The user can browse
            # for them manually after the session settings are restored.
            msg_lines.append(
                "Could not find " + ", ".join(missing) +
                " — you'll need to re-browse for these."
            )

        if notes:
            # Add notes returned by _apply_session_dict(), e.g. manual
            # ellipse or background-selection steps that require new clicks.
            msg_lines += notes

        # Ask whether to immediately rerun all applicable steps with
        # the restored settings.
        msg_lines.append(
            "\nRun the whole pipeline now with these settings?"
        )

        reply = QMessageBox.question(
            self,
            "Session loaded",
            "\n".join(msg_lines),
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No
        )

        if reply == QMessageBox.StandardButton.Yes:
            # Start the automatic run queue using the restored
            # parameters.
            self._on_run_all()

    def _browse_workdir(self):
        """
        Let the user choose a new working directory for future
        EDIFF outputs.

        Changing the working directory only affects where subsequent
        output files are written. It does not reset the pipeline,
        clear loaded input files, or discard already computed results.

        Returns
        -------
        None
            Updates ``self.workdir``, ``self.pipeline.workdir``,
            and the GUI label if the user selects a valid folder.
        """
        # Ask the user to choose an existing directory.
        # The dialog starts from the current working directory.
        path = QFileDialog.getExistingDirectory(
            self,
            "Select working directory",
            str(self.workdir)
        )

        # User cancelled the dialog.
        if not path:
            return

        new_workdir = Path(path)

        # Nothing changed; avoid unnecessary GUI updates.
        if new_workdir == self.workdir:
            return
        try:
            # Ensure the selected directory exists.
            # QFileDialog.getExistingDirectory normally returns
            # an existing folder, but mkdir is still useful if
            # the path is later supplied differently or if parent
            # handling changes.
            new_workdir.mkdir(exist_ok=True, parents=True)

        except Exception as e:
            # The selected path cannot be used, for example because
            # of permissions.
            QMessageBox.warning(
                self,
                "EDIFF",
                f"Could not use that folder as the working directory:\n{e}"
            )
            return

        # Redirect future output files to the new folder.
        # Important:
        # This does NOT recreate the pipeline and does NOT clear
        # progress. Loaded CIF files, ELD image, detected center,
        # generated curves, etc. remain as they are. Only the destination
        # for future saved outputs changes.

        self.workdir = new_workdir
        self.pipeline.workdir = self.workdir

        # Refresh the shortened path label in Step 1.
        self._update_workdir_label()

        # Confirm the change to the user.
        self.statusBar().showMessage(
            f"Working directory changed to: {self.workdir}"
        )

    def _browse_cif(self):
        """
        Let the user select one or more structure files and add
        them to the session.

        The selected files are appended to the existing structure-file
        list rather than replacing it. The GUI currently limits
        the number of selected structure files to three, which is useful
        for keeping multi-phase comparison and calibration manageable.

        Supported file types include CIF/MCIF, VASP/POSCAR/CONTCAR,
        JSON, CSSR, and NC structure files.

        Returns
        -------
        None
            Updates ``self.cif_paths``, refreshes the file list, stores
            the last used folder in the config, and updates workflow
            availability.
        """
        # Start the file dialog in the last directory used for structure
        # files, if it is known.
        last = self._cfg.get("last_cif_dir", "")

        # Allow selecting multiple structure files at once.
        paths, _ = QFileDialog.getOpenFileNames(
            self,
            "Select structure files",
            last,
            "Structure files (*.cif *.mcif *.vasp POSCAR CONTCAR *.json *.cssr *.nc);;"
            "CIF files (*.cif *.mcif);;"
            "VASP files (*.vasp POSCAR CONTCAR);;"
            "All files (*.*)"
        )

        # User cancelled the dialog.
        if not paths:
            return

        # Add only files that are not already selected.
        # This prevents duplicate phases from appearing in the internal
        # list and UI.
        new_paths = [p for p in paths if p not in self.cif_paths]

        # If the user selected only files that were already present,
        # there is nothing to update.
        if not new_paths:
            return

        # The GUI supports at most three structure files.
        # This limit likely corresponds to the number of phase/reference
        # controls supported elsewhere in the application.
        room = 3 - len(self.cif_paths)

        if room <= 0:
            QMessageBox.warning(
                self,
                "EDIFF",
                "Maximum 3 structure files already selected.\n"
                "Remove one before adding another."
            )
            return

        if len(new_paths) > room:
            # Add only as many files as still fit, and tell the user
            # that the remaining selected files were ignored.
            QMessageBox.warning(
                self,
                "EDIFF",
                f"Maximum 3 structure files. Added {room} of the "
                f"{len(new_paths)} new file(s) selected."
            )
            new_paths = new_paths[:room]

        # Append the newly accepted files to the existing selection.
        self.cif_paths = self.cif_paths + new_paths

        # Remember the directory of the first newly added file
        # for the next browse.
        self._cfg["last_cif_dir"] = str(Path(new_paths[0]).parent)
        save_config(self._cfg)

        # Update the visible Step-1 CIF list.
        self._refresh_cif_list()

        # Re-evaluate which workflow steps can now run.
        self._check_files_ready()

    def _image_problems(self, path):
        """Check whether a chosen file can really be used as a diffraction image.

        The file dialog now lets any file be picked, but ediff's image reader
        has rules of its own, and several kinds of file look fine yet fail
        later with a confusing error - or, worse, load as the wrong data. This
        runs the same reader ediff will use and reports, in plain words, what
        is wrong or worth knowing, before the file is accepted.

        Blocking problems (the file would not work):

        - the file type is not an image format the reader supports;
        - the pixel type is not 8-bit or 16-bit unsigned integer (ediff's own
          rule - float and 32-bit images, for example, are rejected);
        - it is read as something other than a single 2-D image (such as a
          GIF, which comes back as a 4-D array).

        Notes (the file works, but the user should know):

        - it holds several pages/frames (a multi-page TIFF): ediff reads only
          the first one, silently;
        - JPEG is lossy, so pixel intensities are altered;
        - a colour image whose channels differ: only the first (red) channel
          is used.

        If the reader cannot be called at all (a different ediff version),
        nothing is reported, and the load itself decides.

        Parameters
        ----------
        path : str
            Path of the chosen file.

        Returns
        -------
        (problems, notes) : (list of str, list of str)
            Plain-language sentences; both empty when the file looks fine.
        """
        import numpy as np
        problems, notes = [], []
        try:
            arr = self.pipeline.ed.io.Data2D.read(path)
        except AttributeError:
            return problems, notes
        except Exception as e:
            text = str(e)
            if "8bit and 16bit" in text:
                problems.append(
                    "its pixel type is not supported - EDIFF accepts only "
                    "8-bit and 16-bit integer images (float and 32-bit "
                    "images, for example, cannot be used)")
            elif "backend" in text:
                problems.append("this file type is not a supported image format")
            else:
                problems.append(f"it could not be read as an image ({type(e).__name__}: {text})")
            return problems, notes
        if arr.ndim != 2:
            problems.append(
                f"it was read as {arr.ndim} dimensions {tuple(arr.shape)} "
                "instead of a single 2-D image (an animated GIF, or another "
                "multi-frame file?)")
        frames, mode = 1, ""
        try:
            from PIL import Image
            with Image.open(path) as im:
                frames, mode = getattr(im, "n_frames", 1), im.mode
            if mode in ("RGB", "RGBA"):
                full = np.asarray(Image.open(path))
                if full.ndim == 3 and not (
                        np.array_equal(full[..., 0], full[..., 1])
                        and np.array_equal(full[..., 0], full[..., 2])):
                    notes.append(
                        "it is a colour image whose channels differ - only "
                        "the first (red) channel is used")
        except Exception:
            pass
        if arr.ndim == 2 and frames > 1:
            notes.append(
                f"the file holds {frames} pages/frames - only the first one "
                "is used")
        if Path(path).suffix.lower() in (".jpg", ".jpeg", ".jpe"):
            notes.append(
                "JPEG uses lossy compression, which changes pixel intensities - "
                "PNG or TIFF is safer for quantitative analysis")
        return problems, notes

    def _browse_eld(self):
        """
        Let the user select an ELD image file and add it
        to the session.

        The selected diffraction image becomes ``self.eld_path``.
        Unlike structure files, only one ELD image is stored
        at a time. After selection, the function pdates the file
        list, main-window title, config, and workflow state.

        Returns
        -------
        None
            Updates ``self.eld_path`` and refreshes the GUI/
            workflow availability.
        """
        # Start from the last directory used for ELD/image files,
        # if available.
        last = self._cfg.get("last_eld_dir", "")

        # Select one diffraction-pattern image.
        # The first filter lists the common image formats ediff's reader
        # handles (each one tested); "All files" lets ANY file be picked.
        # _image_problems() then explains, before the file is accepted, if a
        # chosen file cannot be used (unsupported type or pixel type, GIF,
        # ...) or deserves a warning (JPEG, colour, multi-page TIFF).
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Select ELD image",
            last,
            "Images (*.png *.tif *.tiff *.bmp *.jpg *.jpeg *.webp "
            "*.ppm *.pgm *.tga *.pcx *.jp2);;All files (*)"
        )

        # User cancelled the dialog.
        if not path: return

        # Say so, in plain words, if this file cannot be used (or needs a
        # warning) - and let the user decide. Default is No for a file that
        # would not work, Yes for one that merely needs a note.
        problems, notes = self._image_problems(path)
        if problems or notes:
            text = f"{Path(path).name}\n\n"
            if problems:
                text += ("This file may not be usable as a diffraction image:\n- "
                         + "\n- ".join(problems) + "\n\n")
            if notes:
                text += "".join(f"Note: {n}.\n" for n in notes) + "\n"
            yes, no = (QMessageBox.StandardButton.Yes,
                       QMessageBox.StandardButton.No)
            reply = QMessageBox.question(
                self, "Check the image file",
                text + ("Use it anyway?" if problems else "Use this file?"),
                yes | no, no if problems else yes)
            if reply != yes:
                return

        # Store the selected diffraction image path.
        self.eld_path = path

        # Remember this directory for the next ELD browse action.
        self._cfg["last_eld_dir"] = str(Path(path).parent)
        save_config(self._cfg)

        # Update the visible Step-1 ELD list/label.
        self._refresh_eld_list()

        # Show the current image name in the application title.
        self.setWindowTitle(f"EDIFF - {Path(path).name}")

        # Re-evaluate which workflow steps can now run.
        self._check_files_ready()

    def _merge_into_full_pipeline(self):
        """
        Merge partial CIF-only or ELD-only progress into
        the full workflow.

        This function is called when both structure files
        and an ELD image become available, but not necessarily
        at the same time. It handles two common cases:

        1. CIF files were loaded first
           The user may have already calculated XRD. When the ELD
           image is added later, the image-processing chain should
           unlock without discarding XRD.

        2. ELD image was loaded first
           The user may have already completed ellipse correction,
           center detection, 1D conversion, or background subtraction.
           When CIF files are added later, XRD and calibration should
           unlock without discarding ELD-side progress.

        The function preserves trusted completed steps and activates
        the earliest step that is still required.

        Returns
        -------
        None
            Updates pipeline image loading, StepBlock states, top-bar
            indicators, and the status bar in place.
        """
        # Detect whether the ELD image itself is being replaced with a
        # genuinely different one (not just loaded for the first time),
        # BEFORE reloading it into the pipeline below. This matters
        # because ellipse correction, center detection, 1D conversion,
        # and background subtraction are all specific to the actual
        # image data, not just to "an ELD image existing" -- unlike
        # XRD (step 2), which only depends on the CIF and is already
        # correctly re-verified against self.pipeline.xrd_list further
        # down. Carrying forward steps 3-7's DONE state across a
        # genuinely different image falsely tells the GUI their
        # results are still valid, while the pipeline's own objects
        # for those steps have actually just been reset by the reload
        # below -- producing a frozen workflow: steps stay marked
        # done/locked so they can't be re-run, while anything that
        # unlocks past them ends up using missing or stale data from
        # the previous image.

        eld_image_replaced = (
                self._loaded_eld_path is not None
                and self.eld_path
                and self._loaded_eld_path != self.eld_path
        )

        if self.eld_path and (
                self.pipeline.eld is None
                or self._loaded_eld_path != self.eld_path
        ):
            try:
                self.pipeline.load_eld(self.eld_path)
                self._loaded_eld_path = self.eld_path

            except Exception as e:
                self.statusBar().showMessage(f"Failed to load ELD: {e}")
                return

        # The image was just replaced with a different one: none of
        # the ELD-dependent steps' previous results are valid for it
        # anymore. Reset their GUI state now, before the trusted-chain
        # logic below runs, so it correctly treats them as not-yet-done
        # for this image instead of carrying forward stale progress
        # from the old one.
        if eld_image_replaced:
            for n in (3, 4, 5, 6, 7):
                if self._blocks[n]._state == StepBlock.DONE:
                    self._blocks[n].set_state(StepBlock.LOCKED)
            # Also clear the "raw preview already shown" flag, so the
            # new image's raw pattern gets previewed again in the
            # Ellipse tab, exactly like a brand-new session would.


        # Collect steps currently marked DONE by the GUI.
        # Step 1 is added explicitly because file loading has succeeded.
        done = set(
            n for n in range(2, 8)
            if self._blocks[n]._state == StepBlock.DONE
        )
        done.add(1)

        # Do not trust the visual DONE state of XRD unless there is actual
        # XRD data in the pipeline. This protects against stale GUI state.
        if 2 in done and not self.pipeline.xrd_list:
            done.discard(2)

        # Validate the sequential ELD chain.
        # Steps 3-6 depend on each other:
        # ellipse -> center -> 1D integration -> background subtraction.
        # Therefore, if one earlier step is missing, later ELD-side DONE
        # states are not trusted even if the GUI somehow says they are done.
        trusted = {1}

        if 2 in done:
            trusted.add(2)

        for n in (3, 4, 5, 6):
            if n in done:
                trusted.add(n)
            else:
                # Stop at the first incomplete ELD-side step.
                break

        # Calibration is trusted only if every prerequisite step is trusted.
        if 7 in done and trusted == {1, 2, 3, 4, 5, 6}:
            trusted.add(7)

        # Find the first step in the full workflow that still needs to run.
        # If all steps are trusted, next_step remains None.
        next_step = next(
            (n for n in (2, 3, 4, 5, 6, 7) if n not in trusted),
            None
        )

        # Apply the merged state to all StepBlocks.
        for n in range(1, 8):
            if n in trusted:
                self._blocks[n].set_state(StepBlock.DONE)
            elif n == next_step:
                self._blocks[n].set_state(StepBlock.ACTIVE)
            else:
                self._blocks[n].set_state(StepBlock.LOCKED)

        # Human-readable message for the next recommended action.
        messages = {
            2: "Continue with XRD calculation.",
            3: "Correct ellipse or proceed to center detection.",
            4: "Continue with center detection.",
            5: "Continue by converting to 1D.",
            6: "Continue with background subtraction.",
            7: "Continue with calibration.",
        }
        msg = messages.get(next_step, "All steps already complete.")

        # Update the top progress indicators to reflect the merged state.
        self._update_topbar_indicators_partial(
            done=sorted(trusted),
            active=(next_step or 0),
            locked=[
                n for n in range(1, 8)
                if n not in trusted and n != next_step
            ],
        )

        self.statusBar().showMessage(f"Files loaded — {msg}")

    def _check_files_ready(self):
        """
        Update workflow availability after input files are added
        or removed.

        This function is the central gatekeeper for Step-1 file
        state. It checks whether structure files, an ELD image,
        both, or neither are currently loaded, and then configures
        the workflow accordingly.

        Supported modes
        ---------------
        No files
            Reset the pipeline and lock everything after Step 1.

        CIF only
            Enable XRD calculation. Keep ELD-dependent steps locked.

        ELD only
            Load the ELD image immediately and enable the image
            -processing chain starting at ellipse correction. Keep
            XRD and calibration locked.

        CIF + ELD
            Enable the full workflow. If any progress already exists,
            merge it into the full pipeline instead of resetting from
            the beginning.

        Returns
        -------
        None
            Updates pipeline state, StepBlock states, tabs, result
            cards, top-bar indicators, and status-bar messages in place.
        """
        # Determine current input availability.
        has_cif = bool(self.cif_paths)
        has_eld = bool(self.eld_path)

        if not has_cif and not has_eld:
            # No input files remain.
            # Reset the computational pipeline because later results
            # no longer have valid inputs behind them.
            self.pipeline = EdiffPipeline(self.workdir)
            self._loaded_eld_path = None

            # Clear all result cards because their contents correspond
            # to files that are no longer loaded.
            for card in [
                self._card_xrd,
                self._card_ellipse,
                self._card_eld2d,
                self._card_1d,
                self._card_bkg,
                self._card_comp
            ]:
                card.clear_runs()

            # Disable all result tabs.
            for i in range(self._tabs.count()):
                self._tabs.setTabEnabled(i, False)

            # Return the workflow to Step 1.
            self._go_to_step(1)
            self.statusBar().showMessage("Load files to begin.")
            return

        # At least one input file exists, so Step 1 is complete.
        self._blocks[1].set_state(StepBlock.DONE)

        # Step 1 stays editable even after files are loaded,
        # so users can add/remove files or change the working
        # directory later.
        for w in self._blocks[1]._children: w.setEnabled(True)

        if has_cif and has_eld:
            # Full pipeline mode.
            # Either file type could have been added after partial
            # progress was already made. Preserve existing valid
            # progress if possible.
            any_progress = any(
                self._blocks[n]._state == StepBlock.DONE
                for n in range(2, 8)
            )

            if any_progress:
                # Merge existing CIF-side or ELD-side progress
                # into the full workflow state.
                self._merge_into_full_pipeline()



            else:

                # Fresh full workflow. Ellipse correction (step 3)
                # doesn't depend on XRD at all -- it only needs the ELD
                # image, already loaded here. A plain _go_to_step(2)
                # locks every other step, including step 3, forcing
                # users through XRD first even though nothing about
                # ellipse correction requires it. Steps 2 and 3 are
                # made usable simultaneously instead, exactly like
                # ELD-only mode already treats step 3 as immediately
                # available.

                self._blocks[2].set_state(StepBlock.ACTIVE)
                self._blocks[3].set_state(StepBlock.ACTIVE)

                for n in (4, 5, 6, 7):
                    self._blocks[n].set_state(StepBlock.LOCKED)

                for w in self._blocks[1]._children:
                    w.setEnabled(True)

                self._update_topbar_indicators_partial(
                    done=[1],
                    active=2,
                    locked=[3, 4, 5, 6, 7]
                )

                try:
                    self.pipeline.load_eld(self.eld_path)
                    self._loaded_eld_path = self.eld_path

                except Exception as e:
                    self.statusBar().showMessage(f"Failed to load ELD: {e}")
                    return

                self.statusBar().showMessage(
                    "Files loaded — calculate XRD, or correct ellipse / "
                    "detect center to continue with the diffraction pattern."

                )

        elif has_cif and not has_eld:
            # CIF - only  mode.
            # Only XRD can run because the image-processing chain requires
            # an ELDdiffraction image.
            self._go_to_step(2)

            # Explicitly lock all ELD-dependent and calibration steps.
            for n in range(3, 8):
                self._blocks[n].set_state(StepBlock.LOCKED)

            self._update_topbar_indicators_partial(
                done=[1],
                active=2,
                locked=[3,4,5,6,7]
            )

            self.statusBar().showMessage(
                "CIF loaded - calculate XRD. Load ELD image to unlock center detection."
            )

        elif has_eld and not has_cif:
            # ELD-only mode

            # Steps 3-6 can run because they only need the diffraction image.
            # Step 2 is locked because XRD needs structure files.
            # Step 7 is locked because calibration needs both structure and ELD
            # data. Load the ELD image immediately.
            # In ELD-only mode there is no XRD step that would otherwise trigger
            # image loading.
            try:
                self.pipeline.load_eld(self.eld_path)
                self._loaded_eld_path = self.eld_path

            except Exception as e:
                self.statusBar().showMessage(f"Failed to load ELD: {e}")
                return

            # Step 1 is complete, XRD is unavailable, and ellipse correction
            # is the first active ELD-side step.
            self._blocks[1].set_state(StepBlock.DONE)
            self._blocks[2].set_state(StepBlock.LOCKED)
            self._blocks[3].set_state(StepBlock.ACTIVE)

            # Later image-processing steps and calibration remain locked
            # until the workflow reaches them.
            for n in range(4, 8):
                self._blocks[n].set_state(StepBlock.LOCKED)

            # Keep Step 1 editable.
            for w in self._blocks[1]._children:
                w.setEnabled(True)

            self._update_topbar_indicators_partial(
                done=[1],
                active=3,
                locked=[2, 4, 5, 6, 7]
            )

            self.statusBar().showMessage(
                "ELD loaded. Inspect the 2D pattern — correct ellipse if rings look "
                "oval, otherwise skip to center detection. "
                "Load CIF files to enable XRD and calibration.")


    def _update_topbar_indicators_partial(self, done, active, locked):
        """
        Update the top workflow indicators for a partially available
        pipeline.

        This helper is used when the full linear workflow is not
        available or when only a subset of steps should be marked
        as done/active/locked. For example, it is used in CIF-only
        mode, ELD-only mode, and when merging partial progress
        into the full pipeline.

        Parameters
        ----------
        done : iterable of int
            Step numbers that should be shown as completed.

        active : int
            Step number that should be shown as currently active.
            Use 0 or None when no step should be active.

        locked : iterable of int
            Step numbers that should be shown as locked. The argument
            is currently not used directly because every non-done and
            non-active step is styled as locked, but it documents the
            intended state passed by the caller.

        Returns
        -------
        None
            Updates the badge text/style and label style in place.
        """
        # Iterate over the visual top-bar indicators.
        # Each indicator consists of a badge and a text label.
        for i, (badge, lbl) in enumerate(self._step_indicators):
            # Step numbers are 1-based, while enumerate() is 0-based.
            n = i + 1

            if n in done:
                # Completed step:
                # green filled badge and green label. You currently
                # use "v" as the completion mark. If you want a real
                # checkmark, you can replace it with "✓".
                badge.setText("v")
                badge.setStyleSheet(
                    "background:#34c759; color:#ffffff; border-radius:10px;"
                    "font-size:10px; font-weight:700; border:none;"
                )
                lbl.setStyleSheet(
                    "font-size:12px; color:#34c759; font-weight:600;"
                )

            elif n == active:
                # Active step:
                # blue filled badge showing the step number.
                badge.setText(str(n))
                badge.setStyleSheet(
                    "background:#0071e3; color:#ffffff; border-radius:10px;"
                    "font-size:10px; font-weight:700; border:none;"
                )
                lbl.setStyleSheet(
                    "font-size:12px; color:#0071e3; font-weight:600;"
                )

            else:
                # Locked or inactive step:
                # grey outlined badge and muted label.
                # Note:
                #   The `locked` argument is not explicitly checked here. Anything
                #   that is neither done nor active is styled as locked/inactive.
                badge.setText(str(n))
                badge.setStyleSheet(
                    "background:#f5f5f7; color:#c7c7cc; border-radius:10px;"
                    "border:1.5px solid #e5e5ea; font-size:10px; font-weight:700;"
                )
                lbl.setStyleSheet(
                    "font-size:12px; color:#c7c7cc; font-weight:400;"
                )

    # ── Step 2 ────────────────────────────────────────────────────────────────
    def _on_xrd(self):
        """
        Start theoretical XRD calculation for the selected structure files.

        This function performs the Step-2 workflow:

        1. Check whether a previous XRD calculation is still running.
        2. Warn the user if the chosen wavelength and 2θ range imply
            an extremely small minimum d-spacing, which can make
            reflection enumeration very slow.
        3. Switch the GUI into a running/spinner state.
        4. Capture all Qt widget values on the GUI thread.
        5. Run the XRD calculation in a background worker.
        6. Plot the calculated diffractograms.
        7. Display the result buffers in the XRD result card.
        8. Let the user accept one run, or auto-accept the first run.

        Important threading detail
        --------------------------
        Qt widgets should not be read from the background worker thread.
        Therefore, all widget values are copied into local variables before
        the nested ``do()`` function is executed by ``run_in_thread``.

        Returns
        -------
        None
            Starts an asynchronous XRD calculation and updates the GUI
            state.
        """
        # Prevent the user from launching another XRD calculation while
        # the previous one is still finishing in the background.
        if self._step_busy(2):
            self.statusBar().showMessage(
                "Previous XRD calculation is still finishing — please wait a moment."
            )
            return


        # Sanity check for expensive XRD settings
        # A very short wavelength and a large 2θ range correspond to a very
        # small minimum d-spacing through Bragg's law:
        #     d = λ / (2 sin θ)
        # Very small d_min means many high-index reflections may need to be
        # enumerated. For typical structures, this can become extremely slow.
        wl_check = self._xrd_wl.value()
        tt_max_check = self._xrd_ttmax.value()
        theta_max_rad = math.radians(tt_max_check / 2)

        if theta_max_rad > 0:
            d_min = wl_check / (2 * math.sin(theta_max_rad))

            if d_min < 0.3:
                # Warn the user before starting a potentially huge calculation.
                reply = QMessageBox.warning(
                    self, "This may take a very long time",
                    f"With this wavelength ({wl_check:.4f} Å) and 2θ max "
                    f"({tt_max_check:.1f}°), the smallest accessible "
                    f"d-spacing is about {d_min:.3f} Å.\n\n"
                    "That's extremely fine for typical unit cells — it "
                    "requires checking a very large number of possible "
                    "reflections (potentially millions), and the "
                    "calculation may run for a very long time or appear "
                    "to hang.\n\n"
                    "Continue anyway, or Cancel to adjust the wavelength "
                    "or 2θ max first?",
                    QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.Cancel,
                    QMessageBox.StandardButton.Cancel
                )

                # Cancel means the user wants to edit the parameters first.
                if reply != QMessageBox.StandardButton.Yes:
                    return

        # Put Step 2 into running GUI state

        # Disable the XRD button immediately to avoid double-clicks.
        self._btn_xrd.setEnabled(False)

        # Hide Run, show Stop, enable global Abort,
        # and remember current step.
        self._start_step(2, self._btn_xrd, self._btn_stop2)

        # Show a spinner/progress placeholder in the XRD result card.
        self._card_xrd.show_spinner("Calculating theoretical diffractograms...")

        # Make the XRD tab visible and switch to it so the user
        # sees progress.
        self._tabs.setTabEnabled(self.TAB_XRD, True)
        self._tabs.setCurrentIndex(self.TAB_XRD)

        self.statusBar().showMessage("Calculating XRD...")

        # Capture GUI values on the main thread
        # These values are used inside the background calculation.
        # Reading them here avoids unsafe access to Qt widgets from
        # the worker thread and freezes the calculation settings
        # at the moment the user clicked XRD.
        wavelength   = self._xrd_wl.value()
        tt_range     = (self._xrd_ttmin.value(), self._xrd_ttmax.value())
        peak_sigma   = self._xrd_sigma.value()
        temp_factor  = self._xrd_temp.value()
        peak_profile = self._xrd_profile.currentText()

        # Reuse the previous q/intensity plot limits, if available.
        # XRD and comparison plots use the same physical axis meaning,
        # so preserving the zoom range makes reruns less disruptive.
        recalled_xlim, recalled_ylim = self._recall_axis("q_intensity")

        def do():
            """
            Worker-thread body for XRD calculation and plotting.

            Returns
            -------
            list
                Plot buffers returned by ``self.pipeline.plot_xrd``.
            """
            # Calculate theoretical diffractograms for the currently
            # selected structure files using the captured GUI parameters.
            self.pipeline.calculate_xrd(
                cif_paths=self.cif_paths,
                wavelength=wavelength,
                two_theta_range=tt_range,
                peak_profile_sigma=peak_sigma,
                temp_factors=temp_factor,
                peak_profile=peak_profile,
            )

            # If an ELD image is also selected, make sure it is loaded
            # into the pipeline. This is needed before later ELD-dependent
            # steps run.
            if self.eld_path is not None:
                # Only reload if there is no loaded image, or if the selected
                # raw ELD path changed since the last load.
                # Important:
                # Do NOT compare against pipeline._eld_path directly here.
                # Ellipse correction may change pipeline._eld_path to the
                # corrected image, and reloading based on that would discard
                # the correction.
                if (
                    self.pipeline.eld is None
                    or self._loaded_eld_path != self.eld_path
                ):
                    self.pipeline.load_eld(self.eld_path)

                # Track which original ELD path has been loaded from the GUI side.
                self._loaded_eld_path = self.eld_path

            # Return rendered XRD plot buffers using the remembered axis limits.
            return self.pipeline.plot_xrd(
                xlim=recalled_xlim,
                ylim=recalled_ylim
            )

        def done(bufs):
            """
            GUI-thread completion callback for the XRD worker.

            Parameters
            ----------
            bufs : list
                Plot buffers produced by ``do()`` and returned by
                ``self.pipeline.plot_xrd``.

            Returns
            -------
            None
                Displays XRD plots and either auto-accepts or waits for user
                acceptance of the run.
            """
            # Pair each plot buffer with a phase label.
            # Prefer the description from pipeline.xrd_list; fall back
            # to Phase N.

            items = [
                (self.pipeline.xrd_list[i].description or f"Phase {i+1}", b)
                 for i, b in enumerate(bufs)
            ]

            # Human-readable run label shown in the result card.
            xrd_label = (f"\u03bb={wavelength:.3f} "
                         f"\u03c3={peak_sigma:.3f} "
                         f"{peak_profile}")

            # Display the generated XRD plots in the result card.
            # auto_select=False here because auto-selection is handled
            # explicitly below after the panel snapshot is attached.
            panel = self._card_xrd.show_buffers(
                items,
                run_label=xrd_label,
                auto_select=False
            )

            # Restore any remembered axis limits into the card's axis controls.
            self._card_xrd.set_axis_fields(recalled_xlim, recalled_ylim)

            # Register a refresh callback so the displayed plot can be redrawn
            # with updated axis limits while preserving those limits in axis
            # memory.
            self._card_xrd.set_refresh_callback(
                self._wrap_refresh_with_memory(
                    "q_intensity",
                    lambda xlim, ylim: [
                        (
                            self.pipeline.xrd_list[i].description
                            or f"Phase {i+1}",
                            b,
                        )
                        for i, b in enumerate(
                            self.pipeline.plot_xrd(xlim=xlim, ylim=ylim)
                        )
                    ],
                )
            )

            # Snapshot the current xrd_list for this run.
            # This is important because the user may calculate several XRD
            # runs with different parameters. When they later click "Use this
            # run", the app needs to restore the XRD data corresponding
            # to the selected panel.
            import copy

            panel._xrd_snapshot = copy.copy(self.pipeline.xrd_list)

            # Restore Step-2 buttons from running state to normal state.
            self._finish_step(2, self._btn_xrd, self._btn_stop2)

            # Keep the user on the XRD tab after calculation.
            self._tabs.setCurrentIndex(self.TAB_XRD)

            # Avoid connecting the same signal multiple times after repeated
            # runs.
            try: self._card_xrd.run_selected.disconnect(
                self._on_xrd_run_selected
            )
            except Exception:
                pass

            self._card_xrd.run_selected.connect(self._on_xrd_run_selected)

            # Auto-accept every run, not just the first one. Previously,
            # re-running this step (e.g. to try different XRD
            # parameters) required an explicit "Use this run" click
            # before the workflow would advance -- easy to miss, since
            # that button lives in the plot area rather than the
            # sidebar where the still-locked next step actually is.
            # That produced a confusing "frozen" feeling: later steps
            # stayed locked with no obvious visible reason why.
            # Auto-accepting the newest run keeps the common case
            # (re-running because you want the improved result)
            # frictionless, while "Use this run" remains fully
            # available on every tab -- including older ones -- for
            # anyone who deliberately wants to go back to a previous
            # run instead.
            self._card_xrd.select_run(panel)

            # For later runs, wait for explicit user acceptance.
            self.statusBar().showMessage(
                "XRD calculated — click '✓ Use this run' to accept and continue, "
                "or recalculate with different parameters."
            )

        # Start the worker.
        # busy_key=2 lets _step_busy(2) know that Step 2 is still active,
        # even after Stop/Abort if the underlying calculation has not fully
        # returned yet.
        w = run_in_thread(
            self,
            do,
            on_done=done,
            on_error=lambda e: self._err(
                "XRD failed",
                e,
                card=self._card_xrd
            ),
            busy_key=2,
        )

        # Store the worker so Stop/Abort can cancel or suppress its
        # result later.
        self._active_worker_2 = w

    def _on_xrd_run_selected(self, panel):
        """
        Accept a displayed XRD run and advance the workflow.

        Each XRD result panel stores a snapshot of the ``pipeline.xrd_list``
        that produced it. When the user accepts a panel, this function
        restores that snapshot, marks Step 2 complete, and moves the workflow
        forward.

        If an ELD image is already loaded, the raw ELD preview is shown
        in the ellipse tab so the user can decide whether ellipse correction
        is needed. In CIF-only mode, no ELD preview is shown and
        ``_complete_step`` handles the appropriate workflow state.

        Parameters
        ----------
        panel : object
            Result-card panel selected by the user. Expected to optionally
            contain ``_xrd_snapshot``.

        Returns
        -------
        None
            Restores selected XRD data and advances workflow state.
        """
        # Disconnect this handler to avoid duplicate reactions if
        # the same signal is emitted again or if a new run reconnects
        # it later.
        try:
            self._card_xrd.run_selected.disconnect(
                self._on_xrd_run_selected
            )
        except Exception:
            pass

        # Restore the XRD data corresponding to the selected run.
        # This is necessary because pipeline.xrd_list may have changed
        # after additional recalculations.
        if hasattr(panel, '_xrd_snapshot') and panel._xrd_snapshot:
            self.pipeline.xrd_list = panel._xrd_snapshot

        # Mark Step 2 complete and advance/merge workflow state.
        self._complete_step(2)

        if self.eld_path:
            self.statusBar().showMessage(
                "XRD accepted. Click Preview to inspect the raw pattern, "
                "correct ellipse if rings look oval, or skip to center detection."
            )
            
        # If there is no ELD image, this is CIF-only mode.
        # _complete_step(2) has already set the appropriate message and locked
        # the ELD-dependent steps, so there is nothing else to do here.

    # ── Step 4 ────────────────────────────────────────────────────────────────
    def _on_center(self):
        """
        Start Step 4: detect or load the diffraction-pattern center.

        This method handles center detection for the currently loaded
        ELD image. Depending on the selected GUI option, the center
        can either be detected computationally, refined manually, or
        read from a previously saved center file.

        Workflow
        --------
        1. Check whether a previous center-detection worker is still
           finishing.
        2. Resolve the selected detection and refinement methods
           from the GUI.
        3. If the selected detection method reads from file, verify
           that a file was selected.
        4. If the selected method writes a center file, check whether
           the output file already exists and ask before overwriting
           it.
        5. Capture all GUI values on the main thread.
        6. Optionally show tutorial help for manual detection/refinement.
        7. Run center detection in a background worker.
        8. Display the resulting 2D diffractogram and center-detail
           plot.
        9. Store a snapshot of the detected center so selecting this
           run later restores the correct coordinates.

        Returns
        -------
        None
            Starts an asynchronous center-detection task and updates
            GUI state.
        """
        # Do not start a second center-detection job while the previous
        # one is still finishing in the background.
        if self._step_busy(4):
            self.statusBar().showMessage(
                "Previous center detection is still finishing — please wait a moment.")
            return

        # Translate user-visible combo-box labels into backend method
        # values.
        # DETECTION_OPTIONS and REFINEMENT_OPTIONS are expected to map
        # GUI text to the values accepted by self.pipeline.find_center().
        detection    = DETECTION_OPTIONS[self._ctr_det.currentText()]
        refinement   = REFINEMENT_OPTIONS[self._ctr_ref.currentText()]

        # In your app convention, detection is None when the center
        # should be read from a file rather than detected from the image.
        is_file      = detection is None

        # If the method needs an input file, make sure the user selected
        # one first.
        if is_file and not self._ctr_read_file_path:
            self.statusBar().showMessage(
                "Please browse for a saved center file first.")
            return

        if is_file:
            # Read center coordinates from an existing file.
            # No output file is written by this branch.
            in_file  = self._ctr_read_file_path
            out_file = None

        else:
            # Normal detection branch:
            # write detected center coordinates to an output file
            # in workdir.
            cfile    = self._ctr_file.text().strip() or "center.txt"
            in_file  = None
            out_file = str(self.workdir / cfile)

            # Ask before overwriting an existing center file.
            # This avoids silently replacing data from an earlier
            # run/session.
            if Path(out_file).exists():
                reply = QMessageBox.question(
                    self, "File already exists",
                    f"'{cfile}' already exists in the working directory "
                    "(possibly from an earlier session).\n\n"
                    "Overwrite it?",
                    QMessageBox.StandardButton.Yes
                    | QMessageBox.StandardButton.Cancel,
                )

                if reply != QMessageBox.StandardButton.Yes:
                    return

        # Put Step 4 into running state.
        self._btn_center.setEnabled(False)
        self._start_step(4, self._btn_center, self._btn_stop4)

        # Capture GUI values on the main thread
        # These values are used by the worker. Capturing them here
        # avoids reading Qt widgets from a background thread and keeps
        # the run reproducible even if the user edits GUI fields while
        # the worker is running.
        icut         = self._icut.value()
        cmap         = self._cmap.currentText()

        # Size/scale for the center-detail view.
        csq          = self._ctr_detail_csq.value()

        # Method-specific center-detection parameters.
        cintensity   = self._ctr_cintensity.value()
        csquare      = self._ctr_csq.value()

        # rtype is taken from the first character of the combo text.
        # Example: "1 - something" -> 1.
        rtype        = int(self._ctr_rtype.currentText()[0])

        # Convert string combo-box choices like "True ..." / "False ..."
        # to bools.
        heq          = self._ctr_heq.currentText().startswith("True")
        sobel        = self._ctr_sobel.currentText().startswith("True")
        masking      = self._ctr_masking.currentText().startswith("True")
        normalize_bg = self._ctr_normalize_bg.currentText().startswith("True")
        live_plot = self._ctr_live_plot.currentText().startswith("True")

        # Reuse axis limits from previous 2D image plots.
        # Ellipse and center detection both use pixel-coordinate axes.
        recalled_xlim, recalled_ylim = self._recall_axis("image_2d")

        # Show tutorial panels for interactive modes, if Tutorial
        # mode is enabled
        if detection == "manual":
            self._maybe_show_tutorial("manual_detection")
        elif refinement == "manual":
            self._maybe_show_tutorial("manual_refinement")

        # Show progress in the Center/ELD 2D card.
        self._card_eld2d.show_spinner("Finding center...")
        self._tabs.setTabEnabled(self.TAB_ELD2D, True)
        self._tabs.setCurrentIndex(self.TAB_ELD2D)
        self.statusBar().showMessage("Finding center...")

        def do():
            """
            Worker-thread body for center detection.

            Returns
            -------
            tuple
                ``(buf2d, bufc)``, where ``buf2d`` is the full 2D
                diffractogram plot buffer and ``bufc`` is the center
                -detail plot buffer.
            """
            # Run the backend center-finding routine.
            self.pipeline.find_center(
                detection=detection,
                refinement=refinement,
                icut=icut,
                in_file=in_file,
                out_file=out_file,
                cmap=cmap,
                csquare=csquare,
                cintensity=cintensity,
                rtype=rtype,
                heq=heq,
                sobel=sobel,
                masking=masking,
                normalize_bg=normalize_bg,
                live_plot=live_plot,
            )

            # Return plots for display in the result card.
            return (
                self.pipeline.plot_eld_2d(
                    icut=icut,
                    cmap=cmap,
                    xlim=recalled_xlim,
                    ylim=recalled_ylim
                ),
                self.pipeline.plot_center(
                    csquare=csq,
                    cmap=cmap,
                    xlim=recalled_xlim,
                    ylim=recalled_ylim
                ),
            )

        def done(result):
            """
            GUI-thread completion callback for center detection.

            Parameters
            ----------
            result : tuple
                Tuple returned by ``do()``: ``(buf2d, bufc)``.

            Returns
            -------
            None
                Displays center-detection results and prepares run
                selection.
            """
            buf2d, bufc = result

            # Always show the full 2D diffractogram.
            items = [("2D diffractogram", buf2d)]

            # Show the center-detail plot only if the backend produced one.
            if bufc:
                items.append(("Center detail", bufc))

            # Build a compact run label from the current GUI method names.
            det_short = \
                self._ctr_det.currentText().split("-")[0].strip().split()[0]
            ref_short = \
                self._ctr_ref.currentText().split("-")[0].strip().split()[0]
            ctr_label = \
                f"{det_short} / {ref_short}  icut={self._icut.value()}"

            # If this is the first center run, there is nothing to compare
            # against, so the run is auto-accepted after display.
            is_only_run = len(self._card_eld2d._panels) == 0

            # Display the result buffers without immediately selecting them.
            # Selection is handled explicitly below after snapshots/callbacks
            # exist.
            panel = self._card_eld2d.show_buffers(
                items,
                run_label=ctr_label,
                auto_select=False
            )

            # Restore remembered axis limits into the card's axis controls.
            self._card_eld2d.set_axis_fields(recalled_xlim, recalled_ylim)

            def _refresh_center(xlim, ylim, icut=icut, cmap=cmap, csq=csq):
                """
                Replot the center-detection result for new axis limits.

                This callback is used by the result card when the user changes
                plot limits. Default arguments freeze the parameters from this
                run.

                Parameters
                ----------
                xlim, ylim : tuple or None
                   Axis limits requested by the result card.

                Returns
                -------
                list of tuple
                   List of ``(title, buffer)`` plot items.
                """
                buf2d = self.pipeline.plot_eld_2d(
                    icut=icut,
                    cmap=cmap,
                    xlim=xlim,
                    ylim=ylim
                )

                bufc = self.pipeline.plot_center(
                    csquare=csq,
                    cmap=cmap,
                    xlim=xlim,
                    ylim=ylim
                )

                out = [("2D diffractogram", buf2d)]
                if bufc:
                    out.append(("Center detail", bufc))
                return out

            # Store the refresh callback and remember image-axis
            # limits globally.
            self._card_eld2d.set_refresh_callback(
                self._wrap_refresh_with_memory("image_2d", _refresh_center)
            )

            # Snapshot center coordinates so that accepting this panel
            # later restores the center from this specific run.
            try:
                c = self.pipeline.eld.center
                panel._center_snapshot = (c.x, c.y) if c is not None else None
            except Exception:
                panel._center_snapshot = None

            # Restore Step-4 button state.
            self._finish_step(4, self._btn_center, self._btn_stop4)

            # Avoid duplicate signal connections after repeated runs.
            try:
                self._card_eld2d.run_selected.disconnect(
                    self._on_center_run_selected
                )
            except Exception:
                pass

            self._card_eld2d.run_selected.connect(
                self._on_center_run_selected
            )

            # Auto-accept every run, not just the first one (see the
            # matching comment in the XRD handler for the full
            # rationale). "Use this run" remains available on every
            # tab for anyone who wants to revert to an earlier attempt.
            self._card_eld2d.select_run(panel)

            # Otherwise wait for explicit user acceptance.
            self.statusBar().showMessage(
                "Center detected — click '✓ Use this run' to accept and continue, "
                "or try a different method."
            )

        # Start the worker thread.
        w = run_in_thread(
            self,
            do,
            on_done=done,
            on_error=lambda e: self._err(
                "Center detection failed",
                e,
                card=self._card_eld2d
            ),
            busy_key=4,
        )

        # Keep worker reference for Stop/Abort.
        self._active_worker_4 = w


    def _on_center_run_selected(self, panel):
        """
        Accept a center-detection run and advance the
        workflow.

        When several center-detection runs exist, each panel
        may contain a snapshot of the center coordinates produced
        by that run. Accepting a panel restores those coordinates
        into the pipeline before Step 4 is marked complete.

        Parameters
        ----------
        panel : object
            Result-card panel selected by the user. Expected
            to optionally contain ``_center_snapshot``.

        Returns
        -------
        None
            Restores selected center coordinates and advances
            to Step 5.
        """
        # Disconnect to avoid duplicate calls if the signal is
        # emitted again later.
        try:
            self._card_eld2d.run_selected.disconnect(
                self._on_center_run_selected
            )
        except Exception:
            pass

        # Restore the center coordinates from this specific run
        snap = getattr(panel, '_center_snapshot', None)


        if snap is not None:
            try:
                # Some backend code may use x/y, while other parts
                # may use x2/y2. Set both pairs to keep the accepted
                # center consistent everywhere.
                self.pipeline.eld.center.x  = snap[0]
                self.pipeline.eld.center.y  = snap[1]
                self.pipeline.eld.center.x2 = snap[0]
                self.pipeline.eld.center.y2 = snap[1]
            except Exception:
                # If the center object is missing or shaped unexpectedly,
                # do not crash the GUI. The current pipeline state remains
                # unchanged.
                pass

        # Mark Step 4 complete and activate the next appropriate step.
        self._complete_step(4)

        self.statusBar().showMessage(
            "Center accepted — convert to 1D to continue."
        )

    # ── Step 3 helpers ────────────────────────────────────────────────────────
    def _show_raw_eld_preview(self, switch_tab=True):
        """
        Show the currently loaded raw 2D diffractogram in the dedicated
        Preview tab, which sits before XRD/Ellipse in the tab bar so
        the user can inspect the raw pattern (and adjust icut/colormap,
        via Step 3's icut/colormap fields) immediately after loading an
        image -- before calculating XRD or running ellipse correction.

        This used to write into the Ellipse tab directly, which meant
        the raw pattern and the actual ellipse-correction results
        shared the same space. Separating them means the Ellipse tab
        now only ever shows real correction output (once "Correct
        ellipse" actually runs), while this Preview tab is always the
        one place to check the raw pattern at any time.

        Parameters
        ----------
        switch_tab : bool, optional
            If True, switch the right-side results panel to the
            Preview tab after starting the preview. If False, load the
            preview silently in the background without changing the
            active tab.

        Returns
        -------
        None
            Starts an asynchronous preview plot and displays it in the
            Preview card.
        """
        # No image has been loaded into the pipeline yet, so there is
        # nothing to preview.
        if self.pipeline.eld is None:
            return

        # Capture display settings from the GUI.
        icut = self._icut.value()
        cmap = self._cmap.currentText()

        # Show a temporary progress state in the Preview card.
        self._card_preview.show_spinner("Loading 2D pattern...")

        # Enable the Preview tab; optionally switch to it.
        self._tabs.setTabEnabled(self.TAB_PREVIEW, True)

        if switch_tab:
            self._tabs.setCurrentIndex(self.TAB_PREVIEW)


        def do():
            """
            Worker-thread body for raw ELD preview plotting.

            Returns
            -------
            object
                Plot buffer returned by ``pipeline.plot_eld_2d``.
            """
            return self.pipeline.plot_eld_2d(icut=icut, cmap=cmap)

        def done(buf):
            """
            GUI-thread completion callback for raw ELD preview.

            Parameters
            ----------
            buf : object
                Plot buffer produced by ``do()``.

            Returns
            -------
            None
                Displays the preview in the Ellipse card.
            """
            self._card_preview.show_buffers(
                [("Raw 2D pattern — inspect for elliptical rings", buf)],
                run_label="raw preview",
                auto_select=True,
                show_use_button=False
            )

        # Run preview plotting asynchronously so the GUI remains
        # responsive. If preview plotting fails, simply clear
        # the card instead of raising.
        run_in_thread(
            self,
            do,
            on_done=done,
            on_error=lambda e: self._card_preview.clear_runs())

    # ── Step 3 ────────────────────────────────────────────────────────────────
    def _on_ellipse(self):
        """
        Start Step 3: interactive ellipse-distortion correction.

        This method runs the ellipse-correction workflow for
        the loaded ELD image. The backend opens an interactive
        window where the user clicks four points on a distorted
        diffraction ring. The image is then transformed so those
        points lie on a circular geometry.

        Workflow
        --------
        1. Check whether a previous correction is still finishing.
        2. Put the GUI into running state.
        3. Capture display parameters and remembered axis limits.
        4. Reset the pipeline to the raw ELD image.
        5. Plot the raw image.
        6. Run interactive ellipse correction.
        7. Plot the corrected image.
        8. Display before/after buffers.
        9. Let the user decide whether to use the corrected image
           or skip it.

        Returns
        -------
        None
            Starts an asynchronous interactive ellipse-correction
            task.
        """
        # Do not start a second ellipse-correction job while one
        # is still finishing.
        if self._step_busy(3):
            self.statusBar().showMessage(
                "Previous ellipse correction is still finishing — please wait a moment."
            )
            return

        # Prevent duplicate clicks and hide the difference-view button
        # while a new correction run is being computed.
        self._btn_ellipse.setEnabled(False)
        self._btn_show_diff.setVisible(False)

        # Put Step 3 into running state.
        self._start_step(3, self._btn_ellipse, self._btn_stop3)

        # Capture display settings.
        icut = self._icut.value()
        cmap = self._cmap.currentText()

        # Show tutorial if Tutorial mode is enabled.
        self._maybe_show_tutorial("ellipse_correction")

        # Reuse image-axis limits from earlier 2D plots.
        # Ellipse and center detection both operate in pixel coordinates.
        recalled_xlim, recalled_ylim = self._recall_axis("image_2d")

        # Show progress and guide the user toward the external
        # interactive window.
        self._card_ellipse.show_spinner(
            "Fitting elliptical distortion — follow the window..."
        )
        self._tabs.setTabEnabled(self.TAB_ELLIPSE, True)
        self._tabs.setCurrentIndex(self.TAB_ELLIPSE)

        self.statusBar().showMessage(
            "Ellipse window open — click 4 points on a diffraction ring, then close."
        )


        def do():
            """
            Worker-thread body for ellipse correction.

            Returns
            -------
            tuple
                ``(raw_buf, corrected_buf)`` containing plot buffers
                before and after ellipse correction.
            """
            # Always start ellipse correction from the raw image.
            # This avoids repeatedly correcting an already corrected
            # image.
            self.pipeline.reset_to_raw_eld()
            raw_buf = self.pipeline.plot_eld_2d(
                icut=icut,
                cmap=cmap,
                xlim=recalled_xlim,
                ylim=recalled_ylim
            )

            # Launch the interactive correction routine.
            # This is where the user clicks four points in the external
            # window.
            self.pipeline.correct_ellipse(icut=icut, cmap=cmap)

            # Plot the corrected image for before/after comparison.
            corrected_buf = self.pipeline.plot_eld_2d(
                icut=icut,
                cmap=cmap,
                xlim=recalled_xlim,
                ylim=recalled_ylim
            )

            return raw_buf, corrected_buf

        def done(result):
            """
            GUI-thread completion callback for ellipse correction.

            Parameters
            ----------
            result : tuple
                ``(raw_buf, corrected_buf)`` returned by ``do()``.

            Returns
            -------
            None
                Displays before/after results and enables decision buttons.
            """

            raw_buf, corrected_buf = result

            # Show the before/after comparison.
            # There is no "Use this run" button inside the result card here
            # because Step 3 has dedicated buttons: Use corrected image/Skip.
            self._card_ellipse.show_buffers(
                [("Before correction", raw_buf), ("After correction", corrected_buf)],
                run_label=f"icut={icut}",
                auto_select=False,
                show_use_button=False
            )

            # Restore remembered axis limits into the card controls.
            self._card_ellipse.set_axis_fields(
                recalled_xlim,
                recalled_ylim
            )

            def _refresh_ellipse(xlim, ylim, icut=icut, cmap=cmap):
                """
                Replot before/after ellipse images for new axis limits.

                Parameters
                ----------
                xlim, ylim : tuple or None
                    Axis limits requested by the result card.

                Returns
                -------
                list of tuple
                    ``[("Before correction",raw),("After correction",corrected)]``.
                """
                # plot_raw_eld_2d shows the original uncorrected image.
                raw = self.pipeline.plot_raw_eld_2d(
                    icut=icut,
                    cmap=cmap,
                    xlim=xlim,
                    ylim=ylim
                )

                # plot_eld_2d shows the currently active ELD image, which is the
                # corrected image after correct_ellipse().
                corrected = self.pipeline.plot_eld_2d(
                    icut=icut,
                    cmap=cmap,
                    xlim=xlim,
                    ylim=ylim
                )

                return [
                    ("Before correction", raw),
                    ("After correction", corrected)
                ]

            # Register refresh callback with shared image-axis memory.
            self._card_ellipse.set_refresh_callback(
                self._wrap_refresh_with_memory("image_2d", _refresh_ellipse))

            # Restore normal running controls for Step 3.
            self._finish_step(3, self._btn_ellipse, self._btn_stop3)

            # Do not auto-advance after ellipse correction.
            # The user must explicitly decide whether to use the corrected
            # image or skip and keep the original image.
            self._btn_ellipse.setEnabled(True)
            self._btn_use_ellipse.setVisible(True)
            self._btn_show_diff.setVisible(True)
            self._btn_skip_ellipse.setEnabled(True)

            self.statusBar().showMessage(
                "Ellipse correction done. Click 'Use corrected image' to proceed "
                "with the corrected pattern, or 'Skip' to use the original.")

        def on_err(e):
            """
            Error callback for ellipse correction.

            Parameters
            ----------
            e : Exception
                Exception raised by the worker.

            Returns
            -------
            None
                Restores Step-3 controls and forwards the error
                to ``_err``.
            """

            # Restore Step-3 buttons manually because the worker did
            # not complete through the normal done() path.
            self._btn_stop3.setVisible(False)
            self._btn_ellipse.setVisible(True)
            self._btn_ellipse.setEnabled(True)
            self._btn_show_diff.setVisible(False)

            self._err(
                "Ellipse correction failed",
                e,
                card=self._card_ellipse
            )

        # Start the ellipse-correction worker
        w = run_in_thread(
            self,
            do,
            on_done=done,
            on_error=on_err,
            busy_key=3
        )

        # Keep worker reference so Stop/Abort can cancel/suppress it.
        self._active_worker_3 = w

    def _on_show_ellipse_diff(self):
        """
        Display the absolute difference between the raw and ellipse-
        corrected images.

        The before/after ellipse-correction plots can sometimes look
        very similar, especially when the original diffraction rings
        were already close to circular. This helper computes and
        displays ``|before - after|`` so subtle corrections are easier
        to verify visually.

        The difference image is shown as a new result-card entry
        in the Ellipse tab. Brighter pixels correspond to regions that
        changed more during correction.

        Returns
        -------
        None
            Starts an asynchronous difference-plot calculation and
            updates the Ellipse result card.
        """
        # Capture the current intensity cutoff from the GUI.
        # This controls the display scaling used by the backend plot.
        icut = self._icut.value()

        # Disable the button while the difference image is being
        # computed. This prevents duplicate worker launches from
        # repeated clicks.
        self._btn_show_diff.setEnabled(False)

        def do():
            """
            Worker-thread body for ellipse-difference plotting.

            Returns
            -------
            object
                Plot buffer returned by ``pipeline.plot_ellipse_diff``.
            """
            # Use a high-contrast colormap for difference visualization.
            # The rest of the ellipse UI may use the selected cmap, but
            # "hot" is useful here because the plot represents change
            # magnitude.
            return self.pipeline.plot_ellipse_diff(icut=icut, cmap='hot')

        def done(buf):
            """
            GUI-thread completion callback for ellipse-difference
            plotting.

            Parameters
            ----------
            buf : object
                Plot buffer produced by ``do()``.

            Returns
            -------
            None
                Re-enables the button and displays the difference plot.
            """
            # Re-enable the button after successful computation.
            self._btn_show_diff.setEnabled(True)

            # Show the difference image as its own result-card entry.
            # There is no Use button because this plot is diagnostic only;
            # it does not change which image is used by the pipeline.
            self._card_ellipse.show_buffers(
                [("Absolute difference |Before - After|", buf)],
                run_label="difference",
                auto_select=False,
                show_use_button=False
            )

            self.statusBar().showMessage(
                "Difference shown — brighter areas changed more during correction."
            )


        def on_err(e):
            """
            Error callback for ellipse-difference plotting.

            Parameters
            ----------
            e : Exception
                Exception raised by the worker.

            Returns
            -------
            None
                Re-enables the button and reports the error.
            """
            # Make the button usable again even if plotting failed.
            self._btn_show_diff.setEnabled(True)

            self._err(
                "Could not compute difference",
                e,
                card=self._card_ellipse
            )

        # Start the difference calculation asynchronously so the GUI
        # remains responsive.
        run_in_thread(self, do, on_done=done, on_error=on_err)

    # ── Step 5 ────────────────────────────────────────────────────────────────
    def _on_1d(self):
        """
        Start Step 5: convert the 2D diffraction pattern
        to a 1D profile.

        This method performs radial averaging of the currently
        active ELD diffractogram. The output is a 1D intensity
        profile over pixel radius.

        Workflow
        --------
        1. Check whether a previous 1D conversion is still running.
        2. Resolve the optional output filename.
        3. Ask before overwriting an existing output file.
        4. Put Step 5 into running GUI state.
        5. Capture plotting limits from the GUI and axis-memory
           system.
        6. Run radial averaging in a background worker.
        7. Display the resulting 1D curve.
        8. Mark Step 5 complete and advance to background subtraction.

        Returns
        -------
        None
            Starts an asynchronous 1D-conversion task and updates
            GUI state.
        """
        # Prevent launching a new 1D conversion while the previous one
        # is still finishing in the background.
        if self._step_busy(5):
            self.statusBar().showMessage(
                "Previous 1D conversion is still finishing — please wait a moment."
            )
            return

        # Read optional output filename from the GUI.
        # If the field is empty, no output file path is passed
        # to the backend.
        outfile_name = self._1d_outfile.text().strip()
        outfile = str(self.workdir / outfile_name) if outfile_name else None

        # If an output file would be overwritten, ask the user first.
        if outfile and Path(outfile).exists():
            reply = QMessageBox.question(
                self,
                "File already exists",
                f"'{outfile_name}' already exists in the working directory "
                "(possibly from an earlier session).\n\n"
                "Overwrite it?",
                QMessageBox.StandardButton.Yes
                | QMessageBox.StandardButton.Cancel,
            )

            if reply != QMessageBox.StandardButton.Yes:
                return

        # Put Step 5 into running GUI state
        self._btn_1d.setEnabled(False)
        self._start_step(5, self._btn_1d, self._btn_stop5)

        # Show progress in the 1D result card.
        self._card_1d.show_spinner("Radially averaging...")

        # Enable and activate the 1D result tab.
        self._tabs.setTabEnabled(self.TAB_1D, True)
        self._tabs.setCurrentIndex(self.TAB_1D)

        self.statusBar().showMessage("Radial averaging...")

        # Pixel-radius limits are taken directly from the
        # Step-5 sidebar fields. hese are the intended x-limits
        # for the 1D plot.
        xlim = (self._px_min.value(), self._px_max.value())

        # There is no separate Y-limit control in this sidebar
        # section, so recall the last remembered Y range for
        # pixel/intensity plots.
        _, recalled_ylim = self._recall_axis("pixels_intensity")

        def do():
            """
            Worker-thread body for 1D conversion.

            Returns
            -------
            object
                Plot buffer returned by ``pipeline.plot_eld_1d``.
            """
            # Run radial averaging and optionally save the output
            # file.
            self.pipeline.calculate_1d(
                out_file=outfile
            )

            # Plot the resulting 1D profile using the captured x-range
            # and recalled y-range.
            return self.pipeline.plot_eld_1d(
                xlim=xlim,
                ylim=recalled_ylim
            )

        def done(buf):
            """
            GUI-thread completion callback for 1D conversion.

            Parameters
            ----------
            buf : object
                Plot buffer produced by ``do()``.

            Returns
            -------
            None
                Displays the 1D curve, stores axis memory, and advances workflow.
            """
            # Label the run by the pixel range used for radial
            # averaging display.
            px_label = f"px {xlim[0]}-{xlim[1]}"

            # Show the 1D plot.
            # There is no acceptance step here; the result is accepted
            # immediately.
            self._card_1d.show_buffers(
                [("", buf)],
                run_label=px_label,
                auto_select=True,
                show_use_button=False
            )

            # Populate the result card axis fields with the active limits.
            self._card_1d.set_axis_fields(xlim, recalled_ylim)

            # Remember this step's axis limits for later pixel/intensity
            # plots. This is important for Step 6, because background
            # subtraction has no Pixel min/max fields of its own and would
            # otherwise use a hard-coded fallback range.
            self._remember_axis("pixels_intensity", xlim, recalled_ylim)

            # Register a refresh callback so changing axis limits redraws
            # the 1D plot and stores the new limits in axis memory.
            self._card_1d.set_refresh_callback(
                self._wrap_refresh_with_memory("pixels_intensity",
                    lambda xl, yl: [("", self.pipeline.plot_eld_1d(
                        xlim=xl if xl else xlim, ylim=yl))]))

            # Restore Step-5 controls and advance the workflow.
            self._finish_step(5, self._btn_1d, self._btn_stop5)
            self._complete_step(5)

            self.statusBar().showMessage(
                "1D ready - subtract background to continue."
            )

        # Start the worker thread.
        w = run_in_thread(
            self,
            do,
            on_done=done,
            on_error=lambda e: self._err(
                "Conversion failed", e, card=self._card_1d
            ),
            busy_key=5,
        )

        # Keep worker reference so Stop/Abort can cancel/suppress it.
        self._active_worker_5 = w

    # ── Step 6 ────────────────────────────────────────────────────────────────
    def _on_bkg(self):
        """
        Start Step 6: subtract background from the 1D diffraction
        profile.

        This method applies the selected background-subtraction method
        to the current 1D diffraction curve. It then displays three
        diagnostic plots:
           - original data,
           - background estimate,
           - background-subtracted data.

        Workflow
        --------
        1. Check whether a previous background-subtraction worker is still
           running.
        2. Resolve the selected background method and output/background-file
           path, and - for SimpleFuncs/Baselines - read and validate the
           algorithm, pixel range and parameters (a bad entry is reported
           before anything is started).
        3. For interactive background selection, ask before overwriting
           an existing background file.
        4. Put Step 6 into running GUI state.
        5. Optionally show Tutorial help for interactive background picking.
        6. Recall pixel/intensity axis limits from Step 5.
        7. Run background subtraction in a background worker.
        8. Display original/background/subtracted plots.
        9. Store a deep-copy snapshot of the result so selecting this run
           later restores the correct background-subtracted diffractogram.

        Returns
        -------
        None
            Starts an asynchronous background-subtraction task and updates
            GUI state.
        """
        # Prevent launching a second background-subtraction job while
        # the previous one is still finishing.
        if self._step_busy(6):
            self.statusBar().showMessage(
                "Previous background subtraction is still finishing — please wait a moment."
            )
            return

        # Capture selected method and the background file path.
        method   = self._bkg_method.currentText()
        bkg_file = str(self.workdir / self._bkg_file.text().strip())

        # Read + validate the SimpleFuncs/Baselines options now, on the main
        # thread: do() runs on a worker thread and must not touch widgets,
        # and a bad entry should be reported before the step is marked as
        # running (so a typo doesn't leave the step looking busy).
        try:
            algorithm, bkg_xrange, bkg_params = self._bkg_algorithm_args()
        except ValueError as e:
            QMessageBox.warning(self, "EDIFF", str(e))
            return

        # Tab label identifying exactly what this run used, built from the
        # captured values (not live widgets) so it can't drift if the user
        # edits a field while the calculation is running.
        bkg_label = method
        if algorithm:
            shown = self._read_bkg_params()
            bkg_label = algorithm + "".join(
                f" {k}={v if isinstance(v, bool) else format(v, 'g')}"
                for k, v in shown.items())
            if bkg_params and set(bkg_params) - set(shown):
                bkg_label += " +extra"

        # InteractivePlot writes clicked background points to bkg_file.
        # Ask before overwriting an existing file because this can destroy
        # points from an earlier session/run.
        if method == "InteractivePlot" and Path(bkg_file).exists():
            reply = QMessageBox.question(
                self, "File already exists",
                f"'{self._bkg_file.text().strip()}' already exists in the "
                "working directory (possibly from an earlier session).\n\n"
                "Clicking new points will overwrite it. Continue?",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.Cancel)
            if reply != QMessageBox.StandardButton.Yes:
                return

        # Put Step 6 into running GUI state
        self._btn_bkg.setEnabled(False)
        self._start_step(6, self._btn_bkg, self._btn_stop6)

        # Show progress in the background-subtraction card.
        self._card_bkg.show_spinner("Subtracting background...")

        # Enable and activate the background result tab.
        self._tabs.setTabEnabled(self.TAB_BKG, True)
        self._tabs.setCurrentIndex(self.TAB_BKG)
        self.statusBar().showMessage("Subtracting background...")

        # Show tutorial panel for interactive point-picking background
        # subtraction.
        if method == "InteractivePlot":
            self._maybe_show_tutorial("interactive_bkg")

        # Carry over the same pixel/intensity axis range used for the 1D curve.
        # Background has no Pixel min/max fields of its own, so this keeps the
        # visual range consistent with Step 5.
        recalled_xlim, recalled_ylim = self._recall_axis("pixels_intensity")

        def do():
            """
            Worker-thread body for background subtraction.

            Returns
            -------
            list
                Plot buffers returned by ``pipeline.plot_background``.
            """
            # Apply selected background-subtraction method.
            # This may mutate pipeline.eld.diffractogram in place.
            self.pipeline.subtract_background(
                method=method,
                bkg_file=bkg_file,
                algorithm=algorithm,
                xrange=bkg_xrange,
                params=bkg_params
            )

            # Plot original data, estimated background, and corrected
            # profile.
            return self.pipeline.plot_background(
                xlim=recalled_xlim if recalled_xlim else (0, 500),
                ylim=recalled_ylim
            )

        def done(bufs):
            """
            GUI-thread completion callback for background subtraction.

            Parameters
            ----------
            bufs : list
                Plot buffers for original, background, and background-
                subtracted profiles.

            Returns
            -------
            None
                Displays background-subtraction results and prepares
                run selection.
            """
            # If this is the first background run, there is no alternative
            # to compare against, so auto-accept it after storing its
            # snapshot.
            is_only_run = len(self._card_bkg._panels) == 0

            # Display the three diagnostic plots.
            panel = self._card_bkg.show_buffers(
                list(
                    zip(
                        ["Original","With background","Without background"],
                        bufs)
                ),
                run_label=bkg_label,
                auto_select=False,
                show_use_button=True
            )

            # Display the three diagnostic plots.
            self._card_bkg.set_axis_fields(recalled_xlim, recalled_ylim)

            # Register refresh callback for axis-limit changes.
            self._card_bkg.set_refresh_callback(
                self._wrap_refresh_with_memory(
                    "pixels_intensity",
                    lambda xl, yl: list(
                        zip(
                            ["Original",
                             "With background",
                             "Without background"
                             ],
                            self.pipeline.plot_background(
                                xlim=xl if xl else (0, 500),
                                ylim=yl
                            ),
                        ),
                    ),
                ),
            )
            # Snapshot this run's result so selecting it later restores it.
            # Background subtraction mutates eld.diffractogram in place,
            # so each attempted method/run needs its own independent copy.
            # Otherwise, later runs would overwrite the data backing earlier
            # panels.

            try:
                panel._bkg_snapshot = {
                    'smet': self.pipeline.smet,
                    'diffractogram': self.pipeline.eld.diffractogram.copy(
                        deep=True
                    ),
                }

            except Exception:
                panel._bkg_snapshot = None

            # Restore Step-6 controls.
            self._finish_step(6, self._btn_bkg, self._btn_stop6)

            # Avoid duplicate signal connections after repeated runs.
            try:
                self._card_bkg.run_selected.disconnect(
                    self._on_bkg_run_selected
                )

            except: pass

            self._card_bkg.run_selected.connect(self._on_bkg_run_selected)

            # Auto-accept every run, not just the first one (see the
            # matching comment in the XRD handler for the full
            # rationale). "Use this run" remains available on every
            # tab for anyone who wants to revert to an earlier attempt.
            self._card_bkg.select_run(panel)

        # Start background subtraction asynchronously.
        w = run_in_thread(
            self,
            do,
            on_done=done,
            on_error=lambda e: self._err(
                "Background subtraction failed",
                e,
                card=self._card_bkg
            ),
            busy_key=6,
        )

        # Keep worker reference so Stop/Abort can cancel/suppress it.
        self._active_worker_6 = w

    def _on_bkg_run_selected(self, panel):
        """
        Accept a background-subtraction run and advance the workflow.

        Each background-subtraction result panel may contain a deep-copy
        snapshot of the corrected diffractogram and the selected
        smoothing/background method. When the user accepts a panel, this
        function restores that snapshot into the pipeline before marking
        Step 6 complete.

        Parameters
        ----------
        panel : object
            Result-card panel selected by the user. Expected to optionally
            contain ``_bkg_snapshot``.

        Returns
        -------
        None
            Restores selected background-subtraction result and advances
            to calibration.
        """

        # Disconnect to avoid duplicate callback execution after
        # repeated runs.
        try:
            self._card_bkg.run_selected.disconnect(
                self._on_bkg_run_selected
            )

        except: pass

        # Restore the snapshot associated with the selected background run.
        snap = getattr(panel, '_bkg_snapshot', None)

        if snap is not None:
            try:
                self.pipeline.smet = snap['smet']

                # Copy again when restoring so the snapshot stored on the panel
                # remains independent and reusable.
                self.pipeline.eld.diffractogram = snap['diffractogram'].copy(
                    deep=True
                )

            except Exception:
                # If restore fails, keep whatever pipeline state currently
                # exists rather than crashing the GUI.
                pass

        # Mark Step 6 complete and move to calibration/comparison.
        self._complete_step(6)
        self.statusBar().showMessage(
            "Background subtracted - calibrate to continue."
        )

    # ── Step 7 ────────────────────────────────────────────────────────────────
    def _on_calibrate(self):
        """
        Start Step 7: calibrate the ELD profile and compare
        it with XRD references.

        This method performs the final calibration/comparison step.
        It calibrates the measured ELD 1D profile into q-space, then
        plots it against the theoretical XRD pattern(s) calculated
        from the loaded structure file(s).

        Workflow
        --------
        1. Check whether a previous calibration worker is still
           finishing.
        2. Put Step 7 into running GUI state.
        3. Refresh the reference-phase combo box from ``pipeline.xrd_list``.
        4. Capture calibration settings and q-range on the GUI thread.
        5. Run calibration in a background worker.
        6. Generate comparison plots.
        7. Display comparison result panels.
        8. Store a snapshot of the calibrated diffractogram so accepting
           a run later restores the exact calibrated state.
        9. Let the user accept the run, or auto-accept the first run.

        Important threading detail
        --------------------------
        All values from Qt widgets are read before the worker starts.
        The nested ``do()`` function therefore uses captured Python
        values instead of accessing Qt widgets from a background thread.

        Returns
        -------
        None
            Starts an asynchronous calibration task and updates
            GUI state.
        """
        # Prevent launching a second calibration while the previous one
        # is still finishing in the background.
        if self._step_busy(7):
            self.statusBar().showMessage(
                "Previous calibration is still finishing — please wait a moment."
            )
            return

        # Put Step 7 into running GUI state
        # Disable the button immediately to avoid duplicate clicks.
        self._btn_cal.setEnabled(False)

        # Hide Calibrate, show Stop, enable global Abort, and remember
        # Step 7 as the currently running step.
        self._start_step(7, self._btn_cal, self._btn_stop7)

        # Show progress in the comparison/result card.
        self._card_comp.show_spinner("Calibrating...")

        # Enable and activate the final comparison tab.
        self._tabs.setTabEnabled(self.TAB_COMP, True)
        self._tabs.setCurrentIndex(self.TAB_COMP)
        self.statusBar().showMessage("Calibrating...")

        # Refresh reference-phase selector
        # The reference phase must match the currently accepted XRD result.
        # Rebuilding this combo box here ensures Step 7 uses the current
        # pipeline.xrd_list, including whichever XRD run was accepted earlier.
        self._cal_ref.clear()

        for i, xrd in enumerate(self.pipeline.xrd_list):
            self._cal_ref.addItem(
                xrd.description or f"Phase {i+1}"
            )

        # Capture selected calibration method and reference phase index.
        method  = self._cal_method.currentText()
        ref_idx = self._cal_ref.currentIndex()

        # Collect optional microscope parameters
        # These are needed only by calibration methods that use microscope/
        # camera geometry. For other methods, mic_kwargs remains empty
        # and contributes nothing to pipeline.calibrate().
        mic_kwargs = {}

        if method in ("Microscope", "MicroscopeParams"):
            mic_kwargs = dict(
                voltage_kV          = self._cal_voltage.value(),
                camera_length_mm    = self._cal_cl.value(),
                camera_pixel_size_um = self._cal_pixsize.value(),
                binning             = self._cal_binning.value(),
            )

        # Capture plotting/calibration values on the GUI thread
        # These values are used in the worker and in the result label.
        # Capturing them now ensures the label describes the actual run,
        # even if the user edits the GUI while calibration is still running.
        ft   = self._cal_finetune.value()

        qlo  = self._cal_qmin.value()
        qhi = self._cal_qmax.value()
        xlim = (qlo, qhi)

        # The q-range comes directly from the Step-7 sidebar fields.
        # There is no dedicated Y-range sidebar field here, so reuse
        # the remembered q/intensity Y-axis limits, for example from
        # the XRD tab.
        _, recalled_ylim = self._recall_axis("q_intensity")

        def do():
            """
            Worker-thread body for calibration and comparison plotting.

            Returns
            -------
            list of tuple
                Comparison plot buffers returned by ``_make_comparison_bufs``.
            """
            # Calibrate the ELD diffractogram using the selected method and
            # reference phase. This may mutate pipeline.eld.diffractogram in
            # place by converting its x-axis from pixels to q.
            self.pipeline.calibrate(
                method=method, ref_phase_idx=ref_idx, **mic_kwargs)

            # Generate final comparison plots using the captured fine-tune
            # factor and axis limits.
            return self._make_comparison_bufs(ft, xlim, recalled_ylim)

        def done(bufs):
            """
            GUI-thread completion callback for calibration.

            Parameters
            ----------
            bufs : list of tuple
                Comparison buffers returned by ``do()``.

            Returns
            -------
            None
                Displays comparison plots and prepares run selection.
            """
            # Compact run label shown on the result panel.
            cal_label = (f"{method}  "
                         f"ft={ft:.4f}  q={qlo:.1f}-{qhi:.1f}")

            # If this is the first calibration run, there is no alternative to
            # compare against, so it will be auto-accepted after snapshotting.
            is_only_run = len(self._card_comp._panels) == 0

            # Display comparison plots in the result card.
            # auto_select=False because selection happens explicitly
            # after the snapshot and signal connection are prepared.
            panel = self._card_comp.show_buffers(
                bufs,
                run_label=cal_label,
                auto_select=False,
                show_use_button=True
            )

            # Populate result-card axis fields with the q-range and
            # recalled y-range.
            self._card_comp.set_axis_fields(xlim, recalled_ylim)

            # Remember this step's q-axis range globally.
            # XRD and final comparison share the same q/intensity axis meaning,
            # so if the user recalculates XRD later, it can reuse the q-range
            # chosen during calibration.
            self._remember_axis(
                "q_intensity",
                xlim,
                recalled_ylim
            )

            # Remember this step's q-axis range globally.
            # XRD and final comparison share the same q/intensity axis meaning,
            # so if the user recalculates XRD later, it can reuse the q-range
            # chosen during calibration.
            self._card_comp.set_refresh_callback(
                self._wrap_refresh_with_memory(
                    "q_intensity",
                    lambda xl, yl, ft=ft, xlim=xlim: self._make_comparison_bufs(
                        ft, xl if xl else xlim, yl)
                )
            )

            # Snapshot the calibrated data produced by this specific run.
            # Calibration mutates eld.diffractogram in place, so each attempted
            # method/parameter combination needs its own independent copy.
            # Otherwise, a later calibration run would overwrite the data behind
            # earlier panels.

            try:
                panel._cal_snapshot = {
                    'diffractogram': self.pipeline.eld.diffractogram.copy(deep=True),
                    'calibration_constant': getattr(
                        self.pipeline.eld, 'calibration_constant', None),
                }

            except Exception:
                # If snapshotting fails, the panel can still be shown, but
                # selecting it later will not restore a saved calibrated state.
                panel._cal_snapshot = None

            # Restore Step-7 buttons from running state to normal state.
            self._finish_step(7, self._btn_cal, self._btn_stop7)

            # Avoid duplicate signal connections after repeated calibration
            # runs.
            try:
                self._card_comp.run_selected.disconnect(
                    self._on_calibrate_run_selected
                )
            except: pass

            self._card_comp.run_selected.connect(
                self._on_calibrate_run_selected
            )

            # Auto-accept every run, not just the first one (see the
            # matching comment in the XRD handler for the full
            # rationale). "Use this run" remains available on every
            # tab for anyone who wants to revert to an earlier attempt.
            self._card_comp.select_run(panel)

        # Start calibration asynchronously.
        w = run_in_thread(
            self,
            do,
            on_done=done,
            on_error=lambda e: self._err(
                "Calibration failed",
                e,
                card=self._card_comp
            ),
            busy_key=7,
        )

        # Keep worker reference so Stop/Abort can cancel/suppress it.
        self._active_worker_7 = w

    def _on_calibrate_run_selected(self, panel):
        """
        Accept a calibration run and finalize the workflow.

        Each calibration result panel may contain a snapshot of the
        calibrated diffractogram and calibration constant. When the user
        accepts a panel, this function restores that snapshot into the
        pipeline before marking Step 7 complete.

        After calibration is accepted, the Update plots and Save all
        results buttons  become enabled.

        Parameters
        ----------
        panel : object
            Result-card panel selected by the user. Expected to optionally
            contain ``_cal_snapshot``.

        Returns
        -------
        None
            Restores selected calibration result, marks Step 7 complete, and
             enables final update/save controls.
        """
        # Disconnect to avoid duplicate callback execution
        # after repeated runs.
        try:
            self._card_comp.run_selected.disconnect(
                self._on_calibrate_run_selected
            )
        except: pass

        # Restore the calibrated data produced by this specific run.
        snap = getattr(panel, '_cal_snapshot', None)

        if snap is not None:
            try:
                # Restore the calibrated diffractogram from
                # an independent copy.
                self.pipeline.eld.diffractogram = \
                    snap['diffractogram'].copy(deep=True)

                # Restore calibration constant if it was available
                # in the snapshot.
                if snap['calibration_constant'] is not None:
                    self.pipeline.eld.calibration_constant = \
                        snap['calibration_constant']

            except Exception:
                # If restore fails, keep the current pipeline state rather
                # than crashing the GUI.
                pass

        # Mark Step 7 complete. This also updates the workflow
        # state machine.
        self._complete_step(7)

        # Step 8 (Comparison) isn't a separate sidebar step in this
        # implementation -- accepting a calibration run immediately
        # shows the comparison plot too (see _make_comparison_bufs
        # below), so nothing else marks its topbar badge. _complete_step
        # only knows about steps 1-7, so correct the full picture here
        # explicitly instead of trying to extend that function.
        self._update_topbar_indicators_partial(
            done=[1, 2, 3, 4, 5, 6, 7, 8],
            active=9,
            locked=[]
        )

        # Final controls become available only after a calibration run
        # is accepted.
        self._btn_update.setEnabled(True)
        self._btn_save.setEnabled(True)


        self.statusBar().showMessage(
            "Calibration accepted. Adjust fine-tune if needed, then save."
        )


    def _make_comparison_bufs(self, ft, xlim, ylim=None):
        """
        Generate final comparison plot buffers for calibrated
        ELD vs XRD.

        This helper centralizes the plotting logic used by Step 7.
        It returns a list of named plot buffers suitable for display
        in the comparison result card.

        If multiple XRD phases are loaded, the first plot compares
        the calibrated ELD profile with all phases together. Then one
        individual comparison plot is generated for each phase.

        Parameters
        ----------
        ft : float
            Fine-tune factor applied to the comparison plot. This is
            a multiplicative adjustment used for small manual alignment
            corrections.

        xlim : tuple or None
            q-axis limits used for plotting, usually ``(qmin, qmax)``.

        ylim : tuple or None, optional
            Intensity-axis limits used for plotting. If None, the backend
            plotting function chooses its default y-limits.

        Returns
        -------
        list of tuple
            List of ``(title, buffer)`` pairs. Each buffer is produced
            by the corresponding pipeline plotting function.
        """
        # Accumulate named plot buffers for the result card.
        bufs = []

        # If there is more than one reference phase, add a combined
        # overview plot.
        if len(self.pipeline.xrd_list) > 1:
            bufs.append(
                (
                    "All phases",
                    self.pipeline.plot_all_phases_comparison(
                        fine_tune=ft,
                        xlim=xlim,
                        ylim=ylim,
                    ),
                ),
            )

        # Add one ELD-vs-XRD comparison plot for each phase.
        for i, xrd in enumerate(self.pipeline.xrd_list):
            bufs.append(
                (
                    f"ELD vs {xrd.description or f'Phase {i+1}'}",
                    self.pipeline.plot_final_comparison(
                        fine_tune=ft,
                        xlim=xlim,
                        ylim=ylim,
                        phase_idx=i
                    ),
                ),
            )

        return bufs

    def _on_update_compare(self):
        """
        Re-render the final comparison plots using the current fine-tune
        and q-range.

        This updates the displayed comparison view without rerunning
        calibration. The calibrated data already selected in Step 7 remain
        unchanged; only the visualization is regenerated using the current
        GUI values.

        Unlike long-running pipeline steps, this handler intentionally reads
        widget values live because it runs synchronously on the main GUI
        thread and is meant to reflect the fields exactly as they are
        at the moment the user clicks "Update plots".

        Returns
        -------
        None
            Updates the comparison result card and stored q/intensity
            axis memory.
        """
        # Read current fine-tune and q-range directly from the GUI.
        # This is safe because the handler runs on the main Qt thread.
        ft  = self._cal_finetune.value()

        qlo = self._cal_qmin.value(); qhi = self._cal_qmax.value()
        xlim = (qlo, qhi)

        # Reuse the last remembered y-axis range for q/intensity plots.
        # Step 7 has q min/max controls for x, but no dedicated y-limit
        # fields.
        _, recalled_ylim = self._recall_axis("q_intensity")

        # Recreate comparison plot buffers from the already selected
        # calibrated data.
        bufs = self._make_comparison_bufs(ft, xlim, recalled_ylim)

        # Label this as an updated view rather than a new calibration
        # run.
        label = f"ft={ft:.4f}  q={qlo:.1f}-{qhi:.1f} (updated)"

        # Display the updated plots.
        # show_use_button=False because this is only a different view
        # of the already accepted calibration, not a separate selectable
        # calibration run.
        self._card_comp.show_buffers(
            bufs,
            run_label=label,
            show_use_button=False
        )

        # Update the axis fields in the comparison card.
        self._card_comp.set_axis_fields(xlim, recalled_ylim)

        # Remember these q/intensity limits for other compatible plots,
        # such as XRD.
        self._remember_axis("q_intensity", xlim, recalled_ylim)

        # Register refresh callback so manual axis changes in the card
        # continue to update both the plot and the shared q/intensity
        # axis memory.
        self._card_comp.set_refresh_callback(
            self._wrap_refresh_with_memory("q_intensity",
                lambda xl, yl, ft=ft, xlim=xlim: self._make_comparison_bufs(
                    ft, xl if xl else xlim, yl)))

        # Keep the user on the final comparison tab.
        self._tabs.setCurrentIndex(self.TAB_COMP)

    # ── Save / Download ───────────────────────────────────────────────────────
    def _save_all_step_windows(self):
        """
        Save each step's currently active plot panel as a PNG
        in the working directory.

        This helper is used before saving/downloading results
        so that the exported folder contains the visible result
        from each completed step, not only the final comparison
        plot.

        The method uses each ``PlotCard.active_panel``. Therefore,
        if the user has calculated multiple runs and selected one
        using "Use this run", the saved PNG corresponds to the run
        that is actually active in the pipeline.

        Returns
        -------
        list of str
            Base names of step plots that were successfully saved,
            without the ``.png`` suffix.
        """
        # Map stable output file prefixes to their corresponding
        # result cards.
        steps = [
            ("02_xrd",        self._card_xrd),
            ("03_ellipse",    self._card_ellipse),
            ("04_center",     self._card_eld2d),
            ("05_1d",         self._card_1d),
            ("06_background", self._card_bkg),
            ("07_comparison", self._card_comp),
        ]
        saved = []

        for name, card in steps:
            # Save only cards that currently have an active panel
            # with a matplotlib Figure behind it.
            panel = card.active_panel

            if panel is None or panel.fig is None:
                continue
            try:
                # Save the active figure into the working directory.
                # bbox_inches='tight' removes excessive whitespace.
                # facecolor='#ffffff' prevents transparent/black
                # backgrounds in exported PNGs depending on matplotlib
                # theme settings.

                panel.fig.savefig(str(self.workdir / f"{name}.png"),
                                  dpi=150,
                                  bbox_inches='tight',
                                  facecolor='#ffffff'
                                  )
                saved.append(name)

            except Exception:
                # Exporting one figure should not prevent other figures
                # from being saved. The caller does not currently need
                # detailed per-file errors.
                pass

        return saved

    def _on_save(self):
        """
        Save all final pipeline outputs and current step plot windows.

        This is the handler for "Save all results". It asks the user
        for an output folder, runs the backend ``pipeline.save_all``
        operation in a worker, then saves all active step windows and
        copies PNG/TXT outputs from the working directory into the selected
        destination.

        Returns
        -------
        None
            Starts an asynchronous save operation and reports
            completion to the user.
        """
        # Do nothing if saving is not currently allowed.
        # The Save button is enabled only after accepted calibration.
        if not self._btn_save.isEnabled(): return

        # Prevent overlapping save operations.
        if self._step_busy('save'):
            self.statusBar().showMessage(
                "Already saving — please wait a moment."
            )
            return

        # Start from the last used save directory, if known — otherwise
        # default to the working directory, since that's naturally
        # where the session's inputs and intermediate outputs already
        # live, rather than some arbitrary OS-default folder.
        last = self._cfg.get("last_save_dir", str(self.workdir))

        # Ask the user where final results should be copied.
        out_dir = QFileDialog.getExistingDirectory(
            self,
            "Select output folder",
            last
        )

        # User cancelled the dialog.
        if not out_dir: return

        # Remember this directory for future save/download dialogs.
        self._cfg["last_save_dir"] = out_dir
        save_config(self._cfg)

        self.statusBar().showMessage("Saving...")

        # Capture current fine-tune and q-range from the GUI.
        # These are passed to pipeline.save_all().
        ft   = self._cal_finetune.value()
        xlim = (self._cal_qmin.value(), self._cal_qmax.value())
        out_path = Path(out_dir)

        def do():
            """
            Worker-thread body for backend saving.

            Returns
            -------
            None
                Saves pipeline-level output files into the
                working directory.
            """
            # Let the backend generate its standard output files.
            # These are written to self.pipeline.workdir / self.workdir.
            self.pipeline.save_all(fine_tune=ft, xlim=xlim)

        def done(_):
            """
            GUI-thread completion callback after backend saving.

            Parameters
            ----------
            _ : Any
                Unused result from ``do()``.

            Returns
            -------
            None
                Saves active plot windows, copies outputs, and notifies
                the user.
            """

            # Save current visible/active plot windows on the main
            # thread. This touches Qt-embedded matplotlib figures, so it
            # should not happen in the worker thread.
            self._save_all_step_windows()

            # Copy generated PNG and TXT outputs from the working directory
            # to the user-selected destination folder.
            for f in self._result_files("*.png", "*.txt"):
                shutil.copy(f, out_path / f.name)

            self.statusBar().showMessage(f"Saved to {out_dir}")
            QMessageBox.information(
                self,
                "EDIFF",
                f"Results saved to:\n{out_dir}"
            )

            # Step 9 (Save) is the final step; mark its topbar badge
            # done too, matching Step 8's fix above -- there's no
            # StepBlock 9 to update, only the topbar badge.
            self._update_topbar_indicators_partial(
                done=[1, 2, 3, 4, 5, 6, 7, 8, 9],
                active=0,
                locked=[]
            )

        # Run backend saving asynchronously.
        run_in_thread(
            self,
            do,
            on_done=done,
            on_error=lambda e: self._err("Save failed", e),
            busy_key='save'
        )


    def _result_files(self, *patterns):
        """List result files in the working directory made during this session.

        "Save all results" and "Download results" used to copy every
        ``*.png`` / ``*.txt`` / ``*.bp`` file in the working directory. That
        was fine for a private folder, but the working directory is now the
        user's own directory, which holds many unrelated files, so only files
        created or modified since this session began (``_session_started``)
        are returned. A two-second allowance covers file systems with coarse
        timestamps.

        Parameters
        ----------
        *patterns : str
            Glob patterns relative to the working directory, e.g. ``"*.png"``.

        Returns
        -------
        list of pathlib.Path
            Regular files that match and are new in this session.
        """
        found = []
        for pattern in patterns:
            for f in self.workdir.glob(pattern):
                try:
                    if f.is_file() and f.stat().st_mtime >= self._session_started - 2:
                        found.append(f)
                except OSError:
                    pass
        return found

    def _on_download(self):
        """
        Copy existing result files from the working directory
        to a chosen folder.

        This is similar to saving, but it does not call
        ``pipeline.save_all``. Instead, it first captures the current
        active plot windows, then copies all available PNG, TXT, and BP
        files from the working directory to the selected destination.

        This is useful when the user wants to copy already generated
        outputs without rerunning the backend save step.

        Returns
        -------
        None
            Copies available result files and reports the result
            to the user.
        """
        # Start from the last used output directory, if known.
        last = self._cfg.get("last_save_dir", "")

        # Ask the user where the files should be copied.
        out_dir = QFileDialog.getExistingDirectory(
            self,
            "Select folder to copy results",
            last
        )

        # User cancelled the dialog.
        if not out_dir: return

        # Capture current active step windows first.
        # This ensures Download has PNG plots to copy even if
        # "Save all results" has not been run yet.
        self._save_all_step_windows()

        # Collect supported result files from the working
        # directory.
        files = self._result_files("*.png", "*.txt", "*.bp")

        if not files:
            QMessageBox.warning(
                self,
                "EDIFF",
                "No results yet.\nComplete at least one step before downloading."
            )
            return

        # Copy every available output file to the selected destination
        # folder.
        for f in files:
            shutil.copy(f, Path(out_dir) / f.name)

        self.statusBar().showMessage(
            f"{len(files)} files copied to {out_dir}"
        )

        QMessageBox.information(
            self,
            "EDIFF",
            f"{len(files)} files copied to:\n{out_dir}"
        )

    def _on_change_mode(self):
        """
        Return from the current main window to the startup mode-selection
        dialog.

        The user is warned that unsaved results may be lost. If confirmed,
        active workers are cancelled, a new ``StartupDialog`` is shown, and
        the current main window is hidden.

        Returns
        -------
        None
            Opens the startup dialog if the user confirms mode change.
        """
        # Confirm before leaving the current workflow.
        reply = QMessageBox.question(
            self,
            "Change mode",
            "Return to mode selection?\n\nAny unsaved results will be lost.",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.Cancel,
        )

        if reply != QMessageBox.StandardButton.Yes:
            return

        # Cancel active step workers so they do not continue updating
        # a hidden or soon-to-be-reset window.
        for n in range(2, 8):
            worker = getattr(
                self,
                f"_active_worker_{n}",
                None
            )
            if worker: worker.cancel()

        # Show startup dialog and close this window
        self._startup = StartupDialog()
        self._startup.setStyleSheet(STYLE + "QWidget{background:#ffffff;}")

        # When the user picks a mode, reopen or report unsupported mode.
        self._startup.mode_selected.connect(self._reopen_with_mode)
        self._startup.show()

        # Hide the current main window instead of destroying it
        # immediately.
        self.hide()

    def _reopen_with_mode(self, mode):
        """
        Reopen the main GUI using the selected application mode.

        Currently only polycrystal mode is implemented. Selecting
        polycrystal resets the existing main window to a new session
        and shows it again. Selecting any other mode shows an informational
        message and returns to the existing window.

        Parameters
        ----------
        mode : str
            Mode selected in the startup dialog.

        Returns
        -------
        None
            Resets/shows the main window or reports unsupported mode.
        """
        if mode == "polycrystal":
            # Reset and reshow this window
            self._on_new_session(confirm=False)
            self.show()
        else:
            # Placeholder for future monocrystal or other modes.
            QMessageBox.information(self, "EDIFF",
                "Monocrystal mode is not yet implemented.")
            self.show()

    def _ui_scale_state(self):
        """Describe the interface size: running now, saved, and who set it.

        Returns
        -------
        (applied, saved, external) : (float, float, bool)
            ``applied``: the scale this run is using (from ``QT_SCALE_FACTOR``,
            1.0 if unset). ``saved``: the size stored for the next start
            (defaults to ``DEFAULT_UI_SCALE``, 150%, when nothing is saved yet). ``external``: True when
            ``QT_SCALE_FACTOR`` came from the user's own environment rather
            than from ``apply_saved_ui_scale`` - it then overrides the saved
            size, so changing it here has no effect.
        """
        try:
            applied = float(os.environ.get("QT_SCALE_FACTOR") or 1.0)
        except ValueError:
            applied = 1.0
        saved = float(self._cfg.get("ui_scale", DEFAULT_UI_SCALE))
        external = (bool(os.environ.get("QT_SCALE_FACTOR"))
                    and not os.environ.get("EDIFF_UI_SCALE_FROM_CONFIG"))
        return applied, saved, external

    def _set_ui_scale(self, new):
        """Save `new` as the interface size to use from the next start.

        The size is stored as ``ui_scale`` in the config file and applied by
        ``apply_saved_ui_scale`` on the next launch - Qt cannot change its
        scale factor while the application is running, so nothing on screen
        changes immediately; the status bar says so.

        Parameters
        ----------
        new : float
            The scale factor, e.g. 1.0 for the original size or 1.5 for 150%.
        """
        applied, _saved, external = self._ui_scale_state()
        self._cfg["ui_scale"] = float(new)
        save_config(self._cfg)
        pct = f"{round(new * 100)}%"
        if external:
            msg = (f"Interface size saved ({pct}), but the QT_SCALE_FACTOR "
                   "environment variable is set and overrides it.")
        elif abs(new - applied) < 1e-6:
            msg = f"Interface size: {pct} (this is the current size)."
        else:
            msg = (f"Interface size saved ({pct}) - restart EDIFF to apply it "
                   "(Size menu > Restart EDIFF to apply).")
        self.statusBar().showMessage(msg)

    def _confirm_restart(self, pct):
        """Ask whether to restart EDIFF now; True if the user agrees."""
        reply = QMessageBox.question(
            self, "Interface size",
            f"EDIFF has to restart to use {pct}.\n\n"
            "Restart now? Anything not saved (results, session) will be lost.",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            QMessageBox.StandardButton.No)
        return reply == QMessageBox.StandardButton.Yes

    def _choose_ui_scale(self, new):
        """Handle picking a size in the menu: save it and offer to restart.

        Used by the menu's direct size entries. (The Ctrl+= / Ctrl+- / Ctrl+0
        shortcuts call `_set_ui_scale` / `_change_ui_scale` instead and never
        interrupt with a dialog; the menu's "Restart EDIFF to apply" entry
        applies whatever they saved.) No question is asked when the chosen
        size is the one already running, or when the user's own
        ``QT_SCALE_FACTOR`` would override it.

        Parameters
        ----------
        new : float
            The scale factor chosen.
        """
        self._set_ui_scale(new)
        applied, _saved, external = self._ui_scale_state()
        if external or abs(new - applied) < 1e-6:
            return
        if self._confirm_restart(f"{round(new * 100)}%"):
            self._restart_now()

    def _change_ui_scale(self, direction):
        """Step the saved interface size to the next larger or smaller one.

        Steps through ``UI_SCALES`` starting from the listed size nearest to
        the saved one, and stops at either end of the list. Only saves the
        size (status bar message, no dialog).

        Parameters
        ----------
        direction : int
            +1 for the next larger size, -1 for the next smaller one.
        """
        scales = self.UI_SCALES
        _applied, current, _external = self._ui_scale_state()
        i = min(range(len(scales)), key=lambda k: abs(scales[k] - current))
        self._set_ui_scale(scales[max(0, min(len(scales) - 1, i + direction))])

    def _refresh_size_menu(self):
        """Update the "Size" menu each time it opens.

        Ticks the saved size, shows what is running now versus what is saved
        in the menu's first line, and enables "Restart EDIFF to apply" only
        when a restart would actually change something.
        """
        applied, saved, external = self._ui_scale_state()
        for scale, act in self._size_actions:
            act.setChecked(abs(scale - saved) < 0.005)
        pct = lambda v: f"{round(v * 100)}%"
        if external:
            info = f"Running at {pct(applied)} (set by QT_SCALE_FACTOR)"
        elif abs(saved - applied) < 1e-6:
            info = f"Running at {pct(applied)}"
        else:
            info = f"Running at {pct(applied)}  -  {pct(saved)} after restart"
        self._size_info.setText(info)
        self._size_restart.setEnabled(
            not external and abs(saved - applied) >= 1e-6)

    @staticmethod
    def _relaunch_spec():
        """Build the command line and environment to start EDIFF again.

        Returns
        -------
        (cmd, env) : (list of str, dict) or None
            ``cmd`` repeats how this process was started (``sys.orig_argv``,
            e.g. ``python -m ediff gui``). None if that is unavailable
            (Python older than 3.10).

        Notes
        -----
        ``QT_SCALE_FACTOR`` and its marker are removed from ``env`` when
        ``apply_saved_ui_scale`` set them: they describe THIS run, and
        leaving them would make the new process keep the old size instead
        of reading the newly saved one. A ``QT_SCALE_FACTOR`` from the
        user's own environment is kept - it always wins.
        """
        orig = getattr(sys, "orig_argv", None)
        if not orig or len(orig) < 2:
            return None
        env = dict(os.environ)
        if env.pop("EDIFF_UI_SCALE_FROM_CONFIG", None):
            env.pop("QT_SCALE_FACTOR", None)
        return [sys.executable] + list(orig[1:]), env

    def _restart_from_menu(self):
        """Menu entry "Restart EDIFF to apply...": confirm, then restart."""
        _applied, saved, _external = self._ui_scale_state()
        if self._confirm_restart(f"{round(saved * 100)}%"):
            self._restart_now()

    def _restart_now(self):
        """Start a new EDIFF process, then close this one.

        The new process is started first, with the same command line and
        working directory; only if that succeeds is this window closed. If
        the restart is impossible (no ``sys.orig_argv``, or the process
        cannot be started) a message says so and EDIFF keeps running, so
        nothing is lost.
        """
        import subprocess
        spec = self._relaunch_spec()
        if spec is None:
            QMessageBox.information(
                self, "EDIFF",
                "EDIFF cannot restart itself with this version of Python.\n"
                "Please close it and start it again.")
            return
        cmd, env = spec
        try:
            subprocess.Popen(cmd, env=env, cwd=os.getcwd())
        except Exception as e:
            QMessageBox.warning(
                self, "EDIFF",
                f"Could not restart automatically:\n{e}\n\n"
                "Please close EDIFF and start it again.")
            return
        self.close()
        QApplication.quit()

    def _on_about(self):
        """
        Show the About dialog.

        Returns
        -------
        None
            Opens ``AboutDialog`` modally.
        """
        dlg = AboutDialog(self)

        # block interaction with the main window until the dialog
        # is closed.
        dlg.exec()

    def _remember_axis(self, group, xlim, ylim):
        """
        Store axis limits for a given plot-axis group.

        Axis memory lets different steps share sensible plot limits
        when they use the same physical axis meaning. For example:

        - XRD and final comparison both use q/intensity axes.
        - 1D conversion and background subtraction both use
          pixel/intensity axes.
        - ellipse and center views both use 2D image pixel coordinates.

        Parameters
        ----------
        group : str
            Axis-memory group name, such as ``"q_intensity"``,
            ``"pixels_intensity"``, or ``"image_2d"``.

        xlim : tuple or None
            X-axis limits to remember.

        ylim : tuple or None
            Y-axis limits to remember.

        Returns
        -------
        None
            Updates ``self._axis_memory`` in place.
        """
        # Store exactly what the caller provides.
        # A None value means that no explicit limit is remembered
        # for that axis.
        self._axis_memory[group] = (xlim, ylim)

    def _recall_axis(self, group):
        """
        Return remembered axis limits for a given plot-axis group.

        Parameters
        ----------
        group : str
            Axis-memory group name.

        Returns
        -------
        tuple
            ``(xlim, ylim)`` for the requested group, or ``(None, None)``
            if no limits have been remembered yet.
        """
        # Missing groups default to no explicit axis limits.
        return self._axis_memory.get(group, (None, None))

    def _wrap_refresh_with_memory(self, group, callback):
        """
        Wrap a PlotCard refresh callback so axis changes
        are remembered.

        PlotCards can re-render plots when the user changes axis
        limits. This helper decorates such callbacks so that every
        refresh also updates shared axis memory for the relevant axis
        group.

        Parameters
        ----------
        group : str
            Axis-memory group to update.

        callback : callable
            Original refresh callback. It must accept ``(xlim, ylim)``
            and return the plot buffers expected by the corresponding
            PlotCard.

        Returns
        -------
        callable
            Wrapped callback that first remembers ``(xlim, ylim)`` and
            then calls the original callback.
        """
        def wrapped(xlim, ylim):
            """
            Remember requested axis limits, then redraw the plot.

            Parameters
            ----------
            xlim : tuple or None
                Requested x-axis limits.

            ylim : tuple or None
                Requested y-axis limits.

            Returns
            -------
            object
                Whatever ``callback(xlim, ylim)`` returns.
            """
            self._remember_axis(group, xlim, ylim)
            return callback(xlim, ylim)

        return wrapped

    def _on_new_session(self, confirm=True):
        """
       Reset the application state and start a new EDIFF session.

       This clears loaded files, result cards, step states, cached
       workflow flags, active worker references, axis memory, tutorial
       mode, and the current pipeline object. Files already written
       in the working directory are not deleted.

       Parameters
       ----------
       confirm : bool, optional
           If True, ask the user to confirm before clearing the current
           session. If False, reset immediately. This is useful when
           reopening the app from mode selection.

       Returns
       -------
       None
           Resets GUI and pipeline state in place.
       """
        if confirm:
            # Ask before clearing the current GUI state.
            # Existing files on disk are left untouched.
            reply = QMessageBox.question(
                self, "New session",
                "Clear all results and start from scratch?\n\n"
                "This resets all loaded files and parameters in the app. "
                "It does not delete anything from the working directory — "
                "any files already there are left untouched.",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.Cancel,
            )
            if reply != QMessageBox.StandardButton.Yes:
                return

        # Reset workflow/control flags.
        self._current_step_n = None
        self._btn_abort.setEnabled(False)
        self._busy_keys = set()
        import time
        self._session_started = time.time()   # earlier results are not part of the new session
        self._auto_run_queue = None
        self._axis_memory = {}

        # Reset Tutorial mode and its button appearance.
        self._tutorial_mode = False
        self._btn_tutorial_mode.setText("Tutorial  ○")
        self._btn_tutorial_mode.setStyleSheet(self._tutorial_off_style)

        # Cancel any active step workers.
        # This suppresses their results if they eventually return.
        for n in range(2, 8):
            worker = getattr(self, f"_active_worker_{n}", None)
            if worker:
                worker.cancel()

        # Create a fresh backend pipeline object.
        # Existing files in the working directory are not deleted here.
        # New results  will overwrite EDIFF output files naturally as steps
        # are rerun.
        from ediff.gui.pipeline import EdiffPipeline

        self.pipeline = EdiffPipeline(self.workdir)
        self._loaded_eld_path = None

        # Reset input-file state.
        self.cif_paths = []
        self.eld_path  = None

        # Reset optional center-read-file state and label.
        self._ctr_read_file_path = None
        self._ctr_read_file_lbl.setText("No file selected")
        self._ctr_read_file_lbl.setStyleSheet(
            "color:#86868b; font-size:12px;"
        )

        # Refresh Step-1 file displays.
        self._refresh_cif_list()
        self._refresh_eld_list()

        # Reset main window title.
        self.setWindowTitle("EDIFF")

        # Clear all result cards.
        for card in [
            self._card_preview,
            self._card_xrd,
            self._card_eld2d,
            self._card_ellipse,
            self._card_1d,
            self._card_bkg,
            self._card_comp
        ]:
            card.clear_runs()

        # Disable all result tabs until new data are produced.
        for i in range(self._tabs.count()):
            self._tabs.setTabEnabled(i, False)

        # Hide all Stop buttons.
        for btn in [
            self._btn_stop2,
            self._btn_stop4,
            self._btn_stop3,
            self._btn_stop5,
            self._btn_stop6,
            self._btn_stop7
        ]:
            btn.setVisible(False)

        # Restore all normal Run buttons.
        for btn in [
            self._btn_xrd,
            self._btn_center,
            self._btn_ellipse,
            self._btn_1d,
            self._btn_bkg,
            self._btn_cal
        ]:
            btn.setVisible(True)
            btn.setEnabled(True)

        # Hide the ellipse decision button until a new correction
        # has been run.
        self._btn_use_ellipse.setVisible(False)

        # Disable final controls until a new calibration is accepted.
        self._btn_update.setEnabled(False)
        self._btn_save.setEnabled(False)

        # Return workflow to Step 1.
        self._go_to_step(1)
        self.statusBar().showMessage(
            "New session started - load CIF and ELD files to begin."
        )


    def _err(self, title, msg, card=None):
        """
        Handle errors raised by background workers or GUI-triggered
        pipeline actions.

        This helper centralizes error cleanup for the workflow. It resets
        the active result card, reports the error to the user, closes
        tutorial dialogs that could obscure the modal error message, and
        restores the main step buttons to a usable state.

        A special ``"__CANCELLED__"`` marker is treated as a user-requested
        abort rather than a real error, so no modal error dialog is shown
        in that case.

        Parameters
        ----------
        title : str
           Title used for the modal error dialog.

        msg : Exception or str
           Error message or exception object to display.

        card : PlotCard or None, optional
           Result card associated with the failing step. If provided, its
           spinner is reset to an idle state.

        Returns
        -------
        None
           Updates GUI state and displays an error message if appropriate.
        """
        # Clear the spinner on the card that was running.
        # Without this, a failed computation could leave the GUI looking
        # as if it were still processing forever.
        if card is not None:
            card.reset_idle()

        # "__CANCELLED__" is used as an internal marker for deliberate
        # user aborts. In that case, do not show a scary modal error dialog.
        if '__CANCELLED__' in str(msg):
            self.statusBar().showMessage("Step cancelled.")

        else:
            # "__CANCELLED__" is used as an internal marker for deliberate
            # user aborts. In that case, do not show a scary modal error
            # dialog.
            self.statusBar().showMessage(f"Error: {msg}")

            # Close tutorial dialogs before showing the error dialog.
            # Tutorial dialogs use WindowStaysOnTopHint. If one remains open,
            # the modal QMessageBox may appear behind it. The app would then
            # look frozen because the hidden modal dialog blocks input.
            if hasattr(self, '_tutorial_dialogs'):
                for d in list(self._tutorial_dialogs):
                    try:
                        if d.isVisible():
                            d.close()
                    except RuntimeError:
                        # The Python wrapper may still exist after
                        # the underlying C++ Qt object has already
                        # been deleted.
                        pass

                self._tutorial_dialogs = []

            # Show the actual error dialog.
            QMessageBox.critical(self, title, str(msg))

            # Restore all main Run buttons so the user can recover
            # and try again.
            for btn in [
                self._btn_xrd,
                self._btn_center,
                self._btn_ellipse,
                self._btn_1d,
                self._btn_bkg,
                self._btn_cal,
            ]:
                btn.setEnabled(True)
                btn.setVisible(True)

            # Hide all step-specific Stop buttons because no step should
            # remain in a visible running state after error handling.
            for btn in [
                self._btn_stop2,
                self._btn_stop4,
                self._btn_stop3,
                self._btn_stop5,
                self._btn_stop6,
                self._btn_stop7,
            ]:
                btn.setVisible(False)

            # Hide ellipse-specific decision button after an error.
            # A failed correction should not leave "Use corrected image"
            # available.
            self._btn_use_ellipse.setVisible(False)

            # Disable global Abort because the error handler has already
            # cleaned up the failed/aborted operation.
            self._btn_abort.setEnabled(False)


# ═════════════════════════════════════════════════════════════════════════════
def fix_tooltips(widget):
    """
    Force white tooltip styling recursively for a widget subtree.

    Qt/Fusion styling can sometimes inherit system dark-mode tooltip
    colors, producing unreadable combinations such as dark text
    on a dark background or grey text on black. This helper appends
    a consistent ``QToolTip`` stylesheet to every child widget that
    does not already define tooltip styling.

    Parameters
    ----------
    widget : QWidget
        Root widget whose child widgets should receive the tooltip
        stylesheet.

    Returns
    -------
    None
        Modifies child widget stylesheets in place.
    """

    # Shared tooltip stylesheet used across the application.
    TT = (
        "QToolTip {"
        "background-color: #ffffff; color: #1d1d1f;"
        "border: 1px solid #c7c7cc; border-radius: 5px;"
        "padding: 4px 8px; font-size: 11px;}"
    )

    # Apply to every child widget under the given root widget.
    for w in widget.findChildren(QWidget):
        existing = w.styleSheet()

        # Do not duplicate tooltip styling if the widget already
        # defines it.
        if "QToolTip" not in existing:
            w.setStyleSheet(existing + TT)


# Interface size used when the user has not chosen one yet (1.5 = 150%).
# The "Size" menu overrides it, and an explicitly chosen 100% is respected.
DEFAULT_UI_SCALE = 1.5


def apply_saved_ui_scale():
    """Apply the interface size by setting ``QT_SCALE_FACTOR``.

    The GUI's text and widgets are sized in fixed pixels, so on a large or
    high-resolution screen at 100% display scaling everything looks small.
    Qt can only change its global scale factor when the application starts,
    so this runs first thing in ``main`` and sets ``QT_SCALE_FACTOR`` - Qt
    then scales fonts, widgets and stylesheet pixel sizes uniformly.

    The size is the one saved by the GUI's "Size" menu as ``ui_scale`` in the
    config file, or ``DEFAULT_UI_SCALE`` (150%) when none has been saved or
    the saved value is unusable (not a number, or outside 0.5-4.0). A saved
    1.0 means 100% and applies nothing.

    Must be called BEFORE the ``QApplication`` is created. An explicit
    ``QT_SCALE_FACTOR`` already in the environment always wins and is left
    untouched.

    Returns
    -------
    float or None
        The scale that was applied, or None if nothing was changed (the
        environment already sets one, or the size is 100%).
    """
    if os.environ.get("QT_SCALE_FACTOR"):
        return None
    try:
        scale = float(load_config().get("ui_scale", DEFAULT_UI_SCALE))
    except Exception:
        scale = DEFAULT_UI_SCALE
    if not (0.5 <= scale <= 4.0):
        scale = DEFAULT_UI_SCALE
    if scale == 1.0:
        return None
    os.environ["QT_SCALE_FACTOR"] = repr(scale)
    # Remember that WE set it, so the GUI can tell this apart from a value
    # the user put in their own environment (which would override the menu).
    os.environ["EDIFF_UI_SCALE_FROM_CONFIG"] = "1"
    return scale


def main(prewarmed_pipeline=None):
    """
    Application entry point for the EDIFF GUI.

    This function initializes the Qt application, applies global
    styling and palette settings, shows the startup mode-selection
    dialog, displays a splash screen while the main window loads,
    and finally starts the Qt event loop.

    Current behavior
    ----------------
    Only ``polycrystal`` mode is implemented. If another mode is
    selected, the application exits.

    Returns
    -------
    None
        Starts the GUI event loop and exits the process when
        the app closes.
    """

    # Apply the saved interface size FIRST: Qt reads QT_SCALE_FACTOR only
    # when the QApplication is created.
    apply_saved_ui_scale()

    # Create the Qt application object from command-line a
    # rguments.
    app = QApplication(sys.argv)

    # Pre-warm the interactive worker as early as possible in the app's
    # lifecycle -- before the startup dialog, splash screen, or main
    # window even exist -- to give its cold-start import chain
    # (numpy/scipy/matplotlib/cv2/ediff.center, commonly several
    # seconds on a typical machine) the maximum possible head start
    # before the user could conceivably reach an interactive step. This
    # only launches a subprocess (near-instant); it does not block
    # waiting for that subprocess's own imports to finish.
    if prewarmed_pipeline is None:
        prewarmed_pipeline = EdiffPipeline(default_workdir())
    else:
        # A pipeline made earlier (by ediff.gui.__main__) may still point at
        # an older default folder: start in the user's own directory.
        prewarmed_pipeline.workdir = default_workdir()

    # Use the Fusion style for consistent cross-platform widget
    # rendering.
    app.setStyle("Fusion")

    # Force a fixed numeric locale.
    # Without this, QDoubleSpinBox follows the operating-system locale.
    # Under Czech locale settings, for example, decimal commas may be
    # expected instead of decimal points. If the user types a period but
    # the validator expects a comma, the field can appear focused while
    # silently rejecting input.
    # Using English/United States ensures decimal points behave consistently
    # across platforms.
    QLocale.setDefault(
        QLocale(
            QLocale.Language.English,
            QLocale.Country.UnitedStates,
        )
    )

    # Apply the global application stylesheet.
    app.setStyleSheet(STYLE)

    # (1) Global palette fixes
    # Fusion can still use system palette colors for item views,
    # combo-box dropdowns, and tooltips. Explicitly setting these
    # roles avoids dark-mode mismatches such as grey text on black
    # dropdown backgrounds.
    from qtpy.QtGui import QPalette, QColor

    palette = app.palette()

    # Base/item-view colors.
    palette.setColor(QPalette.ColorRole.Base, QColor("#ffffff"))
    palette.setColor(QPalette.ColorRole.AlternateBase, QColor("#f9f9fb"))
    palette.setColor(QPalette.ColorRole.Text, QColor("#1d1d1f"))

    # General window and button colors.
    palette.setColor(QPalette.ColorRole.Window, QColor("#f5f5f7"))
    palette.setColor(QPalette.ColorRole.WindowText, QColor("#1d1d1f"))
    palette.setColor(QPalette.ColorRole.Button, QColor("#ffffff"))
    palette.setColor(QPalette.ColorRole.ButtonText, QColor("#1d1d1f"))

    # Selection/highlight colors.
    palette.setColor(QPalette.ColorRole.Highlight, QColor("#e8f0fe"))
    palette.setColor(QPalette.ColorRole.HighlightedText, QColor("#0071e3"))

    # Tooltip palette colors.
    palette.setColor(QPalette.ColorRole.ToolTipBase, QColor("#ffffff"))
    palette.setColor(QPalette.ColorRole.ToolTipText, QColor("#1d1d1f"))

    app.setPalette(palette)

    # (2) Application icon
    # Set the app-wide icon as early as possible so it appears
    # on the taskbar even for the startup dialog.
    _icon_px = _pixmap("ediff.png")

    if not _icon_px.isNull():
        app.setWindowIcon(QIcon(_icon_px))

    # (3) Startup mode selection
    # Show the startup dialog before constructing the main window.
    startup = StartupDialog()
    startup.setStyleSheet(STYLE + "QWidget{background:#ffffff;}")

    # Store the selected mode in a mutable list so the nested callback
    # can modify it without needing nonlocal/global state.
    chosen_mode = []

    def on_mode(mode):
        """
        Store the mode selected in the startup dialog.

        Parameters
        ----------
        mode : str
            Mode emitted by ``StartupDialog.mode_selected``.

        Returns
        -------
        None
        """

        chosen_mode.append(mode)

    startup.mode_selected.connect(on_mode)
    startup.show()

    # (4) Manually process events until the startup dialog closes.
    # This is a simple blocking startup-selection loop. The full app
    # event loop starts later with app.exec().
    while startup.isVisible():
        app.processEvents()

    # If the user closed the startup dialog without choosing a mode,
    # quit.
    if not chosen_mode:
        sys.exit(0)

    # (5) Main-window creation
    # Show splash screen while the main window is being constructed.
    splash = make_splash()
    splash.show()
    app.processEvents()

    if chosen_mode[0] == "polycrystal":
        # Main implemented workflow.
        win = EdiffWindow(pipeline=prewarmed_pipeline)

    else:
        # Monocrystal or other modes are not implemented yet.
        # In normal use, the startup dialog probably prevents
        # this path.
        sys.exit(0)

    # Close the splash screen once the main window is ready.
    splash.finish(win)

    # Apply tooltip styling recursively to the main window.
    fix_tooltips(win)

    # Show the main EDIFF GUI.
    win.showMaximized()

    # Start the Qt event loop and exit with its return code.
    sys.exit(app.exec())


if __name__ == "__main__":
    main()