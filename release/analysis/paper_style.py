"""Shared plotting style of the paper's figures (PRD, REVTeX two-column).

Every script that draws a figure of the paper calls `setup()` before drawing and `save()`
to write it (fig_flux_identity.py, a quick-look diagnostic outside the paper, does not).
The figures
are designed at their printed size (figsize in inches = size on the page), so the font
sizes below are the printed sizes in points:

  single-column figure  COLUMN = 3.4 in  (\\columnwidth)
  two-column figure*    TEXT   = 7.0 in  (\\textwidth)

Text is typeset by LaTeX (Computer Modern); if no working LaTeX is found the mathtext
Computer Modern fonts are used instead (same look, slightly different spacing).

Colours (colour-blind safe, Paul Tol's "bright" set) are tied to meanings, the same in
every figure:
  RES[0..3]     grid resolutions, coarse -> finest (RES_MARK: their markers);
  CAT_A, CAT_B  the two members of a pair (fixed grid / moving mesh, uniform / refined
                grids, aether / scalar flux, v0 = 1 / v0 = 1.267);
  sphere_colors(Thetas)  spheres of fixed areal radius: "plasma" from dark violet
                       (innermost) to orange (outermost) in every figure;
  black         extrapolations to infinity and exact relations; grey dotted: references.

save(fig, out) writes out.png (300 dpi) and out.pdf (vector).  save_data(out, **arrays)
writes the plotted arrays to <dir of out>/data/<stem>.npz, and `replot()` redraws a figure
from that file (`--replot FILE.npz [--out NAME.png]`), without the run directories.
"""
import argparse
import os
import shutil
import subprocess
import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

COLUMN = 3.4          # inches, \columnwidth of REVTeX PRD
TEXT = 7.0            # inches, \textwidth of REVTeX PRD

FS_LABEL = 10         # axis labels
FS_TICK = 9           # tick labels
FS_LEGEND = 9         # legends
FS_PANEL = 10         # panel tags (a), (b), ...
FS_SMALL = 8          # small text (crowded legends and annotations)
FS_TINY = FS_SMALL - 1   # smallest text: the legends of the narrow panels of fig:sharp only

# Paul Tol "bright" (colour-blind safe)
BLUE, CYAN, GREEN, YELLOW, RED, PURPLE, GREY = ("#4477AA", "#66CCEE", "#228833", "#CCBB44",
                                                "#EE6677", "#AA3377", "#BBBBBB")
RES = (BLUE, GREEN, RED, PURPLE)            # coarse, mid, fine, finest
RES_MARK = ("s", "o", "^", "D")
CAT_A, CAT_B = BLUE, RED

_usetex = None


def _no_type1_subsetting():
    """Work around a matplotlib (3.10/3.11) bug: with text.usetex the PDF backend subsets
    the embedded Type 1 TeX fonts and drops the glyph at code 0 of cmsy, i.e. the minus
    sign (the PNG is fine, the PDF loses every minus).  Embed the full fonts instead
    (a few tens of kB per font)."""
    try:
        import matplotlib._type1font as t1
        if hasattr(t1.Type1Font, "subset"):
            t1.Type1Font.subset = lambda self, chars, prefix: self
    except Exception:          # private API: if it changes, just leave it alone
        pass


def latex_ok():
    """True if latex + dvipng + ghostscript are available and work."""
    global _usetex
    if _usetex is None:
        _usetex = all(shutil.which(x) for x in ("latex", "dvipng")) and \
            any(shutil.which(x) for x in ("gs", "gswin64c"))
        if _usetex:
            try:
                subprocess.run(["latex", "--version"], capture_output=True, timeout=30,
                               check=True)
            except Exception:
                _usetex = False
    return _usetex


def setup():
    """Set the rcParams of the paper style (call before creating figures)."""
    matplotlib.rcdefaults()
    rc = {
        "font.size": FS_LABEL,
        "axes.labelsize": FS_LABEL,
        "axes.titlesize": FS_LABEL,
        "xtick.labelsize": FS_TICK,
        "ytick.labelsize": FS_TICK,
        "legend.fontsize": FS_LEGEND,
        "legend.frameon": False,
        "legend.handlelength": 1.8,
        "legend.handletextpad": 0.5,
        "legend.columnspacing": 1.2,
        "legend.labelspacing": 0.3,
        "legend.borderaxespad": 0.4,
        "axes.linewidth": 0.6,
        "lines.linewidth": 1.2,
        "lines.markersize": 4,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "xtick.minor.visible": True,
        "ytick.minor.visible": True,
        "xtick.major.size": 3.5,
        "ytick.major.size": 3.5,
        "xtick.minor.size": 2,
        "ytick.minor.size": 2,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.minor.width": 0.5,
        "ytick.minor.width": 0.5,
        "axes.grid": False,
        "axes.prop_cycle": matplotlib.cycler(color=[BLUE, RED, GREEN, PURPLE, CYAN, YELLOW]),
        "axes.formatter.use_mathtext": True,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.02,
        "figure.dpi": 100,
        "pdf.fonttype": 42,
    }
    if latex_ok():
        _no_type1_subsetting()
        rc.update({"text.usetex": True, "font.family": "serif",
                   "font.serif": ["Computer Modern Roman"],
                   "text.latex.preamble": r"\usepackage{amsmath}"})
    else:
        print("[paper_style] LaTeX not available: mathtext Computer Modern fallback",
              file=sys.stderr)
        rc.update({"text.usetex": False, "mathtext.fontset": "cm", "font.family": "serif",
                   "font.serif": ["cmr10"], "axes.unicode_minus": False})
    plt.rcParams.update(rc)


def AE():
    """The ae ligature (LaTeX or unicode)."""
    return r"\text{\ae}" if plt.rcParams["text.usetex"] else "æ"


def sphere_colors(thetas):
    """Colours of a family of spheres of fixed areal radius, innermost to outermost:
    "plasma" from dark violet (smallest Theta) to orange (largest Theta)."""
    n = len(thetas)
    order = np.argsort(np.argsort(thetas))
    return [plt.cm.plasma(0.85 * k / max(n - 1, 1)) for k in order]


def panel_tag(ax, tag, loc="upper left", **kw):
    """Panel label '(a)' etc. inside the axes corner."""
    x, ha = ((0.04, "left") if "left" in loc else (0.96, "right") if "right" in loc
             else (0.5, "center"))
    y, va = (0.95, "top") if "upper" in loc else (0.05, "bottom")
    kw.setdefault("fontsize", FS_PANEL)
    return ax.text(x, y, tag, transform=ax.transAxes, ha=ha, va=va, **kw)


def headroom(ax, frac=0.2):
    """Raise the upper y limit by `frac` of the axis height (room for a tag or legend)."""
    lo, hi = ax.get_ylim()
    if ax.get_yscale() == "log":
        ax.set_ylim(lo, hi * (hi / lo) ** frac)
    else:
        ax.set_ylim(lo, hi + frac * (hi - lo))


def rescaled_norm_curves(ax, t, C, ns, scal):
    """Constraint norms C[k] (one row per grid N_r = ns[k]) rescaled by scal[k] (N_r^2
    relative to the coarsest grid), with the resolution colours and markers: the curves
    collapse for second-order convergence (fig:conv-collapse, fig:conv-constraint)."""
    for k, n in enumerate(ns):
        s = scal[k]
        lab = rf"${s}\times$, $N_r={n}$" if s > 1 else rf"$N_r={n}$"
        ax.plot(t, C[k] * s, "-", marker=RES_MARK[k], color=RES[k], ms=3.5,
                mec="white", mew=0.4, label=lab, zorder=3 + k)


def _sci(v):
    """Tick label m x 10^e of a positive value v (10^e when m = 1)."""
    if v <= 0:
        return ""
    e = int(np.floor(np.log10(v) + 1e-9))
    m = v / 10.0 ** e
    if abs(m - 1) < 1e-6:
        return rf"$10^{{{e}}}$"
    return rf"${m:.0f}\times10^{{{e}}}$" if abs(m - round(m)) < 1e-6 else rf"${m:.1f}\times10^{{{e}}}$"


def log_ticks(ax, axis="y", subs=(1, 2, 5), plain=False):
    """Readable log-axis labels: decades only if the axis spans >= 2 decades, otherwise
    the m x 10^e ticks with m in `subs` (or m = 1..9 if that gives fewer than two ticks
    inside the range), all written as m x 10^e (plain=True: as plain numbers, %g); minor
    ticks unlabelled.  Call after the data are drawn and the limits are final."""
    a = ax.yaxis if axis == "y" else ax.xaxis
    lo, hi = (ax.get_ylim() if axis == "y" else ax.get_xlim())
    if np.log10(hi / lo) >= 2.0:
        a.set_major_locator(mticker.LogLocator(base=10, numticks=12))
    else:
        def cand(ms):
            es = range(int(np.floor(np.log10(lo))) - 1, int(np.ceil(np.log10(hi))) + 1)
            v = [m * 10.0 ** e for e in es for m in ms]
            return [x for x in v if lo * (1 - 1e-9) <= x <= hi * (1 + 1e-9)]
        ticks = cand(subs)
        if len(ticks) < 2:
            ticks = cand(range(1, 10))
        a.set_major_locator(mticker.FixedLocator(ticks))
    a.set_major_formatter(mticker.FuncFormatter(lambda v, _: rf"${v:g}$" if plain else _sci(v)))
    a.set_minor_locator(mticker.LogLocator(base=10, subs=np.arange(1, 10), numticks=20))
    a.set_minor_formatter(mticker.NullFormatter())


def save(fig, out):
    """Write <stem>.png (300 dpi) and <stem>.pdf; returns the two paths."""
    stem = os.path.splitext(out)[0]
    d = os.path.dirname(stem)
    if d:
        os.makedirs(d, exist_ok=True)
    fig.savefig(stem + ".png", dpi=300)
    fig.savefig(stem + ".pdf")
    plt.close(fig)
    w, h = fig.get_size_inches()
    print(f"wrote {stem}.png and {stem}.pdf (designed at {w:.2f} x {h:.2f} in)")
    return stem + ".png", stem + ".pdf"


def data_path(out):
    """Path of the data file of figure `out`: <dir of out>/data/<stem>.npz."""
    stem = os.path.splitext(out)[0]
    return os.path.join(os.path.dirname(stem), "data", os.path.basename(stem) + ".npz")


def save_data(out, figure, **arrays):
    """Save the plotted arrays of figure `figure` next to `out` (data/<stem>.npz)."""
    p = data_path(out)
    os.makedirs(os.path.dirname(p), exist_ok=True)
    np.savez_compressed(p, _figure=np.array(figure), **{k: np.asarray(v)
                                                       for k, v in arrays.items()})
    print(f"wrote {p}")
    return p


def load_data(path):
    """Arrays of a save_data file as a dict (0-d arrays as Python scalars)."""
    with np.load(path, allow_pickle=False) as z:
        D = {k: z[k] for k in z.files}
    for k, v in D.items():
        if v.ndim == 0:
            D[k] = v.item()
    return D


def add_replot_arg(ap):
    ap.add_argument("--replot", metavar="NPZ", default=None,
                    help="redraw the figure from the .npz saved by a previous run "
                         "(data/<figure>.npz next to the figure); no run directory needed")
    return ap


def replot(plotters, default_out):
    """If --replot NPZ is on the command line, redraw from it and return True.

    plotters: a function plot(D, out) or a dict {figure name: plot function};
    default_out: output name when --out is not given (for a dict: a directory, used when
    --outdir is not given either; the figure is written as <dir>/<npz stem>.png)."""
    pre = argparse.ArgumentParser(add_help=False)
    pre.add_argument("--replot", default=None)
    pre.add_argument("--out", default=None)
    pre.add_argument("--outdir", default=None)
    k, _ = pre.parse_known_args()
    if not k.replot:
        return False
    setup()
    D = load_data(k.replot)
    name = D.pop("_figure")
    if isinstance(plotters, dict):
        stem = os.path.splitext(os.path.basename(k.replot))[0]
        out = k.out or os.path.join(k.outdir or default_out, stem + ".png")
        plotters[name](D, out)
    else:
        plotters(D, k.out or default_out)
    return True
