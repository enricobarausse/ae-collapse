#!/usr/bin/env python3
"""Scattering convergence ladder: paper fig:conv-fields-A001, fig:conv-fields-A006 and
fig:conv-constraint (conv_scatter_fields_A001.pdf, conv_scatter_fields_A006.pdf,
conv_scatter_constraint.pdf) and the scattering rows of tab:conv-scatter (sec:conv,
paragraph "Scattering", and Appendix sec:app-conv).

Runs (README §5, "Scattering convergence ladder"): fixed grid (--ale=0), r_max = 150,
N = 6401 / 12801 / 25601 (dr halves exactly, so the coarse points coincide with every 2nd /
4th point of the finer grids and the differences are taken pointwise, no interpolation),
A0 = 0.01 / 0.04 / 0.06, dumps at t = 8, 9, 9.5, ..., 14.  One directory per run, named
<workdir>/A<A0>_N<N> (--pattern).

Conventions (standard NR self-convergence):
  fields       ||f_h - f_h/2||_2 against 4 ||f_h/2 - f_h/4||_2 over r in [--rlo, --rhi]
               (default 0..20): the curves overlap iff the scheme is second order;
  constraint   ||Theta C||_2 rescaled by N_r^2 (1x, 4x, 16x): the curves collapse iff
               ||C|| -> 0 at second order.  Theta C is the code's regularised Hamiltonian
               constraint (column ThetaC); never monitor C = Theta C / Theta itself, it is
               0/0 at the origin;
  order        p = log2(||f_h - f_h/2|| / ||f_h/2 - f_h/4||) (for Theta C: log2 of the
               successive norm ratios), averaged over --fit_times (default 12, 13, 14, after
               the origin bounce; t = 11 sits on the bounce, where a norm can pass through
               an accidental minimum).

Outputs (in --outdir, each as .png and .pdf, with the plotted arrays in --outdir/data/):
conv_scatter_fields_A001.png and conv_scatter_fields_A006.png (one per --field_amps,
2 x 3 panels), conv_scatter_constraint.png, conv_scatter_orders.png (the order table as a
colour grid, not in the paper) and the printed order table
(paper tab:conv-scatter: A0 = 0.01: 2.35 2.52 2.17 2.16 2.92 2.16 1.95; 0.04: 2.62 2.46
2.20 2.14 3.36 2.12 1.90; 0.06: 2.81 1.99 1.45 1.47 2.42 1.42 1.49 for psi, P, K, a_r,
Theta, Q, ||Theta C||).  Also printed: per amplitude, the minimum of the central lapse
over the dumps of the finest grid (paper: N_min ~ 0.98, 0.70, 0.4) and the constraint-norm
ratios coarse/mid and mid/fine at the last fit time (paper, A0 = 0.06 at t = 14: 1.8 and
2.5).

Usage:  fig_conv_scatter.py WORKDIR [--amps 0.01,0.04,0.06] [--grids 6401,12801,25601]
                            [--field_amps 0.01,0.06] [--outdir .]
        fig_conv_scatter.py --replot OUTDIR/data/<figure>.npz [--outdir DIR | --out F.png]
"""
import argparse
import functools
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
import spin0_common as sc
import paper_style as ps

FIELDS = list(sc.CONV_FIELDS)
LATEX = {"Psi": r"$\psi$", "P": r"$P$", "trK": r"$K$", "Acr": r"$a_r$",
         "Theta": r"$\Theta$", "Q": r"$Q$"}


@functools.lru_cache(maxsize=None)
def load(work, pattern, A, n, t):
    """Profile of run (A, n) at time t, or None (with a warning on stderr) if the dump is
    missing; a missing dump leaves a gap in the figures and is skipped in the orders.
    Cached: every dump is read once (callers must not modify the profiles)."""
    f = sc.dump_path(os.path.join(work, pattern.format(A=A, N=n)), t)
    if not os.path.exists(f):
        print(f"warning: missing dump {f}", file=sys.stderr)
        return None
    return sc.load_profile(f)


def nrm(v, r, lo, hi):
    """L2 norm over lo <= r <= hi."""
    return sc.norm_l2(v, r, lo, hi)


class Ladder:
    def __init__(self, a):
        self.work, self.pattern = a.workdir, a.pattern
        self.ns = [int(x) for x in a.grids.split(",")]
        self.stride = {n: (n - 1) // (self.ns[0] - 1) for n in self.ns}
        for n in self.ns:
            if (n - 1) != self.stride[n] * (self.ns[0] - 1):
                raise SystemExit(f"grid N={n} is not nested in N={self.ns[0]}")
        self.lo, self.hi = a.rlo, a.rhi

    def diffs(self, A, t, f):
        """(||coarse - mid||, ||mid - fine||) on the coarse grid's points."""
        D = {n: load(self.work, self.pattern, A, n, t) for n in self.ns}
        if any(v is None for v in D.values()):
            return np.nan, np.nan
        base = D[self.ns[0]]["r"]
        v = [D[n][f][::self.stride[n]] for n in self.ns]
        k = min(min(map(len, v)), len(base))
        v = [x[:k] for x in v]
        r = base[:k]
        return nrm(v[0] - v[1], r, self.lo, self.hi), nrm(v[1] - v[2], r, self.lo, self.hi)

    def cnorms(self, A, t):
        """||Theta C|| on each grid (NaN where the dump is missing)."""
        out = []
        for n in self.ns:
            D = load(self.work, self.pattern, A, n, t)
            out.append(nrm(D["ThetaC"], D["r"], self.lo, self.hi) if D is not None else np.nan)
        return out

    def order(self, A, f, fit_times):
        """Self-convergence order of field f ('C': ||Theta C||) averaged over fit_times."""
        orders = []
        for t in fit_times:
            if f == "C":
                a, b, c = self.cnorms(A, t)
                if b > 0 and c > 0:
                    orders += [np.log2(a / b), np.log2(b / c)]
            else:
                a, b = self.diffs(A, t, f)
                if b > 0:
                    orders.append(np.log2(a / b))
        return np.mean(orders) if orders else np.nan


def plot_orders(D, out):
    """Order table as a colour grid (not in the paper)."""
    M, amps = D["M"], [str(x) for x in D["amps"]]
    cmap = LinearSegmentedColormap.from_list("div", [ps.RED, "#f4f4f4", ps.BLUE])
    fig, ax = plt.subplots(figsize=(ps.TEXT, 0.5 + 0.45 * len(amps)))
    # diverging colours centred on the target order 2, from 1 to 2.9 (a colour scale only:
    # larger orders, e.g. Theta at A0 = 0.01 and 0.04, saturate; the numbers are printed)
    ax.imshow(M, cmap=cmap, norm=TwoSlopeNorm(vmin=1.0, vcenter=2.0, vmax=2.9), aspect="auto")
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            ax.text(j, i, f"{M[i, j]:.2f}", ha="center", va="center")
    ax.set_xticks(range(7))
    ax.set_xticklabels([LATEX[f] for f in FIELDS] + [r"$\|\Theta C\|$"])
    ax.set_yticks(range(len(amps)))
    ax.set_yticklabels([f"$A_0={A}$" for A in amps])
    ax.tick_params(which="both", length=0, top=False, right=False)
    ax.minorticks_off()
    for s in ax.spines.values():
        s.set_visible(False)
    ps.save(fig, out)


def plot_constraint(D, out):
    """||Theta C||_2 rescaled by N_r^2, one panel per amplitude."""
    amps, ns, scal, t = [str(x) for x in D["amps"]], D["ns"], D["scal"], D["t"]
    fig, axes = plt.subplots(1, len(amps), figsize=(ps.TEXT, 2.6), squeeze=False)
    for i, (ax, A) in enumerate(zip(axes[0], amps)):
        ps.rescaled_norm_curves(ax, t, D["C"][i], ns, scal)
        ax.set_yscale("log")
        ps.headroom(ax, 0.2)
        ps.log_ticks(ax)
        ax.set_xlabel("$t$")
        ps.panel_tag(ax, f"({'abcdef'[i]}) $A_0={A}$", loc="upper center")
    axes[0][0].set_ylabel(r"$N_r^2$-rescaled $\|\Theta C\|_2$")
    h, lab = axes[0][0].get_legend_handles_labels()
    fig.legend(h, lab, loc="lower center", ncol=len(ns), bbox_to_anchor=(0.5, 0.97))
    fig.tight_layout(w_pad=0.8)
    ps.save(fig, out)


def plot_fields(D, out):
    """Self-convergence of the six fields at one amplitude, 2 x 3 panels."""
    t = D["t"]
    fig, axes = plt.subplots(2, 3, figsize=(ps.TEXT, 3.5), sharex=True)
    for k, (ax, f) in enumerate(zip(axes.flat, FIELDS)):
        ax.plot(t, D["dc"][k], "-", marker="o", color=ps.RES[0], ms=3.5, mec="white",
                mew=0.4, zorder=3)
        ax.plot(t, 4 * D["df"][k], "--", marker="s", color=ps.RES[2], ms=3.5, mec="white",
                mew=0.4, zorder=4)
        ax.set_yscale("log")
        ax.margins(y=0.12)
        ps.log_ticks(ax)
        ps.panel_tag(ax, LATEX[f], loc="upper left")
    for ax in axes[1]:
        ax.set_xlabel("$t$")
    h = [plt.Line2D([], [], color=ps.RES[0], marker="o", ms=3.5, mec="white", mew=0.4,
                    label=r"$\|f_h-f_{h/2}\|_2$"),
         plt.Line2D([], [], color=ps.RES[2], ls="--", marker="s", ms=3.5, mec="white", mew=0.4,
                    label=r"$4\,\|f_{h/2}-f_{h/4}\|_2$")]
    fig.legend(handles=h, loc="lower center", ncol=2, bbox_to_anchor=(0.5, 0.975))
    fig.tight_layout(w_pad=0.6, h_pad=0.4)
    ps.save(fig, out)


PLOTTERS = {"conv_scatter_orders": plot_orders, "conv_scatter_constraint": plot_constraint,
            "conv_scatter_fields": plot_fields}


def main():
    if ps.replot(PLOTTERS, "."):
        return
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("workdir", help="directory containing the run directories A<A0>_N<N>")
    ap.add_argument("--pattern", default="A{A}_N{N}",
                    help="run-directory name pattern with {A} and {N} (default A{A}_N{N})")
    ap.add_argument("--amps", default="0.01,0.04,0.06",
                    help="amplitudes, as written in the directory names (comma list)")
    ap.add_argument("--grids", default="6401,12801,25601",
                    help="nested grids N, coarse to fine (comma list)")
    ap.add_argument("--times", default="8,9,9.5,10,10.5,11,11.5,12,13,14",
                    help="dump times plotted (comma list)")
    ap.add_argument("--fit_times", default="12,13,14",
                    help="dump times averaged for the order table (comma list)")
    ap.add_argument("--rlo", type=float, default=0.0, help="norm window: inner radius")
    ap.add_argument("--rhi", type=float, default=20.0, help="norm window: outer radius")
    ap.add_argument("--field_amps", default="0.01,0.06",
                    help="amplitudes with a field-convergence figure, one figure each "
                         "(conv_scatter_fields_A001.png, ...; default 0.01,0.06)")
    ap.add_argument("--outdir", default=".", help="directory for the figures (png + pdf; "
                    "the plotted arrays go to OUTDIR/data/*.npz)")
    ps.add_replot_arg(ap)
    a = ap.parse_args()
    amps = a.amps.split(",")
    times = [float(x) for x in a.times.split(",")]
    fit_times = [float(x) for x in a.fit_times.split(",")]
    L = Ladder(a)
    ns = L.ns
    os.makedirs(a.outdir, exist_ok=True)
    ps.setup()

    # ---- order table (colour grid; not a paper figure, the paper prints the table) ----
    M = np.array([[L.order(A, f, fit_times) for f in FIELDS + ["C"]] for A in amps])
    out = os.path.join(a.outdir, "conv_scatter_orders.png")
    D = dict(M=M, amps=np.array(amps))
    ps.save_data(out, "conv_scatter_orders", **D)
    plot_orders(D, out)

    # ---- Hamiltonian constraint, rescaled by N_r^2 ----
    C = np.array([[[L.cnorms(A, t)[j] for t in times] for j in range(len(ns))] for A in amps])
    out = os.path.join(a.outdir, "conv_scatter_constraint.png")
    D = dict(t=np.array(times), C=C, amps=np.array(amps), ns=np.array(ns),
             scal=np.array([L.stride[n] ** 2 for n in ns]))
    ps.save_data(out, "conv_scatter_constraint", **D)
    plot_constraint(D, out)

    # ---- field self-convergence overlays: one figure per amplitude ----
    for A in a.field_amps.split(","):
        dd = np.array([[L.diffs(A, t, f) for t in times] for f in FIELDS])   # (6, nt, 2)
        out = os.path.join(a.outdir, f"conv_scatter_fields_A{A.replace('.', '')}.png")
        D = dict(t=np.array(times), dc=dd[:, :, 0], df=dd[:, :, 1], A=A, ns=np.array(ns))
        ps.save_data(out, "conv_scatter_fields", **D)
        plot_fields(D, out)

    print(f"\n{'A0':>6} " + "".join(f"{f:>8}" for f in FIELDS) + f"{'||ThC||':>10}")
    for i, A in enumerate(amps):
        print(f"{A:>6} " + "".join(f"{M[i, j]:8.2f}" for j in range(6)) + f"{M[i, 6]:10.2f}")
    print(f"\n(target 2.00; orders averaged over t = {a.fit_times}, norm over "
          f"{a.rlo:g} < r < {a.rhi:g})")
    tl = fit_times[-1]
    print(f"\nminimum central lapse over the dumps (N = {ns[-1]}) and ||ThC|| ratios at t = {tl:g}"
          "   [paper: N_min ~ 0.98, 0.70, 0.4; A0 = 0.06, t = 14: 1.8 and 2.5]")
    for A in amps:
        nmin = min(D["NN"][0] for D in (load(L.work, L.pattern, A, ns[-1], t) for t in times)
                   if D is not None)
        c = L.cnorms(A, tl)
        print(f"  A0 = {A}: N_min = {nmin:.3f}   ratios " +
              " ".join(f"{c[k] / c[k + 1]:.2f}" for k in range(len(c) - 1)))


if __name__ == "__main__":
    main()
