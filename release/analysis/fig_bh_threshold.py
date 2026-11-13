#!/usr/bin/env python3
"""Black-hole formation threshold of the time-symmetric shell: paper fig:threshold (bh_threshold.pdf)
(Sec. "Convergence"; paper: the smallest amplitude with an early horizon is 0.135 on the
coarsest grid and 0.125 on the two finest, so A* lies between 0.120 and 0.125).

Reads the bisection runs of README §5 ("Black-hole threshold": r_max = 300, t <= 60, dumps
every 2, directories bis_N<N>_A<A0>) and plots, against the amplitude A0 and for each grid,
the largest areal radius Theta of a trapped surface (outermost zero crossing of
theta_+ = Theta' + (3Q - 2K) Theta/6 at proper radius r > 0.5, spin0_common.apparent_horizon)
found at t <= --t_early (default 20), i.e. the apparent
horizon that forms when the shell first passes through the centre (Theta = 0 where none
forms).  The script also prints, per run, the trapped-surface radius at the last dump
(t = 60), which below the threshold picks up the small central region that becomes trapped
later while the central lapse collapses (a resolution effect).

The amplitude of each run is read by spin0_common.run_amplitude (cmd.txt, otherwise the
directory name).

Usage:  fig_bh_threshold.py RUNDIR [RUNDIR ...] [--out bh_threshold.png]
        fig_bh_threshold.py --replot data/bh_threshold.npz [--out bh_threshold.png]
"""
import argparse
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import spin0_common as sc
import paper_style as ps


def plot(D, out):
    """Draw the figure from the arrays of save_data (also used by --replot): the areal radius
    of the apparent horizon that forms when the shell reaches the centre (zero where none
    forms) against A0, one curve per grid; single column."""
    fig, ax = plt.subplots(1, 1, figsize=(ps.COLUMN, 2.9))
    Ns = sorted(set(D["N"].astype(int)))
    for i, N in enumerate(Ns):
        m = D["N"].astype(int) == N
        o = np.argsort(D["A0"][m])
        A, early = (D[k][m][o] for k in ("A0", "early"))
        dr = D["rmax"] / (N - 1)
        lab = rf"$N_r={N}$, $dr={dr:.4f}$"
        kw = dict(marker=ps.RES_MARK[i % 4], color=ps.RES[i % 4], ms=4.5, lw=1.0,
                  mec="white", mew=0.4)
        ax.plot(A, early, label=lab, **kw)
    ax.set_ylim(-0.08, 2.1)
    ax.set_ylabel(rf"$\Theta_{{\rm AH}}$ at $t\leq{D['t_early']:g}$")
    ax.set_xlabel(r"$A_0$")
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=2, frameon=False,
              fontsize=ps.FS_SMALL, handlelength=1.6, columnspacing=1.0)
    fig.tight_layout()
    ps.save(fig, out)


def main():
    if ps.replot(plot, "bh_threshold.png"):
        return
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("rundirs", nargs="+", help="run directories (one per N, A0)")
    ap.add_argument("--out", default="bh_threshold.png",
                    help="output figure (png; a pdf is written next to it)")
    ps.add_replot_arg(ap)
    ap.add_argument("--t_early", type=float, default=20.0,
                    help="horizons found at t <= this count as 'early' (default 20)")
    ap.add_argument("--rmax", type=float, default=300.0, help="r_max of the runs, for dr")
    a = ap.parse_args()

    rows = []
    for d in a.rundirs:
        name = os.path.basename(d.rstrip("/"))
        if not os.path.isdir(d):
            print(f"{name:>22}  (not a directory, skipped)")
            continue
        try:
            t, paths = sc.list_dumps(d, exact=False)
        except FileNotFoundError:
            print(f"{name:>22}  (no dumps yet, skipped)")
            continue
        early, last, N = 0.0, 0.0, 0
        for tt, f in zip(t, paths):
            p = sc.load_profile(f)
            N = len(p["r"])
            th = sc.apparent_horizon(p)[2]          # zeros at r > 0.5
            th = th if np.isfinite(th) else 0.0
            if tt <= a.t_early:
                early = max(early, th)
            last = th
        rows.append((N, sc.run_amplitude(d), early, last, t[-1]))
        print(f"{name:>22} N={N:6d} A0={rows[-1][1]:.4f} "
              f"early={early:.3f} last={last:.3f} (t={t[-1]:g})")
    rows = np.array(rows)
    ps.setup()
    D = dict(N=rows[:, 0], A0=rows[:, 1], early=rows[:, 2], last=rows[:, 3], t_end=rows[:, 4],
             rmax=a.rmax, t_early=a.t_early)
    ps.save_data(a.out, "bh_threshold", **D)
    plot(D, a.out)


if __name__ == "__main__":
    main()
