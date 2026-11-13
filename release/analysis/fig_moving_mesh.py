#!/usr/bin/env python3
"""Moving-mesh validation: paper fig:mesh (moving_mesh.pdf) (Sec. "Formulation and numerical method",
paragraph "Discretization"; Theta_AH(t) also in Sec. "The mass function and the horizon").

Runs (README §5, "Moving-mesh validation"): A0 = 0.35, N = 25600, r_max = 300, t <= 150,
dumps every 2, once on the fixed grid (--ale=0) and once with the moving mesh (--ale=1).

Left panel: areal radius Theta(r_max) of the outer boundary against t on both grids (the
fixed grid falls inward with the gamma_rr = 1 slicing; the mesh holds it) and the boundary
coordinate r(r_max) of the moving mesh.  Right panel: areal radius of the apparent horizon
(outermost zero of theta_+, linear interpolation) on both grids.

Prints Theta_AH at --t (28) on both grids, the maximum difference over the common dump
times, the range of Theta_AH on the moving mesh and Theta(r_max) first -> last.
Paper: the fixed-grid Theta(r_max) shrinks and the moving mesh holds it; Theta_AH agrees
on the two grids, reaches its maximum at t ~ 28 and decreases from 12.6 to 11.8 by
t = 150.  Public-code runs: Theta(r_max) 284 -> 124 by t = 150 on the fixed grid (not
quoted in the paper); Theta_AH equal to 8.0e-4; 12.60 (paper 12.6) at t = 28 -> 11.8.
The axis limits (Theta <= 500, Theta_AH <= 13.1, t <= 152) are those of the paper's runs
(r_max = 300, t <= 150).

Usage:  fig_moving_mesh.py FIXED_DIR ALE_DIR [--out moving_mesh.png]
        fig_moving_mesh.py --replot data/moving_mesh.npz [--out moving_mesh.png]
"""
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import spin0_common as sc
import paper_style as ps


def series(rundir):
    t, paths = sc.list_dumps(rundir, exact=False)
    Tb, Rb, AH = [], [], []
    for f in paths:
        p = sc.load_profile(f)
        Tb.append(p["Theta"][-1])
        Rb.append(p["r"][-1])
        AH.append(sc.apparent_horizon(p)[2])
    return np.array(t), np.array(Tb), np.array(Rb), np.array(AH)


def plot(D, out):
    """Draw the figure from the arrays of save_data (also used by --replot)."""
    fig, ax = plt.subplots(1, 2, figsize=(ps.TEXT, 2.6))
    ax[0].plot(D["t0"], D["Tb0"], "-", color=ps.CAT_B, label=r"$\Theta(r_{\max})$, fixed grid")
    ax[0].plot(D["t1"], D["Tb1"], "-", color=ps.CAT_A, label=r"$\Theta(r_{\max})$, moving mesh")
    ax[0].plot(D["t1"], D["Rb1"], "--", color=ps.CAT_A, label=r"boundary coordinate $r$, moving mesh")
    ax[0].set_xlabel("$t$")
    ax[0].set_ylabel("outer boundary radius")
    ax[0].set_ylim(0, 500)               # axis limits: the paper's runs (r_max = 300, t <= 150)
    ax[0].legend(loc="lower left")
    ax[1].plot(D["t0"], D["AH0"], "o", color=ps.CAT_B, ms=4.5, mfc="none", mew=0.9,
               label="fixed grid")
    ax[1].plot(D["t1"], D["AH1"], ".", color=ps.CAT_A, ms=3.5, label="moving mesh")
    ax[1].set_xlabel("$t$")
    ax[1].set_ylabel(r"$\Theta_{\rm AH}$")
    ax[1].set_ylim(top=13.1)
    ax[1].legend(loc="lower right")
    for a_, tag in zip(ax, "ab"):
        a_.set_xlim(0, 152)
        ps.panel_tag(a_, f"({tag})")
    fig.tight_layout(w_pad=1.5)
    ps.save(fig, out)


def main():
    if ps.replot(plot, "moving_mesh.png"):
        return
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("fixed", help="run directory of the fixed-grid run (--ale=0)")
    ap.add_argument("ale", help="run directory of the moving-mesh run (--ale=1)")
    ap.add_argument("--t", type=float, default=28.0,
                    help="dump time at which Theta_AH of the two runs is printed (default 28)")
    ap.add_argument("--out", default="moving_mesh.png", help="output figure (png; a pdf is "
                    "written next to it)")
    ps.add_replot_arg(ap)
    a = ap.parse_args()

    t0, Tb0, _, AH0 = series(a.fixed)
    t1, Tb1, Rb1, AH1 = series(a.ale)
    D = dict(t0=t0, Tb0=Tb0, AH0=AH0, t1=t1, Tb1=Tb1, Rb1=Rb1, AH1=AH1)
    ps.setup()
    ps.save_data(a.out, "moving_mesh", **D)
    plot(D, a.out)

    common = sorted(set(t0) & set(t1))
    d = [abs(AH0[t0 == t][0] - AH1[t1 == t][0]) for t in common
         if np.isfinite(AH0[t0 == t][0]) and np.isfinite(AH1[t1 == t][0])]
    if a.t not in common:
        raise SystemExit(f"--t {a.t:g} is not a dump time of both runs")
    print(f"Theta_AH at t = {a.t:g}: fixed {float(AH0[t0 == a.t][0]):.4f}, "
          f"moving mesh {float(AH1[t1 == a.t][0]):.4f}"
          f";  max |difference| over {len(d)} common times: {max(d):.2e}   [paper: the two agree]")
    g = np.isfinite(AH1)
    print(f"Theta_AH on the moving mesh: first at t = {t1[g][0]:g}, max {np.nanmax(AH1):.3f} at "
          f"t = {t1[np.nanargmax(AH1)]:g}, last {AH1[g][-1]:.3f} at t = {t1[g][-1]:g}"
          "   [paper: maximum at t ~ 28, 12.6 -> 11.8 by t = 150]")
    print(f"Theta(r_max): fixed {Tb0[0]:.1f} -> {Tb0[-1]:.1f};  moving mesh {Tb1[0]:.1f} -> "
          f"{Tb1[-1]:.1f} (r_boundary {Rb1[0]:.1f} -> {Rb1[-1]:.1f})   [public-code run; not in the paper]")


if __name__ == "__main__":
    main()
