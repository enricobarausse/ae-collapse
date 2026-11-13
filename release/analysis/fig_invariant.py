#!/usr/bin/env python3
"""The spin-0 wave in invariant variables, figure invariant_wout_dm; v0 = 1.267, c2 = 0.

NOT IN THE PAPER: the script is kept as a supplementary check of the 1/Theta
falloff of w_out and of the growth of delta m with radius (the Coulomb term of eq.
dm-coulomb, which the paper shows with fig:coulomb instead).

Definitions:
  w_out = -(b_N/v0) trK + a_r                                   [eq:wout]
  u = t - Theta_*/v0,  Theta_* = Theta + 2 M_ADM ln(Theta/(2 M_ADM) - 1)   [eq:retarded]
  m = Misner-Sharp mass [eq:ms];  delta m(Theta, t) = m(Theta, t) - m(Theta, t_first),
  the change since the first dump on the same sphere.
Left: Theta w_out on spheres of fixed areal radius (default 120, 180, 240, 300, 360)
against u: the profiles overlap (asymptotic waveform).  Right: delta m against u: same
retarded-time profile, amplitude growing with radius (the aether Coulomb term, eq:dm-coulomb).

Prints per sphere: the main maximum and minimum of Theta w_out and their u, and the maximum
and minimum of delta m (public-code run: Theta w_out peaks 0.41 -> 0.36 and minima -0.14 ->
-0.19 from Theta = 120 to 360; delta m maxima ~9 -> ~23.5).

Usage:  fig_invariant.py RUNDIR --c2 0 [--madm M]
        fig_invariant.py --replot data/invariant_wout_dm.npz [--out invariant_wout_dm.png]
"""
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import spin0_common as sc
import paper_style as ps


def plot(D, out):
    """(a) Theta w_out and (b) delta m on fixed-Theta spheres against u."""
    fig, (axl, axr) = plt.subplots(1, 2, figsize=(ps.TEXT, 2.7))
    cols = ps.sphere_colors(D["thetas"])
    for k, R in enumerate(D["thetas"]):
        c = cols[k]
        axl.plot(D["u"][:, k], D["Thw"][:, k], color=c, label=rf"$\Theta={R:.0f}$")
        axr.plot(D["u"][:, k], D["dm"][:, k], color=c)
    for ax, tag in ((axl, "(a)"), (axr, "(b)")):
        ax.axhline(0, color="k", lw=0.5)
        ax.set_xlim(D["umin"], D["umax"])
        ax.set_xlabel(r"$u=t-\Theta_*/v_0$")
        ps.panel_tag(ax, tag, loc="upper right")
    axl.set_ylabel(r"$\Theta\,w_{\rm out}$")
    axr.set_ylabel(r"$\delta m$")
    h, lab = axl.get_legend_handles_labels()
    fig.legend(h, lab, loc="lower center", ncol=len(lab), bbox_to_anchor=(0.5, 0.97))
    fig.tight_layout(w_pad=1.5)
    ps.save(fig, out)


def main():
    if ps.replot(plot, "invariant_wout_dm.png"):
        return
    ap = sc.add_common_args(argparse.ArgumentParser(description=__doc__,
                            formatter_class=argparse.RawDescriptionHelpFormatter),
                            out="invariant_wout_dm.png")
    ap.set_defaults(c2=0.0)
    ap.add_argument("--thetas", default="120,180,240,300,360",
                    help="areal radii of the spheres (comma list)")
    ap.add_argument("--umin", type=float, default=-30.0, help="retarded-time window: start")
    ap.add_argument("--umax", type=float, default=120.0, help="retarded-time window: end")
    ps.add_replot_arg(ap)
    a = ap.parse_args()
    cp = sc.couplings(a.c1, a.c2)
    M = sc.resolve_madm(a)
    v0 = cp["v0"]
    thetas = np.array([float(x) for x in a.thetas.split(",")])
    print(f"v0 = {v0:.4f}, M_ADM = {M:.4f}")

    t, S = sc.sphere_series(a.rundir, thetas, ["wout", "m"], cp=cp, nproc=a.nproc)
    U, Y, DM = (np.full((len(t), len(thetas)), np.nan) for _ in range(3))
    print(f"\n{'Theta':>6} {'max Th*w':>9} {'u':>7} {'min Th*w':>9} {'u':>7} "
          f"{'max dm':>8} {'min dm':>8}")
    for k, R in enumerate(thetas):
        u = sc.retarded_time(t, R, M, v0)
        y = R * S["wout"][:, k]
        dm = S["m"][:, k] - S["m"][0, k]
        s = (u >= a.umin) & (u <= a.umax) & np.isfinite(y)
        if not s.any():
            print(f"{R:6.0f}  (no data in the u window)")
            continue
        jx, jn = np.argmax(y[s]), np.argmin(y[s])
        print(f"{R:6.0f} {y[s][jx]:9.4f} {u[s][jx]:7.1f} {y[s][jn]:9.4f} {u[s][jn]:7.1f} "
              f"{np.nanmax(dm[s]):8.3f} {np.nanmin(dm[s]):8.3f}")
        U[:, k], Y[s, k], DM[s, k] = u, y[s], dm[s]
    out = a.out
    ps.setup()
    D = dict(thetas=thetas, u=U, Thw=Y, dm=DM, umin=a.umin, umax=a.umax)
    ps.save_data(out, "invariant_wout_dm", **D)
    plot(D, out)


if __name__ == "__main__":
    main()
