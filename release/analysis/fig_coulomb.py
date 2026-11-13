#!/usr/bin/env python3
"""The aether Coulomb field left behind by the spin-0 pulse, paper fig:coulomb (aether_coulomb.pdf)
(Sec. "Mass conservation and radiated energy"; v0 = 1.267).

Definitions:
  Theta^2 a_r on spheres of fixed areal radius (default Theta = 75, 100, 150, 200), against
  u = t - Theta_*/v0 [eq:retarded];
  m = Misner-Sharp mass [eq:ms];  m_c = m - (c14/2) Theta^2 a_r  [eq:mc].
(a) Theta^2 a_r vs u.  (b) m - M_ADM (thin) and m_c - M_ADM (thick) vs u.

Prints per sphere:
  * the maximum of Theta^2 a_r during the pulse, and at the end the power-law exponent of
    this maximum against Theta (paper: the amplitude grows roughly as Theta^0.7;
    public-code run: 19.7 at Theta = 75 -> 40.5 at 200, exponent 0.73);
  * its value behind the pulse, averaged over u in [--late] (default 200..250)
    (paper: all spheres relax to the same value; public-code run 4.6-5.0);
  * the swing (max - min over u in [--umin, --umax]) of m and of m_c, and the peak
    excursions max(m - M_ADM), max|m_c - M_ADM| (paper: the swing of m grows with radius
    and removing the Coulomb term reduces it; public-code run: 7.1 -> 18.8 for m, 0.9 ->
    0.7 for m_c on Theta = 75 -> 200);
  * m - M_ADM and m_c - M_ADM averaged over the late window (paper: m settles above
    M_ADM, so that the integral of L du is negative on every sphere, and m_c settles below
    M_ADM).

Usage:  fig_coulomb.py RUNDIR --c2 0 [--madm M]
        fig_coulomb.py --replot data/aether_coulomb.npz [--out aether_coulomb.png]
"""
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import spin0_common as sc
import paper_style as ps


def plot(D, out):
    """(a) Theta^2 a_r; (b) m - M_ADM (thin) and m_c - M_ADM (thick), against u."""
    fig, (axl, axr) = plt.subplots(1, 2, figsize=(ps.TEXT, 2.7))
    cols = ps.sphere_colors(D["thetas"])
    for k, R in enumerate(D["thetas"]):
        c = cols[k]
        axl.plot(D["u"][:, k], D["Th2a"][:, k], color=c, label=rf"$\Theta={R:.0f}$")
        axr.plot(D["u"][:, k], D["dm"][:, k], color=c, lw=0.7, alpha=0.6)
        axr.plot(D["u"][:, k], D["dmc"][:, k], color=c, lw=1.6)
    for ax, tag in ((axl, "(a)"), (axr, "(b)")):
        ax.axhline(0, color="k", lw=0.5)
        ax.set_xlim(D["umin"], D["umax"])
        ax.set_xlabel(r"$u=t-\Theta_*/v_0$")
        ps.panel_tag(ax, tag, loc="upper right")
    axl.set_ylabel(r"$\Theta^2 a_r$")
    axr.set_ylabel(r"mass $-\,M_{\rm ADM}$")
    axr.set_ylim(-1.5, 3.0)
    kh = [plt.Line2D([], [], color="0.35", lw=0.7, alpha=0.6, label=r"$m-M_{\rm ADM}$"),
          plt.Line2D([], [], color="0.35", lw=1.6, label=r"$m_c-M_{\rm ADM}$")]
    axr.legend(handles=kh, loc="lower right")
    h, lab = axl.get_legend_handles_labels()
    fig.legend(h, lab, loc="lower center", ncol=len(lab), bbox_to_anchor=(0.5, 0.97))
    fig.tight_layout(w_pad=1.5)
    ps.save(fig, out)


def main():
    if ps.replot(plot, "aether_coulomb.png"):
        return
    ap = sc.add_common_args(argparse.ArgumentParser(description=__doc__,
                            formatter_class=argparse.RawDescriptionHelpFormatter),
                            out="aether_coulomb.png")
    ap.set_defaults(c2=0.0)
    ap.add_argument("--thetas", default="75,100,150,200",
                    help="areal radii of the spheres (comma list)")
    ap.add_argument("--umin", type=float, default=-60.0, help="retarded-time window: start")
    ap.add_argument("--umax", type=float, default=260.0, help="retarded-time window: end")
    ap.add_argument("--late", default="200,250",
                    help="retarded-time window u0,u1 behind the pulse, for the late averages "
                         "(default 200,250)")
    ps.add_replot_arg(ap)
    a = ap.parse_args()
    cp = sc.couplings(a.c1, a.c2)
    M = sc.resolve_madm(a)
    v0 = cp["v0"]
    thetas = np.array([float(x) for x in a.thetas.split(",")])
    ul0, ul1 = [float(x) for x in a.late.split(",")]
    print(f"v0 = {v0:.4f}, c14/2 = {cp['c14'] / 2:.3f}, M_ADM = {M:.4f}")

    t, S = sc.sphere_series(a.rundir, thetas, ["m", "Acr"], cp=cp, nproc=a.nproc)
    # per sphere (columns): u, Theta^2 a_r, m - M_ADM, m_c - M_ADM in the u window
    U, Th2a_u, dm_u, dmc_u = (np.full((len(t), len(thetas)), np.nan) for _ in range(4))
    cmax = np.full(len(thetas), np.nan)
    print("swing = max - min over the u window;  peak = max (m - M_ADM), max |m_c - M_ADM|")
    print(f"\n{'Theta':>6} {'max Th2a':>9} {'late Th2a':>10} {'swing m':>8} {'peak m-M':>9} "
          f"{'swing m_c':>10} {'peak|m_c-M|':>12} {'late m-M':>9} {'late m_c-M':>11} {'u_last':>7}")
    for k, R in enumerate(thetas):
        u = sc.retarded_time(t, R, M, v0)
        Th2a = R ** 2 * S["Acr"][:, k]
        m = S["m"][:, k]
        mc = m - 0.5 * cp["c14"] * Th2a          # m_c of eq:mc (sc.m_coulomb)
        s = (u >= a.umin) & (u <= a.umax) & np.isfinite(Th2a)
        if not s.any():
            print(f"{R:6.0f}  (no data in the u window)")
            continue
        late = (u >= ul0) & (u <= ul1) & np.isfinite(Th2a)

        def late_mean(y):
            return np.mean(y[late]) if late.any() else np.nan

        print(f"{R:6.0f} {np.max(Th2a[s]):9.3f} {late_mean(Th2a):10.3f} {np.ptp(m[s]):8.3f} "
              f"{np.max(m[s]) - M:9.3f} {np.ptp(mc[s]):10.3f} {np.max(np.abs(mc[s] - M)):12.3f} "
              f"{late_mean(m) - M:9.3f} {late_mean(mc) - M:11.3f} {u[s][-1]:7.1f}")
        U[:, k], Th2a_u[s, k], dm_u[s, k], dmc_u[s, k] = u, Th2a[s], m[s] - M, mc[s] - M
        cmax[k] = np.max(Th2a[s])
    g = np.isfinite(cmax)
    if g.sum() >= 2:
        print(f"max Theta^2 a_r ~ Theta^{np.polyfit(np.log(thetas[g]), np.log(cmax[g]), 1)[0]:.2f}"
              f" over Theta = {thetas[g][0]:g}..{thetas[g][-1]:g}   [paper: roughly Theta^0.7]")
    out = a.out
    ps.setup()
    D = dict(thetas=thetas, u=U, Th2a=Th2a_u, dm=dm_u, dmc=dmc_u, umin=a.umin, umax=a.umax)
    ps.save_data(out, "aether_coulomb", **D)
    plot(D, out)


if __name__ == "__main__":
    main()
