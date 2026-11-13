#!/usr/bin/env python3
"""Misner-Sharp mass function at t = 28: paper fig:mass (mass_function.pdf) and the numbers of
Sec. "The mass function and the horizon" (M_hor = Theta_AH/2, M_ADM, M_hor/M_ADM, the
exterior energy M_ADM - M_hor and the sign of the aether energy density).

Run (README §5, "Mass function and horizon"): A0 = 0.35, N = 51201, r_max = 300, moving
mesh, profile at t = 28.

Definitions:  m(Theta) = (Theta/2)(1 - h), h = Theta'^2 - (3Q - 2K)^2 Theta^2/36 [eq:ms],
with Theta' = np.gradient(Theta, r) on the (non-uniform) r column; Theta_AH = outermost
zero of theta_+ = Theta' + (3Q - 2K) Theta/6 (linear interpolation between grid points);
M_hor = Theta_AH/2 (the horizon mass of the gamma_rr = 1 slicing, since h = 0 there);
M_ADM = m at Theta = --th_adm (default 250, far in the asymptotic region);
rho_ae = aether energy density measured by the aether (spin0_common.aether_set:
c14 (a_r' + a_r^2/2 + 2 a_r Theta'/Theta) - (c13 + 3 c2) K^2/6 - (3/4) c13 Q^2, couplings
--c1, --c2); printed: the Theta intervals outside the horizon where rho_ae > 0.

Paper: M_hor = 6.35 (from the Richardson limit 12.70 of Theta_AH, Sec. "Convergence"),
M_ADM = 3.7111, M_ADM - M_hor = -2.64, M_hor/M_ADM = 1.71; rho_ae < 0 everywhere outside
the horizon except across the outgoing scalar shell (15 <~ Theta <~ 30).  Public-code
run: Theta_AH = 12.677 on this grid (M_hor = 6.338), M_ADM = 3.7111, ratio 1.708,
M_ADM - M_hor = -2.639 with the Richardson Theta_AH; rho_ae > 0 only at Theta = 18.3-20.3,
inside the scalar shell (rho_psi > 1e-2 of its maximum outside the horizon at
Theta = 12.7-20.6).

Usage:  fig_mass_function.py RUNDIR [--t 28] [--th_adm 250] [--th_max 290] [--th_rich 12.70]
                             [--c1 0.7] [--c2 -0.5119403] [--out mass_function.png]
        fig_mass_function.py --replot data/mass_function.npz [--out mass_function.png]
"""
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import spin0_common as sc
import paper_style as ps


def plot(D, out):
    """m(Theta) with Theta_AH, M_hor and M_ADM marked (single column)."""
    fig, ax = plt.subplots(figsize=(ps.COLUMN, 2.5))
    ax.plot(D["Theta"], D["m"], color=ps.BLUE, lw=1.3, zorder=3)
    ax.axvline(D["th_ah"], color=ps.RED, lw=1.0, ls="--")
    ax.axhline(D["M_adm"], color=ps.GREEN, lw=1.2, ls=":")
    ax.plot([D["th_ah"]], [D["M_hor"]], "o", color=ps.RED, ms=5, mec="white", mew=0.6,
            zorder=5)
    ax.set_xlim(0, D["Theta"].max() + 5)
    ax.set_ylim(-0.2, 7.4)          # room above M_hor ~ 6.3 for its label (the paper's run)
    ax.text(D["th_ah"] + 8, D["M_hor"] + 0.3,
            r"$M_{\rm hor}=\Theta_{\rm AH}/2$", va="bottom")
    ax.text(D["Theta"].max() - 5, D["M_adm"] + 0.25, rf"$M_{{\rm ADM}}={D['M_adm']:.4f}$",
            ha="right", va="bottom")
    ax.text(D["th_ah"] + 4, 0.3, r"$\Theta_{\rm AH}$", color=ps.RED, ha="left")
    ax.set_xlabel(r"areal radius $\Theta$")
    ax.set_ylabel(r"Misner--Sharp mass $m$")
    fig.tight_layout()
    ps.save(fig, out)


def main():
    if ps.replot(plot, "mass_function.png"):
        return
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("rundir", help="run directory with prof_t*.csv dumps")
    ap.add_argument("--t", type=float, default=28.0, help="dump time (default 28)")
    ap.add_argument("--th_adm", type=float, default=250.0,
                    help="areal radius at which M_ADM is read (default 250)")
    ap.add_argument("--th_max", type=float, default=290.0,
                    help="largest areal radius plotted (default 290)")
    ap.add_argument("--c1", type=float, default=sc.C1_DEFAULT, help="coupling c1 (default 0.7)")
    ap.add_argument("--c2", type=float, default=sc.C2_BASELINE,
                    help="coupling c2 (default -0.5119403)")
    ap.add_argument("--th_rich", type=float, default=12.70,
                    help="Richardson limit of Theta_AH from fig_conv_collapse.py (default 12.70)")
    ap.add_argument("--out", default="mass_function.png",
                    help="output figure (png; a pdf is written next to it)")
    ps.add_replot_arg(ap)
    a = ap.parse_args()

    p = sc.load_dump(a.rundir, a.t)
    Th = p["Theta"]
    m = sc.misner_sharp(p, np.gradient(Th, p["r"]))
    _, r_ah, th_ah = sc.apparent_horizon(p)
    M_adm = m[np.argmin(np.abs(Th - a.th_adm))]
    M_hor = th_ah / 2.0

    sel = (Th > 1) & (Th < a.th_max)
    ps.setup()
    D = dict(Theta=Th[sel], m=m[sel], th_ah=th_ah, M_hor=M_hor, M_adm=M_adm, t=a.t)
    ps.save_data(a.out, "mass_function", **D)
    plot(D, a.out)
    print(f"t = {a.t:g}: Theta_AH = {th_ah:.3f} (r_AH = {r_ah:.3f})  M_hor = {M_hor:.3f}  "
          f"M_ADM(Theta = {a.th_adm:g}) = {M_adm:.4f}  M_hor/M_ADM = {M_hor / M_adm:.3f}"
          "   [paper: M_hor 6.35 (Richardson 12.70), 3.7111, 1.71]")
    print(f"M_ADM - M_hor = {M_adm - M_hor:.3f} on this grid, {M_adm - a.th_rich / 2:.3f} with "
          f"M_hor = {a.th_rich:g}/2   [paper: -2.64]")
    z = np.zeros_like(Th)
    rho = sc.aether_set(p, sc.couplings(a.c1, a.c2), z, z, z)[0]   # rho needs no time derivatives
    outside = (p["r"] > r_ah) & (Th < a.th_max)        # outside the horizon

    def spans(mask):
        e = np.flatnonzero(np.diff(np.concatenate([[0], mask.astype(int), [0]])))
        return ", ".join(f"{Th[i]:.1f}-{Th[k - 1]:.1f}" for i, k in zip(e[::2], e[1::2])) or "nowhere"
    rho_psi = sc.scalar_set(p)[0]
    print("aether energy density rho_ae > 0 outside the horizon at Theta = "
          + spans((rho > 0) & outside) + ";  scalar shell (rho_psi > 1e-2 of its maximum outside "
          "the horizon) at Theta = " + spans(outside & (rho_psi > 1e-2 * rho_psi[outside].max()))
          + "   [paper: rho_ae < 0 except across the outgoing scalar shell, 15 <~ Theta <~ 30]")


if __name__ == "__main__":
    main()
