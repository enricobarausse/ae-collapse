#!/usr/bin/env python3
"""Steepening of the inner shell of the white-hole runs at t = 30 across the seven-grid
ladder: the numbers of Sec. "Large c1: transient white holes" (in each group of runs the
shell width shrinks roughly as dr^0.8, |P| grows as dr^-0.5 and rho as dr^-1; the growth
slows on the finest grids, and on the finest run the width levels off once the refinement
is centred on the shell) and the finest run's t = 29 ... 31 time series.  Also printed, not
quoted in the paper: 10-16 points across the peak, |Theta C| = 28 / 49 / 85 on the uniform
grids, growing with resolution.

Runs (README §5, "Inner shell: the seven-grid ladder at t = 30"): positional run
directories, coarse to fine (the three uniform grids, then the zone-B refined ones);
--series is the finest run (F = 16) whose quarter-time dumps 29, 29.25, ..., 31 are listed.

Per run at t = --t: the shell = peak of |P| at r > 3, Theta < 35; spacing r[k+1] - r[k] at
the peak; FWHM of |P| (and points across it); |P|; within |r - r_k| < 0.2 the maxima of
rho = (P^2 + psi'^2)/2, T_in = (P + psi')^2, |psi'^2 - P^2| (the scalar field's
invariant d_a psi d^a psi), |m|, |Theta''|, |K|, |a_r|,
|Q|, and the integrals of rho and T_in.  Then the power-law exponents of each quantity
against the spacing (all grids, and the three finest) and, per family (the first
--n_uniform runs uniform, the others refined), of the width, |P| and rho (public-code
ladder: width dr^0.84 / dr^0.79, |P| dr^-0.50 / dr^-0.47, rho dr^-0.99 / dr^-0.95).

Usage:  wh_formation_ladder.py D1 D2 ... D7 [--t 30] [--series D7]
"""
import argparse
import os
import numpy as np
import spin0_common as sc


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("rundirs", nargs="+", help="run directories, coarse to fine")
    ap.add_argument("--t", type=float, default=30.0, help="dump time of the ladder (default 30)")
    ap.add_argument("--series", default=None,
                    help="run directory whose dumps at --series_times are listed")
    ap.add_argument("--series_times", default="29,29.25,29.5,29.75,30,30.25,30.5,30.75,31",
                    help="dump times of the time series (comma list)")
    ap.add_argument("--n_uniform", type=int, default=3,
                    help="the first n run directories are the uniform grids, the rest the "
                         "refined ones, for the per-family fits (default 3)")
    ap.add_argument("--window", type=float, default=0.2,
                    help="half-width in r of the window around the peak (default 0.2)")
    a = ap.parse_args()

    rows, labels = [], []
    print(f"t = {a.t:g}:  run  spacing@shell  width(pts)  |P|   rho   T_in   |psi'^2-P^2|   m_max"
          "   |Th''|   |K|  |a|  |Q|   |ThC|   int rho   int T_in")
    for d in a.rundirs:
        f = sc.dump_path(d, a.t)
        lab = os.path.basename(d.rstrip("/"))
        if not os.path.exists(f):
            print(f"  {lab}: missing dump at t = {a.t:g}")
            continue
        p = sc.load_profile(f)
        r, psi, P, K, A, Th, Q, C = (p[k] for k in ("r", "Psi", "P", "trK", "Acr", "Theta", "Q",
                                                     "ThetaC"))
        k = sc.peak_index(p, "P", rmin=3.0, thmax=35.0)
        dx = r[k + 1] - r[k]
        W, n = sc.fwhm(r, P, k)
        w = np.abs(r - r[k]) < a.window
        Phi = np.gradient(psi, r)
        dTh = np.gradient(Th, r)
        ddTh = np.gradient(dTh, r)
        m = sc.misner_sharp(p, dTh)
        rho = 0.5 * (P ** 2 + Phi ** 2)
        Tin = (P + Phi) ** 2
        Rs = np.abs(Phi ** 2 - P ** 2)

        def I(g, w=w):
            """integral over the window around the peak"""
            return np.trapezoid(g[w], r[w])

        rows.append((dx, W, abs(P[k]), rho[w].max(), Tin[w].max(), Rs[w].max(), np.abs(m[w]).max(),
                     np.abs(ddTh[w]).max(), I(rho), I(Tin)))
        labels.append(lab)
        print(f"  {lab:>20s} {dx:.5f}  {W:.5f}({n:3d})  {abs(P[k]):.3f}  {rho[w].max():.3g}  "
              f"{Tin[w].max():.3g}  {Rs[w].max():.4g}  {np.abs(m[w]).max():.1f}  "
              f"{np.abs(ddTh[w]).max():.3g}  {np.abs(K[w]).max():.2f} {np.abs(A[w]).max():.2f} "
              f"{np.abs(Q[w]).max():.2f}  {np.abs(C[w]).max():.3g}   {I(rho):.4g}  {I(Tin):.4g}")
    if len(rows) >= 2:
        R = np.array(rows)
        lx = np.log(R[:, 0])
        print(f"power-law exponents against the spacing ({len(rows)} grids; last three in brackets)"
              "   [paper: width ~ dr^0.8, |P| ~ dr^-0.5, rho ~ dr^-1]")
        for j, nm in [(1, "width"), (2, "|P|"), (3, "rho"), (4, "T_in"), (5, "|psi'^2-P^2|"),
                      (6, "m_max"), (7, "|Th''|"), (8, "int rho"), (9, "int T_in")]:
            pfit = np.polyfit(lx, np.log(R[:, j]), 1)[0]
            p3 = np.polyfit(lx[-3:], np.log(R[-3:, j]), 1)[0] if len(rows) >= 3 else np.nan
            print(f"  {nm:12s} ~ dr^{pfit:+.2f}  (dr^{p3:+.2f})")
        nu = a.n_uniform
        if 2 <= nu <= len(rows) - 2:
            print(f"per family (first {nu} runs uniform, the others refined)"
                  "   [paper: in each group width ~ dr^0.8, |P| ~ dr^-0.5, rho ~ dr^-1]")
            for j, nm in [(1, "width"), (2, "|P|"), (3, "rho")]:
                fu = np.polyfit(lx[:nu], np.log(R[:nu, j]), 1)[0]
                fr = np.polyfit(lx[nu:], np.log(R[nu:, j]), 1)[0]
                print(f"  {nm:12s} uniform dr^{fu:+.2f}   refined dr^{fr:+.2f}")
    if a.series:
        print(f"\ntime series of {a.series}:")
        for t in [float(x) for x in a.series_times.split(",")]:
            f = sc.dump_path(a.series, t)          # opened by the requested time
            if not os.path.exists(f):
                continue
            p = sc.load_profile(f)
            r, psi, P, Th, C = p["r"], p["Psi"], p["P"], p["Theta"], p["ThetaC"]
            k = sc.peak_index(p, "P", rmin=3.0, thmax=35.0)
            W, n = sc.fwhm(r, P, k)
            w = np.abs(r - r[k]) < a.window
            rho = 0.5 * (P ** 2 + np.gradient(psi, r) ** 2)
            print(f"  t = {t:5.2f}: shell r = {r[k]:.4f} Theta = {Th[k]:.3f} |P| = {abs(P[k]):.3f} "
                  f"width = {W:.5f} ({n} pts) rho_max = {rho[w].max():.3g} |ThC| = {np.abs(C[w]).max():.0f}")


if __name__ == "__main__":
    main()
