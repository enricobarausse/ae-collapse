#!/usr/bin/env python3
"""Inner-shell steepening in the other shell families (Sec. "Large c1: transient white
holes": intermediate data form a transient white hole whose inner shell steepens as for
the baseline shell; public-code ladders, not quoted in the paper: when the shell forms,
t = 20-22 for (r0, sigma) = (10, 3), A0 = 0.06 and t = 24-26 for (15, 4), A0 = 0.020, it
spans 13/14/17 and 14/15/18 points on the three grids, and |P| at the shell grows by ~1.3
per halving of dr).

Runs (README §5, "Other initial data"): the uniform ladders N = 19200 / 38400 / 76800 on
r_max = 60 (dx = 0.0031 / 0.0016 / 0.0008), dumps every 1 from t = 20 to 45.  Positional run
directories, coarse to fine.

For every dump time present in any run and for each run: the areal radius of the |P| peak
in 2 < Theta < 35 and r > 1 (--theta_window, --rmin: r > 1 excludes the centre, Theta < 35
the outgoing front), |P|, its FWHM (and points), rho = (P^2 + psi'^2)/2 and |Theta C| at
the peak, and the minimum of the lapse.

Usage:  wh_family_ladder.py D1 D2 D3
"""
import argparse
import os
import numpy as np
import spin0_common as sc


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("rundirs", nargs="+", help="run directories, coarse to fine")
    ap.add_argument("--theta_window", default="2,35",
                    help="Theta range in which the shell is searched (default 2,35)")
    ap.add_argument("--rmin", type=float, default=1.0, help="exclude r < rmin (default 1)")
    a = ap.parse_args()
    th_lo, th_hi = (float(x) for x in a.theta_window.split(","))
    dirs, dump_t = [], []
    for d in a.rundirs:
        try:
            dump_t.append(set(sc.list_dumps(d, exact=False)[0]))
            dirs.append(d)
        except FileNotFoundError:
            print(f"(skipping {d}: no dumps yet)")
    a.rundirs = dirs
    names = [os.path.basename(d.rstrip("/")) for d in a.rundirs]
    times = sorted(set().union(*dump_t))
    print("runs:", names)
    print("   t  | per run: Theta_pk  |P|   FWHM(pts)   rho   |ThC|@pk   N_min")
    for t in times:
        row = []
        for d in a.rundirs:
            f = sc.dump_path(d, t)
            if not os.path.exists(f):
                row.append("--")
                continue
            p = sc.load_profile(f)
            r, psi, P, Th, N, C = (p[k] for k in ("r", "Psi", "P", "Theta", "NN", "ThetaC"))
            j = sc.peak_index(p, "P", rmin=a.rmin, thmin=th_lo, thmax=th_hi)
            w, npts = sc.fwhm(r, P, j)
            rho = 0.5 * (P[j] ** 2 + np.gradient(psi, r)[j] ** 2)
            row.append(f"{Th[j]:6.2f} {abs(P[j]):.3f} {w:.4f}({npts:3d}) {rho:.3f} {abs(C[j]):7.1f} {N.min():.2f}")
        print(f"{t:5.1f} | " + " | ".join(row))


if __name__ == "__main__":
    main()
