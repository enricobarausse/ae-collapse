#!/usr/bin/env python3
"""Trapped / anti-trapped classification of the white-hole family runs (Sec. "Large c1:
transient white holes": amplitude dependence, the A0 = 0.10 anti-trapped history and the
other shell families (r0, sigma) = (10, 3) and (15, 4)): for each run, the dumps at
which the classification changes, with the areal-radius ranges of the trapped
(theta_+ < 0 < theta_-) and anti-trapped (theta_- < 0 < theta_+) regions and the central
lapse.  A region made of several disjoint pieces is printed as a comma list of ranges
(e.g. the A0 = 0.10 anti-trapped region after it splits).

Runs (README §5): the A0 = 0.10 / 0.12 runs at c1 = 0.8 and the (r0, sigma) = (10, 3) and
(15, 4) scans and ladders.  Positional run directories (shell globs are fine).  If the
run's stderr was saved as --err (err.txt) inside the directory, a FATAL stop is flagged.
With --extent the ranges are printed at EVERY dump (paper, A0 = 0.12: theta_+ < 0 by
t ~ 2.5; A0 = 0.10: theta_+ dips below zero at t ~ 6-17, the anti-trapped region appears
at t ~ 20, splits at t ~ 36, its inner part shrinks after t ~ 40 and is gone by t ~ 44;
public-code run wh_wh_N25600: trapped 6.0-6.1 at t = 6 to 4.3 at t = 16, anti-trapped
4.2-4.6 at t = 20, 2.7-17.6 at t = 36, split into 2.3-13.9 and 18.7-19.9 at t = 38,
4.1-12.9 and 21.3-22.2 at t = 40, only the layer 23.7-24.4 at t = 42).  Paper, other families: sigma = 3
disperses at A0 = 0.01, 0.015, 0.03; A0 = 0.08 trapped then untrapped at dr ~ 0.009 but a
black hole at dr ~ 0.0047; r0 = 15, A0 = 0.028 trapped region gone by t ~ 20, anti-trapped
from t ~ 24 on three grids; the sigma = 3 central collapse at A0 = 0.03 occurs at
dr = 0.0031 and 0.0008 but not 0.0016, and at A0 = 0.04 later on each finer grid.

Usage:  family_transitions.py RUNDIR [RUNDIR ...] [--extent]
"""
import argparse
import os
import numpy as np
import spin0_common as sc


def classify(path, rmin):
    p = sc.load_profile(path)
    r, Th = p["r"], p["Theta"]
    tp, tm = sc.theta_pm(p, np.gradient(Th, r))
    m = r > rmin
    tp, tm, Th = tp[m], tm[m], Th[m]
    trap = (tp < 0) & (tm > 0)
    anti = (tm < 0) & (tp > 0)
    def rng(k):
        """areal-radius range of each contiguous run of grid points with k True."""
        e = np.flatnonzero(np.diff(np.concatenate([[0], k.astype(int), [0]])))
        return ",".join(f"{Th[i:j].min():.1f}-{Th[i:j].max():.1f}"
                        for i, j in zip(e[::2], e[1::2])) or "-"
    return rng(trap), rng(anti), p["NN"][0]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("rundirs", nargs="+", help="run directories")
    ap.add_argument("--rmin", type=float, default=0.3,
                    help="exclude r < rmin (the central throat) from the classification")
    ap.add_argument("--err", default="err.txt",
                    help="name of the saved stderr inside each directory (optional)")
    ap.add_argument("--extent", action="store_true",
                    help="print the trapped / anti-trapped extent at every dump")
    a = ap.parse_args()
    for d in a.rundirs:
        if not os.path.isdir(d):
            print(f"== {os.path.basename(d.rstrip('/'))} (not a directory, skipped)")
            continue
        t, paths = sc.list_dumps(d, exact=False)
        ef = os.path.join(d, a.err)
        fatal = " FATAL" if os.path.exists(ef) and "FATAL" in open(ef).read() else ""
        prev, out = None, []
        for tt, f in zip(t, paths):
            tr, an, n0 = classify(f, a.rmin)
            key = (tr != "-", an != "-")
            if a.extent or key != prev:
                out.append(f"t={tt:g}: trapped={tr} anti={an} N0={n0:.2g}")
                prev = key
        print(f"== {os.path.basename(d.rstrip('/'))}{fatal} (last dump t = {t[-1]:g})")
        for o in out:
            print("    " + o)


if __name__ == "__main__":
    main()
