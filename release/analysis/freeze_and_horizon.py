#!/usr/bin/env python3
"""Central-lapse freezing times and apparent-horizon histories of the white-hole /
black-hole runs at c1 = 0.8 (Sec. "Large c1: transient white holes": the freezing time of
the frozen centre moves later on every finer grid, roughly as dr^-1/3, with no sign of
converging; at A0 = 0.12 a black hole forms.  Public-code runs, not quoted in the paper:
freezing at t ~ 53, 63, 77, 101, 135 for dr = 0.0625, 0.031, 0.016, 0.0078, 0.0039;
Theta_AH = 4.85-5.07 for t > 30 at A0 = 0.12).

Runs (README §5, "Frozen centre vs resolution" and "White hole vs black hole"):
  --freeze DIR[:dr] ...   A0 = 0.10 runs with dumps every 2 to t = 150: prints the central
                          lapse N(r = 0) at every dump and the last dump time at which it
                          exceeds 0.1, 0.01, 0.001 (the centre "freezes for good" between
                          that time and the next dump; dr is only a label);
  --horizon DIR ...       A0 = 0.12 (or any) runs: prints the areal radius of the outermost
                          apparent horizon (theta_+ = 0, linear interpolation) at every
                          dump, its minimum and maximum for t > --t_from (30).

Usage:  freeze_and_horizon.py --freeze F1:0.0625 F2:0.031 ... --horizon H1 H2
"""
import argparse
import os
import numpy as np
import spin0_common as sc


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--freeze", nargs="*", default=[],
                    help="run directories of the freeze ladder, optionally DIR:dr")
    ap.add_argument("--horizon", nargs="*", default=[],
                    help="run directories whose Theta_AH(t) is printed")
    ap.add_argument("--thresholds", default="0.1,0.01,0.001",
                    help="central-lapse thresholds (comma list)")
    ap.add_argument("--t_from", type=float, default=30.0,
                    help="Theta_AH range quoted for t > this (default 30)")
    a = ap.parse_args()
    thr = [float(x) for x in a.thresholds.split(",")]

    if a.freeze:
        print("== central lapse N0(t) of the freeze ladder ==   [paper: later on every finer grid, "
              "~ dr^-1/3; public-code runs: t ~ 53, 63, 77, 101, 135 for dr = 0.0625 ... 0.0039]")
    for spec in a.freeze:
        d, _, dr = spec.partition(":")
        t, paths = sc.list_dumps(d, exact=False)
        n0 = np.array([sc.load_profile(f)["NN"][0] for f in paths])
        print(f"{os.path.basename(d.rstrip('/'))}" + (f"  dr = {dr}" if dr else "")
              + f"  (last dump t = {t[-1]:g})")
        print("   " + " ".join(f"{tt:.0f}:{n:.2g}" for tt, n in zip(t, n0)))
        for th in thr:
            above = np.where(n0 > th)[0]
            if len(above) == 0:
                s = "never"
            else:
                k = above[-1]                   # last dump with N0 above the threshold
                s = f"{t[k]:g}"
                if k + 1 < len(t):
                    s += f" (N0 < {th:g} from t = {t[k + 1]:g} on)"
                else:
                    s += " (still above at the last dump: no freeze)"
            print(f"   last t with N0 > {th:g}: {s}")

    if a.horizon:
        print("== apparent horizon Theta_AH(t) (outermost theta_+ = 0) ==   [not in the paper; "
              "public-code runs, A0 = 0.12: 4.85-5.07 for t > 30, ~4% in time, ~2% between "
              "N = 25600 and 51200]")
    for d in a.horizon:
        t, paths = sc.list_dumps(d, exact=False)
        out = [(tt, sc.apparent_horizon(sc.load_profile(f))[2]) for tt, f in zip(t, paths)]
        print(os.path.basename(d.rstrip("/")) + "  " + " ".join(f"{tt:g}:{x:.3f}" for tt, x in out
                                                                 if np.isfinite(x)))
        late = [(tt, x) for tt, x in out if tt > a.t_from and np.isfinite(x)]
        if late:
            lo = min(late, key=lambda q: q[1])
            hi = max(late, key=lambda q: q[1])
            print(f"   t > {a.t_from:g}: min {lo[1]:.3f} at t = {lo[0]:g}, max {hi[1]:.3f} at t = {hi[0]:g}, "
                  f"last {late[-1][1]:.3f} at t = {late[-1][0]:g}")


if __name__ == "__main__":
    main()
