#!/usr/bin/env python3
"""Amplitude scan of the time-symmetric shell at N = 6400 (dr = 0.023): the passage of
Sec. "Convergence" on what happens above A0 = 0.06.  On this grid A0 = 0.05-0.06 disperse,
A0 = 0.07-0.13 collapse the centre into a cell-sized frozen throat (Theta_AH ~ 0.1-1, the
resolution artefact of the white-hole section; 0.12 disperses), and only A0 = 0.14-0.15
form a genuine detached horizon.  The threshold itself is located by bisection on finer
grids with bisect_bh_threshold.py and fig_bh_threshold.py (early horizon from A0 = 0.135 /
0.130 / 0.125 / 0.125 at N = 12801 / 25601 / 51201 / 102401: A* in 0.120-0.125 on the two
finest grids).

Runs (README §5, "Amplitude scan"): N = 6400, r_max = 150, t <= 60, dumps every 1,
A0 = 0.05 ... 0.15 (one directory per amplitude, given as positional arguments).

For each run: the minimum over all dumps of theta_+ = Theta' + (3Q - 2K) Theta/6 at
r > --rmin (0.5), the central lapse at the last dump, the dump times at which theta_+ < 0
at r > --rmin, the areal radius Theta_AH of the outermost zero of theta_+ at r > 0.5 at the
last dump (spin0_common.apparent_horizon: the cut is on the proper radius r, the criterion
of every script, no cut on Theta is applied), and, if the run's stdout was saved as
--log (log.txt), the first/last time of the code's "Horizon detected" flag.  A frozen central throat (Theta_AH ~ 0.1) is
NOT a horizon: read Theta_AH(end) > 0.5 as a detached horizon, the criterion of
bisect_bh_threshold.py (fig:threshold).

Usage:  scan_bh_threshold.py RUNDIR [RUNDIR ...] [--log log.txt]
"""
import argparse
import os
import re
import numpy as np
import spin0_common as sc


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("rundirs", nargs="+", help="run directories (one per amplitude)")
    ap.add_argument("--rmin", type=float, default=0.5,
                    help="exclude r < rmin from the min of theta_+ (default 0.5, the cut of "
                         "spin0_common.apparent_horizon)")
    ap.add_argument("--log", default="log.txt",
                    help="name of the saved stdout of the run inside each directory (optional)")
    a = ap.parse_args()
    print(f"{'run':>14} {'min theta+':>11} {'N0(end)':>8} {'Theta_AH(end)':>14}  "
          "theta+<0 at t   |  code's 'Horizon detected' first/last t")
    for d in a.rundirs:
        t, paths = sc.list_dumps(d, exact=False)
        rows = []
        for tt, f in zip(t, paths):
            p = sc.load_profile(f)
            tp, _ = sc.theta_pm(p, np.gradient(p["Theta"], p["r"]))
            m = p["r"] > a.rmin
            rows.append((tt, tp[m].min(), p["NN"][0]))
        th_ah = sc.apparent_horizon(p)[2]            # last dump
        neg = [f"{x[0]:g}" for x in rows if x[1] < 0]
        neg_s = (", ".join(neg[:4]) + (f", ... {neg[-1]}" if len(neg) > 4 else "")) if neg else "-"
        hz = None
        lf = os.path.join(d, a.log)
        if os.path.exists(lf):
            hz = re.findall(r"Horizon detected.*t=([\d.]+)", open(lf).read())
        hz_s = f"{hz[0]} / {hz[-1]}" if hz else ("-" if hz is not None else "(no log)")
        print(f"{os.path.basename(d.rstrip('/')):>14} {min(x[1] for x in rows):11.3f} "
              f"{rows[-1][2]:8.3f} {th_ah:14.3f}  {neg_s:<16} |  {hz_s}")


if __name__ == "__main__":
    main()
