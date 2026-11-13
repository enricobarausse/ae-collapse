#!/usr/bin/env python3
"""Black-hole formation threshold of the time-symmetric shell by bisection in A0
(Sec. "Convergence", fig:threshold: the early horizon appears from A0 = 0.135 / 0.130 / 0.125 /
0.125 at N = 12801 / 25601 / 51201 / 102401, so A* lies in 0.120-0.125 on the two finest grids;
the figure itself is made by fig_bh_threshold.py).

Runs (README §5, "Black-hole threshold"): r_max = 300, t <= 60, dumps every 2, one
directory per (N, A0), e.g. bis_N12801_A0.130 (the amplitude is read from --A0= in the
run's cmd.txt, an optional copy of the command line, if present, otherwise from the
directory name, "A<value>"; spin0_common.run_amplitude).

Classification: a run counts as "black hole" if at some dump the outermost zero crossing of
theta_+ = Theta' + (3Q - 2K) Theta/6 at proper radius r > 0.5 (spin0_common.apparent_horizon,
the criterion of every script, analysis/README.md "Horizon criteria") lies at
Theta_AH > --th_min (0.5), i.e. a trapped surface detached from the centre (the Theta cut
is the black-hole verdict, not part of the horizon finder: a frozen central throat can sit
at r > 0.5 with Theta ~ 0.1).  A frozen central throat and dispersal both count as "no
black hole".  Near the switch this single criterion is ambiguous (the paper's
fig:threshold): the central region that becomes trapped at t ~ 26-28 when the lapse
collapses can reach Theta ~ 0.5-1 on the finer grids, so --t_early restricts the search
to the early horizon (t <= 20) when asked.  The code's stdout flag "Horizon detected" is NOT
used (it latches onto transient focusing artefacts and onto sign changes near the r = 45
search edge).

Prints, per run: N, A0, the largest Theta_AH found (and when), the central lapse at the
last dump searched (N0(t<=t_early)) and at the last dump of the run (N0(end)), the time
of the last dump of the run and the verdict; then, per N, the bracketing amplitudes
(largest A0 without a black hole, smallest with one), and whether the runs have reached
t_final.

Usage:  bisect_bh_threshold.py RUNDIR [RUNDIR ...] [--t_early 20]
"""
import argparse
import os
import numpy as np
import spin0_common as sc


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("rundirs", nargs="+", help="run directories (one per N, A0)")
    ap.add_argument("--th_min", type=float, default=0.5,
                    help="a horizon counts only if Theta_AH > this (default 0.5)")
    ap.add_argument("--t_early", type=float, default=np.inf,
                    help="only dumps at t <= this are searched (e.g. 20 for the early horizon, "
                         "as in fig_bh_threshold.py)")
    a = ap.parse_args()
    rows = []
    print(f"{'run':>22} {'N':>6} {'A0':>7} {'max Theta_AH':>12} {'at t':>6} {'N0(t<=t_early)':>14} "
          f"{'N0(end)':>8} {'t_end':>6}  verdict")
    for d in a.rundirs:
        if not os.path.isdir(d):
            print(f"{os.path.basename(d.rstrip('/')):>22}  (not a directory, skipped)")
            continue
        try:
            t, paths = sc.list_dumps(d, exact=False)
        except FileNotFoundError:
            print(f"{os.path.basename(d.rstrip('/')):>22}  (no dumps yet)")
            continue
        best, t_best, n0, N = np.nan, np.nan, np.nan, 0
        for tt, f in zip(t, paths):
            if tt > a.t_early:
                continue
            p = sc.load_profile(f)
            N = len(p["r"])
            n0 = p["NN"][0]
            th = sc.apparent_horizon(p)[2]          # zeros at r > 0.5
            if np.isfinite(th) and th > a.th_min and (np.isnan(best) or th > best):
                best, t_best = th, tt
        n0_end = sc.load_profile(paths[-1])["NN"][0]      # central lapse at the last dump
        A0 = sc.run_amplitude(d)
        bh = np.isfinite(best)
        rows.append((N, A0, bh, t[-1]))
        print(f"{os.path.basename(d.rstrip('/')):>22} {N:6d} {A0:7.4f} {best:12.3f} {t_best:6.4g} "
              f"{n0:14.3f} {n0_end:8.3g} {t[-1]:6g}  {'black hole' if bh else 'no black hole (throat/dispersal)'}")
    print("\nbracketing amplitudes per N (largest A0 without a black hole < A* < smallest with one):")
    for N in sorted(set(r[0] for r in rows)):
        R = [r for r in rows if r[0] == N]
        no = [r[1] for r in R if not r[2]]
        yes = [r[1] for r in R if r[2]]
        lo = f"{max(no):.4f}" if no else "-"
        hi = f"{min(yes):.4f}" if yes else "-"
        tmin = min(r[3] for r in R)
        print(f"  N = {N}: {lo} < A* < {hi}   ({len(R)} runs, earliest last dump t = {tmin:g})")


if __name__ == "__main__":
    main()
