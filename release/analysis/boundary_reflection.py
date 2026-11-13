#!/usr/bin/env python3
"""Outer-boundary reflection tests of sec:infinity ("The spin-0 wave far from the source";
paper: "These conditions reflect very little").

Each test compares a run with the boundary at R against a reference run with the
boundary at 2R and the same dr: until the reference's own boundary signal returns,
any difference is what the boundary at R sends back.

  spin0    linear spin-0 pulse in flat space (--A0=0 --probe_amp=1e-6 --probe_r=50, fixed
           grid). Reflected w_in = (b_N/v0) trK + a_r at r = R - 10, relative to the
           incident w_out = -(b_N/v0) trK + a_r at the same station.  The coefficients
           a_N, b_N, v0 follow from --c1, --c2 (spin0_common.couplings, eq:aNbN; give
           the couplings of the runs).  If the run's stderr was saved as RUN_R/err.txt
           (./collapse ... 2> err.txt), its line "spin-0 speed v0 = ..." is a
           cross-check (both printed; a warning if they disagree beyond its 6 decimals).
  scalar   weak scalar shell in flat space (--A0=1e-4, fixed grid). |Theta dpsi/dTheta|
           along the ingoing light ray from the edge after the pulse hits it, relative to
           the incident peak at r = R - 1.
  late     collapse background (A0 = 0.35, moving mesh). Signal left after the scalar
           pulse has crossed the boundary at R: max |Theta dpsi/dTheta| difference against
           the reference in the Theta window (--theta_window) and the time window
           (--t_window), relative to the reference's maximum in the same window. Keep the
           Theta window outside the horizon.  The moving mesh coarsens toward the
           boundary, so the two runs do not share the spacing at large Theta: the result
           bounds the reflection from above.

Usage:
  boundary_reflection.py spin0  RUN_R RUN_2R [--R 150] [--c1 0.7] [--c2 -0.5119403]
  boundary_reflection.py scalar RUN_R RUN_2R [--R 150]
  boundary_reflection.py late   RUN_R RUN_2R [--theta_window 15,110] [--t_window 236,300]
"""
import argparse
import os
import sys
import numpy as np
import spin0_common as sc

SPIN0_STATION = 10.0     # spin0: w_in, w_out measured at r = R - 10
SCALAR_STATION = 1.0     # scalar: incident peak measured at r = R - 1
SCALAR_STRIDE = 3        # scalar: every 3rd integer dump after the hit
RAY_HALF_WIDTH = 6.0     # scalar: half-width in r of the strip around the ingoing ray
THETA_MARGIN = 5.0       # late: the reference is interpolated from a Theta window 5 wider


def dump_times(d):
    """Nominal dump times of run d (file names)."""
    return list(sc.list_dumps(d, exact=False)[0])


def err_coeffs(d):
    """(a_N, b_N, v0) as printed by the code on stderr (6 decimals), if the run's stderr was
    saved as <d>/err.txt; None otherwise."""
    path = os.path.join(d, "err.txt")
    if not os.path.exists(path):
        return None
    with open(path) as f:
        line = next((x for x in f if x.startswith("spin-0 speed v0")), None)
    if line is None:
        return None
    v0 = float(line.split("=")[1].split()[0])
    aN = float(line.split("aN=")[1].split()[0])
    bN = float(line.split("bN=")[1].rstrip(")\n "))
    return aN, bN, v0


def spin0_coeffs(a):
    """(a_N, b_N, v0) from the couplings --c1, --c2 (spin0_common.couplings, eq:aNbN).
    Cross-check against <run>/err.txt if present: both are printed if they differ by more
    than 1e-12 relative (always the case, since the code prints 6 decimals), with a warning
    if they differ beyond the printed precision (wrong --c1/--c2 for this run)."""
    cp = sc.couplings(a.c1, a.c2)
    mine = (cp["aN"], cp["bN"], cp["v0"])
    log = err_coeffs(a.run)
    if log is not None:
        rel = max(abs(x - y) / abs(x) for x, y in zip(mine, log))
        if rel > 1e-12:
            print(f"couplings c1={a.c1:.7g} c2={a.c2:.7g}: aN={mine[0]:.12g} bN={mine[1]:.12g} "
                  f"v0={mine[2]:.12g};  err.txt: aN={log[0]:.6f} bN={log[1]:.6f} v0={log[2]:.6f}")
        if max(abs(x - y) for x, y in zip(mine, log)) > 6e-7:
            print("WARNING: the couplings differ from err.txt beyond its 6 printed decimals: "
                  "check --c1/--c2 against the run", file=sys.stderr)
    return mine


def theta_dpsi(prof, n):
    """Theta dpsi/dTheta on the first n grid points (n = the grid of the run at R: the
    reference at 2R with the same dr is compared on those points)."""
    th, psi = prof["Theta"][:n], prof["Psi"][:n]
    return th * np.gradient(psi, th)


def spin0(a):
    """Reflected w_in against incident w_out at r = R - 10, before the reference's echo."""
    _, bN, v0 = spin0_coeffs(a)
    R, tK = a.R, bN / v0
    st = R - SPIN0_STATION
    # before the reference's own boundary signal comes back to the station
    ts = [t for t in dump_times(a.run) if t <= 2 * a.R - SPIN0_STATION]
    inc, refl = 0.0, 0.0
    for t in ts:
        ref, run = sc.load_dump(a.ref, t), sc.load_dump(a.run, t)
        i = np.argmin(abs(run["r"] - st))
        inc = max(inc, abs(-tK * ref["trK"][i] + ref["Acr"][i]))
        refl = max(refl, abs((tK * run["trK"][i] + run["Acr"][i])
                             - (tK * ref["trK"][i] + ref["Acr"][i])))
    print(f"spin-0 v0={v0:.4f} R={R:g}: incident max|w_out| at r={st:g}: {inc:.3e}; "
          f"reflected max|w_in|: {refl:.3e}; ratio {refl / inc:.2e}")


def scalar(a):
    """Signal along the ingoing light ray from the edge, after the pulse has hit it."""
    R = a.R
    ts = [int(t) for t in dump_times(a.ref) if t == int(t)]
    n = len(sc.load_dump(a.run, ts[0])["r"])
    inc, t_hit = 0.0, None
    for t in [t for t in ts if t < R + 10]:              # the pulse front reaches R at t ~ R
        ref = sc.load_dump(a.ref, t)
        i = np.argmin(abs(ref["r"][:n] - (R - SCALAR_STATION)))
        w = abs(theta_dpsi(ref, n))[i]
        if w > inc:
            inc, t_hit = w, t
    out = []
    for t in range(t_hit + SCALAR_STRIDE, ts[-1] + 1, SCALAR_STRIDE):
        ref, run = sc.load_dump(a.ref, t), sc.load_dump(a.run, t)
        m = abs(run["r"] - (R - (t - t_hit))) < RAY_HALF_WIDTH      # along the ingoing ray
        out.append(np.abs(theta_dpsi(run, n) - theta_dpsi(ref, n))[m].max() / inc)
    print(f"scalar R={R:g}, dr={R / (n - 1):.5f}: incident max|Theta dpsi/dTheta| at "
          f"r={R - SCALAR_STATION:g}: {inc:.3e} (t={t_hit})")
    print("  reflected/incident along the ingoing ray: " + " ".join(f"{x:.1e}" for x in out))
    print(f"  median {np.median(out):.2e}, max {max(out):.2e}")


def late(a):
    """Difference against the reference in a (Theta, t) window after the pulse has left."""
    th_lo, th_hi = a.theta_window
    t_lo, t_hi = a.t_window
    t_ref = set(dump_times(a.ref))
    ts = [t for t in dump_times(a.run) if t_lo <= t <= t_hi and t in t_ref]
    dmax, wmax, worst = 0.0, 0.0, None
    for t in ts:
        run, ref = sc.load_dump(a.run, t), sc.load_dump(a.ref, t)
        th = run["Theta"]
        keep = (th > th_lo) & (th < th_hi)
        w_run = theta_dpsi(run, len(th))[keep]
        thr, wr = ref["Theta"], theta_dpsi(ref, len(ref["Theta"]))
        kr = (thr > th_lo - THETA_MARGIN) & (thr < th_hi + THETA_MARGIN)
        w_ref = np.interp(th[keep], thr[kr], wr[kr])
        d = np.abs(w_run - w_ref).max()
        wmax = max(wmax, np.abs(w_ref).max())
        if d > dmax:
            dmax, worst = d, t
    print(f"late signal, {t_lo:g} <= t <= {t_hi:g}, {len(ts)} dumps, {th_lo:g} < Theta < "
          f"{th_hi:g}: max|difference| {dmax:.3e} (t={worst}), max|reference| {wmax:.3e}; "
          f"relative {dmax / wmax:.1%}")


def pair(s):
    """'lo,hi' -> (lo, hi) as floats."""
    lo, hi = (float(x) for x in s.split(","))
    return lo, hi


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("test", choices=["spin0", "scalar", "late"])
    ap.add_argument("run", help="run with the boundary at R")
    ap.add_argument("ref", help="reference run with the boundary at 2R, same dr")
    ap.add_argument("--R", type=float, default=150.0,
                    help="spin0, scalar: boundary radius of RUN_R (default 150)")
    ap.add_argument("--c1", type=float, default=sc.C1_DEFAULT,
                    help=f"spin0: aether coupling c1 of the runs (default {sc.C1_DEFAULT:g})")
    ap.add_argument("--c2", type=float, default=sc.C2_BASELINE,
                    help=f"spin0: aether coupling c2 of the runs (default {sc.C2_BASELINE:.7g}, "
                         "v0 = 1; 0 for v0 = 1.267)")
    ap.add_argument("--theta_window", type=pair, default=(15.0, 110.0),
                    help="late: areal-radius window lo,hi (default 15,110)")
    ap.add_argument("--t_window", type=pair, default=(236.0, 300.0),
                    help="late: time window lo,hi (default 236,300)")
    a = ap.parse_args()
    {"spin0": spin0, "scalar": scalar, "late": late}[a.test](a)


if __name__ == "__main__":
    main()
