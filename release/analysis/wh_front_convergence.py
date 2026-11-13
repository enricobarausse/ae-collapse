#!/usr/bin/env python3
"""Convergence of the escaping white-hole front with zone-A refinement locked on before
the compression: the escaping-pulse numbers of Sec. "Large c1: transient white holes"
(paper: the pulse width shrinks roughly as exp(-0.15 t); the thin anti-trapped layer on
the pulse narrows from about 0.9 at t = 40 to 0.06 at t = 80 on these refined grids, which
agree on its width to within 20 %, while its depth min theta_- grows with resolution.
Public-code runs, not in the paper: front width and position converge at about first
order to t ~ 70, the peak |P| more slowly; the constraint residual at the front decreases
with resolution up to t ~ 70).  These are the runs of panel (b) of fig:sharp.

Runs (README §5, "Escaping pulse with zone A locked on"): c1 = 0.8, c2 = -0.7111111,
ingoing, A0 = 0.10, r_max = 150, base grids N = 12800 / 25600 / 51200 with --ref_F=8
--ref_ton=40 --ref_w=3, dumps every 0.5 from t = 40 to 80.  Three run directories, coarse
to fine, as positional arguments.

Front tracking: at the first dump (t = 40) the front is the maximum of |P| with
Theta in --theta_window (default 15..140, i.e. outside the inner shell); at every later
dump it is the LARGEST local maximum of |P| within a window of half-width 1 + 0.5 dt in
Theta around the position predicted from the previous dump, Theta_prev + v dt with
v = --v_front (1, the front moves outward at about unit speed in Theta).  Fallbacks: if no
grid point lies in that window, the whole --theta_window is searched; if the window holds
no local maximum, the largest |P| in it is taken.  This keeps the track on the front when
the inner shell, which has a larger |P| on the coarse grids after t ~ 45, would otherwise
be picked by a global search.

At each time of --times (default 40, 45, ..., 80): the grid spacing at the front, its
areal radius Theta_f, |P|, FWHM W (and points across it), the anti-trapped layer on the
front (lay = areal-radius extent of the points with theta_- < 0 within |Theta - Theta_f|
< 1, i.e. between the outermost zeros of theta_- around the peak; on the finest grids
the layer can contain a small theta_- > 0 island near the peak, which this definition
spans), and within |r - r_front| < 3: max |K|, min theta_-, max m, the scalar energy
integral E = int Theta^2 (P^2 + psi'^2) dr and max |Theta C|.  The spread of the layer
width across the three grids, (max - min)/min, is printed at each time.  For the three
grids the successive differences of Theta_f, W, |P|,
max|K| and E, the observed order p that solves (h1^p - h2^p)/(h2^p - h3^p) = d12/d23
with the actual spacings h at the front (log2 of the ratio when the spacings halve
exactly), and the Richardson value z + (z - y)/((h2/h3)^p - 1).  At the end, per grid, the
least-squares rate k of W ~ exp(-k t) over the times >= --fit_from (default 50; before
that the FWHM measures the broad hump the pulse sits on), with and without the last time.

Public-code runs (whconv_N12800 / N25600 / N51200): layer 0.89 / 0.90 / 0.84 at t = 40 and
0.059 / 0.062 / 0.067 at t = 80, spread 7-18 % (<= 20 %) at every time; min theta_- -0.21 /
-0.25 / -0.28 at t = 40, -0.85 / -1.12 / -1.36 at t = 80; k = 0.12 / 0.13 / 0.14 over t = 50-80 and
0.14 / 0.15 / 0.16 over t = 50-75 (paper: roughly exp(-0.15 t)).

Usage:  wh_front_convergence.py COARSE MID FINE [--times 40,45,50,55,60,65,70,75,80]
                                [--fit_from 50]
"""
import argparse
import os
from collections import namedtuple
import numpy as np
import spin0_common as sc

# one grid at one time: N, spacing dr at the front, Theta_f, |P|, FWHM W (and points),
# anti-trapped layer width, then within the r window: max|K|, min theta_-, max m, E, max|ThC|
Row = namedtuple("Row", "N dr Theta_f P W npts layer maxK min_thm max_m E maxThC")


def track_front(rundir, th_lo, th_hi, v_front):
    """dict t -> index of the front at every dump of rundir (see the module docstring)."""
    t, paths = sc.list_dumps(rundir, exact=False)
    idx, th_prev, t_prev = {}, None, None
    for tt, f in zip(t, paths):
        p = sc.load_profile(f)
        Th, P = p["Theta"], np.abs(p["P"])
        if th_prev is None:
            j = sc.peak_index(p, "P", rmin=-np.inf, thmin=th_lo, thmax=th_hi)
        else:
            dt = tt - t_prev
            pred = th_prev + v_front * dt
            hw = 1.0 + 0.5 * dt
            cand = np.where((Th > pred - hw) & (Th < pred + hw))[0]
            if len(cand) == 0:
                cand = np.where((Th > th_lo) & (Th < th_hi))[0]
            # local maxima of |P| among the candidates
            loc = [k for k in cand if 0 < k < len(P) - 1 and P[k] >= P[k - 1] and P[k] >= P[k + 1]]
            if not loc:
                loc = list(cand)
            j = max(loc, key=lambda k: P[k])
        idx[tt] = j
        th_prev, t_prev = Th[j], tt
    return idx


def observed_order(h, f):
    """Order p from three values f on spacings h (h1 > h2 > h3): solves
    (h1^p - h2^p)/(h2^p - h3^p) = (f1 - f2)/(f2 - f3) by bisection; nan if no solution."""
    d12, d23 = f[0] - f[1], f[1] - f[2]
    if d23 == 0 or d12 / d23 <= 0:
        return np.nan
    target = d12 / d23

    def g(p):
        return (h[0] ** p - h[1] ** p) / (h[1] ** p - h[2] ** p) - target

    lo, hi = 0.02, 12.0
    if g(lo) * g(hi) > 0:
        return np.nan
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        if g(lo) * g(mid) <= 0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("rundirs", nargs=3, help="run directories, coarse, medium, fine")
    ap.add_argument("--times", default="40,45,50,55,60,65,70,75,80", help="dump times (comma list)")
    ap.add_argument("--fit_from", type=float, default=50.0,
                    help="first time of the exp(-k t) fit of the front width (default 50)")
    ap.add_argument("--theta_window", default="15,140",
                    help="Theta range of the initial front search (default 15,140)")
    ap.add_argument("--v_front", type=float, default=1.0,
                    help="predicted outward speed of the front in Theta per unit t (default 1)")
    ap.add_argument("--window", type=float, default=3.0,
                    help="half-width in r of the window around the front (default 3)")
    a = ap.parse_args()
    th_lo, th_hi = (float(x) for x in a.theta_window.split(","))
    tracks = [track_front(d, th_lo, th_hi, a.v_front) for d in a.rundirs]
    widths = {d: [] for d in a.rundirs}
    for t in [float(x) for x in a.times.split(",")]:
        out = []
        for d, tr in zip(a.rundirs, tracks):
            f = sc.dump_path(d, t)
            if not os.path.exists(f) or t not in tr:
                continue
            p = sc.load_profile(f)
            r, psi, P, K, Th, C = (p[k] for k in ("r", "Psi", "P", "trK", "Theta", "ThetaC"))
            dTh = np.gradient(Th, r)
            m = sc.misner_sharp(p, dTh)
            _, thm = sc.theta_pm(p, dTh)
            j = tr[t]
            W, n = sc.fwhm(r, P, j)
            w = np.abs(r - r[j]) < a.window
            E = np.trapezoid((Th ** 2 * (P ** 2 + np.gradient(psi, r) ** 2))[w], r[w])
            lay = np.where((np.abs(Th - Th[j]) < 1.0) & (thm < 0))[0]
            lw = Th[lay].max() - Th[lay].min() if len(lay) else 0.0
            widths[d].append((t, W))
            out.append(Row(len(r), r[j + 1] - r[j], Th[j], abs(P[j]), W, n, lw,   # dr at the front
                           np.abs(K[w]).max(), thm[w].min(), m[w].max(), E, np.abs(C[w]).max()))
        print(f"\nt = {t:g}")
        for o in out:
            print("  N=%6d dr=%.5f Th_f=%6.3f |P|=%.4f W=%.4f (%4d pts) lay=%.4f max|K|=%.4f "
                  "min th-=%+.3f max m=%6.2f E=%.4f |ThC|=%.3f" % tuple(o))
        lws = [o.layer for o in out]
        if len(lws) > 1 and min(lws) > 0:
            print(f"    anti-trapped layer width: spread across grids "
                  f"{100 * (max(lws) - min(lws)) / min(lws):.0f} %   "
                  "[paper: about 0.9 at t = 40, 0.06 at t = 80, within 20 %]")
        if len(out) == 3:
            h = [o.dr for o in out]
            for k, nm in [("Theta_f", "Theta_f"), ("W", "W"), ("P", "|P|"), ("maxK", "max|K|"),
                          ("E", "E")]:
                x, y, z = (getattr(o, k) for o in out)
                pp = observed_order(h, (x, y, z))
                rich = z + (z - y) / ((h[1] / h[2]) ** pp - 1) if np.isfinite(pp) else float("nan")
                print(f"    {nm:8s}: differences {x - y:+.4g}, {y - z:+.4g}  (spacings "
                      f"{h[0]:.5f}/{h[1]:.5f}/{h[2]:.5f})  -> order p = {pp:.2f}  Richardson {rich:.4g}")
    print(f"\nfront width W ~ exp(-k t), least squares over t >= {a.fit_from:g}   [paper: exp(-0.15 t)]")
    for d in a.rundirs:
        tw = np.array([x for x in widths[d] if x[0] >= a.fit_from - 1e-9])
        line = []
        for n in sorted({len(tw), len(tw) - 1}):    # all fitted times, and without the last
            if n >= 2:
                k = -np.polyfit(tw[:n, 0], np.log(tw[:n, 1]), 1)[0]
                line.append(f"k = {k:.3f} (t = {tw[0, 0]:g} .. {tw[n - 1, 0]:g})")
        print(f"  {os.path.basename(d.rstrip('/')):>20}: " + ",  ".join(line))


if __name__ == "__main__":
    main()
