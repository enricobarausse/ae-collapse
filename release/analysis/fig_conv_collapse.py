#!/usr/bin/env python3
"""Collapse convergence ladder: paper fig:conv-collapse (conv_collapse.pdf), the collapse
row of tab:conv-scatter (the orders p), tab:conv-ratio (constraint-norm ratios),
Theta_AH(t = 28) with its Richardson limit, and the ADM-mass drift (sec:conv, paragraph
"Collapse", and Appendix sec:app-conv).

Runs (README §5, "Collapse convergence ladder"): fixed grid (--ale=0), r_max = 300,
A0 = 0.35, nested N = 12801 / 25601 / 51201, dumps every 2 up to t = 28.  Three run
directories, coarse to fine, as positional arguments.

Definitions:
  norms        ||.||_2 over r in [r_AH + 2, --rmax_norm] (default 150), with r_AH the largest
               of the three apparent-horizon radii at that time (outermost zero of
               theta_+ = Theta' + (3Q - 2K) Theta/6); i.e. OUTSIDE the horizon, where the
               slicing has not frozen;
  order        p = log2(||f_h - f_h/2|| / ||f_h/2 - f_h/4||), differences taken pointwise
               on the coarse grid (nested grids), averaged over the dump times >= --tmin
               (default 4) present on all three grids; for Theta C the log2 of the two
               successive norm ratios;
  ratios       ||Theta C||_coarse / ||Theta C||_mid and mid / fine at --ratio_times
               (4 for second order);
  Theta_AH     linear-interpolated zero of theta_+ at t = --t_ah (28), Richardson limit
               assuming second order, and the observed order from the three values;
  M_ADM        Misner-Sharp mass m(Theta) at the grid point nearest Theta = --th_adm (250)
               at t = 2 and t_ah.  On the fixed grid the outer edge falls inward (Theta(r_max)
               = 284 at t = 2, ~235 at t = 28), so from t ~ 20 the sphere Theta = 250 is off
               the grid and the value is read at the GRID EDGE (the printout gives the Theta
               actually used).  The paper's M_ADM drift is read at Theta = 250 on the moving
               mesh: pass that run with --mesh_run (README: mass_hi, N = 51201, --ale=1).
               With --mesh_run the script also prints: the first dump at which the sphere
               --th_adm has left the fixed grid; m at t_ah on the fixed-grid edge sphere on
               the three fixed grids and on the moving mesh (their relative differences);
               the change of m between t = 2 and t_ah on the spheres --drift_spheres
               (default 100, 150, 200, 250) of the moving mesh, and the fit
               A + B/Theta^3 of that change, A being its extrapolation to Theta -> infinity.
  history      on the finest grid at t = 2 and t_ah: Theta_AH, r_AH, the central lapse and
               the areal radius at r = 150 (Appendix "Convergence tests").

Paper (public-code runs): p = (2.62, 2.55, 2.34, 2.38, 2.62, 2.42, 2.24) for psi, P, K,
a_r, Theta, Q, ||Theta C||; ratios 3.69-4.69 (coarse pair) and 3.95-6.26 (fine pair);
Theta_AH = 12.293 / 12.604 / 12.677, differences shrinking by 4.3, Richardson 12.70;
M_ADM(Theta = 250, moving mesh) 3.70836 -> 3.71107 (+0.07 %); Theta = 250 leaves the fixed
grid at t ~ 20 (public-code run: Theta(r_max) = 253.9 / 249.8 at t = 18 / 20); at t = 28
the fixed-grid edge is at Theta ~ 235, where the three fixed
grids agree to 1e-4 and the fixed grid and the moving mesh to 2e-5 (the fixed-grid ladder
gives 3.70836 -> 3.71143 / 3.71164 / 3.71168 there, +0.08-0.09 %); the change of m falls off as
Theta^-3 on Theta = 100-250 (1.2 % at Theta = 100) and its extrapolation to infinity is
about 0.01 % (public-code run: 1.217 % at Theta = 100, A = -0.0102 % for p = 3 and
+0.0037 % with the power free, p = 3.12); Theta_AH ~ 8.5 at t = 2, central lapse below
1e-6 by t = 28, r_AH ~ 9 -> ~ 65 and Theta(r = 150) ~ 137 -> ~ 90 between t = 2 and 28.

The Misner-Sharp mass of m_at uses Theta' = spin0_common.ddr(Theta, r) (second order also
at the two end points, which the grid-edge sphere of the --mesh_run printout uses).

Figure: ||Theta C||_2 rescaled by N_r^2 against t on the three grids (single column); the
orders, formerly the left panel, are printed (and tabulated in the paper).

Usage:  fig_conv_collapse.py COARSE MID FINE [--mesh_run mass_hi] [--out conv_collapse.png]
        fig_conv_collapse.py --replot data/conv_collapse.npz [--out conv_collapse.png]
"""
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import spin0_common as sc
import paper_style as ps

FIELDS = list(sc.CONV_FIELDS)


def nrm(v, r, lo, hi):
    """L2 norm over lo <= r <= hi; NaN with 10 points or fewer in the window."""
    return sc.norm_l2(v, r, lo, hi, min_points=11)


def m_at(p, th):
    """Misner-Sharp mass at the grid point nearest the areal radius th -> (m, Theta there)."""
    Th = p["Theta"]
    m = sc.misner_sharp(p, sc.ddr(Th, p["r"]))
    i = np.argmin(np.abs(Th - th))
    return m[i], Th[i]


def plot(D, out):
    """||Theta C||_2 outside the horizon, rescaled by N_r^2 (single column)."""
    fig, ax = plt.subplots(figsize=(ps.COLUMN, 2.5))
    ps.rescaled_norm_curves(ax, D["t"], D["C"], D["N"], D["scal"])
    ax.set_yscale("log")
    ps.headroom(ax, 0.2)
    ps.log_ticks(ax)
    ax.set_xlabel("$t$")
    ax.set_ylabel(r"$N_r^2$-rescaled $\|\Theta C\|_2$")
    ax.legend(loc="upper right")
    fig.tight_layout()
    ps.save(fig, out)


def main():
    if ps.replot(plot, "conv_collapse.png"):
        return
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("rundirs", nargs=3, help="run directories, coarse, medium, fine")
    ap.add_argument("--tmin", type=float, default=4.0,
                    help="first dump time used for the orders and the figure (default 4)")
    ap.add_argument("--rmax_norm", type=float, default=150.0,
                    help="outer edge of the norm window (default 150)")
    ap.add_argument("--ratio_times", default="4,8,12,16,20,24,28",
                    help="dump times of the constraint-ratio table (comma list)")
    ap.add_argument("--t_ah", type=float, default=28.0,
                    help="time at which Theta_AH and M_ADM are quoted (default 28)")
    ap.add_argument("--th_adm", type=float, default=250.0,
                    help="areal radius at which M_ADM is read (default 250; on the fixed grid "
                         "the nearest point, i.e. the grid edge once Theta(r_max) < 250)")
    ap.add_argument("--mesh_run", default=None,
                    help="optional moving-mesh run (README: mass_hi) on which M_ADM is also "
                         "read at Theta = --th_adm, which stays on that grid")
    ap.add_argument("--drift_spheres", default="100,150,200,250",
                    help="spheres of the M_ADM-drift falloff on --mesh_run (default 100,150,200,250)")
    ap.add_argument("--out", default="conv_collapse.png",
                    help="output figure (png; a pdf is written next to it)")
    ps.add_replot_arg(ap)
    a = ap.parse_args()

    dirs = a.rundirs
    T = [sc.list_dumps(d, exact=False) for d in dirs]
    common = sorted(set(T[0][0]) & set(T[1][0]) & set(T[2][0]))
    times = [t for t in common if t >= a.tmin]
    print("times used (all three grids present):", [f"{t:g}" for t in times])

    def load(k, t):
        return sc.load_dump(dirs[k], t)

    cache = {}

    def trio(t):
        if t not in cache:
            D = [load(k, t) for k in range(3)]
            n = [len(d["r"]) for d in D]
            stride = [(n[k] - 1) // (n[0] - 1) for k in range(3)]
            for k in range(3):
                if n[k] - 1 != stride[k] * (n[0] - 1):
                    raise SystemExit(f"grids are not nested: N = {n}")
            lo = max(sc.apparent_horizon(d)[1] for d in D) + 2.0
            cache[t] = (D, stride, lo)
        return cache[t]

    def norms(t, f):
        """norms on the coarse points: f = field -> (||c-m||, ||m-f||); 'C' -> the 3 norms."""
        D, stride, lo = trio(t)
        base = D[0]["r"]
        key = "ThetaC" if f == "C" else f
        v = [D[k][key][::stride[k]] for k in range(3)]
        kk = min(min(map(len, v)), len(base))
        v = [x[:kk] for x in v]
        r = base[:kk]
        if f == "C":
            return tuple(nrm(x, r, lo, a.rmax_norm) for x in v)
        return nrm(v[0] - v[1], r, lo, a.rmax_norm), nrm(v[1] - v[2], r, lo, a.rmax_norm)

    def order(f):
        """Self-convergence order of field f ('C': ||Theta C||) averaged over the times."""
        orders = []
        for t in times:
            if f == "C":
                x, y, z = norms(t, "C")
                if y > 0 and z > 0:
                    orders += [np.log2(x / y), np.log2(y / z)]
            else:
                x, y = norms(t, f)
                if y > 0:
                    orders.append(np.log2(x / y))
        return np.nanmean(orders) if orders else np.nan

    M = np.array([order(f) for f in FIELDS + ["C"]])
    N = [len(load(k, times[0])["r"]) for k in range(3)]
    scal = [s ** 2 for s in trio(times[0])[1]]      # N_r^2 rescaling: 1, 4, 16

    # ---- figure: rescaled constraint norms (the orders are printed, tab:conv-scatter) ----
    ps.setup()
    D = dict(t=np.array(times), C=np.array([norms(t, "C") for t in times]).T,
             N=np.array(N), scal=np.array(scal))
    ps.save_data(a.out, "conv_collapse", **D)
    plot(D, a.out)

    # ---- printed numbers ----
    print(f"\nmeasured collapse convergence orders (outside the horizon, norm window "
          f"[r_AH+2, {a.rmax_norm:g}], averaged over t = {times[0]:g}..{times[-1]:g}):")
    print("  " + "  ".join(f"{f}: {M[j]:.2f}" for j, f in enumerate(FIELDS + ["||ThC||"])))
    print("  [paper: 2.62 2.55 2.34 2.38 2.62 2.42 2.24]")
    print(f"\nconstraint-norm ratios   ||ThC||_{N[0]}/||ThC||_{N[1]}   ||ThC||_{N[1]}/||ThC||_{N[2]}"
          "   (4 = second order)")
    for t in [float(x) for x in a.ratio_times.split(",")]:
        if t not in common:
            continue
        x, y, z = norms(t, "C")
        print(f"  t = {t:4g}   {x / y:.2f}   {y / z:.2f}")
    print("  [paper: 3.69-4.69 coarse pair, 3.95-6.26 fine pair]")
    th = [sc.apparent_horizon(load(k, a.t_ah))[2] for k in range(3)]
    pah = np.log2((th[1] - th[0]) / (th[2] - th[1]))
    print(f"\nTheta_AH at t = {a.t_ah:g}: " + ", ".join(f"{x:.3f}" for x in th)
          + f"   Richardson (2nd order) {th[2] + (th[2] - th[1]) / 3:.3f}, observed order {pah:.2f}"
          f" (differences shrink by {2 ** pah:.1f})   [paper: 12.293 / 12.604 / 12.677, 4.3, 12.70]")
    for t in (2.0, a.t_ah):
        if t in common:
            q = load(2, t)
            _, rah, thah = sc.apparent_horizon(q)
            print(f"finest grid, t = {t:g}: Theta_AH = {thah:.2f}, r_AH = {rah:.1f}, central lapse "
                  f"{q['NN'][0]:.1e}, Theta(r = 150) = {np.interp(150.0, q['r'], q['Theta']):.1f}")
    print("  [paper: Theta_AH ~ 8.5 at t = 2; lapse < 1e-6 by t = 28; r_AH ~ 9 -> ~ 65; "
          "Theta(r = 150) ~ 137 -> ~ 90]")
    for k in range(3):
        m2, th2 = m_at(load(k, 2.0), a.th_adm) if 2.0 in common else (np.nan, np.nan)
        m1, th1 = m_at(load(k, a.t_ah), a.th_adm)
        edge = "  (grid edge)" if abs(th1 - a.th_adm) > 1.0 else ""
        print(f"M_ADM fixed grid N = {N[k]}: t = 2: {m2:.5f} at Theta = {th2:.1f}   "
              f"t = {a.t_ah:g}: {m1:.5f} at Theta = {th1:.1f}{edge}"
              f"   drift {100 * (m1 / m2 - 1):+.2f} %")
    print("  [fixed grid: 3.70836 -> 3.71143 / 3.71164 / 3.71168 at the grid edge Theta ~ 235, +0.08-0.09 %]")
    if a.mesh_run:
        def mload(t):
            return sc.load_dump(a.mesh_run, t)
        m2, th2 = m_at(mload(2.0), a.th_adm)
        m1, th1 = m_at(mload(a.t_ah), a.th_adm)
        print(f"M_ADM moving mesh ({a.mesh_run}): t = 2: {m2:.5f} at Theta = {th2:.1f}   "
              f"t = {a.t_ah:g}: {m1:.5f} at Theta = {th1:.1f}   drift {100 * (m1 / m2 - 1):+.2f} %")
        print("  [paper: 3.70836 -> 3.71107 at Theta = 250, +0.07 %]")
        t_off = [t for t in T[2][0] if load(2, t)["Theta"][-1] < a.th_adm]
        if t_off:
            print(f"Theta = {a.th_adm:g} off the fixed grid from the dump at t = {t_off[0]:g}"
                  "   [paper: t ~ 20]")
        me = [m_at(load(k, a.t_ah), 1e9) for k in range(3)]          # grid-edge spheres
        th_e = me[2][1]
        mm = m_at(mload(a.t_ah), th_e)
        spread = (max(x[0] for x in me) - min(x[0] for x in me)) / me[2][0]
        print(f"t = {a.t_ah:g}, fixed-grid edge Theta = {th_e:.1f}: three fixed grids agree to "
              f"{spread:.1e}, finest fixed grid vs moving mesh (Theta = {mm[1]:.1f}) "
              f"{abs(mm[0] / me[2][0] - 1):.1e}   [paper: 1e-4, 2e-5]")
        S = np.array([float(x) for x in a.drift_spheres.split(",")])
        dm = np.array([m_at(mload(a.t_ah), x)[0] - m_at(mload(2.0), x)[0] for x in S])
        m0 = m_at(mload(2.0), S[-1])[0]
        X = np.vstack([np.ones_like(S), S ** -3]).T
        A, B = np.linalg.lstsq(X, dm, rcond=None)[0]
        print("change of m between t = 2 and " f"{a.t_ah:g} on the moving mesh: "
              + ", ".join(f"Theta = {x:g}: {100 * d / m0:+.3f} %" for x, d in zip(S, dm)))
        print(f"  fit A + B/Theta^3: A = {100 * A / m0:+.4f} % of M_ADM (Theta -> infinity), "
              f"max residual {100 * np.abs(X @ [A, B] - dm).max() / m0:.4f} %"
              "   [paper: Theta^-3 falloff, 1.2 % at Theta = 100, extrapolation about 0.01 %]")
        best = None
        for q in np.arange(2.5, 4.001, 0.01):       # the same fit with the power left free
            Xq = np.vstack([np.ones_like(S), S ** -q]).T
            cq = np.linalg.lstsq(Xq, dm, rcond=None)[0]
            rq = np.abs(Xq @ cq - dm).max()
            if best is None or rq < best[0]:
                best = (rq, q, cq[0])
        print(f"  fit A + B/Theta^p, p free: p = {best[1]:.2f}, A = {100 * best[2] / m0:+.4f} % "
              f"of M_ADM, max residual {100 * best[0] / m0:.4f} %")


if __name__ == "__main__":
    main()
