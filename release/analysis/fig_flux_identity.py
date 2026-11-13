#!/usr/bin/env python3
"""Flux identity of the paper (Sec. "Mass conservation and radiated energy") from the profile
dumps.  NOT the paper's fig:fluxid (flux_identity.pdf): that figure and its 0.02 % come from the
per-step sphere records, flux_identity_rec.py; this script writes flux_identity_dumps.png
from the Delta t = 0.5 dumps, where the burst is under-sampled.  A quick-look diagnostic,
drawn outside the paper's figure style (paper_style.py).

Checks, on spheres of fixed areal radius Theta = 40 ... 200, the exact Misner-Sharp
identity [paper eq:lum]

    L = -dm/dt|_Theta = F_ae + F_psi,
    F_X = (N Theta^2/2) [ A (rho_X + p_X) + ((A^2 + B^2)/B) j_X ],
    A = Theta (Q/2 - K/3),  B = dTheta/dr,

with rho, j, p the aether-frame projections of the stress tensors:
  * aether: spin0_common.aether_set (covariant SET of arXiv:1501.07274 eq. 2.34 reduced to
    the code variables; the time derivatives of N, beta, Q at fixed proper radius r enter
    as W = d_t ln N, B = d_t beta, Qt = d_t Q and are taken across dumps, after
    interpolating the neighbouring dumps onto the r grid of the central one, which matters
    on the moving mesh);
  * scalar: rho = p = (P^2 + psi'^2)/2, j = -P psi'.
-dm/dt|_Theta: m (eq:ms) interpolated to the sphere on every dump, then differentiated in
time across the dumps.  Both time derivatives use the same centred stencil on the exact
dump times (phi.dat): order 8 (9 dumps) by default, and the residual is also printed for
the lower orders of --orders.  With the standard dumps (Delta t = 0.5) the residual at the
burst is set by the time sampling of the burst (it falls with the order of the stencil);
the paper's figure and its 0.02 % residual therefore come from flux_identity_rec.py, which
does the same check on the per-step records of the code's sphere recorder
(--rec_spheres); this script, on Delta t = 0.5 dumps of the dr = 0.023 run, gives 0.6-3 %
at the burst (rms 0.03-0.17 %), a time-sampling error, not a formula error.

Prints, per sphere: the peak of |dm/dt|, the peaks of F_ae and |F_psi|, and the maximum
residual |F_ae + F_psi + dm/dt| as a fraction of the peak flux.  Paper (v0 = 1.267,
dr = 0.023, r_max = 600, t <= 600, per-step records): agreement to within 0.02 % of the
peak flux on every sphere from 40 to 200; |F_psi| <~ 0.02 against aether-flux peaks of 0.9-4
(the scalar-flux and peak numbers are reproduced here; the residual is not, see above).

Usage:  fig_flux_identity.py RUNDIR --c2 0 [--thetas 40,60,...,200] [--plot 40,100,200]
                              [--out flux_identity_dumps.png]
"""
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import spin0_common as sc


def dump_series(rundir, thetas, cp, orders=(4, 8), nproc=4, tmin=-np.inf, tmax=np.inf,
                verbose=True):
    """Flux-identity series from the profile dumps of rundir.

    Returns (times, S, dmdt, res, ok, order): S holds m, F_psi, F_ae<o> for every order o
    and F_ae (= the last order); dmdt and res = F_ae + F_psi + dm/dt use the last order;
    ok marks the dumps with a centred stencil of every order.  Also used by
    flux_identity_rec.py for the comparison with the per-step sphere recorder."""
    thetas = np.asarray(thetas, float)
    orders = list(orders)
    assert all(o in (2, 4, 6, 8, 10) for o in orders)
    order = orders[-1]
    win = max(orders) // 2
    times, paths = sc.list_dumps(rundir, tmin, tmax)
    if verbose:
        print(f"{len(paths)} dumps, t = {times[0]:.2f} .. {times[-1]:.2f};  v0 = {cp['v0']:.4f}")

    def per_dump(i, cache, t):
        p = cache[i]
        r, Th = p["r"], p["Theta"]
        ks = sorted(cache)
        d = {k: sc.ddr(p[k], r) for k in ("Theta", "Acr", "trK", "Q", "Psi")}
        m = sc.misner_sharp(p, d["Theta"])
        out = dict(m=m)
        beta = sc.beta_shift(p)
        for o in orders:
            h = o // 2
            if i - h < 0 or i + h >= len(t):   # first/last dumps: no centred stencil
                out[f"Fae{o}"] = np.full_like(Th, np.nan)
                continue
            st = [k for k in ks if abs(k - i) <= h]
            dt = sc.fixed_r_time_derivs([cache[k] for k in st], t[st], st.index(i))
            W, B, Qt = dt["lnN"], dt["beta"], dt["Q"]
            rho, j, pr = sc.aether_set(p, cp, W, B, Qt, beta=beta, d=d)
            out[f"Fae{o}"] = sc.sphere_flux(p, rho, j, pr, d["Theta"])
        rps, jps, pps = sc.scalar_set(p, d["Psi"])
        out["Fps"] = sc.sphere_flux(p, rps, jps, pps, d["Theta"])
        return sc.at_theta(Th, thetas, out)

    rows = sc.map_dumps(per_dump, times, paths, nproc=nproc, window=win)
    S = {k: np.array([row[k] for row in rows]) for k in rows[0]}
    ok = np.zeros(len(times), bool)
    ok[win:len(times) - win] = True          # centred stencils on both sides, all orders
    if verbose:
        print("\nmax residual / peak flux vs order of the time differences (time sampling):")
        for o in orders:
            dm_o = sc.ddt_series(times, S["m"], order=o)
            r_o = np.abs(S[f"Fae{o}"] + S["Fps"] + dm_o)[ok]
            pk_o = np.nanmax(np.abs(dm_o[ok]), axis=0)
            print(f"   order {o}: " + "  ".join(f"{R:.0f}:{100 * x:.2f}%" for R, x in
                                                zip(thetas, np.nanmax(r_o, axis=0) / pk_o)))
    dmdt = sc.ddt_series(times, S["m"], order=order)
    S["Fae"] = S[f"Fae{order}"]
    res = S["Fae"] + S["Fps"] + dmdt
    return times, S, dmdt, res, ok, order


def main():
    ap = sc.add_common_args(argparse.ArgumentParser(description=__doc__,
                            formatter_class=argparse.RawDescriptionHelpFormatter), madm=False,
                            out="flux_identity_dumps.png")
    ap.add_argument("--thetas", default="40,60,80,100,120,140,160,180,200",
                    help="areal radii of the spheres (comma list)")
    ap.add_argument("--plot", default="40,100,200", help="spheres shown in the figure")
    ap.add_argument("--orders", default="4,8",
                    help="orders of the centred time differences across dumps; the LAST "
                         "one is plotted and quoted")
    ap.add_argument("--tmin", type=float, default=-np.inf, help="first dump time used")
    ap.add_argument("--tmax", type=float, default=np.inf, help="last dump time used")
    a = ap.parse_args()
    cp = sc.couplings(a.c1, a.c2)
    thetas = np.array([float(x) for x in a.thetas.split(",")])
    orders = [int(x) for x in a.orders.split(",")]
    times, S, dmdt, res, ok, order = dump_series(a.rundir, thetas, cp, orders, a.nproc,
                                                 a.tmin, a.tmax)

    print(f"\n{'Theta':>6} {'max|dm/dt|':>11} {'F_ae max':>9} {'F_ae min':>9} "
          f"{'max|F_psi|':>11} {'max|res|':>10} {'res/peak':>9} {'rms res/peak':>12}")
    worst = 0.0
    for k, R in enumerate(thetas):
        pk = np.nanmax(np.abs(dmdt[ok, k]))
        rr = np.abs(res[ok, k])
        frac = np.nanmax(rr) / pk
        worst = max(worst, frac)
        print(f"{R:6.0f} {pk:11.4f} {np.nanmax(S['Fae'][ok, k]):9.4f} "
              f"{np.nanmin(S['Fae'][ok, k]):9.4f} {np.nanmax(np.abs(S['Fps'][ok, k])):11.2e} "
              f"{np.nanmax(rr):10.2e} {100 * frac:8.3f}% {100 * np.sqrt(np.nanmean(rr ** 2)) / pk:11.4f}%")
    print(f"\n[order {order}] F_ae + F_psi and -dm/dt|_Theta agree to within {100 * worst:.2f}% of the peak "
          f"flux on spheres Theta = {thetas[0]:.0f}..{thetas[-1]:.0f}  [paper, from the per-step records: 0.02%]")
    print(f"max |F_psi| over all spheres = {np.nanmax(np.abs(S['Fps'][ok])):.3g}  "
          f"[paper: <~0.02];  aether-flux peaks max|F_ae| = "
          f"{np.nanmin(np.nanmax(np.abs(S['Fae'][ok]), axis=0)):.2f} .. "
          f"{np.nanmax(np.nanmax(np.abs(S['Fae'][ok]), axis=0)):.2f}  [paper: 0.9-4]")

    # ---------------- figure ----------------
    show = [float(x) for x in a.plot.split(",")]
    fig, axs = plt.subplots(2, len(show), figsize=(4.6 * len(show), 6.2), sharex="col",
                            gridspec_kw=dict(height_ratios=[2.2, 1]), squeeze=False)
    for c, R in enumerate(show):
        k = int(np.argmin(np.abs(thetas - R)))
        ax, axr = axs[0, c], axs[1, c]
        tt = times[ok]
        ax.plot(tt, (S["Fae"] + S["Fps"])[ok, k], color="0.75", lw=4,
                label=r"$F_{\rm æ}+F_\psi$")
        ax.plot(tt, -dmdt[ok, k], "k--", lw=1.2, label=r"$-dm/dt|_\Theta$")
        ax.plot(tt, S["Fae"][ok, k], color="C0", lw=1, label=r"$F_{\rm æ}$ (aether)")
        ax.plot(tt, S["Fps"][ok, k], color="C1", lw=1, label=r"$F_\psi$ (scalar)")
        ax.axhline(0, color="k", lw=0.6)
        ax.set_title(rf"$\Theta={thetas[k]:.0f}$")
        ax.grid(alpha=0.3)
        pk = np.nanmax(np.abs(dmdt[ok, k]))
        axr.plot(tt, res[ok, k], color="C3", lw=1)
        axr.axhline(0, color="k", lw=0.6)
        axr.set_ylim(-0.01 * pk, 0.01 * pk)
        axr.text(0.03, 0.85, r"$\pm1\%$ of max $|dm/dt|$", transform=axr.transAxes,
                 fontsize=8)
        axr.set_xlabel(r"$t$")
        axr.grid(alpha=0.3)
        if c == 0:
            ax.set_ylabel("energy flux through the sphere")
            axr.set_ylabel(r"$F_{\rm æ}+F_\psi+dm/dt$")
            ax.legend(fontsize=8)
    fig.tight_layout()
    out = a.out
    fig.savefig(out, dpi=150)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
