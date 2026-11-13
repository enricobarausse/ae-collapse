#!/usr/bin/env python3
"""Misner-Sharp luminosity and the burst as a clock, paper fig:lum (luminosity_clock.pdf)
(Sec. "The spin-0 waveform and its speed"; v0 = 1.267, c2 = 0).

Definitions:
  m = Misner-Sharp mass [eq:ms] on spheres of fixed areal radius (default Theta = 120,
  160, 200, 240, 300, 350, 400);  L = -dm/dt|_Theta from centred differences across the
  dumps (5-point stencil on the exact dump times);  u = t - Theta_*/v0 [eq:retarded].
  Burst arrival time t_arr(Theta): the time of the deepest minimum of L (the sharp negative
  spike of the burst; --arrival max uses the sharp positive spike instead) within
  BURST_HALF_WINDOW (15) of the largest |L| on the sphere, refined by a parabola through
  the three samples around it.
  Clock: lag = t_arr(Theta) - t_arr(Theta_0) against the tortoise distance
  Delta Theta_* = Theta_* - Theta_*0; a straight-line fit lag = c + Delta Theta_*/v gives the
  speed v of the burst.
(a) L against u.  (b) lag vs Delta Theta_* (points coloured by sphere) with the fit (solid),
v0 (dashed) and the scalar speed 1 (dotted).

Prints the arrival times of the two spikes and their retarded times u, and the fitted
speed (paper: the burst at u ~ 10-15 on every sphere; v = 1.241, 2 % below v0 = 1.267,
attributed to the Schwarzschild tortoise), and the slope 1 that the scalar speed 1 would
give (not in the paper).  Public-code run (sl600_N25601): deepest minimum at u = 9.7-14.7
on Theta = 120-400 (the positive spike just before it at u = 7.2-12.0), v = 1.2406.

Usage:  fig_luminosity_clock.py RUNDIR --c2 0 [--madm M]
        fig_luminosity_clock.py --replot data/luminosity_clock.npz [--out luminosity_clock.png]
"""
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import spin0_common as sc
import paper_style as ps

BURST_HALF_WINDOW = 15.0   # the burst spikes lie within +-15 in t of the largest |L|


def plot(D, out):
    """(a) L = -dm/dt on fixed-Theta spheres vs u; (b) burst arrival lag vs Delta Theta_*."""
    fig, (axl, axr) = plt.subplots(1, 2, figsize=(ps.TEXT, 2.7))
    cols = ps.sphere_colors(D["thetas"])
    ccl = dict(zip(D["thetas"], cols))
    for k, R in enumerate(D["thetas"]):
        axl.plot(D["u"][:, k], D["L"][:, k], color=cols[k], lw=1.0,
                 label=rf"$\Theta={R:.0f}$")
    axl.axhline(0, color="k", lw=0.5)
    axl.set_xlim(D["umin"], D["umax"])
    axl.set_xlabel(r"$u=t-\Theta_*/v_0$")
    axl.set_ylabel(r"$L=-dm/dt|_\Theta$")
    xx = np.linspace(0, D["dT"].max() * 1.05, 50)
    axr.plot(xx, xx, ":", color="0.45", lw=1.2, label=r"$v=1$ (scalar)")
    axr.plot(xx, xx / D["v0"], "--", color="k", lw=1.0, label=rf"$v_0={D['v0']:.3f}$")
    axr.plot(xx, D["icpt"] + D["slope"] * xx, "-", color=ps.BLUE, lw=1.2,
             label=rf"fit, $v={D['v']:.3f}$", zorder=2)
    for x, y, R in zip(D["dT"], D["lag"], D["th_clock"]):
        axr.plot(x, y, "o", color=ccl[R], ms=5, mec="k", mew=0.4, zorder=3)
    axr.set_xlabel(r"$\Delta\Theta_*$")
    axr.set_ylabel("arrival lag of the burst")
    axr.set_xlim(left=-10)
    axr.set_ylim(bottom=-10)
    axr.legend(loc="upper left", bbox_to_anchor=(0.0, 0.88))
    for ax, tag in ((axl, "(a)"), (axr, "(b)")):
        ps.panel_tag(ax, tag, loc="upper left")
    h, lab = axl.get_legend_handles_labels()
    fig.legend(h, lab, loc="lower center", ncol=len(lab), bbox_to_anchor=(0.5, 0.97),
               handlelength=1.4, columnspacing=1.0)
    fig.tight_layout(w_pad=1.5)
    ps.save(fig, out)


def main():
    if ps.replot(plot, "luminosity_clock.png"):
        return
    ap = sc.add_common_args(argparse.ArgumentParser(description=__doc__,
                            formatter_class=argparse.RawDescriptionHelpFormatter),
                            out="luminosity_clock.png")
    ap.set_defaults(c2=0.0)
    ap.add_argument("--thetas", default="120,160,200,240,300,350,400",
                    help="areal radii of the spheres (comma list)")
    ap.add_argument("--arrival", choices=("min", "max"), default="min",
                    help="burst arrival = deepest minimum (default) or highest maximum of L")
    ap.add_argument("--umin", type=float, default=-40.0, help="retarded-time window: start")
    ap.add_argument("--umax", type=float, default=220.0, help="retarded-time window: end")
    ps.add_replot_arg(ap)
    a = ap.parse_args()
    cp = sc.couplings(a.c1, a.c2)
    M = sc.resolve_madm(a)
    v0 = cp["v0"]
    thetas = np.array([float(x) for x in a.thetas.split(",")])
    print(f"v0 = {v0:.4f}, M_ADM = {M:.4f}")

    t, S = sc.sphere_series(a.rundir, thetas, ["m"], cp=cp, nproc=a.nproc)
    L = -sc.ddt_series(t, S["m"], order=4)
    Ts = sc.tortoise(thetas, M)

    U, LL = (np.full((len(t), len(thetas)), np.nan) for _ in range(2))
    tarr = np.full(len(thetas), np.nan)
    tmin_arr = np.full(len(thetas), np.nan)
    tmax_arr = np.full(len(thetas), np.nan)
    print(f"\n{'Theta':>6} {'Theta_*':>8} {'t_arr(min L)':>12} {'u':>6} {'L_min':>7} "
          f"{'t_arr(max L)':>12} {'u':>6} {'L_max':>7}")
    for k, R in enumerate(thetas):
        y = L[:, k]
        u = sc.retarded_time(t, R, M, v0)
        ok = np.isfinite(y)
        if not ok.any():
            print(f"{R:6.0f}  (sphere outside the grid)")
            continue
        # burst = the sharp spike: search near the largest |L| on this sphere
        jb = int(np.nanargmax(np.abs(y)))
        win = ok & (np.abs(t - t[jb]) < BURST_HALF_WINDOW)
        idx = np.nonzero(win)[0]
        jn = idx[np.argmin(y[idx])]
        jx = idx[np.argmax(y[idx])]
        tn, yn = sc.parabolic_peak(t, y, jn)
        tx, yx = sc.parabolic_peak(t, y, jx)
        tmin_arr[k], tmax_arr[k] = tn, tx
        tarr[k] = tn if a.arrival == "min" else tx
        print(f"{R:6.0f} {Ts[k]:8.2f} {tn:12.3f} {tn - Ts[k] / v0:6.2f} {yn:7.3f} "
              f"{tx:12.3f} {tx - Ts[k] / v0:6.2f} {yx:7.3f}")
        s = ok & (u >= a.umin) & (u <= a.umax)
        U[:, k], LL[s, k] = u, y[s]

    g = np.isfinite(tarr)
    ub = tmin_arr[g] - Ts[g] / v0
    print(f"burst (deepest minimum of L) at u = {ub.min():.1f} .. {ub.max():.1f}   "
          "[paper: u ~ 10-15]")
    dT = Ts[g] - Ts[g][0]
    lag = tarr[g] - tarr[g][0]
    slope, icpt = np.polyfit(dT, lag, 1)
    v = 1.0 / slope
    resid = lag - (icpt + slope * dT)
    print(f"\nburst speed from the lag fit ({a.arrival} of L): v = {v:.4f}  "
          f"({100 * (v / v0 - 1):+.2f}% vs v0 = {v0:.4f}); rms residual {resid.std():.3f}"
          f"   [paper: 1.241]")
    # the other spike as arrival marker: a check of the definition
    alt = tmax_arr if a.arrival == "min" else tmin_arr
    sl2 = np.polyfit(dT, alt[g] - alt[g][0], 1)[0]
    print(f"   (other spike as marker: v = {1.0 / sl2:.4f})")
    print("scalar speed 1 would give lag = Delta Theta_*: slope 1 vs measured "
          f"{slope:.4f}   [not in the paper]")

    out = a.out
    ps.setup()
    D = dict(thetas=thetas, u=U, L=LL, umin=a.umin, umax=a.umax, th_clock=thetas[g],
             dT=dT, lag=lag, slope=slope, icpt=icpt, v=v, v0=v0)
    ps.save_data(out, "luminosity_clock", **D)
    plot(D, out)


if __name__ == "__main__":
    main()
