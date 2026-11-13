#!/usr/bin/env python3
"""Spin-0 pulse extraction in the v0 = 1 baseline, paper fig:wout (wout_background.pdf)
(Sec. "The spin-0 waveform and its speed").

Definitions (paper eq:wout, eq:retarded, eq:wbar):
  w_out = -(b_N/v0) trK + a_r     (= 0.8955 trK + a_r for the baseline couplings)
  on spheres of fixed areal radius Theta = 50, 70, 90, 110, 130, 150, 180;
  background  w_bar(t) = Gaussian running mean of w_out(t) of width sigma = 10 at fixed
  Theta [eq:wbar] (kernel renormalised at the ends; plotted/used only where the window is
  inside the data at the late end);
  pulse  Theta * delta w_out = Theta (w_out - w_bar), plotted against u = t - Theta_*/v0,
  Theta_* = Theta + 2 M_ADM ln(Theta/(2 M_ADM) - 1).

Prints, per sphere, the pulse maximum and minimum in the window u in [--umin, --umax]
(default -30..80) and the retarded times of the extrema; the constancy of the amplitude
across spheres ((max - min)/mean of the extremum over the spheres) for the widths of
--sigma_scan (default sigma = 5, 7.5, 10, 12.5, 15);
and the drift in u of the pulse minimum between the innermost and outermost sphere.
Paper: with sigma = 10 the pulse amplitude is the same on every sphere, and it depends on
sigma but is the same on every sphere for any sigma in [5, 15]; the narrow pulse drifts
slowly in u.  Public-code run (not quoted in the paper): amplitude +0.0606/-0.0718,
constant to a few percent from Theta = 50 to 180 (the constancy of the trough holds for any
sigma in 5-15, that of the peak is more sensitive to the filter and to the innermost
sphere); drift 10.7 / 11.3 in u.

Usage:  fig_wout_background.py RUNDIR [--madm M] [--sigma 10]
        fig_wout_background.py --replot data/wout_background.npz [--out wout_background.png]
"""
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import spin0_common as sc
import paper_style as ps


def pulse_stats(t, W, thetas, M, v0, sigma, umin, umax):
    """Per sphere: (max, u at max, min, u at min) of Theta*delta w_out in the u window,
    as an array (n_spheres, 4)."""
    rows = []
    for k, R in enumerate(thetas):
        bg, full = sc.gaussian_running_mean(t, W[:, k], sigma)
        dw = R * (W[:, k] - bg)
        u = sc.retarded_time(t, R, M, v0)
        s = full & (u >= umin) & (u <= umax) & np.isfinite(dw)
        uu, yy = u[s], dw[s]
        jmax, jmin = int(np.argmax(yy)), int(np.argmin(yy))
        umx, ymx = sc.parabolic_peak(uu, yy, jmax)
        umn, ymn = sc.parabolic_peak(uu, yy, jmin)
        rows.append((ymx, umx, ymn, umn))
    return np.array(rows)


def plot(D, out):
    """(a) Theta w_out and its background at fixed Theta vs t; (b) the pulse vs u."""
    fig, (axl, axr) = plt.subplots(1, 2, figsize=(ps.TEXT, 2.7))
    cols = ps.sphere_colors(D["thetas"])
    for k, R in enumerate(D["thetas"]):
        c = cols[k]
        axl.plot(D["t"], D["ThW"][:, k], color=c, label=rf"$\Theta={R:.0f}$")
        axl.plot(D["t"], D["ThBG"][:, k], "--", color=c, lw=0.8)
        axr.plot(D["u"][:, k], D["ThdW"][:, k], color=c)
    axl.set_xlabel(r"$t$")
    axl.set_ylabel(r"$\Theta w_{\rm out}$")
    axl.set_xlim(0, D["t"].max())
    axr.axhline(0, color="k", lw=0.5)
    ulab = r"$u=t-\Theta_*$" if abs(D["v0"] - 1) < 1e-3 else r"$u=t-\Theta_*/v_0$"
    axr.set_xlabel(ulab)
    axr.set_ylabel(r"$\Theta\,\delta w_{\rm out}$")
    ok = np.isfinite(D["ThdW"])
    axr.set_xlim(np.min(D["u"][ok]), np.max(D["u"][ok]))
    for ax, tag in ((axl, "(a)"), (axr, "(b)")):
        ps.panel_tag(ax, tag, loc="upper left")
    h, lab = axl.get_legend_handles_labels()
    fig.legend(h, lab, loc="lower center", ncol=len(lab), bbox_to_anchor=(0.5, 0.97),
               handlelength=1.4, columnspacing=1.0)
    fig.tight_layout(w_pad=1.5)
    ps.save(fig, out)


def main():
    if ps.replot(plot, "wout_background.png"):
        return
    ap = sc.add_common_args(argparse.ArgumentParser(description=__doc__,
                            formatter_class=argparse.RawDescriptionHelpFormatter),
                            out="wout_background.png")
    ap.add_argument("--thetas", default="50,70,90,110,130,150,180",
                    help="areal radii of the spheres (comma list)")
    ap.add_argument("--sigma", type=float, default=10.0,
                    help="width of the Gaussian running mean (background)")
    ap.add_argument("--sigma_scan", default="5,7.5,10,12.5,15",
                    help="widths for the constancy scan (comma list)")
    ap.add_argument("--umin", type=float, default=-30.0, help="retarded-time window: start")
    ap.add_argument("--umax", type=float, default=80.0, help="retarded-time window: end")
    ps.add_replot_arg(ap)
    a = ap.parse_args()
    cp = sc.couplings(a.c1, a.c2)
    M = sc.resolve_madm(a)
    v0 = cp["v0"]
    thetas = np.array([float(x) for x in a.thetas.split(",")])
    print(f"v0 = {v0:.4f}, -b_N/v0 = {-cp['bN'] / v0:.4f}, M_ADM = {M:.4f}")

    t, S = sc.sphere_series(a.rundir, thetas, ["wout"], cp=cp, nproc=a.nproc)
    W = S["wout"]

    st = pulse_stats(t, W, thetas, M, v0, a.sigma, a.umin, a.umax)
    print(f"\npulse Theta*delta w_out (sigma = {a.sigma:g}), window u in [{a.umin:g},{a.umax:g}]")
    print(f"{'Theta':>6} {'max':>8} {'u(max)':>8} {'min':>8} {'u(min)':>8}")
    for R, (ymx, umx, ymn, umn) in zip(thetas, st):
        print(f"{R:6.0f} {ymx:8.4f} {umx:8.2f} {ymn:8.4f} {umn:8.2f}")

    def sp(x):
        """relative spread (max - min)/|mean|"""
        return (x.max() - x.min()) / abs(x.mean())

    print(f"amplitude: max = {st[:, 0].mean():+.4f} (spread {100 * sp(st[:, 0]):.1f}%), "
          f"min = {st[:, 2].mean():+.4f} (spread {100 * sp(st[:, 2]):.1f}%)   "
          f"[paper: the same on every sphere]")
    # spread excluding the innermost sphere, which is still close to the source
    print(f"   excluding Theta = {thetas[0]:.0f}: spread of max {100 * sp(st[1:, 0]):.1f}%, "
          f"of min {100 * sp(st[1:, 2]):.1f}%")
    print(f"drift of the pulse in u (minimum): {st[-1, 3] - st[0, 3]:+.2f} between "
          f"Theta = {thetas[0]:.0f} and {thetas[-1]:.0f};  (maximum): "
          f"{st[-1, 1] - st[0, 1]:+.2f}   [paper: drifts slowly in u]")

    print("\nsigma scan (amplitude depends on sigma, constancy should not):")
    print(f"{'sigma':>6} {'<max>':>8} {'spread':>7} {'<min>':>8} {'spread':>7}")
    for sg in [float(x) for x in a.sigma_scan.split(",")]:
        s2 = pulse_stats(t, W, thetas, M, v0, sg, a.umin, a.umax)
        print(f"{sg:6.1f} {s2[:, 0].mean():8.4f} {100 * sp(s2[:, 0]):6.1f}% "
              f"{s2[:, 2].mean():8.4f} {100 * sp(s2[:, 2]):6.1f}%")

    # ---------------- figure ----------------
    nt, nk = W.shape
    ThBG, U, DW = (np.full((nt, nk), np.nan) for _ in range(3))
    for k, R in enumerate(thetas):
        bg, full = sc.gaussian_running_mean(t, W[:, k], a.sigma)
        ThBG[full, k] = R * bg[full]
        u = sc.retarded_time(t, R, M, v0)
        s = full & (u >= a.umin) & (u <= a.umax)
        U[:, k] = u
        DW[s, k] = R * (W[s, k] - bg[s])
    out = a.out
    ps.setup()
    D = dict(t=t, thetas=thetas, ThW=thetas[None, :] * W, ThBG=ThBG, u=U, ThdW=DW,
             v0=v0, sigma=a.sigma)
    ps.save_data(out, "wout_background", **D)
    plot(D, out)


if __name__ == "__main__":
    main()
