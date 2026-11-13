#!/usr/bin/env python3
"""The transient white hole: paper fig:wh (white_hole.pdf) (Sec. "Large c1: transient white holes").

Runs (README §5, "Inner shell: the seven-grid ladder at t = 30"): c1 = 0.8,
c2 = -0.7111111, ingoing shell, A0 = 0.10, r_max = 40; three uniform grids
N = 12800 / 25600 / 51200 (dx = 0.0031 / 0.0016 / 0.0008) and four grids refined around the
shell (zone B pre-positioned at x0 = 9.9, F = 8 on the same three N, and F = 16 on
N = 51200), all with a dump at t = 30; the finest uniform run (dumps every 2 to t = 54)
also gives panel (a).

Panel (a), from --panel_a (default: the last of --uniform), at the dumps requested at
--times (the panel-(a) run's --times list of the README recipe; requested times without a
file are skipped).  The dumps are opened by their requested time and placed at their
exact times (phi.dat): the file names keep one decimal, so e.g. the dump requested at
29.25 is prof_t029.2.csv.  At every dump: the areal radii of the zeros of
theta_- = Theta' - (3Q - 2K) Theta/6 at r > 1 (edges of the anti-trapped region), the areal
radius of the inner shell (peak of |P| at r > 3, Theta < 35) from t >= 30, of the escaping
pulse (largest local maximum of |P| beyond the shell) from t >= --t_pulse, of the inner
front (steepest |dP/dTheta| at 1.5 < Theta < 5) up to t = --t_front (these three time cuts
compare the REQUESTED times), a dotted line at T_STEEPEN = 29.5 (the inner shell steepens,
caption of fig:wh), a grey line at --t_slice (the time of fig:sharp(a)), and the outgoing
light ray dr/dt = N - beta started at r = --r_ray (9.966, the shell at t = 30) and
integrated across the dumps (velocity interpolated linearly in time, NSUB = 400 sub-steps
per dump interval).  beta is rebuilt from beta' = N (Q + K/3).  The left edge of panel (a)
is T_LEFT = 19.
Panels (b, c): FWHM and height of the |P| peak of the shell at t = --t (30) against the
grid spacing r[k+1] - r[k] at the peak, for the uniform (blue) and refined (red) groups,
each with a power-law fit (family_fits).  Paper: the width shrinks roughly as dr^0.8 and the height
grows as dr^-0.5 (per-family fits).

Printed: spacing, FWHM and |P| of the shell per grid and the power laws (public-code runs:
width dr^0.84 uniform / dr^0.79 refined, |P| dr^-0.50 / dr^-0.47); the distance of the
light ray behind the escaping pulse at the times of --ray_report (paper: 3.2 at t = 30,
0.3 at t = 34, 0.1 at t = 40); and the inner-front track of panel (a): areal radius,
steepest slope, jump of P across it, lapse there and its coordinate speed against the
ingoing light speed -beta - N (paper: a jump of about 0.3 in P, along ingoing light rays,
N ~ 0.13 where it travels, nearly stalling at Theta ~ 3.2, falling to the centre at
t ~ 44).  The anti-trapped history (paper: from t ~ 20, split at t ~ 36, inner part gone
by t ~ 44) is printed by family_transitions.py --extent; the layer on the escaping pulse
by wh_front_convergence.py.

Usage:  fig_white_hole.py --uniform D1,D2,D3 --refined D4,D5,D6,D7 [--panel_a D3]
                          [--times T1,T2,...] [--out white_hole.png]
        fig_white_hole.py --replot data/white_hole.npz [--out white_hole.png]
"""
import argparse
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import spin0_common as sc
import paper_style as ps

T_STEEPEN = 29.5    # the inner shell steepens (caption of fig:wh): dotted line in panel (a)
T_LEFT = 19.0       # left edge of panel (a)
NSUB = 400          # sub-steps of the light ray per dump interval
# requested dump times of the panel-(a) run (README recipe: T(20,54,2) and the quarter-time
# dumps around the formation of the inner shell)
TIMES_A = ",".join([f"{t:g}" for t in range(20, 55, 2)]
                   + ["29", "29.25", "29.5", "29.75", "30.25", "30.5", "30.75", "31"])


def shell_point(rundir, t):
    """Spacing r[k+1] - r[k], FWHM and height of the |P| peak of the inner shell at t."""
    p = sc.load_dump(rundir, t)
    r, P = p["r"], p["P"]
    k = sc.peak_index(p, "P", rmin=3.0, thmax=35.0)
    return r[k + 1] - r[k], sc.fwhm(r, P, k)[0], abs(P[k])


def family_fits(dx, y, refined):
    """Power-law fits y ~ dx^s per family -> [(mask, s, intercept)] for (uniform, refined),
    s = nan for a family with fewer than two grids."""
    out = []
    for m in (~refined, refined):
        if m.sum() >= 2:
            s = np.polyfit(np.log(dx[m]), np.log(y[m]), 1)
            out.append((m, s[0], s[1]))
        else:
            out.append((m, np.nan, np.nan))
    return out


def plot(D, out):
    """(a) anti-trapped edges, inner shell and light ray; (b, c) width and height of the
    |P| peak against the spacing at the shell, per family, with power-law fits."""
    fig = plt.figure(figsize=(ps.TEXT, 2.55))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.35, 1, 1])
    ax = fig.add_subplot(gs[0])
    ax.plot(D["zt"], D["zth"], "o", color=ps.PURPLE, ms=2.2, mew=0,
            label=r"$\theta_-=0$")
    ax.plot(D["shell_t"], D["shell_th"], "k-", lw=1.3, label="inner shell")
    if "pulse_t" in D:
        ax.plot(D["pulse_t"], D["pulse_th"], "-", color=ps.GREEN, lw=1.0, label="escaping pulse")
        ax.plot(D["front_t"], D["front_th"], "-", color=ps.YELLOW, lw=1.4, label="inner front")
    ax.plot(D["ray_t"], D["ray_th"], "--", color=ps.BLUE, lw=1.3, label="light ray from the shell")
    ax.axvline(T_STEEPEN, color="0.5", lw=.7, ls=":")
    if "t_slice" in D:
        ax.axvline(float(D["t_slice"]), color="0.75", lw=.7, ls="-", zorder=0)
    ax.set_xlabel("$t$")
    ax.set_ylabel(r"areal radius $\Theta$")
    ax.set_ylim(0, 40)
    ax.set_xlim(T_LEFT, D["tmax"])
    ax.legend(loc="upper left", bbox_to_anchor=(0.0, 0.92), markerscale=1.8,
              fontsize=ps.FS_SMALL, handlelength=1.5, labelspacing=0.25)
    ps.panel_tag(ax, "(a)")
    dx, ref = 1e3 * D["dx"], D["refined"].astype(bool)     # spacing in units of 1e-3
    for k, (key, yl, tag) in enumerate((("width", "width of the $|P|$ peak", "(b)"),
                                        ("height", r"max $|P|$", "(c)"))):
        ax = fig.add_subplot(gs[k + 1])
        y = D[key]
        fams = ((ps.CAT_A, "o", "uniform"), (ps.CAT_B, "s", "refined"))
        fits = family_fits(dx, y, ref)
        gm = [np.exp(np.mean(np.log(y[m]))) if m.any() else np.nan for m, _, _ in fits]
        notes = []
        for f_, ((m, s0, s1), (c, mk, lab2)) in enumerate(zip(fits, fams)):
            if m.sum() < 2:
                continue
            ax.loglog(dx[m], y[m], mk, color=c, ms=4.5, mec="white", mew=0.4, zorder=3,
                      label=lab2)
            xx = np.logspace(np.log10(dx[m].min()), np.log10(dx[m].max()), 10)
            ax.loglog(xx, np.exp(s1) * xx ** s0, color=c, lw=0.9, ls="--")
            notes.append((np.sqrt(dx[m].min() * dx[m].max()), (s0, s1), c, gm[f_] > gm[1 - f_],
                          lab2))
        ax.set_xlabel(r"$10^3\,dr$ at the shell")
        ax.set_ylabel(yl)
        ps.log_ticks(ax, "x", subs=(1, 2, 3, 5), plain=True)
        ps.log_ticks(ax, "y", subs=(1, 2, 3, 5), plain=True)
        ps.panel_tag(ax, tag, loc="upper left" if key == "width" else "upper right")
        fig.canvas.draw()
        # exponent beside the fit, on the side away from the other family
        for xm, s, c, above, lab2 in notes:
            p0 = ax.transData.transform((xm, np.exp(s[1]) * xm ** s[0]))
            p1 = ax.transData.transform((1.5 * xm, np.exp(s[1]) * (1.5 * xm) ** s[0]))
            d = (p1 - p0) / np.hypot(*(p1 - p0))
            n = np.array([-d[1], d[0]]) * (1 if d[0] > 0 else -1)   # normal pointing up
            off = 22 * (n if above else -n)
            ax.annotate(rf"{lab2}" + "\n" + rf"$\propto dr^{{{s[0]:.2f}}}$", (xm, np.exp(s[1]) * xm ** s[0]),
                        xytext=tuple(off), textcoords="offset points", color=c,
                        ha="center", va="center", fontsize=ps.FS_LEGEND, linespacing=1.0)
    fig.tight_layout(w_pad=0.6)
    ps.save(fig, out)


def track_panel_a(rundir, times, t_shell, t_pulse, t_front):
    """Features of panel (a) at the dumps requested at `times` (see the module docstring).

    Returns a dict: rows [(t, [Theta of the theta_- zeros], Theta_shell or nan)], pulse
    [(t, Theta)], front [(t, Theta, r, N, beta, jump, slope)], Z [(t, r, N - beta, Theta)]
    and keys [(t_req, t)], with t the exact time of each dump."""
    t_req, t_all, paths = sc.requested_dumps(rundir, times)
    rows, pulse, front, Z = [], [], [], []
    for tq, t, f in zip(t_req, t_all, paths):
        p = sc.load_profile(f)
        r, Th, N = p["r"], p["Theta"], p["NN"]
        _, thm = sc.theta_pm(p, sc.ddr(Th, r))
        beta = sc.beta_shift(p)
        s = np.sign(thm)
        zm = [Th[i] for i in np.where(s[:-1] * s[1:] < 0)[0] if r[i] > 1]
        ks = sc.peak_index(p, "P", rmin=3.0, thmax=35.0)      # the inner shell
        rows.append((t, zm, Th[ks] if tq >= t_shell - 1e-9 else np.nan))
        Pa = np.abs(p["P"])
        if tq >= t_pulse - 1e-9:   # escaping pulse: largest local max of |P| beyond the shell
            lm = np.where((Pa[1:-1] > Pa[:-2]) & (Pa[1:-1] >= Pa[2:]))[0] + 1
            lm = lm[(Th[lm] > Th[ks] + 0.5) & (Th[lm] < Th[-1] - 1.0)]
            if len(lm):
                pulse.append((t, Th[lm[np.argmax(Pa[lm])]]))
        if tq <= t_front + 1e-9:   # inner front: steepest |dP/dTheta| at 1.5 < Theta < 5
            jf, slope = sc.steepest_slope(p, 1.5, 5.0)
            jump = np.interp(Th[jf] + 0.05, Th, p["P"]) - np.interp(Th[jf] - 0.05, Th, p["P"])
            front.append((t, Th[jf], r[jf], N[jf], beta[jf], jump, slope))
        Z.append((t, r, N - beta, Th))
    return dict(rows=rows, pulse=pulse, front=front, Z=Z, keys=list(zip(t_req, t_all)))


def light_ray(Z, i0, r0, nsub=NSUB):
    """Outgoing light ray dr/dt = N - beta from r0 at the dump i0 of Z [(t, r, N - beta,
    Theta)], the velocity interpolated linearly in time between dumps -> [(t, Theta)]."""
    rr = r0
    ray = [(Z[i0][0], np.interp(rr, Z[i0][1], Z[i0][3]))]
    for k in range(i0, len(Z) - 1):
        h = (Z[k + 1][0] - Z[k][0]) / nsub
        for s_ in range(nsub):
            al = s_ / nsub
            v = (1 - al) * np.interp(rr, Z[k][1], Z[k][2]) + al * np.interp(rr, Z[k + 1][1], Z[k + 1][2])
            rr += h * v
        ray.append((Z[k + 1][0], np.interp(rr, Z[k + 1][1], Z[k + 1][3])))
    return ray


def main():
    if ps.replot(plot, "white_hole.png"):
        return
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--uniform", required=True,
                    help="comma list of the uniform-grid run directories (coarse to fine)")
    ap.add_argument("--refined", required=True,
                    help="comma list of the refined (zone B) run directories")
    ap.add_argument("--panel_a", default=None,
                    help="run directory for panel (a) (default: the last uniform one)")
    ap.add_argument("--times", default=TIMES_A,
                    help="requested dump times of the panel-(a) run (its --times; default: "
                         "the README recipe, 20,22,...,54 and 29,29.25,...,31)")
    ap.add_argument("--t", type=float, default=30.0,
                    help="dump time of panels (b, c) and start of the light ray (default 30)")
    ap.add_argument("--r_ray", type=float, default=9.966,
                    help="starting radius r of the outgoing light ray at --t (default 9.966)")
    ap.add_argument("--tmax", type=float, default=55.0, help="right edge of panel (a)")
    ap.add_argument("--t_pulse", type=float, default=29.25,
                    help="first dump at which the escaping pulse is traced (default 29.25)")
    ap.add_argument("--t_front", type=float, default=42.0,
                    help="last dump at which the inner front is traced (default 42)")
    ap.add_argument("--ray_report", default="30,34,40",
                    help="times at which the distance of the light ray behind the escaping "
                         "pulse is printed (default 30,34,40)")
    ap.add_argument("--t_slice", type=float, default=40.0,
                    help="time of the profile of fig:sharp(a), marked in panel (a)")
    ap.add_argument("--out", default="white_hole.png",
                    help="output figure (png; a pdf is written next to it)")
    ps.add_replot_arg(ap)
    a = ap.parse_args()
    uni = a.uniform.split(",")
    ref = a.refined.split(",")
    run_a = a.panel_a or uni[-1]

    # ---- panels (b, c): the inner shell at t on every grid ----
    pts = [shell_point(d, a.t) for d in uni + ref]
    refined = np.array([False] * len(uni) + [True] * len(ref))
    dx = np.array([q[0] for q in pts])
    width, height = np.array([q[1] for q in pts]), np.array([q[2] for q in pts])

    # ---- panel (a): anti-trapped edges, shell, pulse, front, light ray ----
    A = track_panel_a(run_a, [float(x) for x in a.times.split(",")], a.t, a.t_pulse,
                      a.t_front)
    rows, pulse, front, Z = A["rows"], A["pulse"], A["front"], A["Z"]
    t_of = dict(A["keys"])                     # requested time -> exact time
    if a.t not in t_of:
        raise SystemExit(f"{run_a}: no dump at t = {a.t:g}")
    i0 = [k for k, (tq, _) in enumerate(A["keys"]) if tq == a.t][0]
    ray = light_ray(Z, i0, a.r_ray)

    zt = np.concatenate([[q[0]] * len(q[1]) for q in rows]) if rows else np.array([])
    zth = np.concatenate([q[1] for q in rows]) if rows else np.array([])
    shell = np.array([(q[0], q[2]) for q in rows if np.isfinite(q[2])])
    ps.setup()
    D = dict(zt=zt, zth=zth, shell_t=shell[:, 0], shell_th=shell[:, 1],
             ray_t=np.array([x[0] for x in ray]), ray_th=np.array([x[1] for x in ray]),
             dx=dx, width=width, height=height, refined=refined, tmax=a.tmax,
             pulse_t=np.array([q[0] for q in pulse]), pulse_th=np.array([q[1] for q in pulse]),
             front_t=np.array([q[0] for q in front]), front_th=np.array([q[1] for q in front]),
             t_slice=a.t_slice)
    ps.save_data(a.out, "white_hole", **D)
    plot(D, a.out)

    print(f"t = {a.t:g}: spacing at the shell, FWHM, |P| (uniform / refined):")
    for d, q, rf in zip(uni + ref, pts, refined):
        print(f"  {os.path.basename(d.rstrip('/')):>22} {q[0]:.5f}  {q[1]:.5f}  {q[2]:.3f}  "
              f"{'refined' if rf else 'uniform'}")
    s_all = {j: np.polyfit(np.log(dx), np.log(y), 1)[0] for j, y in ((1, width), (2, height))}
    print(f"power laws, all grids: width ~ dr^{s_all[1]:+.2f}, |P| ~ dr^{s_all[2]:+.2f}   "
          "[paper, per family: ~ dr^0.8, dr^-0.5]")
    for yl, y in (("width (FWHM)", width), (r"max $|P|$", height)):
        for (m, s, _), fam in zip(family_fits(dx, y, refined),
                                  ("uniform grid", "refined near the shell")):
            if m.sum() >= 2:
                print(f"   {yl:12s} {fam:24s} dr^{s:+.2f}")
    print(f"light ray from r = {a.r_ray} at t = {a.t:g}: Theta = {ray[0][1]:.2f} -> {ray[-1][1]:.2f} "
          f"at t = {ray[-1][0]:g}")
    print("light ray behind the escaping pulse (Theta_pulse - Theta_ray)   "
          "[paper: 3.2 at t = 30, 0.3 at t = 34, 0.1 at t = 40]")
    rt, pt = dict((x[0], x[1]) for x in ray), dict(pulse)
    for tr in [float(x) for x in a.ray_report.split(",")]:
        te = t_of.get(tr)                      # exact time of the dump requested at tr
        if te is not None and te in rt and te in pt:
            print(f"  t = {tr:g}: pulse {pt[te]:.2f}, ray {rt[te]:.2f}, behind by {pt[te] - rt[te]:.2f}")
    print("inner front (steepest |dP/dTheta| at 1.5 < Theta < 5): Theta, steepest slope, "
          "jump P(Theta+0.05) - P(Theta-0.05), lapse N, dr/dt between dumps vs -beta - N   "
          "[paper: jump ~0.3, N ~ 0.13, along rdot = -beta - N, stalls at Theta ~ 3.2, "
          "falls to the centre at t ~ 44]")
    for k, (tf, thf, rf, nf, bf, jp, gm) in enumerate(front):
        v = (rf - front[k - 1][2]) / (tf - front[k - 1][0]) if k else np.nan
        print(f"  t = {tf:5.2f}: Theta = {thf:.3f}  slope {gm:6.1f}  jump {jp:+.3f}  N = {nf:.3f}"
              f"  dr/dt = {v:+.3f}  -beta-N = {-bf - nf:+.3f}")


if __name__ == "__main__":
    main()
