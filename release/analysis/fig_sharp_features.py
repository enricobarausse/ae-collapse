#!/usr/bin/env python3
"""The three sharp features of the white-hole runs: paper fig:sharp (sharp_features.pdf) (Sec.
"Large c1: transient white holes").

Runs (README §5): c1 = 0.8, c2 = -0.7111111, ingoing shell, A0 = 0.10.
  --uniform   the three uniform grids of the inner-shell ladder (r_max = 40,
              N = 12800 / 25600 / 51200, dx = 0.0031 / 0.0016 / 0.0008)
  --refined   the inner-shell grids with zone B pre-positioned at x0 = 9.9 (F = 8 on
              N = 51200, F = 16 on N = 51200)
  --escaping  the run with zone A locked on the escaping pulse (F = 8, N = 51200, r_max = 150)

Panels, all showing P:
  (a) |P| against Theta at t = --t_over (40) on the finest uniform grid, with the three
      features marked;
  (b) the escaping pulse on the zone-A run at --t_esc (40, 50, 64, 72, 80), against
      Theta - Theta_peak (peak of |P| at Theta > 15);
  (c) the inner shell at t = 30 on the uniform and refined grids, against Theta - Theta_peak
      (peak of |P| at 5 < Theta < 12);
  (d) the inner front at t = --t_front (40) on the three uniform grids, against Theta.
Printed: the steepest |dP/dTheta| of the inner front (1.5 < Theta < 5) on the uniform
grids at t = --t_steep (24) and --t_front (40) (paper: little change at t = 24, doubling
per halving of dr at t = 40) and the peak of the escaping pulse at each --t_esc.
Legends: the uniform runs of panels (c, d) are labelled by the grid spacing at the
feature (read from the data), the refined runs by their refinement factors --refined_F
(8,16: the code's --refB_F of each run, which the profiles do not record); line styles go
from dotted (coarsest) to solid (finest) in each family.

Usage:  fig_sharp_features.py --uniform D1,D2,D3 --refined D4,D5 --escaping D6
                              [--out sharp_features.png]
        fig_sharp_features.py --replot data/sharp_features.npz [--out sharp_features.png]
"""
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import spin0_common as sc
import paper_style as ps

LINESTYLES = (":", "--", "-")     # coarse -> fine within a family (the finest is solid)


def peak(p, lo, hi):
    """Index of the |P| peak at lo < Theta < hi."""
    return sc.peak_index(p, "P", rmin=-np.inf, thmin=lo, thmax=hi)


def steepest(d, t, lo=1.5, hi=5.0):
    """Steepest |dP/dTheta| of the inner front of run d at t -> (slope, Theta there)."""
    p = sc.load_dump(d, t)
    j, g = sc.steepest_slope(p, lo, hi)
    return g, p["Theta"][j]


def family_styles(refined):
    """Colour and line style of each run of panel (c): uniform CAT_A, refined CAT_B, line
    styles from LINESTYLES ending with the solid one in each family."""
    out = []
    for fam, c in ((False, ps.CAT_A), (True, ps.CAT_B)):
        n = int(np.sum(refined == fam))
        out += [(c, LINESTYLES[len(LINESTYLES) - n + k]) for k in range(n)]
    return out


def plot(D, out):
    """(a) |P| against Theta with the three features; (b) the escaping pulse; (c) the inner
    shell on every grid; (d) the inner front on the uniform grids."""
    fig, ax = plt.subplots(1, 4, figsize=(ps.TEXT, 2.6),
                           gridspec_kw=dict(width_ratios=[1.3, 1, 1, 1]))
    a = ax[0]
    a.semilogy(D["o_th"], np.abs(D["o_P"]), "k-", lw=0.7)
    a.set_xlim(0, 30)
    a.set_ylim(1e-5, 1e3)
    a.set_xlabel(r"$\Theta$")
    a.set_ylabel(r"$|P|$")
    # labels at |P| = 8, 40, 1 (text) with arrows down to |P| = 0.012, 0.9, 0.012
    for x, s, y in zip(D["o_mark_th"], ("inner\nfront", "inner\nshell", "escaping\npulse"),
                       (8, 40, 1.0)):
        a.annotate(s, (x, 0.9 if y > 2 else 0.012), (x, y), ha="center", va="bottom",
                   fontsize=ps.FS_SMALL, arrowprops=dict(arrowstyle="->", lw=0.6))
    ps.panel_tag(a, "(a)", loc="lower left")
    a.text(0.97, 0.04, rf"$t={D['t_over']:g}$", transform=a.transAxes, ha="right",
           fontsize=ps.FS_SMALL)

    a = ax[1]
    cols = plt.cm.viridis(np.linspace(0, 0.85, len(D["e_t"])))
    for k, (t, c) in enumerate(zip(D["e_t"], cols)):
        n = int(D["e_n"][k])
        a.plot(D["e_x"][k, :n], D["e_P"][k, :n], color=c, lw=0.8, label=rf"$t={t:g}$")
    a.set_xlim(-0.3, 0.3)
    a.set_xlabel(r"$\Theta-\Theta_{\rm peak}$")
    a.set_ylabel(r"$P$")
    a.legend(fontsize=ps.FS_TINY, loc="center right", handlelength=1.0)
    ps.panel_tag(a, "(b)")

    a = ax[2]
    styles = family_styles(D["s_refined"].astype(bool))
    for k, ((c, ls), lab) in enumerate(zip(styles, D["s_lab"])):
        n = int(D["s_n"][k])
        a.plot(D["s_x"][k, :n], D["s_P"][k, :n], color=c, ls=ls, lw=0.8, label=str(lab))
    a.set_xlim(-0.05, 0.05)
    a.set_ylim(-1.35, 1.6)
    a.set_xlabel(r"$\Theta-\Theta_{\rm peak}$")
    a.legend(fontsize=ps.FS_TINY, loc="upper right", handlelength=1.6)
    ps.panel_tag(a, "(c)", loc="lower left")

    a = ax[3]
    nf = len(D["f_n"])
    for k in range(nf):
        n = int(D["f_n"][k])
        a.plot(D["f_th"][k, :n], D["f_P"][k, :n], "k", ls=LINESTYLES[len(LINESTYLES) - nf + k],
               lw=0.8, label=str(D["f_lab"][k]))
    a.set_xlim(3.42, 3.68)       # the front of the paper's runs at t = 40 (Theta ~ 3.5)
    a.set_ylim(-0.3, 0.55)
    a.set_xlabel(r"$\Theta$")
    a.legend(fontsize=ps.FS_TINY, loc="upper right", handlelength=1.6)
    ps.panel_tag(a, "(d)", loc="lower left")
    fig.tight_layout(w_pad=0.6)
    ps.save(fig, out)


def pack(rows):
    """List of 1-d arrays -> (NaN-padded 2-d array, lengths)."""
    n = max(len(x) for x in rows)
    A = np.full((len(rows), n), np.nan)
    for k, x in enumerate(rows):
        A[k, :len(x)] = x
    return A, np.array([len(x) for x in rows])


def main():
    if ps.replot(plot, "sharp_features.png"):
        return
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--uniform", required=True, help="three uniform runs, coarse to fine")
    ap.add_argument("--refined", required=True, help="refined inner-shell runs (F = 8, F = 16)")
    ap.add_argument("--escaping", required=True, help="zone-A run on the escaping pulse")
    ap.add_argument("--t_over", type=float, default=40.0, help="time of panel (a)")
    ap.add_argument("--t_shell", type=float, default=30.0, help="time of panel (c)")
    ap.add_argument("--t_front", type=float, default=40.0, help="time of panel (d)")
    ap.add_argument("--t_steep", type=float, default=24.0,
                    help="earlier time of the steepest-slope printout (default 24)")
    ap.add_argument("--t_esc", default="40,50,64,72,80", help="times of panel (b)")
    ap.add_argument("--refined_F", default="8,16",
                    help="refinement factor F (--refB_F) of each --refined run, for the legend "
                         "(default 8,16)")
    ap.add_argument("--out", default="sharp_features.png",
                    help="output figure (png; a pdf is written next to it)")
    ps.add_replot_arg(ap)
    a = ap.parse_args()
    uni, ref = a.uniform.split(","), a.refined.split(",")
    ref_F = a.refined_F.split(",")
    if len(ref_F) != len(ref):
        raise SystemExit("--refined_F needs one factor per --refined run")

    p = sc.load_dump(uni[-1], a.t_over)
    Th, P = p["Theta"], p["P"]
    marks = [steepest(uni[-1], a.t_over)[1], Th[peak(p, 5.0, 15.0)], Th[peak(p, 15.0, 35.0)]]

    e_t = [float(x) for x in a.t_esc.split(",")]
    ex, eP = [], []
    for t in e_t:
        q = sc.load_dump(a.escaping, t)
        k = peak(q, 15.0, 1e9)
        m = np.abs(q["Theta"] - q["Theta"][k]) < 0.5
        ex.append(q["Theta"][m] - q["Theta"][k])
        eP.append(q["P"][m])
        print(f"escaping pulse t = {t:g}: Theta_peak = {q['Theta'][k]:.2f}, P = {q['P'][k]:.4f}")

    sx, sP, slab = [], [], []
    for d, F in [(d, None) for d in uni] + list(zip(ref, ref_F)):
        q = sc.load_dump(d, a.t_shell)
        k = peak(q, 5.0, 12.0)
        m = np.abs(q["Theta"] - q["Theta"][k]) < 0.1
        sx.append(q["Theta"][m] - q["Theta"][k])
        sP.append(q["P"][m])
        if F is not None:
            slab.append(rf"refined, $F={F}$")
        else:
            slab.append(f"uniform {q['r'][k + 1] - q['r'][k]:.4f}")    # spacing at the shell

    fth, fP, flab = [], [], []
    for d in uni:
        q = sc.load_dump(d, a.t_front)
        m = (q["Theta"] > 3.0) & (q["Theta"] < 4.2)
        j, _ = sc.steepest_slope(q)
        fth.append(q["Theta"][m])
        fP.append(q["P"][m])
        flab.append(rf"$dr={q['r'][j + 1] - q['r'][j]:.4f}$")   # spacing at the front

    print("inner front, steepest |dP/dTheta| (uniform grids, coarse to fine):")
    for t in (a.t_steep, a.t_front):
        g = [steepest(d, t) for d in uni]
        print(f"  t = {t:g}: " + "  ".join(f"{x[0]:.1f} (Theta {x[1]:.2f})" for x in g)
              + "   ratios " + " ".join(f"{g[i+1][0]/g[i][0]:.2f}" for i in range(2)))

    e_x, e_n = pack(ex)
    e_P, _ = pack(eP)
    s_x, s_n = pack(sx)
    s_P, _ = pack(sP)
    f_th, f_n = pack(fth)
    f_P, _ = pack(fP)
    ps.setup()
    D = dict(o_th=Th, o_P=P, o_mark_th=np.array(marks), t_over=a.t_over,
             e_t=np.array(e_t), e_x=e_x, e_P=e_P, e_n=e_n,
             s_x=s_x, s_P=s_P, s_n=s_n, s_lab=np.array(slab),
             s_refined=np.array([False] * len(uni) + [True] * len(ref)),
             f_th=f_th, f_P=f_P, f_n=f_n, f_lab=np.array(flab))
    ps.save_data(a.out, "sharp_features", **D)
    plot(D, a.out)


if __name__ == "__main__":
    main()
