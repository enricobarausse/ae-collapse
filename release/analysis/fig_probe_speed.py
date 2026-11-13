#!/usr/bin/env python3
"""Physical speed of a linear spin-0 probe on the collapse background, paper fig:probe (probe_speed.pdf)
(Sec. "The spin-0 waveform and its speed").

Runs (release/README.md, section 5, "Probe runs"): for each c2, a background run and the
same run with a small trK pulse added to the initial data (--probe_amp=1e-6
--probe_r=150; the pulse is linear, so the difference of the two runs is the probe
propagating on the dynamical background).  Dumps every 1 (or 0.5) time unit.

Definitions:
  delta X = X(probe run) - X(background run), dump by dump (on the moving mesh the probe
  run is first interpolated onto the r grid of the background run if the grids differ);
  delta w_out = -(b_N/v0) delta trK + delta a_r  (outgoing spin-0 characteristic, eq:wout);
  r_pk(t): position of the maximum of Theta |delta w_out| on the OUTGOING pulse, tracked
  from one dump to the next inside a window around the position predicted by the
  characteristic dr/dt = -beta + N v0 (parabolic refinement through 3 samples);
  physical speed  v(t) = (dr_pk/dt + beta(r_pk, t)) / N(r_pk, t), with dr_pk/dt centred
  differences in time, N from the background run, and beta rebuilt from the background
  run by integrating beta' = N (Q + trK/3) outward from beta(0) = 0.
Points with t < --tstart (the probe still overlaps the collapsing shell / ingoing half),
or with the peak within --edge of the outer boundary, or with the maximum on the edge of
the search window, are excluded.

Prints, per pair: median and MAD of v(t), the relative deviation of the median from v0,
mean and standard deviation, the largest |v - v0|/v0 after --t_late (default 40, the end
of the early scatter from the splitting of the initial pulse), and the constancy of
Theta |delta w_out| at the peak.
Paper: the deviation stays below 3e-4 after t ~ 40; the medians agree with v0 to 0.01 %
(v0 = 1) and 0.02 % (v0 = 1.267) (public-code runs: medians 1.00012 and 1.26706, MAD
4e-5 and 8e-5; largest deviation for t >= 40: 3.0e-4 and 1.9e-4).

Usage:
  fig_probe_speed.py --pair BG_DIR,PROBE_DIR,-0.5119403 --pair BG_DIR2,PROBE_DIR2,0
                     [--probe_r 150] [--tstart 15] [--out probe_speed.png]
  fig_probe_speed.py --replot data/probe_speed.npz [--out probe_speed.png]

Figure (single column): (v - v0)/v0 against t for every pair (in units of 1e-3), with the
medians (dotted) and v = v0 (dashed).
"""
import argparse
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import spin0_common as sc
import paper_style as ps


def track(bg, pr, cp, probe_r, tstart, edge, hw, verbose=True):
    """Track the outgoing probe peak -> dict of arrays: t, r (= r_pk), amp (Theta |delta
    w_out| at the peak), N, beta (there), ok (sample usable), v (physical speed) and good
    (ok with both neighbours ok, i.e. a valid centred dr_pk/dt)."""
    tA, pA = sc.list_dumps(bg)
    tB, pB = sc.list_dumps(pr)
    nA = {round(sc.nominal_time(p), 6): (t, p) for t, p in zip(tA, pA)}
    nB = {round(sc.nominal_time(p), 6): (t, p) for t, p in zip(tB, pB)}
    common = sorted(set(nA) & set(nB))
    if len(common) < 5:
        raise RuntimeError(f"only {len(common)} common dump times in {bg} and {pr}")
    v0, k_bn = cp["v0"], cp["bN"] / cp["v0"]
    rows = []
    r_prev, t_prev = probe_r, 0.0
    checked = False
    for tn in common:
        t, fa = nA[tn]
        tb, fb = nB[tn]
        if abs(t - tb) > 1e-9:
            print(f"   warning: dump {tn}: times differ ({t} vs {tb})", file=sys.stderr)
        A = sc.load_profile(fa)
        B = sc.load_profile(fb)
        r = A["r"]
        if len(B["r"]) == len(r) and np.array_equal(B["r"], r):
            dK, da = B["trK"] - A["trK"], B["Acr"] - A["Acr"]
        else:
            Ki, ai = sc.lagrange_interp(B["r"], [B["trK"], B["Acr"]], r)
            dK, da = Ki - A["trK"], ai - A["Acr"]
        dw = -k_bn * dK + da
        y = A["Theta"] * np.abs(dw)
        N = A["NN"]
        beta = sc.beta_shift(A)
        if not checked and verbose:
            # sanity: the difference must vanish outside the probe's causal cone (initial
            # half-width 40 ~ 10 probe widths; speed bound 3 > max(v0, 1) + |beta|)
            far = np.abs(r - probe_r) > 40.0 + 3.0 * t
            if far.any():
                print(f"   background cancellation at t = {t:.2f}: max |delta trK| outside "
                      f"the cone = {np.max(np.abs(dK[far])):.2e} "
                      f"(probe max {np.max(np.abs(dK)):.2e})")
            checked = True
        # predicted position: one step along the outgoing characteristic from the last peak
        Np, bp = sc.lagrange_interp(r, [N, beta], [r_prev])
        r_pred = r_prev + (t - t_prev) * (-bp[0] + Np[0] * v0)
        win = np.nonzero(np.abs(r - r_pred) <= hw)[0]
        if len(win) < 3:
            break
        j = win[np.argmax(y[win])]
        rpk, apk = sc.parabolic_peak(r, y, j)
        at_edge = j in (win[0], win[-1])
        Nk, bk = sc.lagrange_interp(r, [N, beta], [rpk])
        ok = (t >= tstart) and (rpk <= r[-1] - edge) and not at_edge
        rows.append((t, rpk, apk, Nk[0], bk[0], ok))
        r_prev, t_prev = rpk, t
        if rpk > r[-1] - 1.0:            # the peak has (nearly) reached the outer boundary
            break
    R = np.array(rows)
    out = dict(t=R[:, 0], r=R[:, 1], amp=R[:, 2], N=R[:, 3], beta=R[:, 4],
               ok=R[:, 5].astype(bool))
    drdt = np.gradient(out["r"], out["t"])
    out["v"] = (drdt + out["beta"]) / out["N"]
    # a centred difference needs both neighbours to be valid peaks
    good = out["ok"].copy()
    good[1:-1] &= out["ok"][:-2] & out["ok"][2:]
    good[[0, -1]] = False
    out["good"] = good
    return out


def plot(D, out):
    """Relative deviation (v - v0)/v0 of the probe speed for every pair, one panel
    (single column); dotted: the medians; dashed: zero (v = v0)."""
    fig, ax = plt.subplots(figsize=(ps.COLUMN, 2.4))
    cols, marks = (ps.CAT_A, ps.CAT_B, ps.GREEN, ps.PURPLE), ("o", "s", "^", "D")
    ax.axhline(0, color="k", ls="--", lw=0.9)
    for ip in range(int(D["npair"])):
        v0, med = D[f"v0{ip}"], D[f"med{ip}"]
        lab = D[f"label{ip}"] or rf"$v_0={v0:.3f}$"
        c = cols[ip % 4]
        ax.plot(D[f"t{ip}"], 1e3 * (D[f"v{ip}"] / v0 - 1), marks[ip % 4], color=c, ms=2.6,
                mew=0, alpha=0.9, label=lab)
        ax.axhline(1e3 * (med / v0 - 1), color=c, ls=":", lw=1.1)
    ax.set_xlabel("$t$")
    ax.set_ylabel(r"$10^3\,(v-v_0)/v_0$")
    ps.headroom(ax, 0.25)
    ax.legend(loc="upper right", ncol=2, markerscale=1.6, handletextpad=0.2)
    fig.tight_layout()
    ps.save(fig, out)


def main():
    if ps.replot(plot, "probe_speed.png"):
        return
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pair", action="append", required=True,
                    help="BG_DIR,PROBE_DIR,c2[,label]  (repeat: one curve per pair)")
    ap.add_argument("--c1", type=float, default=sc.C1_DEFAULT, help="coupling c1 (default 0.7)")
    ap.add_argument("--probe_r", type=float, default=150.0,
                    help="centre of the probe at t = 0 (the code's --probe_r)")
    ap.add_argument("--tstart", type=float, default=15.0,
                    help="drop samples before this time (probe still overlapping the shell)")
    ap.add_argument("--t_late", type=float, default=40.0,
                    help="print the largest deviation from v0 after this time (default 40)")
    ap.add_argument("--edge", type=float, default=20.0,
                    help="exclude peaks closer than this to the outer boundary")
    ap.add_argument("--hw", type=float, default=8.0,
                    help="half width of the search window around the predicted position")
    ap.add_argument("--out", default="probe_speed.png",
                    help="output figure (png; a pdf is written next to it)")
    ps.add_replot_arg(ap)
    a = ap.parse_args()

    D = dict(npair=len(a.pair))
    for ip, spec in enumerate(a.pair):
        f = spec.split(",")
        bg, pr, c2 = f[0], f[1], float(f[2])
        cp = sc.couplings(a.c1, c2)
        v0 = cp["v0"]
        print(f"\n{bg} / {pr}:  c2 = {c2}, v0 = {v0:.6f}")
        T = track(bg, pr, cp, a.probe_r, a.tstart, a.edge, a.hw)
        g = T["good"]
        v = T["v"][g]
        med = np.median(v)
        mad = np.median(np.abs(v - med))
        print(f"   {g.sum()} speed samples, t = {T['t'][g].min():.1f} .. {T['t'][g].max():.1f}")
        print(f"   median v = {med:.5f}, MAD = {mad:.5f} (1.4826 MAD = {1.4826 * mad:.5f});"
              f"  (median - v0)/v0 = {100 * (med / v0 - 1):+.3f}%   [paper: 0.01% (v0=1), "
              f"0.02% (v0=1.267)]")
        print(f"   mean v = {v.mean():.5f}, std = {v.std():.5f}")
        late = T["t"][g] >= a.t_late
        if late.any():
            print(f"   max |v - v0|/v0 for t >= {a.t_late:g}: {np.abs(v[late] / v0 - 1).max():.2e}"
                  f"   [paper: below 3e-4 after t ~ 40]")
        amp = T["amp"][g]
        print(f"   Theta|delta w_out| at the peak: {amp.min():.3e} .. {amp.max():.3e} "
              f"(spread {100 * (amp.max() - amp.min()) / amp.mean():.1f}%)")
        print(f"   metric on the track: N = {T['N'][g].min():.3f}..{T['N'][g].max():.3f}, "
              f"beta = {T['beta'][g].min():.3f}..{T['beta'][g].max():.3f}")
        D.update({f"t{ip}": T["t"][g], f"v{ip}": v, f"v0{ip}": v0, f"med{ip}": med,
                  f"mad{ip}": mad, f"label{ip}": f[3] if len(f) > 3 else ""})
    ps.setup()
    ps.save_data(a.out, "probe_speed", **D)
    plot(D, a.out)


if __name__ == "__main__":
    main()
