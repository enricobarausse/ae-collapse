#!/usr/bin/env python3
"""Flux identity of the paper (fig:fluxid (flux_identity.pdf); Sec. "Mass conservation and radiated
energy") from the per-step SPHERE RECORDER.

Same check as fig_flux_identity.py,

    L = -dm/dt|_Theta = F_ae + F_psi,
    F_X = (N Theta^2/2) [ A (rho_X + p_X) + ((A^2 + B^2)/B) j_X ],
    A = Theta (Q/2 - K/3),  B = dTheta/dr                       [paper eq:lum],

but on the records written by `collapse --rec_spheres=Th1,... [--rec_every=k]`
(RUNDIR/spheres.bin, layout in RUNDIR/spheres.txt): at every recorded step, for every
sphere, the 13 grid points i-6..i+6 around the outermost crossing Theta[i] <= Theta_s <
Theta[i+1] of r, Psi, P, trK, Acr, Theta, N, Q and the shift beta.  So the time
derivatives are taken across steps (dt = CFL dr, or k dt) instead of across the
Delta t = 0.5 profile dumps, and the burst is resolved in time.

Per record (vectorised over records, one sphere at a time):
  * r-derivatives: second-order centred on the 13-point stencil (the same non-uniform
    formula as np.gradient(edge_order=2), applied row by row);
  * W = d_t ln N|_r, B = d_t beta|_r, Qt = d_t Q|_r: the neighbouring records are
    interpolated (cubic Lagrange in r) onto the r points of the central record -- this
    handles the moving mesh and the change of the crossing index i from step to step --
    and differentiated with centred Fornberg weights on the exact record times;
  * rho, j, p of the aether (spin0_common.aether_set) and of the scalar (scalar_set), the
    fluxes (sphere_flux) and m (misner_sharp, eq:ms) on the stencil, then cubic Lagrange
    interpolation in Theta to Theta_s;
  * -dm/dt|_Theta: centred Fornberg differences of m(Theta_s, t) over the records.
Prints per sphere the peak |dm/dt|, max and rms of |F_ae + F_psi + dm/dt| relative to the
peak, for the time-stencil orders --orders; the effect of the dump cadence on the dm/dt side
alone (m subsampled at the dump times, order 8); and, with --dumps RUNDIR (profile dumps
of the same run), the dump-based numbers of fig_flux_identity.py over the same time window.

Figure (paper fig:fluxid, the spheres of --plot): top, F_ae + F_psi, F_ae, F_psi and
-dm/dt|_Theta; bottom, the residual eps = (F_ae + F_psi + dm/dt)/max|dm/dt| in units of 1e-4
(grey dots: the dump-based residual, with --dumps).

Paper (dr = 0.023, v0 = 1.267): the two sides agree to within 0.02 % of the peak flux on
every sphere Theta = 40-200; |F_psi| <~ 0.02 against aether-flux peaks 0.9-4 (the F_ae
max/min columns: 0.86 at Theta = 40 to 3.87 at 200).  Public-code runs, not quoted in the
paper: 0.020 % at Theta = 40 falling to 0.002 % at 200; 0.098 % / 0.006 % at dr = 0.047.

Usage:  flux_identity_rec.py RUNDIR --c2 0 [--plot 40,100,200] [--dumps RUNDIR]
                             [--out flux_identity.png]
        flux_identity_rec.py --replot data/flux_identity.npz [--out flux_identity.png]
"""
import argparse
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import spin0_common as sc
import paper_style as ps

# the recorded fields, in the order of spheres.bin (NR = the shift beta)
FIELDS = ("r", "Psi", "P", "trK", "Acr", "Theta", "NN", "Q", "NR")
NF = len(FIELDS)


# ----------------------------------------------------------------------------------------
# reading the recorder output
# ----------------------------------------------------------------------------------------
def read_header(rundir):
    """Layout and run parameters of the recorder (RUNDIR/spheres.txt) as a dict."""
    h = {}
    with open(os.path.join(rundir, "spheres.txt")) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            k, *v = line.split()
            h[k] = v
    assert h["format"] == ["float64_le"], h["format"]
    assert h["fields"] == list(FIELDS), h["fields"]
    out = dict(ns=int(h["n_spheres"][0]), nst=int(h["n_stencil"][0]),
               L=int(h["record_doubles"][0]), N=int(h["N"][0]),
               thetas=np.array([float(x) for x in h["theta_s"]]),
               rec_every=int(h["rec_every"][0]))
    for k in ("c1", "c2", "A0", "rmax"):
        out[k] = float(h[k][0])
    assert out["L"] == 3 + out["ns"] * (1 + NF * out["nst"])
    return out


def read_records(rundir, tmin=-np.inf, tmax=np.inf):
    """(header, step, t, dt, raw) with raw a (nrec, L) memmap; a partial last record (run
    killed while writing) is dropped."""
    h = read_header(rundir)
    path = os.path.join(rundir, "spheres.bin")
    nrec = os.path.getsize(path) // (8 * h["L"])
    raw = np.memmap(path, dtype="<f8", mode="r", shape=(nrec, h["L"]))
    t = np.array(raw[:, 1])
    sel = np.nonzero((t >= tmin) & (t <= tmax))[0]
    raw = raw[sel[0]:sel[-1] + 1]
    return h, np.array(raw[:, 0]), np.array(raw[:, 1]), np.array(raw[:, 2]), raw


def sphere_block(h, raw, s):
    """i_lo (nrec,) and dict field -> (nrec, nst) for sphere number s."""
    nst = h["nst"]
    o = 3 + s * (1 + NF * nst)
    i_lo = np.array(raw[:, o]).astype(int)
    blk = np.array(raw[:, o + 1:o + 1 + NF * nst]).reshape(-1, NF, nst)
    return i_lo, {k: blk[:, f, :] for f, k in enumerate(FIELDS)}


# ----------------------------------------------------------------------------------------
# row-wise numerics on the (nrec, nst) stencils
# ----------------------------------------------------------------------------------------
def ddr_rows(y, r):
    """d/dr along the last axis on a non-uniform grid (r may differ row by row): the
    np.gradient(edge_order=2) formulas -- second-order centred inside, second-order
    one-sided at the two ends."""
    y = np.asarray(y, float)
    r = np.asarray(r, float)
    if y.ndim == 1 and r.ndim == 1:
        return np.gradient(y, r, edge_order=2)
    r = np.broadcast_to(r, y.shape)
    out = np.empty_like(y)
    h1 = r[..., 1:-1] - r[..., :-2]
    h2 = r[..., 2:] - r[..., 1:-1]
    out[..., 1:-1] = (-h2 / (h1 * (h1 + h2)) * y[..., :-2]
                      + (h2 - h1) / (h1 * h2) * y[..., 1:-1]
                      + h1 / (h2 * (h1 + h2)) * y[..., 2:])
    a, b = r[..., 1] - r[..., 0], r[..., 2] - r[..., 1]
    out[..., 0] = (-(2 * a + b) / (a * (a + b)) * y[..., 0] + (a + b) / (a * b) * y[..., 1]
                   - a / (b * (a + b)) * y[..., 2])
    a, b = r[..., -2] - r[..., -3], r[..., -1] - r[..., -2]
    out[..., -1] = (b / (a * (a + b)) * y[..., -3] - (a + b) / (a * b) * y[..., -2]
                    + (a + 2 * b) / (b * (a + b)) * y[..., -1])
    return out


def interp_rows(x, y, xq, npts=4):
    """Row-wise piecewise Lagrange interpolation (default cubic): x, y (n, m) with x
    increasing along each row, xq (n, q) -> (n, q); NaN outside [x[:,0], x[:,-1]]."""
    n, m = x.shape
    j = np.sum(x[:, None, :] < xq[:, :, None], axis=2) - 1
    lo = np.clip(j - (npts // 2 - 1), 0, m - npts)
    idx = lo[:, :, None] + np.arange(npts)[None, None, :]
    X = np.take_along_axis(x[:, None, :].repeat(xq.shape[1], 1), idx, axis=2)
    Y = np.take_along_axis(y[:, None, :].repeat(xq.shape[1], 1), idx, axis=2)
    w = np.ones_like(X)
    for a in range(npts):
        for b in range(npts):
            if a != b:
                w[..., a] *= (xq - X[..., b]) / (X[..., a] - X[..., b])
    v = (w * Y).sum(axis=2)
    v[(xq < x[:, :1]) | (xq > x[:, -1:])] = np.nan
    return v


def time_weights(t, order):
    """Centred first-derivative Fornberg weights (n, order+1) on the record times; NaN
    rows where the stencil would leave the series.  Uniform spacing reuses one set."""
    h = order // 2
    n = len(t)
    W = np.full((n, 2 * h + 1), np.nan)
    cache = {}
    for i in range(h, n - h):
        x = t[i - h:i + h + 1] - t[i]
        key = tuple(np.round(x / x[-1], 12))
        if key not in cache:
            cache[key] = sc.fd_weights(0.0, x / x[-1], 1)
        W[i] = cache[key] / x[-1]
    return W


def shifted(a, k):
    """a[n + k] at row n (NaN beyond the ends)."""
    out = np.full_like(a, np.nan)
    if k > 0:
        out[:-k] = a[k:]
    elif k < 0:
        out[-k:] = a[:k]
    else:
        out[:] = a
    return out


def ddt_fixed_r(F, Wt):
    """d_t at fixed r of lnN, beta, Q on the stencil points of each record."""
    h = Wt.shape[1] // 2
    r = F["r"]
    g = {"lnN": np.log(np.maximum(F["NN"], 1e-300)), "beta": F["NR"], "Q": F["Q"]}
    out = {k: np.zeros_like(r) for k in g}
    for c, k in enumerate(range(-h, h + 1)):
        rk = shifted(r, k)
        ok = np.isfinite(rk[:, 0])
        for name, v in g.items():
            G = np.full_like(r, np.nan)
            if k == 0:
                G = v
            else:
                G[ok] = interp_rows(rk[ok], shifted(v, k)[ok], r[ok])
            out[name] = out[name] + Wt[:, c:c + 1] * G
    return out


def sphere_series_rec(h, raw, s, cp, Wt):
    """m, F_ae, F_psi at Theta_s for sphere s at every record (NaN where invalid); Wt =
    time_weights of the records (for the time derivatives at fixed r)."""
    i_lo, F = sphere_block(h, raw, s)
    ths = h["thetas"][s]
    r = F["r"]
    d = {k: ddr_rows(F[k], r) for k in ("Theta", "Acr", "trK", "Q", "Psi")}
    dt = ddt_fixed_r(F, Wt)
    rho, j, pr = sc.aether_set(F, cp, dt["lnN"], dt["beta"], dt["Q"], beta=F["NR"], d=d,
                               deriv=ddr_rows)
    Fae = sc.sphere_flux(F, rho, j, pr, d["Theta"])
    rps, jps, pps = sc.scalar_set(F, d["Psi"])
    Fps = sc.sphere_flux(F, rps, jps, pps, d["Theta"])
    m = sc.misner_sharp(F, d["Theta"])
    valid = (i_lo >= h["nst"] // 2) & (i_lo + h["nst"] // 2 <= h["N"] - 1)
    q = np.full((len(i_lo), 1), ths)
    out = {}
    for name, v in (("m", m), ("Fae", Fae), ("Fps", Fps)):
        x = interp_rows(F["Theta"], v, q)[:, 0]
        x[~valid] = np.nan
        out[name] = x
    return out


def ddt(y, Wt):
    """d/dt of the record series y with the weights Wt of time_weights."""
    h = Wt.shape[1] // 2
    out = np.zeros_like(y)
    for c, k in enumerate(range(-h, h + 1)):
        out = out + Wt[:, c] * shifted(y, k)
    return out


def stats(res, dmdt):
    """(peak |dm/dt|, max |res| / peak, rms res / peak) over the finite entries."""
    ok = np.isfinite(res) & np.isfinite(dmdt)
    pk = np.max(np.abs(dmdt[ok]))
    return pk, np.max(np.abs(res[ok])) / pk, np.sqrt(np.mean(res[ok] ** 2)) / pk


# ----------------------------------------------------------------------------------------
def plot(D, out):
    """Top: F_ae + F_psi, -dm/dt, F_ae, F_psi on each sphere; bottom: the residual
    eps = (F_ae + F_psi + dm/dt)/max|dm/dt| in units of 1e-4."""
    t = D["t"]
    n = len(D["thetas"])
    fig, axs = plt.subplots(2, n, figsize=(ps.TEXT, 3.7), sharex=True,
                            gridspec_kw=dict(height_ratios=[1.7, 1]), squeeze=False)
    ae = ps.AE()
    for c in range(n):
        ax, axr = axs[0, c], axs[1, c]
        ax.plot(t, D["Fae"][c] + D["Fps"][c], color="0.75", lw=3.2,
                label=rf"$F_{ae}+F_\psi$")
        ax.plot(t, D["Fae"][c], color=ps.CAT_A, lw=0.9, label=rf"$F_{ae}$")
        ax.plot(t, D["Fps"][c], color=ps.CAT_B, lw=0.9, label=r"$F_\psi$")
        ax.plot(t, -D["dmdt"][c], "k--", lw=0.8, dashes=(3, 2), label=r"$-dm/dt|_\Theta$")
        ax.axhline(0, color="k", lw=0.5)
        ps.panel_tag(ax, rf"$\Theta={D['thetas'][c]:.0f}$", loc="upper right")
        pk = np.nanmax(np.abs(D["dmdt"][c]))
        axr.plot(t, 1e4 * D["res"][c] / pk, color="k", lw=0.8)
        if "t_dump" in D:
            axr.plot(D["t_dump"], 1e4 * D["res_dump"][c] / pk, ".", color="0.5", ms=2)
        axr.axhline(0, color="k", lw=0.5)
        axr.set_xlabel(r"$t$")
        axr.set_xlim(t[0], t[-1])
    axs[0, 0].set_ylabel("energy flux")
    axs[1, 0].set_ylabel(r"$10^4\,\epsilon$")
    axs[0, 0].legend(loc="lower right", fontsize=ps.FS_SMALL + 0.5, handlelength=1.6)
    fig.tight_layout(w_pad=0.8, h_pad=0.3)
    ps.save(fig, out)


def main():
    if ps.replot(plot, "flux_identity.png"):
        return
    ap = sc.add_common_args(argparse.ArgumentParser(description=__doc__,
                            formatter_class=argparse.RawDescriptionHelpFormatter), madm=False,
                            out="flux_identity.png")
    ap.add_argument("--orders", default="2,4",
                    help="orders of the centred time differences across records; the LAST "
                         "one is plotted and quoted")
    ap.add_argument("--plot", default=None, help="spheres shown in the figure (default: all, "
                    "at most 3)")
    ap.add_argument("--dumps", default=None,
                    help="run directory with the profile dumps of the same run: also print "
                         "the dump-based residuals of fig_flux_identity.py")
    ap.add_argument("--dump_dt", type=float, default=0.5,
                    help="dump cadence emulated on the m series (default 0.5)")
    ap.add_argument("--tmin", type=float, default=-np.inf, help="first record time used")
    ap.add_argument("--tmax", type=float, default=np.inf, help="last record time used")
    ps.add_replot_arg(ap)
    a = ap.parse_args()
    cp = sc.couplings(a.c1, a.c2)
    h, _, t, _, raw = read_records(a.rundir, a.tmin, a.tmax)
    if abs(h["c1"] - a.c1) > 1e-12 or abs(h["c2"] - a.c2) > 1e-6:
        print(f"WARNING: couplings in spheres.txt (c1={h['c1']}, c2={h['c2']}) differ from "
              f"--c1/--c2 ({a.c1}, {a.c2})", file=sys.stderr)
    thetas = h["thetas"]
    orders = [int(x) for x in a.orders.split(",")]
    order = orders[-1]
    print(f"{len(t)} records (rec_every={h['rec_every']}), t = {t[0]:.4f} .. {t[-1]:.4f}, "
          f"record spacing {np.median(np.diff(t)):.5g};  N={h['N']} rmax={h['rmax']:g}  "
          f"dr={h['rmax'] / (h['N'] - 1):.5g};  v0 = {cp['v0']:.4f}")

    S = {}
    print("\nmax / rms residual |F_ae+F_psi+dm/dt| relative to the peak, vs order of the "
          "time differences across records:")
    for o in orders:
        Wt = time_weights(t, o)
        line = []
        for s, R in enumerate(thetas):
            v = sphere_series_rec(h, raw, s, cp, Wt)
            v["dmdt"] = ddt(v["m"], Wt)
            v["res"] = v["Fae"] + v["Fps"] + v["dmdt"]
            if o == order:
                S[s] = v
            _, mx, rms = stats(v["res"], v["dmdt"])
            line.append(f"{R:.0f}: {100 * mx:.3f}% / {100 * rms:.4f}%")
        print(f"   order {o}:  " + "   ".join(line))

    print(f"\n[recorder, order {order}]")
    print(f"{'Theta':>6} {'max|dm/dt|':>11} {'F_ae max':>9} {'F_ae min':>9} "
          f"{'max|F_psi|':>11} {'max|res|':>10} {'res/peak':>9} {'rms res/peak':>12}")
    for s, R in enumerate(thetas):
        v = S[s]
        pk, mx, rms = stats(v["res"], v["dmdt"])
        print(f"{R:6.0f} {pk:11.4f} {np.nanmax(v['Fae']):9.4f} {np.nanmin(v['Fae']):9.4f} "
              f"{np.nanmax(np.abs(v['Fps'])):11.2e} {mx * pk:10.2e} {100 * mx:8.3f}% "
              f"{100 * rms:11.4f}%")

    # the dump cadence on the dm/dt side alone: m sampled at the dump times (first record
    # with t >= k*dump_dt, as the code writes its profiles), order-8 centred differences,
    # against the recorder fluxes at those records
    if a.dump_dt > 0:
        treq = np.arange(np.ceil(t[0] / a.dump_dt), np.floor(t[-1] / a.dump_dt) + 1) * a.dump_dt
        k = np.unique(np.searchsorted(t, treq - 1e-12))
        k = k[k < len(t)]
        print(f"\n[dm/dt from m sampled every {a.dump_dt:g} (order 8), fluxes from the recorder]")
        for s, R in enumerate(thetas):
            v = S[s]
            dm8 = sc.ddt_series(t[k], v["m"][k], order=8)
            ok = np.zeros(len(k), bool)
            ok[4:-4] = True               # centred order-8 stencil: half-width 4
            res = np.where(ok, v["Fae"][k] + v["Fps"][k] + dm8, np.nan)
            pk, mx, rms = stats(res, np.where(ok, dm8, np.nan))
            print(f"   Theta={R:.0f}: max {100 * mx:.3f}%   rms {100 * rms:.4f}%")

    D = None
    if a.dumps:
        import fig_flux_identity as ffi
        print(f"\n[dumps of {a.dumps}: fig_flux_identity.py, same spheres]")
        tt, SD, dmD, resD, okD, _ = ffi.dump_series(a.dumps, thetas, cp, (4, 8), a.nproc,
                                                    a.tmin, a.tmax)
        D = (tt, SD, dmD, resD, okD)
        lo, hi = max(t[0], tt[okD][0]), min(t[-1], tt[okD][-1])
        print(f"\ncomparison over the common window t = {lo:.2f} .. {hi:.2f} "
              f"(residual / peak |dm/dt|):")
        print(f"{'Theta':>6} {'rec max':>9} {'rec rms':>9} {'dump max':>9} {'dump rms':>9}")
        wr = (t >= lo) & (t <= hi)
        wd = okD & (tt >= lo) & (tt <= hi)
        for s, R in enumerate(thetas):
            v = S[s]
            _, mx, rms = stats(np.where(wr, v["res"], np.nan), np.where(wr, v["dmdt"], np.nan))
            _, mxd, rmsd = stats(np.where(wd, resD[:, s], np.nan), np.where(wd, dmD[:, s], np.nan))
            print(f"{R:6.0f} {100 * mx:8.3f}% {100 * rms:8.4f}% {100 * mxd:8.3f}% {100 * rmsd:8.4f}%")

    # ---------------- figure (paper fig:fluxid) ----------------
    show = ([float(x) for x in a.plot.split(",")] if a.plot else list(thetas[:3]))
    idx = [int(np.argmin(np.abs(thetas - R))) for R in show]
    Dp = dict(t=t, thetas=thetas[idx],
              Fae=np.array([S[s]["Fae"] for s in idx]), Fps=np.array([S[s]["Fps"] for s in idx]),
              dmdt=np.array([S[s]["dmdt"] for s in idx]), res=np.array([S[s]["res"] for s in idx]))
    if D is not None:
        tt, _, _, resD, okD = D
        Dp.update(t_dump=tt[okD], res_dump=resD[okD][:, idx].T)
    out = a.out
    ps.setup()
    ps.save_data(out, "flux_identity", **Dp)
    plot(Dp, out)


if __name__ == "__main__":
    main()
