#!/usr/bin/env python3
"""Energy carried by the spin-0 wave, paper fig:fluxinf (spin0_flux_infinity.pdf) and eq:Einf
(Sec. "Mass conservation and radiated energy"; v0 = 1.267, c2 = 0, A0 = 0.35).

Definitions:
  m = Misner-Sharp mass [eq:ms];  m_c = m - (c14/2) Theta^2 a_r  [eq:mc];
  E(Theta, u) = M_ADM - m_c(Theta, u): energy that has crossed the sphere of areal radius
  Theta up to retarded time u = t - Theta_*/v0 [eq:retarded];
  L = -d_u m_c|_Theta = -d_t m_c|_Theta (centred 5-point differences across the dumps);
  E_inf(u), L_inf(u): least-squares fits E = E_inf + b/Theta (same for L) over the spheres
  (default Theta = 75, 100, 150, 200) at each u on a common u grid (each sphere's series is
  interpolated linearly in u), shown for u >= --u_fit (default 95; before that the pulse
  is crossing the spheres and there is no common 1/Theta form);
  radiated energy: E averaged over u in [--avg] (default 200..250) on each sphere, then
  fitted as E_inf + b/Theta over the spheres -> E_inf;
  uncertainty: half the peak-to-peak oscillation of E_inf(u) over u in [--osc]
  (default 150..250), the residual oscillation of m_c;
  scalar contribution: E_psi(Theta, u) = integral of F_psi dt on the sphere since the first
  dump, F_psi the scalar flux of eq:lum (rho = p = (P^2 + psi'^2)/2, j = -P psi'),
  averaged and extrapolated in the same way;
  corrected vs uncorrected rate: per sphere, the peak of |L| = |d_u m_c| against the peak
  of the Misner-Sharp luminosity |d_u m| (the quantity of fig:lum(a)) on
  the same sphere, and their ratio.
Optional resolution study (--ladder coarse,medium,fine): E_inf on each grid (dr = spacing
of the uniform inner grid), the observed order for a factor-2 ladder, and the Richardson
extrapolation with that order (and with order 2).

Paper (M_ADM = 3.708): E_inf = 0.532 +- 0.023 = (0.144 +- 0.006) M_ADM, scalar part
E_psi = 0.197, i.e. E_ae = 0.335 = 0.090 M_ADM (9 %; the conclusion rounds E_psi to 5 %);
the corrected rate is more than an order of magnitude smaller than the uncorrected
luminosity of fig:lum(a).  Also printed, not quoted in the paper
(public-code runs): E_inf = 0.6231 / 0.5443 / 0.5323 at dr = 0.0469 / 0.0234 / 0.0117,
observed order 2.71, Richardson 0.5302 (the runs before the scheme corrections of
VALIDATION.md gave 0.6219 / 0.5435 / 0.5317, order 2.73, 14.3 %); E_psi 0.197-0.202 on
Theta = 75-200; the 1/Theta fit of E(u) accurate to better than 1e-3 for u >~ 120; peak
|d_u m| / peak |d_u m_c| = 5.3 / 8.1 / 15 / 25 on Theta = 75 / 100 / 150 / 200 (the
spheres Theta = 120-400 of fig:lum(a) have peaks 2.3-7.5 against <= 0.27
here, ratio >~ 10).

Usage:  fig_energy.py RUNDIR --c2 0 [--madm M] [--ladder DIR047,DIR023,DIR0117]
        fig_energy.py --replot data/spin0_flux_infinity.npz [--out spin0_flux_infinity.png]
(without --madm, M_ADM is estimated from m(Theta) at large radius on the first dump)
"""
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import spin0_common as sc
import paper_style as ps


def fit_inv(thetas, Y):
    """Least-squares Y = Y_inf + b/Theta over the spheres; Y has shape (..., n_theta).
    Returns (Y_inf, b, max |residual|), NaN where any sphere is missing."""
    X = np.vstack([np.ones_like(thetas), 1.0 / thetas]).T
    Y = np.atleast_2d(Y)
    out = np.full((Y.shape[0], 3), np.nan)
    for i, y in enumerate(Y):
        if np.all(np.isfinite(y)):
            c, *_ = np.linalg.lstsq(X, y, rcond=None)
            out[i] = c[0], c[1], np.max(np.abs(X @ c - y))
    return out.T


def analyse(rundir, cp, M, thetas, nproc, ugrid, avg, osc):
    """Sphere series -> dict with u grid, E, L, E_psi per sphere, fits and E_inf.
    ugrid: u range of the common grid; avg: u window of the E average; osc: u window of
    the uncertainty (each a pair lo, hi)."""
    t, S = sc.sphere_series(rundir, thetas, ["m", "Acr", "Fpsi"], cp=cp, nproc=nproc)
    v0 = cp["v0"]
    mc = sc.m_coulomb(S["m"], thetas[None, :], S["Acr"], cp)
    L = -sc.ddt_series(t, mc, order=4)
    Lm = -sc.ddt_series(t, S["m"], order=4)          # uncorrected Misner-Sharp luminosity
    Fp = np.nan_to_num(S["Fpsi"])
    Epsi = np.concatenate([np.zeros((1, len(thetas))),
                           np.cumsum(0.5 * np.diff(t)[:, None] * (Fp[1:] + Fp[:-1]), axis=0)])
    ug = np.arange(ugrid[0], ugrid[1] + 1e-9, 0.5)
    E_u, L_u, Ep_u = (np.full((len(ug), len(thetas)), np.nan) for _ in range(3))
    for k, R in enumerate(thetas):
        u = sc.retarded_time(t, R, M, v0)
        ok = np.isfinite(mc[:, k])
        inside = (ug >= u[ok][0]) & (ug <= u[ok][-1])
        E_u[inside, k] = np.interp(ug[inside], u[ok], M - mc[ok, k])
        L_u[inside, k] = np.interp(ug[inside], u[ok], L[ok, k])
        Ep_u[inside, k] = np.interp(ug[inside], u[ok], Epsi[ok, k])
    Einf_u, _, resE_u = fit_inv(thetas, E_u)
    Linf_u, _, _ = fit_inv(thetas, L_u)
    w = (ug >= avg[0]) & (ug <= avg[1])
    Eavg = np.mean(E_u[w], axis=0)                 # NaN if a sphere does not cover the window
    Epavg = np.mean(Ep_u[w], axis=0)
    Einf, b, _ = (x[0] for x in fit_inv(thetas, Eavg))
    Epinf, _, _ = (x[0] for x in fit_inv(thetas, Epavg))
    o = (ug >= osc[0]) & (ug <= osc[1]) & np.isfinite(Einf_u)
    err = 0.5 * np.ptp(Einf_u[o]) if o.any() else np.nan
    p0 = sc.load_profile(sc.list_dumps(rundir)[1][0])
    return dict(t=t, ug=ug, E_u=E_u, L_u=L_u, Ep_u=Ep_u, Einf_u=Einf_u, Linf_u=Linf_u,
                resE_u=resE_u, Eavg=Eavg, Epavg=Epavg, Einf=Einf, b=b, Epinf=Epinf, err=err,
                dr=p0["r"][1] - p0["r"][0], Ep_end=Epsi[-1],
                Lpk=np.nanmax(np.abs(L), axis=0), Lmpk=np.nanmax(np.abs(Lm), axis=0))


def plot(D, out):
    """(a) L = -d_u m_c and (b) E = M_ADM - m_c on the spheres, with their 1/Theta
    extrapolations (black) for u >= u_fit; single column."""
    ug, M = D["u"], D["M"]
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(ps.COLUMN, 4.3), sharex=True)
    cols = ps.sphere_colors(D["thetas"])
    show = ug >= D["u_fit"]
    for k, th in enumerate(D["thetas"]):
        ax1.plot(ug, D["L"][:, k], color=cols[k], lw=0.9, label=rf"$\Theta={th:.0f}$")
        ax2.plot(ug, D["E"][:, k], color=cols[k], lw=0.9)
    ax1.plot(ug[show], D["Linf"][show], "k", lw=1.8, label=r"$\Theta\to\infty$")
    ax2.plot(ug[show], D["Einf_u"][show], "k", lw=1.8)
    ax2.axhspan(D["Einf"] - D["err"], D["Einf"] + D["err"], color="0.6", alpha=0.35, lw=0)
    for ax in (ax1, ax2):
        ax.axhline(0, color="k", lw=0.5)
        ax.axvspan(ug[0], D["u_fit"], color="0.92", lw=0, zorder=0)
    ax1.set_ylim(-0.08, 0.08)      # fixed y range of the paper's figure (chosen for its runs)
    ax1.yaxis.set_major_locator(plt.MultipleLocator(0.04))
    ax1.set_ylabel(r"$L=-\partial_u m_c|_\Theta$")
    ax1.legend(ncol=1, loc="lower right", fontsize=ps.FS_SMALL + 0.5, handlelength=1.4,
               labelspacing=0.2)
    ax2.set_ylabel(r"$E=M_{\rm ADM}-m_c$")
    ax2.set_xlabel(r"$u=t-\Theta_*/v_0$")
    ax2.set_xlim(ug[0], ug[-1])
    sec = ax2.secondary_yaxis("right", functions=(lambda e: e / M, lambda f: f * M))
    sec.set_ylabel(r"$E/M_{\rm ADM}$")
    sec.tick_params(which="both", direction="in")
    ax2.tick_params(which="both", right=False)
    ps.panel_tag(ax1, "(a)", loc="upper right")
    ps.panel_tag(ax2, "(b)", loc="lower right")
    fig.tight_layout(h_pad=0.3)
    ps.save(fig, out)


def main():
    if ps.replot(plot, "spin0_flux_infinity.png"):
        return
    ap = sc.add_common_args(argparse.ArgumentParser(description=__doc__,
                            formatter_class=argparse.RawDescriptionHelpFormatter),
                            out="spin0_flux_infinity.png")
    ap.set_defaults(c2=0.0)
    ap.add_argument("--thetas", default="75,100,150,200",
                    help="areal radii of the spheres (comma list)")
    ap.add_argument("--avg", default="200,250", help="u window for E")
    ap.add_argument("--osc", default="150,250", help="u window for the uncertainty")
    ap.add_argument("--u_fit", type=float, default=95.0,
                    help="show E_inf(u), L_inf(u) for u >= this")
    ap.add_argument("--ugrid", default="-80,255", help="u range of the figure")
    ap.add_argument("--ladder", default=None,
                    help="comma list of run dirs, coarse -> fine (factor 2 in dr)")
    ps.add_replot_arg(ap)
    a = ap.parse_args()
    a.avg = [float(x) for x in a.avg.split(",")]
    a.osc = [float(x) for x in a.osc.split(",")]
    a.ugrid = [float(x) for x in a.ugrid.split(",")]
    cp = sc.couplings(a.c1, a.c2)
    M = sc.resolve_madm(a)
    thetas = np.array([float(x) for x in a.thetas.split(",")])
    print(f"v0 = {cp['v0']:.4f}, c14 = {cp['c14']}, M_ADM = {M:.4f}")

    R = analyse(a.rundir, cp, M, thetas, a.nproc, a.ugrid, a.avg, a.osc)
    print(f"\nrun {a.rundir}  (inner dr = {R['dr']:.4f}, t <= {R['t'][-1]:.1f})")
    print(f"{'Theta':>6} {'<E>_u':>8} {'<E_psi>_u':>10} {'E_psi(t_end)':>13}")
    for k, th in enumerate(thetas):
        print(f"{th:6.0f} {R['Eavg'][k]:8.4f} {R['Epavg'][k]:10.2e} {R['Ep_end'][k]:13.2e}")
    print(f"E averaged over u in [{a.avg[0]:g},{a.avg[1]:g}], fit E_inf + b/Theta over "
          f"Theta = {thetas[0]:g}..{thetas[-1]:g}:  E_inf = {R['Einf']:.4f}, b = {R['b']:.3f}")
    print(f"uncertainty (half peak-to-peak of E_inf(u), u in [{a.osc[0]:g},{a.osc[1]:g}]): "
          f"+-{R['err']:.3f}")
    print(f"E_inf = {R['Einf']:.3f} +- {R['err']:.3f} = ({100 * R['Einf'] / M:.1f} +- "
          f"{100 * R['err'] / M:.1f})% of M_ADM = {M:.3f}   [paper: 0.532 +- 0.023 = (14.4+-0.6)%]")
    print(f"scalar contribution: E_psi extrapolated = {R['Epinf']:.2e} "
          f"(per sphere {np.nanmin(R['Epavg']):.2e}..{np.nanmax(R['Epavg']):.2e})"
          f"   [paper: 0.197]")
    for u0 in (100, 110, 120, 140):
        g = (R["ug"] >= u0) & np.isfinite(R["resE_u"])
        if g.any():
            print(f"1/Theta fit of E(u): max residual for u >= {u0}: "
                  f"{np.nanmax(R['resE_u'][g]):.2e}   [not in the paper: < 1e-3 for u >~ 120]")
    print("peak |L| = |d_u m_c| vs peak Misner-Sharp luminosity |d_u m| on each sphere   "
          "[paper: more than an order of magnitude smaller]")
    for k, th in enumerate(thetas):
        print(f"  Theta = {th:5.0f}: {R['Lpk'][k]:.3f} vs {R['Lmpk'][k]:.3f}, "
              f"ratio {R['Lmpk'][k] / R['Lpk'][k]:.1f}")

    if a.ladder:
        dirs = a.ladder.split(",")
        rows = []
        for d in dirs:
            Rd = R if d.rstrip("/") == a.rundir.rstrip("/") else analyse(
                d, cp, M, thetas, a.nproc, a.ugrid, a.avg, a.osc)
            rows.append((Rd["dr"], Rd["Einf"], Rd["err"]))
            print(f"   ladder {d}: dr = {Rd['dr']:.4f}  E_inf = {Rd['Einf']:.4f}  "
                  f"(+-{Rd['err']:.3f})")
        print("\nresolution table   [not in the paper; public-code runs: 0.6231 / 0.5443 / 0.5323, "
              "order 2.71, Richardson 0.5302]")
        for dr, e, _ in rows:
            print(f"   dr = {dr:.4f}   E_inf = {e:.4f}")
        if len(rows) >= 3:
            e1, e2, e3 = (x[1] for x in rows[-3:])
            q = (e1 - e2) / (e2 - e3)
            p = np.log2(q) if q > 0 else np.nan
            print(f"   observed order p = log2[(E1-E2)/(E2-E3)] = {p:.2f}")
            if np.isfinite(p):
                print(f"   Richardson (order p):  {e3 + (e3 - e2) / (2 ** p - 1):.4f}")
            print(f"   Richardson (order 2):  {e3 + (e3 - e2) / 3:.4f}")

    # ---------------- figure ----------------
    out = a.out
    ps.setup()
    D = dict(thetas=thetas, u=R["ug"], L=R["L_u"], E=R["E_u"], Linf=R["Linf_u"],
             Einf_u=R["Einf_u"], Einf=R["Einf"], err=R["err"], u_fit=a.u_fit, M=M)
    ps.save_data(out, "spin0_flux_infinity", **D)
    plot(D, out)


if __name__ == "__main__":
    main()
