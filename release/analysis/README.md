# Analysis scripts

Scripts that produce every figure and number of the paper from the profile dumps (and the
sphere records) of `collapse`: the convergence, mass-function, moving-mesh and white-hole
material (second table below) and the spin-0 material of Sec. "The spin-0 mode" (first
table). Where the paper rounds a number or describes it in words, the tables give the
paper's statement and the value the script prints on the final public-code runs
(`VALIDATION.md`, "Final rerun"); values marked "not in the paper" are printed as
additional checks.

## Figures of the paper

The paper's figures and tables are cited by their LaTeX labels (`fig:…`, `tab:…`, `eq:…`,
`sec:…` in the paper source), here and in the scripts, since their numbers move between
versions.

| paper | file | script | runs (README §5) |
|---|---|---|---|
| fig:mesh | `moving_mesh.pdf` | `fig_moving_mesh.py` | Moving-mesh validation |
| fig:threshold | `bh_threshold.pdf` | `fig_bh_threshold.py` | Black-hole threshold |
| fig:mass | `mass_function.pdf` | `fig_mass_function.py` | Mass function and horizon |
| fig:wout | `wout_background.pdf` | `fig_wout_background.py` | Spin-0 wave, v0 = 1 baseline |
| fig:probe | `probe_speed.pdf` | `fig_probe_speed.py` | Probe runs |
| fig:lum | `luminosity_clock.pdf` | `fig_luminosity_clock.py` | Superluminal spin-0, dr = 0.023 |
| fig:fluxid | `flux_identity.pdf` | `flux_identity_rec.py --plot 40,100,200 --out flux_identity.png` | Superluminal spin-0, recorder run at dr = 0.023 |
| fig:coulomb | `aether_coulomb.pdf` | `fig_coulomb.py` | Superluminal spin-0, dr = 0.0117 |
| fig:fluxinf | `spin0_flux_infinity.pdf` | `fig_energy.py` | Superluminal spin-0, dr = 0.0117 (ladder: all three) |
| fig:sharp | `sharp_features.pdf` | `fig_sharp_features.py` | Inner-shell ladder and escaping-pulse runs |
| fig:wh | `white_hole.pdf` | `fig_white_hole.py` | Inner-shell ladder |
| fig:conv-fields-A001, fig:conv-fields-A006 | `conv_scatter_fields_A001.pdf`, `conv_scatter_fields_A006.pdf` | `fig_conv_scatter.py` | Scattering convergence ladder |
| fig:conv-constraint | `conv_scatter_constraint.pdf` | `fig_conv_scatter.py` | Scattering convergence ladder |
| fig:conv-collapse | `conv_collapse.pdf` | `fig_conv_collapse.py` | Collapse convergence ladder (+ mass_hi) |
| — | `invariant_wout_dm.pdf` | `fig_invariant.py` | **not in the paper**; kept as a supplementary check |
| — | `flux_identity_dumps.png` of `fig_flux_identity.py`, `conv_scatter_orders.png` | | **not in the paper** (dump-based flux identity; order table as a colour grid) |

Tables tab:conv-scatter and tab:conv-ratio are printed by `fig_conv_scatter.py` and
`fig_conv_collapse.py`. Plain numpy + matplotlib (no scipy/pandas/sympy): Python ≥ 3.9, numpy ≥ 2.0
(`np.trapezoid`; tested with 2.5) and matplotlib ≥ 3.5 (tested with 3.11). numpy is pinned
to one thread, and the spin-0 scripts read the dumps with `--nproc` forked worker
processes (default 4; the fork start method, i.e. Linux or macOS). Every script takes
`--help`, and most print the paper's value next to their own (as `[paper: …]`; the
classification scripts `scan_bh_threshold.py`, `bisect_bh_threshold.py`,
`fig_bh_threshold.py`, `family_transitions.py` and `wh_family_ladder.py` print tables,
whose paper values are given in the tables below).

## Input: run directories

One directory per run, as written by `collapse` (README section 5 gives the exact
command lines):

    RUNDIR/prof_tTTT.T.csv   one per requested time, columns
                             r, Psi, P, trK, Acr, Theta, NN, Q, ThetaC
    RUNDIR/phi.dat           (optional) step, t, dt, phi at every step

* The spin-0 analyses need dumps every Δt = 0.5 (`--times=T(0.5,600,0.5)`); the probe
  runs use Δt = 1.
* With the moving mesh (`--ale=1`, default) the `r` column is non-uniform and moves from
  dump to dump. All radial derivatives are `np.gradient(f, r)` (second order at the end
  points, `edge_order=2`, in `spin0_common.ddr` and `apparent_horizon`); quantities "at
  fixed Θ" are cubic Lagrange interpolations in the areal radius over the monotonic exterior;
  time derivatives "at fixed r" interpolate the neighbouring dumps onto the r grid of
  the central one.
* A dump is written at the first step with t ≥ t_req and named after t_req with one
  decimal (29.25 → `prof_t029.2.csv`); if `phi.dat` is present, the exact dump times are
  read from it (otherwise the nominal times of the file names are used). Scripts that use
  quarter-time dumps open them by their requested time (`spin0_common.requested_dumps`,
  `load_dump`; the module docstring says which call to use when).
* Log files the scripts read if present: `log.txt` (the run's stdout: the code's
  "Horizon detected" lines, `scan_bh_threshold.py --log`), `err.txt` (the run's stderr:
  FATAL stops for `family_transitions.py --err`; `boundary_reflection.py spin0` checks its
  `--c1/--c2` against the line "spin-0 speed v0 = …") and `cmd.txt` (an
  optional copy of the command line, from which `fig_bh_threshold.py` and
  `bisect_bh_threshold.py` read `--A0=`; otherwise the amplitude is taken from the directory
  name). `collapse` writes none of them: redirect stdout and stderr when you launch it
  (`./collapse … > log.txt 2> err.txt`).
* Optional: `export SPIN0_CACHE=/some/dir` caches the sphere time series (.npz), so that
  re-running a script, or a second script with the same spheres, does not re-read the
  dumps.

## Figure style and replotting

`paper_style.py` is the common style of every figure: LaTeX text in Computer Modern (falls
back to the mathtext Computer Modern fonts if LaTeX is missing), figures designed at their
printed size (3.4 in for a single-column figure, 7.0 in for a two-column figure*, so the
9–10 pt fonts are the printed sizes), inward ticks, no grid, no titles (panel tags (a), (b),
... inside the axes), and colours tied to meanings: grid resolutions coarse → finest in
blue, green, red, purple (markers square, circle, triangle, diamond); pairs (fixed grid /
moving mesh, uniform / refined grids, æther / scalar flux, v0 = 1 / 1.267) in blue / red;
families of spheres of fixed Θ in the "plasma" map from dark violet (innermost) to orange
(outermost); extrapolations to infinity in black.

Every figure is written as `NAME.png` (300 dpi) and `NAME.pdf` (vector, for the paper),
and the plotted arrays are saved to `data/NAME.npz` next to the figure. Every figure
script takes `--replot data/NAME.npz [--out NAME.png]` (or `--outdir DIR` for
`fig_conv_scatter.py`) and redraws the figure from that file alone, without the run
directories (the numbers are printed only by a full run).

## Shared module

`spin0_common.py` (shared by all scripts): couplings (a_N, b_N, v0, eq:aNbN), dump
listing/loading (`list_dumps`, `requested_dumps`, `dump_path`, `load_dump`,
`load_profile`), ∂_r, Fornberg time-difference weights, interpolation at fixed Θ, θ±,
Misner–Sharp m (eq:ms), w_out / w_in (eq:wout), m_c (eq:mc), tortoise radius and retarded
time (eq:retarded), shift β (β' = N(Q + K/3), β(0) = 0), the aether stress tensor in the
aether frame (`aether_set`, returning ρ_æ, j_æ, p_æ; ported from the validated symbolic
reduction of arXiv:1501.07274 eq. 2.34, c3 = c4 = 0), the scalar one (`scalar_set`), the
flux of eq:lum (`sphere_flux`), the Gaussian running mean, `map_dumps` and
`sphere_series`; for the convergence, threshold and white-hole scripts also
`apparent_horizon`, `peak_index`, `steepest_slope`, `fwhm`, `norm_l2` and `run_amplitude`.

### Horizon criteria

All scripts locate horizons with one criterion: the outermost zero of θ₊ at proper radius
r > 0.5 (`spin0_common.apparent_horizon`; the cut excludes the few cells at the centre,
where a frozen central throat makes θ₊ change sign). No cut on Θ is applied by the finder.
On top of it:

* `bisect_bh_threshold.py` counts a run as a black hole only if Θ_AH > `--th_min` (0.5),
  i.e. a horizon detached from the frozen central throat (the throat can sit at r > 0.5
  with Θ ≈ 0.1);
* `fig_bh_threshold.py` plots Θ_AH as found (a throat shows as a small value; 0 where
  there is no zero);
* `scan_bh_threshold.py` takes the minimum of θ₊ over r > `--rmin` (0.5, the same cut).

## Spin-0 scripts, figures and printed numbers

All take the run directory as first argument (except `fig_probe_speed.py`), `--c1`
(default 0.7), `--c2`, `--madm` (default: estimated from m(Θ) at large radius on the first
dump), `--nproc`, `--out`. `--help` prints the docstring with the exact definitions.

| script | figure | run | prints (paper values) |
|---|---|---|---|
| `fig_wout_background.py RUN [--madm 3.711]` | `wout_background.png` (fig:wout) | v0 = 1 baseline (c2 = −0.5119403), dr = 0.0117 (N = 25600, r_max = 300), t ≤ 300 | pulse Θδw_out extrema per sphere Θ = 50…180 and their u; amplitude spread across spheres; drift in u; σ scan 5–15 (paper: the amplitude depends on σ but is the same on every sphere for any σ in [5, 15]; the pulse drifts slowly in u. Not in the paper: +0.061/−0.072, trough constant to a few %, drift 10.7/11.3) |
| `fig_invariant.py RUN --c2 0 [--madm M]` | `invariant_wout_dm.png` (**not in the paper**) | v0 = 1.267 | Θw_out extrema and δm extrema on Θ = 120…360 (profiles overlap in u; δm grows with Θ) |
| `fig_flux_identity.py RUN --c2 0` | `flux_identity_dumps.png` (**not** the paper's fig:fluxid, which is made by `flux_identity_rec.py`; a quick-look diagnostic outside the paper's figure style) | v0 = 1.267, t ≤ 600, Δt = 0.5 dumps | per sphere Θ = 40…200: peak flux, F_æ range, max ∣F_ψ∣, max residual ∣F_æ + F_ψ + dm/dt∣ / peak; residual vs order of the time stencil (∣F_ψ∣ ≲ 0.02 vs peaks 0.9–4; the residual from the dumps is 0.6–3 % at the burst, a time-sampling error — the paper's figure and number come from `flux_identity_rec.py`) |
| `flux_identity_rec.py RUN --c2 0 [--plot 40,100,200] [--dumps RUN] [--out flux_identity.png]` | `flux_identity.png` (= paper fig:fluxid, `flux_identity.pdf`, with `--plot 40,100,200`; without `--plot` the first three spheres are drawn) | a run with the sphere recorder (`--rec_spheres=40,60,…,200 --rec_every=5`), dr = 0.023 | same identity from the per-step records (`spheres.bin`): per sphere peak \|dm/dt\|, F_æ range, max \|F_ψ\|, max and rms residual / peak for the record-stencil orders 2, 4 (paper: within 0.02 % of the peak on every sphere Θ = 40–200, \|F_ψ\| ≲ 0.02 against æther-flux peaks 0.9–4; printed 0.020 % at Θ = 40 → 0.002 % at Θ = 200, \|F_ψ\| ≤ 0.0197, peaks 0.86–3.87; not in the paper: 0.098 % / 0.006 % at dr = 0.047); the residual with dm/dt from m resampled at the dump cadence; with `--dumps` the dump-based numbers of `fig_flux_identity.py` on the same window |
| `fig_coulomb.py RUN --c2 0 [--madm M]` | `aether_coulomb.png` (fig:coulomb) | v0 = 1.267 | per sphere Θ = 75…200: max Θ²a_r (19.7 → 40.5) and its power law in Θ (paper: roughly Θ^0.7; printed Θ^0.73), late Θ²a_r (paper: the same on every sphere; printed 4.6–5.0), swing / peak excursion of m (growing with radius, 7.1 → 18.8) and of m_c (0.9 → 0.7), late m − M_ADM (> 0, so ∫L du < 0) and m_c − M_ADM (< 0) |
| `fig_luminosity_clock.py RUN --c2 0 [--madm M]` | `luminosity_clock.png` (fig:lum) | v0 = 1.267, needs Θ up to 400 | burst arrival times on Θ = 120…400 and their retarded times (paper: burst at u ≃ 10–15; printed u = 9.7–14.7), fitted speed (paper 1.241 vs v0 = 1.267; printed 1.2406), and the slope 1 that the scalar speed would give (not in the paper) |
| `fig_energy.py RUN --c2 0 [--madm M] [--ladder R047,R023,R0117]` | `spin0_flux_infinity.png` (fig:fluxinf) | v0 = 1.267, finest grid for the figure | ⟨E⟩ over u ∈ [200,250] per sphere, E_∞ from E_∞ + b/Θ, uncertainty (half peak-to-peak of E_∞(u), u ∈ [150,250]), E_∞/M_ADM, scalar contribution E_ψ = ∫F_ψ dt, 1/Θ-fit residual, peak of the corrected rate \|∂_u m_c\| against the peak Misner–Sharp luminosity \|∂_u m\| on each sphere; with `--ladder` the resolution table, observed order and Richardson value (paper, M_ADM = 3.708: E_∞ = 0.532 ± 0.023 = (0.144 ± 0.006) M_ADM; E_ψ = 0.197, i.e. E_æ = 0.335 = 0.090 M_ADM; corrected rate more than an order of magnitude below the luminosity of fig:lum(a): printed ratio 25 at Θ = 200, 15 at 150, but only 8.1 / 5.3 at Θ = 100 / 75. Not in the paper: E_∞ = 0.6231/0.5443/0.5323 at dr = 0.0469/0.0234/0.0117, order 2.71, Richardson 0.5302; E_ψ 0.197–0.202 on Θ = 75–200; the runs before the scheme corrections of VALIDATION.md gave 0.6219/0.5435/0.5317, order 2.73, 14.3%) |
| `fig_probe_speed.py --pair BG,PROBE,c2[,label] ...` | `probe_speed.png` (fig:probe; single column: (v − v0)/v0 against t for all pairs) | background and probe runs (`--probe_amp=1e-6 --probe_r=150`, r_max = 350, Δt = 1) | per pair: median and MAD of (ṙ_pk + β)/N, (median − v0)/v0 and the largest \|v − v0\|/v0 after `--t_late` (40) (paper: 0.01 % for v0 = 1 and 0.02 % for v0 = 1.267, deviation below 3e-4 after t ≃ 40; printed medians 1.00012 and 1.26706, MAD 4e-5 / 8e-5, largest deviation for t ≥ 40 3.0e-4 / 1.9e-4) |

Typical calls (run directories as in README §5: basef = the v0 = 1 baseline run,
sl600_N12801 / sl600_N25601 / sl600_N51201 = the c2 = 0 runs at dr = 0.047 / 0.023 /
0.0117, sl600_N25601_rec / sl600_N12801_rec = the recorder runs, probe_*_bg / probe_*_pr =
the background / probe pairs):

    python3 fig_wout_background.py  basef            --madm 3.711
    python3 flux_identity_rec.py    sl600_N25601_rec --c2 0 --plot 40,100,200   # fig:fluxid
    python3 flux_identity_rec.py    sl600_N12801_rec --c2 0                     # dr = 0.047 residuals
    python3 fig_flux_identity.py    sl600_N25601     --c2 0
    python3 fig_invariant.py        sl600_N51201     --c2 0
    python3 fig_coulomb.py          sl600_N51201     --c2 0
    python3 fig_luminosity_clock.py sl600_N25601     --c2 0
    python3 fig_energy.py           sl600_N51201     --c2 0 \
                                    --ladder sl600_N12801,sl600_N25601,sl600_N51201
    python3 fig_probe_speed.py --pair probe_v1_bg,probe_v1_pr,-0.5119403 \
                               --pair probe_v1267_bg,probe_v1267_pr,0

On M_ADM: without `--madm` the scripts estimate it from m(Θ) at large radius on the first
dump, which gives 3.7084 for every A0 = 0.35 run (the initial data do not depend on c2);
the paper quotes M_ADM = 3.7111 at t = 28 on the finest v0 = 1 moving-mesh run (+0.07 %
drift since t = 2 at Θ = 250) and 3.708 for the c2 = 0 energy runs.

## Convergence, mass function, moving mesh, amplitude scan and white-hole scripts

All use the shared helpers of `spin0_common.py` (`load_profile`, `list_dumps`,
`load_dump`, `apparent_horizon`, `misner_sharp`, `theta_pm`, `beta_shift`, `fwhm`,
`peak_index`, …; `boundary_reflection.py` takes the couplings as `--c1/--c2`, checked
against the run's `err.txt` if present); run directories are positional
arguments except for `fig_white_hole.py`, `fig_sharp_features.py` and
`freeze_and_horizon.py` (options), and `boundary_reflection.py` takes the test name
first; `--out` (or `--outdir` for `fig_conv_scatter.py`, which writes several figures)
sets the figure path.

| script | figure | runs (README §5) | prints (paper values) |
|---|---|---|---|
| `fig_conv_scatter.py WORKDIR [--field_amps 0.01,0.06]` | `conv_scatter_fields_A001.png`, `conv_scatter_fields_A006.png` (fig:conv-fields-A001, fig:conv-fields-A006; one 2 × 3 figure per amplitude of `--field_amps`), `conv_scatter_constraint.png` (fig:conv-constraint), `conv_scatter_orders.png` (not in the paper) | scattering ladder, `WORKDIR/A<A0>_N<N>` | order table (tab:conv-scatter: 2.35 2.52 2.17 2.16 2.92 2.16 1.95 / 2.62 2.46 2.20 2.14 3.36 2.12 1.90 / 2.81 1.99 1.45 1.47 2.42 1.42 1.49); minimum central lapse (paper N_min ≃ 0.98, 0.70, 0.4; printed 0.979 / 0.703 / 0.383 over the dumps); constraint-norm ratios at t = 14 (paper, A0 = 0.06: 1.8 and 2.5; printed 1.79 / 2.46) |
| `fig_conv_collapse.py COARSE MID FINE --mesh_run mass_hi` | `conv_collapse.png` (fig:conv-collapse, the rescaled constraint norms; the orders are printed) | collapse ladder col_N12801/25601/51201, moving-mesh run mass_hi | p = (2.62, 2.55, 2.34, 2.38, 2.62, 2.42, 2.24) (tab:conv-scatter); constraint ratios 3.69–4.69 / 3.95–6.26 (tab:conv-ratio); Θ_AH(28) = 12.293/12.604/12.677, differences shrinking by 4.3, Richardson 12.70; Θ_AH ≃ 8.5, r_AH ≃ 9, Θ(r = 150) ≃ 137 at t = 2 and central lapse 7e-7, r_AH ≃ 65, Θ(r = 150) ≃ 88 at t = 28 (paper ≃ 90); M_ADM 3.70836 → 3.71107 at Θ = 250 on the moving mesh (+0.07 %); Θ = 250 off the fixed grid from the dump at t = 20 (Θ(r_max) = 249.8; paper "t ≃ 20"); at the fixed-grid edge Θ = 234.5 at t = 28 the fixed grids agree to 6.7e-5 and the fixed grid and moving mesh to 6.9e-6 (paper 1e-4, 2e-5); change of m between t = 2 and 28: 1.217 % at Θ = 100 … 0.073 % at 250, extrapolation A + B/Θ³ → −0.0102 % (power free: p = 3.12, +0.0037 %; paper: Θ^−3, 1.2 %, about 0.01 %) |
| `fig_mass_function.py RUN` | `mass_function.png` (fig:mass) | mass_hi (N = 51201, ALE, t = 28) | paper: M_hor = 6.35 (Richardson Θ_AH = 12.70), M_ADM = 3.7111, M_ADM − M_hor = −2.64, M_hor/M_ADM = 1.71, ρ_æ < 0 outside the horizon except across the outgoing scalar shell (15 ≲ Θ ≲ 30). Printed: Θ_AH = 12.677 on this grid, M_hor = 6.338, 3.7111, 1.708, −2.639 with the Richardson Θ_AH; ρ_æ > 0 only at Θ = 18.3–20.3, inside the scalar shell at Θ = 12.7–20.6 (the paper's upper bound 30 is loose). The public-code run before the moving-mesh corrector correction gave M_ADM = 3.7109 |
| `fig_moving_mesh.py FIXED ALE` | `moving_mesh.png` (fig:mesh) | mesh_ale0, mesh_ale1 | Θ_AH equal to < 1e-3 on the two grids (8.0e-4); max 12.60 at t ≈ 28 → 11.8 at t = 150 (paper 12.6 → 11.8); Θ(r_max) 284 → 124 on the fixed grid (not in the paper) |
| `scan_bh_threshold.py RUN ...` | — | astar_A0.05 … 0.15 (N = 6400) | per amplitude: min θ₊, final central lapse, Θ_AH at the end, times with θ₊ < 0 (0.05–0.06 disperse, 0.07–0.11 and 0.13 central throat only, 0.12 disperses, 0.14–0.15 genuine horizons) |
| `fig_white_hole.py --uniform D1,D2,D3 --refined D4,…,D7` | `white_hole.png` (fig:wh) | seven-grid ladder at t = 30 (panel (a): finest uniform run, `box40_N51200_r2`; its dumps are opened by the requested times `--times` and placed at their exact times) | spacing, FWHM and \|P\| of the shell per grid; per-family power laws (final rerun: width ∝ dr^0.84 uniform / dr^0.79 refined, \|P\| ∝ dr^−0.50 / dr^−0.47; paper: dr^0.8, dr^−0.5); the light ray behind the escaping pulse (paper 3.2 / 0.3 / 0.1 at t = 30 / 34 / 40; printed 3.16 / 0.30 / 0.09); the inner front of panel (a) at each dump: Θ, steepest slope, jump of P, lapse, dr/dt against −β − N (paper: jump ≃ 0.3, N ≃ 0.13, along ingoing light rays, stalling at Θ ≃ 3.2; printed jump 0.28–0.32 and N = 0.125–0.154 for t = 20–30, Θ = 3.54 → 3.15, dr/dt within 0.005 of −β − N up to t = 36) |
| `fig_sharp_features.py --uniform D1,D2,D3 --refined D4,D5 --escaping D6` | `sharp_features.png` (fig:sharp) | the three sharp features in P: overview at t = 40, escaping pulse t = 40–80 (zone A), inner shell at t = 30 (uniform + refined), inner front at t = 40 | peak of the escaping pulse per time (paper: Θ ≃ 21 → 62; printed 21.17 → 61.51); steepest \|dP/dΘ\| of the inner front on the uniform grids at t = 24 (4.9 / 6.2 / 6.9, changing little) and t = 40 (7.1 / 14.8 / 30.7, ×2.1 per halving of dr: paper "doubles") |
| `wh_formation_ladder.py D1 … D7 --series D7` | — | same ladder | per grid: spacing, width (points), \|P\|, ρ, …, \|ΘC\| (28 / 49 / 85 on the uniform grids, not in the paper); exponents over all grids and the three finest (paper: ρ ∝ dr^−1; printed dr^−0.99 / dr^−0.95); finest run t = 29 … 31 (paper: the width levels off once the refinement is centred on the shell) |
| `wh_front_convergence.py COARSE MID FINE` | — | whconv_N12800/25600/51200 (zone A, F = 8; the escaping-pulse runs of fig:sharp(b)) | the front tracked from t = 40 (largest local \|P\| maximum in a window around the position predicted from the previous dump): spacing, Θ_f, width, \|P\|, the anti-trapped layer on the front, min θ₋, max\|K\|, E per grid at t = 40, 45, …, 80, observed orders with the actual spacing ratios, Richardson values; \|ΘC\| at the front; the rate k of width ∝ e^(−kt) (paper: the pulse width shrinks roughly as e^(−0.15t); the layer narrows from about 0.9 at t = 40 to 0.06 at t = 80 on these grids, which agree on its width to within 20 %, while its depth min θ₋ grows with resolution. Printed: layer 0.89/0.90/0.84 at t = 40, 0.059/0.062/0.067 at t = 80, spread 7–18 %; min θ₋ −0.21/−0.25/−0.28 at t = 40, −0.85/−1.12/−1.36 at t = 80; k = 0.14/0.15/0.16 over t = 50–75, 0.12/0.13/0.14 over 50–80. Not in the paper: at dr = 0.0015/0.0007/0.0004 at the front Θ_f converges with p = 0.6–0.9 and the width with p = 0.7–1.1 for t = 60–80, the peak \|P\| with p ≈ 0.2–0.3 up to t = 75 (0.9 at t = 80); \|ΘC\| at the front 0.19/0.08/0.05 at t = 60 and 2.1/1.4/0.7 at t = 70, decreasing with resolution up to t = 70 and growing beyond) |
| `fig_bh_threshold.py bis_N* …` | `bh_threshold.png` (fig:threshold) | threshold bisection bis_N<N>_A<A0> (r_max = 300, t ≤ 60) | per grid: Θ_AH of the outermost trapped surface at t ≤ 20 (early horizon) and at t = 60; the early horizon appears from A0 = 0.135 / 0.130 / 0.125 / 0.125 at N = 12801 / 25601 / 51201 / 102401 (A* in 0.120–0.125 on the two finest grids); Θ_AH ≃ 1.25–1.44 at the switch, 1.4–1.9 at formation, 0.8–1.4 by t = 60 |
| `bisect_bh_threshold.py bis_N* … [--t_early 20]` | — | same runs | per run: largest detached (Θ > 0.5) Θ_AH and when, central lapse N₀ at the last dump searched and at the last dump of the run, verdict; per N the bracketing amplitudes (with `--t_early 20`: 0.134/0.135, 0.125/0.130, 0.120/0.125, 0.120/0.125; without it the late central trapped region makes the N ≥ 25601 brackets non-monotonic) |
| `freeze_and_horizon.py --freeze D:dr … --horizon D …` | — | freeze_N1601 … 38401_r150; wh_bh_N25600/51200 | central lapse history and last time N₀ > 0.1/0.01/0.001 (paper: the freezing time moves later on every finer grid, roughly as dr^−1/3, with no sign of converging; printed t ≃ 53, 63, 77, 101, 135 for dr = 0.0625, 0.031, 0.016, 0.0078, 0.0039, not quoted in the paper); Θ_AH(t) (not in the paper: 4.85–5.07 for t > 30) |
| `family_transitions.py RUN … [--extent]` | — | wh_bh/wh_wh runs, fam_* scans and ladders | dumps at which the trapped / anti-trapped classification changes (every dump with `--extent`), with Θ ranges (disjoint pieces as a comma list) and N₀ (paper, A0 = 0.10: θ₊ dips below zero at t ≃ 6–17, anti-trapped from t ≃ 20, split at t ≃ 36, inner part gone by t ≃ 44; printed: trapped at t = 6–16, anti-trapped 4.2–4.6 at t = 20, 2.7–17.6 at t = 36, split at t = 38, only the layer on the escaping pulse from t = 42; the σ = 3 central collapse — a trapped region a few cells across with N₀ → 0 — at A0 = 0.04 at t = 22 / 24 / 27 / 32.5 for dr = 0.0094 / 0.0031 / 0.0016 / 0.0008, and at A0 = 0.03 at t = 17 for dr = 0.0031, not up to t = 35 for dr = 0.0016 and again at t = 22 for dr = 0.0008, i.e. non-monotonic in the resolution) |
| `wh_family_ladder.py D1 D2 D3` | — | fam_s3_A0.06 and fam_r15_A0.020 ladders | per dump and grid: Θ, \|P\|, FWHM (points), ρ, \|ΘC\| of the inner shell, min N (paper: the inner shell steepens as for the baseline shell; public-code ladders, not quoted in the paper: when the shell forms, t = 20–22 for σ = 3 and t = 24–26 for r0 = 15, it spans 13/14/17 and 14/15/18 points, \|P\| grows by 1.3 per halving and \|ΘC\| by 1.5) |
| `boundary_reflection.py {spin0,scalar,late} RUN_R RUN_2R` | — | flat-space spin-0 pulse and weak scalar shell (boundary at R and 2R); reflbg_R150 vs basef | spin-0 (`--c1/--c2` = the couplings of the runs, default 0.7 / −0.5119403; checked against `RUN_R/err.txt` if present): reflected w_in / incident w_out at r = R − 10 (paper < 1e-5); scalar: reflected / incident along the ingoing ray (paper 0.03 %); late (not in the paper): relative difference of Θ∂_Θψ in Θ = 15–110, t = 236–300, 7.3 %, an upper bound that includes the different resolution of the two moving meshes |

Typical calls (run directories as in README §5; four runs of the final rerun lost dumps
when the disk filled and were repeated as `*_r2`, VALIDATION.md "Final rerun": the numbers
and fig:sharp and fig:wh come from `box40_N51200_r2` and `whconv_N51200_r2`; the paper's
fig:wh was drawn from the incomplete `box40_N51200`, whose panel (a) lacks four dumps at
t = 46–54, and with the one-decimal nominal times of its quarter-time dumps; the call
below uses the complete run and the exact dump times):

    python3 fig_conv_scatter.py   REPRO --outdir figs        # REPRO contains A0.01_N6401, ...
    python3 fig_conv_collapse.py  col_N12801 col_N25601 col_N51201 --mesh_run mass_hi
    python3 fig_mass_function.py  mass_hi
    python3 fig_moving_mesh.py    mesh_ale0 mesh_ale1
    python3 scan_bh_threshold.py  astar_A0.*
    python3 fig_white_hole.py --uniform box40_N12800,box40_N25600,box40_N51200_r2 \
                              --refined box40pp_N12800,box40pp_N25600,box40pp_N51200,box40pp_N51200_F16
    python3 wh_formation_ladder.py box40_N12800 box40_N25600 box40_N51200 box40pp_N12800 \
                              box40pp_N25600 box40pp_N51200 box40pp_N51200_F16 --series box40pp_N51200_F16
    python3 fig_sharp_features.py --uniform box40_N12800,box40_N25600,box40_N51200 \
                              --refined box40pp_N51200,box40pp_N51200_F16 --escaping whconv_N51200
    python3 wh_front_convergence.py whconv_N12800 whconv_N25600 whconv_N51200
    python3 freeze_and_horizon.py --freeze freeze_N1601:0.0625 freeze_N3201:0.031 freeze_N6401:0.016 \
                              freeze_N12801:0.0078 freeze_N38401_r150:0.0039 --horizon wh_bh_N25600 wh_bh_N51200
    python3 family_transitions.py wh_bh_N25600 wh_wh_N25600 fam_*_N6400 [--extent]
    python3 fig_bh_threshold.py   bis_N12801_A0.1* bis_N25601_A0.1* bis_N51201_A0.1* bis_N102401_A0.1*
    python3 bisect_bh_threshold.py bis_N12801_A0.1* bis_N25601_A0.1* bis_N51201_A0.1* bis_N102401_A0.1* --t_early 20
    python3 wh_family_ladder.py   fam_s3_A0.06_N19200 fam_s3_A0.06_N38400 fam_s3_A0.06_N76800
    python3 boundary_reflection.py spin0  RUN_R100 RUN_R200 --R 100 --c2 0   # v0 = 1.267 runs
    python3 boundary_reflection.py scalar screfl_R150 screfl_R300
    python3 boundary_reflection.py late   reflbg_R150 basef --theta_window 15,110

## Notes

* **Flux identity and the dump cadence.** `fig_flux_identity.py` needs time derivatives
  on both sides (−dm/dt|_Θ, and d_t ln N, d_t β, d_t Q inside T^æ_ss). From Δt = 0.5
  dumps the sharp burst is sampled by only a few points, and the residual at the burst
  falls with the order of the time stencil (it is a sampling error, not a formula
  error); the default is the 9-dump (eighth-order) stencil. The paper's number (0.02 % of
  the peak at dr = 0.023, converging to zero with resolution) comes from the per-step
  records: `collapse --rec_spheres=…` writes them (`spheres.bin`, README §4) and
  `flux_identity_rec.py` evaluates the identity on them (time differences across steps,
  neighbouring records interpolated in r onto the central record's stencil, so the moving
  mesh is handled). See VALIDATION.md, "Sphere recorder".
* **E_ψ.** `fig_energy.py` prints the scalar energy that has crossed each sphere,
  ∫F_ψ dt; it is part of E = M_ADM − m_c (0.197 of the 0.532 radiated, i.e. the scalar
  carries 5.3 % and the æther 9.0 % of M_ADM; the scalar contribution is not negligible).
* **Probe runs.** The probe is tracked from the first dump along the outgoing
  characteristic; samples before `--tstart` (default 15), near the outer boundary
  (`--edge`) or with the maximum at the edge of the search window are dropped. The
  tracker was tested on a synthetic pair (static metric with a shift, moving non-uniform
  grid, a 1/r pulse following ṙ = −β + N v0): it returns v0 to 1e-6.
