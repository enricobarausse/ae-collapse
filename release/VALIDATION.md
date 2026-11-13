# Validation: the public collapse.cpp against the development code

This file records, in chronological order, how the public source was checked against the
development version it was extracted from and what was changed afterwards. At extraction
the public code was bit-identical to the development code for time-symmetric data (every
case of the first sections). Then the ingoing initial data (`--ingoing=1`) were corrected
in two steps (sections "Initial-data fix" and "Ingoing initial data, stage evaluation"),
and after the second review six points of the scheme were corrected for every kind of data
(section "Scheme corrections"); every production run of the paper was redone with the
corrected code (section "Full rerun"). Finally the outer boundary condition for the incoming
spin-0 field was made exact and every run was redone once more with the code shipped here
(last section, "Final rerun"). The development source, the driver scripts and the run
directories named here are not part of the release; the numbers are those obtained at
validation time.

Reference: the development `collapse.cpp` (not distributed; md5 038e70a67fe65d9645f32c578877102c).
Release: `release/collapse.cpp` (md5 762b40e901f66ea2854596e183b27ece at validation time;
the sphere recorder, the initial-data iteration and the header were added afterwards, see
below).
Compiler: g++ 8.5.0 (Red Hat 8.5.0-22). Date: 2026-10-09.

## Build

```
g++ -O2 -std=c++11 -o bin/orig  collapse_dev.cpp           # reference (development source)
g++ -O2 -std=c++11 -o bin/relP  release/collapse.cpp       # release, same flags
g++ -O2 -std=c++11 -Wall -Wextra -o bin/rel release/collapse.cpp   # release with warnings
```

`-Wall -Wextra` gives **no warnings** for the release code (`make warn`).

## Cases

Each case runs in its own directory. The reference and `rel`/`relP` use the original
positional command line (the release binary accepts it); `relN` runs the release binary
with the equivalent named options (check of the new CLI).

| case | what it exercises | positional command (reference) |
|---|---|---|
| a | fixed-grid scattering, + auxiliary dumps (`COLLAPSE_DUMP_EVERY=500`) | `3201 0.04 14 "8,10,12,14" 150 0.7 -0.5119403 0 0` |
| b | ALE collapse to a black hole, horizon finder | `3200 0.35 30 "5,10,15,20,25,30" 300` |
| c | white-hole ID (c1=0.8, ingoing), zone B pre-positioned (F=8, w=1.5, x0=9.9, on at 20, tracking from 30.5) | `3200 0.10 35 "20,25,29,29.5,30,31,33,35" 40 0.8 -0.7111111 1 1 0 0 50 5 8 20 1.5 0 0 9.9 30.5` |
| d | zone A tracking the outgoing front (F=4, t_on=40, w=3), acquires at t=40.0 | `6400 0.10 50 "30,40,42,45,48,50" 150 0.8 -0.7111111 1 1 0 4 40 3` |
| e | both zones: B acquires by search (t=31.19), A at t=40.0 | `3200 0.10 45 "25,30,35,40,45" 40 0.8 -0.7111111 1 1 0 4 40 3 4 20 1.5` |
| f | superluminal couplings c2=0 (v0=1.267, spin-0 projector with v0≠1), ALE | `3200 0.35 20 "10,20" 300 0.7 0` |
| g | r0/sigma arguments (sigma=3), white-hole ID | `3200 0.06 30 "10,20,30" 60 0.8 -0.7111111 1 1 0 0 50 5 0 20 2 0 0 0 0 10 3` |

Named form of the same cases (`relN`), e.g. case c:
`--N=3200 --A0=0.10 --t_final=35 --times=20,25,29,29.5,30,31,33,35 --rmax=40 --c1=0.8
--c2=-0.7111111 --ingoing=1 --refB_F=8 --refB_ton=20 --refB_w=1.5 --refB_x0=9.9
--refB_ttrack=30.5`. (Driver scripts `cases.sh`, `cases_named.sh`, `compare.sh`, not
included.)

## Comparison

`compare.sh` checks that the reference and test directories contain the same set of output
files and runs `cmp` on every one of them: all `prof_t*.csv`, `phi.dat` (every step),
`psi_times.dat` and, for case a, the 54 auxiliary dumps (`psi_dump_*`, `*_dump_*`).

```
### orig vs rel        (release, -O2 -Wall -Wextra)
case a: 66 files compared (4 profiles, t_end=1.40015624999993040e+01): 0 differ
case b: 8 files compared (6 profiles, t_end=3.00000000000013678e+01): 0 differ
case c: 10 files compared (8 profiles, t_end=3.50000244812937140e+01): 0 differ
case d: 8 files compared (6 profiles, t_end=5.00004239370112060e+01): 0 differ
case e: 7 files compared (5 profiles, t_end=4.50000567025613023e+01): 0 differ
case f: 4 files compared (2 profiles, t_end=2.00031259768684784e+01): 0 differ
case g: 5 files compared (3 profiles, t_end=3.00018755861122131e+01): 0 differ
### orig vs relN       (release, named options)
  identical result: every file of every case identical (0 differ)
### orig vs relP       (release, exactly -O2 -std=c++11)
  identical result: every file of every case identical (0 differ)
```

**All 108 output files of each release run are byte-identical to the reference.**
The stderr event lines (outer-BC gate switches, refinement pre-positioning / acquisition /
tracking with times and positions) are identical in cases c, d, e, and so are the
apparent-horizon lines on stdout in case b.

## Other checks

* Probe (new, off by default so not covered above): `--N=3200 --rmax=350 --t_final=3
  --probe_amp=1e-6` vs the same without the probe: max |ΔK| = 5.8e-7 at r = 149.1, and
  < 4.2e-12 for |r − 150| > 30, i.e. a linear, localised perturbation.
* Removed features are refused rather than silently ignored: positional slot 10 = 3 →
  `FATAL: the inner taper (slot 10, ale_inner) was removed`; slot 17 = 16 →
  `FATAL: the adaptive zone-B factor ... was removed`; an unknown named option →
  `FATAL: unknown option`.

## Not covered

* The removed features (inner taper, `qform`, `refB_Fmax`) by construction.
* Long runs and large N: the cases above go through every code path the paper uses
  (fixed grid, ALE, gate, both refinement zones including acquisition and tracking,
  pre-positioning, v0 ≠ 1, ingoing ID, r0/sigma); the arithmetic is the same at any N.

## Initial-data fix for the ingoing shell (2026-10-09, after the equation audit)

A symbolic (SymPy) audit of the equations found that `set_initial_data`
solved Theta from the Hamiltonian constraint with Q = 0 and computed Q only afterwards.
For the ingoing shell P != 0 sources Q through the momentum constraint, and Q^2 enters the
Hamiltonian constraint, so the data violated the constraint at t = 0 by an amount that did
not shrink with resolution.  The two solves are now iterated to a fixed point (6 iterations,
|dTheta| < 1e-13) when `--ingoing != 0`; time-symmetric data (Q = 0 identically) skip the
loop and remain bit-identical to the development code.

Check (c1=0.8, c2=-0.7111111, ingoing, A0=0.10, rmax=40, dump at t=0.01), max |Theta C| for r<35:

| N    | before   | after    |
|------|----------|----------|
| 1600 | 2.93e-2  | 1.17e-2  |
| 3200 | 3.01e-2  | 5.92e-3  |
| 6400 | 3.08e-2  | 2.97e-3  |

`cmp` of the time-symmetric run (N=1601, A0=0.35, t=6) before/after: identical.
Consequence: the white-hole results of the paper were recomputed with the corrected data
(the numbers in README §5 — Θ_AH = 4.85–5.07 at A0 = 0.12, freeze times 55/63/77/101/135 —
are the recomputed ones, 53 instead of 55 for dr = 0.0625 after the final rerun at the end
of this file; the earlier 4.678/4.676 and 60/65/78/102 came from the uncorrected data).

## Sphere recorder (`--rec_spheres`, `--rec_every`; 2026-10-10)

Purpose: the flux identity F_æ + F_ψ = −dm/dt|_Θ (paper Fig. flux_identity) needs time
derivatives on both sides; from the Δt = 0.5 profile dumps the burst is undersampled
(0.6–3 % of the peak at dr = 0.023, whereas fine-step records of the development code gave
0.3 %). The recorder writes, every k steps, 13-point
stencils of r, ψ, P, K, a_r, Θ, N, Q, β around each sphere (README §4);
`analysis/flux_identity_rec.py` evaluates the identity on them.

**Bit identity.** Built with `-O2 -std=c++11` (`make warn` clean) and compared with `cmp`
against the binary built from the previous release source, every output
file including stdout/stderr, with the recorder OFF and ON (`--rec_spheres=20,30,100
--rec_every=3`):

| case | options | files | differ (off) | differ (on) |
|---|---|---|---|---|
| time-symmetric, ALE | `--N=1601 --A0=0.35 --t_final=8 --times=2,4,8` | 7 | 0 | 0 |
| ingoing, v0 ≠ 1 | `--N=1600 --A0=0.10 --t_final=10 --times=5,10 --rmax=40 --c1=0.8 --c2=-0.7111111 --ingoing=1` (Θ_s = 100 > Θ(r_max): i = −1 path) | 6 | 0 | 0 |
| zone-A refinement, c2 = 0 | `--N=1601 --A0=0.35 --c2=0 --t_final=6 --times=3,6 --ref_F=4 --ref_ton=2 --ref_w=5` | 6 | 0 | 0 |

The recorded r stencils coincide with the profile r column at the same grid indices. The refactor of
`fig_flux_identity.py` (computation moved into `dump_series()`, reused by the new script)
gives identical stdout and an identical PNG on a test run.

**Validation runs** (c1 = 0.7, c2 = 0 i.e. v0 = 1.267, A0 = 0.35, time-symmetric, dumps
every 0.5, `--rec_every=1`, order-4 time stencils across records; residual =
|F_æ + F_ψ + dm/dt| / max|dm/dt| on the sphere; "dumps" = `fig_flux_identity.py`, order 8,
on the profile dumps of the SAME run over the same window):

| run | dr | Θ_s | peak \|dm/dt\| | recorder max / rms | dumps max / rms |
|---|---|---|---|---|---|
| N=1601, r_max=300, t≤160 | 0.1875 | 40 | 2.07 | 2.75 % / 0.33 % | 2.97 % / 0.34 % |
| | | 70 | 3.85 | 1.02 % / 0.16 % | 1.06 % / 0.16 % |
| | | 100 | 5.65 | 1.16 % / 0.19 % | 1.34 % / 0.19 % |
| N=3201, r_max=300, t≤160 | 0.094 | 40 | 0.86 | 0.31 % / 0.17 % | 0.30 % / 0.17 % |
| | | 70 | 1.43 | 0.13 % / 0.06 % | 0.18 % / 0.06 % |
| | | 100 | 2.01 | 0.06 % / 0.03 % | 0.75 % / 0.10 % |
| N=1601, r_max=75, t≤69 | 0.047 | 40 | 0.71 | 0.20 % / 0.07 % | 0.23 % / 0.08 % |
| | | 60 | 1.01 | 0.04 % / 0.03 % | 0.07 % / 0.03 % |
| N=3201, r_max=75, t≤69 | **0.023** | 40 | 0.84 | **0.06 %** / 0.014 % | **1.52 %** / 0.25 % |
| | | 60 | 1.22 | **0.02 %** / 0.005 % | **0.35 %** / 0.06 % |

The r_max = 75 runs have the paper's dr = 0.023 (and dr = 0.047) at N ≤ 3201; the box is
small (outer boundary at Θ ≈ 75), which changes the solution but not the local identity;
Θ = 40, 60 lie in the ALE taper (x > 0.4 r_max = 30), so they also test the moving-mesh
handling (neighbouring records interpolated in r onto the central stencil).

Findings:
* At dr = 0.023 the per-step records give 0.02–0.06 % (max) where the dumps of the same
  run give 0.35–1.5 %: the dump-based excess is time sampling of the burst, and per-step
  sampling brings the identity well below the 0.3 % of the development-code records. At dr ≥ 0.094 the recorder
  and the dumps agree and the residual is a smooth, slowly varying spatial-truncation term
  rather than sampling (its absolute size falls with dr, but this was not traced further).
* The dm/dt side alone is not the issue: m(Θ_s, t) resampled at the dump times with the
  order-8 stencil gives nearly the recorder value (0.07 % vs 0.06 % at Θ = 40, dr = 0.023);
  the sampling error of the dumps comes from the fixed-r time derivatives of ln N, β, Q
  inside T^æ_ss (W, B, Q_t), which change on the burst time scale.
* Order 2 vs 4 across records and `--rec_every=5` (emulated by keeping every 5th record)
  change the residual by < 0.05 % of the peak at dr = 0.023: `--rec_every=5` is enough.
* Caveat: dr ≥ 0.094 (N ≤ 3201 at r_max = 300) is pre-asymptotic for this ID — even
  m(Θ = 20, t = 9.5) is 3.88 / 4.11 / 4.23 / 4.26 / 4.26 for dr = 0.19 / 0.094 / 0.047 /
  0.023 / 0.0117, and the burst peaks differ by 2.4× between dr = 0.19 and 0.094 — so only
  the dr ≤ 0.047 rows say anything quantitative about the paper's run.
* Side observation (v0 = 1 couplings, t < 10, Θ = 20–40): fluxes are small (|dm/dt| ~
  0.01) there and the residual relative to that local peak is O(10–100 %) at dr ≥ 0.094;
  this was not investigated further (the dump analysis would see the same).

**Production runs** (README §5, r_max = 600, t ≤ 600, `--rec_every=5`, spheres Θ = 40 … 200;
the paper quotes the bound "within 0.02 % of the peak flux on every sphere"): max residual / peak 0.020 % at Θ = 40 falling to
0.002 % at Θ = 200 for dr = 0.023, and 0.098 % / 0.006 % for dr = 0.047 — a factor 3–5 per
halving of dr, i.e. the residual converges to zero with resolution. The radiated energy from
the same three profile-dump runs (`analysis/fig_energy.py`, M_ADM = 3.708): E_∞ = 0.6219 /
0.5435 / 0.5317, observed order 2.73, Richardson 0.530, E_∞ = 0.532 ± 0.023; scalar part
E_ψ = 0.197. (These are the values before the scheme corrections below; the rerun with the
corrected scheme gives 0.6231 / 0.5443 / 0.5323, order 2.71, Richardson 0.5302, the same
E_∞ = 0.532 ± 0.023 and E_ψ = 0.197, and the same flux-identity residuals — see "Full
rerun" at the end of this file.)

## Ingoing initial data, stage evaluation (2026-10-10)

After the iteration above, the RK4 solve of the Hamiltonian constraint still used P, Q
(and a_r, K) from the grid point i at the stages r_i + dx/2 and r_i + dx, which left the
ingoing data first order in dx (the time-symmetric data have P = Q = 0 there and were
unaffected). The public code now evaluates P = ψ' + ψ/r analytically at the stage radii
and interpolates Q, a_r, K linearly to them (`integrate_HamiltonianConstraint_initial`,
stage weight w = 0, 1/2, 1). Check: `--A0=0.10 --rmax=40 --c1=0.8 --c2=-0.7111111
--ingoing=1 --times=0.01 --t_final=0.02`, max |Θ C| for r < 35 at the first dump
(t = 0.01), before/after this change (the "before" column is the "after" column of the
previous table):

| N    | before (iterated, grid-point stages) | after (stage evaluation) | ratio |
|------|------|------|------|
| 1600 | 1.17e-2 | 6.48e-4 | — |
| 3200 | 5.92e-3 | 1.57e-4 | 4.13 |
| 6400 | 2.97e-3 | 3.86e-5 | 4.07 |

i.e. second order, like the time-symmetric data (3.27e-4 / 8.10e-5 / 2.02e-5 on the same
grids, unchanged to the bit). The Θ/Q iteration still converges in 6 iterations.

Size of the change on an ingoing run (`--N=1600 --A0=0.10 --rmax=40 --c1=0.8
--c2=-0.7111111 --ingoing=1`, max |new − old| over the grid relative to max |old|):
at t = 5: ψ 2.8e-5, P 1.2e-4, K 1.2e-4, a_r 9.0e-5, Θ 2.8e-5, N 1.6e-5, Q 1.2e-4; at
t = 10: ψ 5.0e-5, P 4.0e-4, K 2.4e-4, a_r 1.9e-4, Θ 2.8e-5, N 3.3e-5, Q 3.1e-4. The
evolved constraint residual at t = 10 (max |Θ C|, r < 35: 7.84e-2) is unchanged to three
digits, i.e. it is dominated by the evolution, not by the initial data.

Bit identity after this change, binary built from the previous source vs the new one,
`cmp` on every output file plus stdout/stderr: time-symmetric ALE (N = 1601, t = 8), zone-A
refinement with c2 = 0, the legacy positional case a with auxiliary dumps (68 files), the
probe case, and the three time-symmetric constraint runs above: **0 files differ**. The
ingoing cases (N = 1600 r_max = 40 to t = 10 with the recorder; the zone-B pre-positioned
case c; the three ingoing constraint runs) differ, as they must.

## Header, comments and packaging (2026-10-10)

The citation/licence/build lines of the file header, the `.gitignore` and the
documentation edits of the final polish touch no code: the binary built from the polished
source is byte-identical (same md5) to the one built from the previous source, and the
four quick cases `--N=1601 --A0=0.35 --t_final=8 --times=2,4,8` (time-symmetric, ALE),
`--N=1600 --A0=0.10 --t_final=10 --times=5,10 --rmax=40 --c1=0.8 --c2=-0.7111111 --ingoing=1
--rec_spheres=20,30,100 --rec_every=3` (ingoing, recorder), `--N=1601 --A0=0.35 --c2=0
--t_final=6 --times=3,6 --ref_F=4 --ref_ton=2 --ref_w=5` (zone A, v0 ≠ 1) and the legacy
positional case a with `COLLAPSE_DUMP_EVERY=500` give identical output files (89 in all),
stdout and stderr with both binaries.

## Second review pass (2026-10-12)

Code changes of the second readability pass and their check. Changes to the arithmetic:
none. Changes to the behaviour: (i) `phi.dat` and `psi_times.dat` are now static streams
so that a FATAL stop (`exit()`, which does not destroy automatic objects) flushes them —
before, a run stopped by a guard lost the buffered tail of `phi.dat` (checked with a
deliberate FATAL at step 401, N = 1601: 325 of 403 lines were written, `psi_times.dat` was
empty; now both are complete and the old file is a byte-identical prefix of the new one);
(ii) non-numeric entries of `--times` / `--rec_spheres` are rejected instead of being read
as 0 by `atof`. Readability: the three copies of the shell formulas and the two copies of
the quintic taper were merged into one helper each (same floating-point expressions).

Bit identity, binary built from the previous source vs the new one (`-O2 -std=c++11`,
`make warn` clean), md5 of every output file (profiles, auxiliary dumps, `phi.dat`,
`psi_times.dat`, `spheres.bin/.txt`, stdout, stderr), 89 files in six cases: **0 differ**.
Cases: time-symmetric fixed grid with auxiliary dumps (`--N=3201 --A0=0.35 --t_final=6
--times=2,4,6 --rmax=300 --ale=0 --dump_every=200`); ingoing white-hole data with zone B
pre-positioned, held, tracking and released, recorder on including the i = −1 path
(`--N=1601 --ingoing=1 --c1=0.8 --c2=-0.7111111 --A0=0.10 --rmax=40 --t_final=4
--times=1,2,3,4 --refB_F=8 --refB_ton=1 --refB_w=1.5 --refB_x0=9.9 --refB_ttrack=2.5
--rec_spheres=5,10,20,100 --rec_every=3`); moving mesh with recorder (`--N=3201 --A0=0.35
--t_final=6 --times=2,4,6 --rmax=300 --ale=1 --rec_spheres=20,40 --rec_every=5`); zone A
with c2 = 0 (`--N=1601 --A0=0.35 --c2=0 --t_final=6 --times=3,6 --rmax=300 --ref_F=4
--ref_ton=2 --ref_w=5`); both zones by search in the legacy positional form with sigma = 3
(`1601 0.06 4 "2,4" 60 0.8 -0.7111111 1 1 0 4 3 3 4 1 1.5 0 0 0 0 10 3`); the probe
(`--N=1601 --A0=0.35 --rmax=350 --t_final=3 --times=3 --probe_amp=1e-6 --probe_r=150`).
The 21 analysis scripts were run on the paper's run directories before and after the pass:
every printed number and all 16 figure files are identical (the one removed line is an empty
header that `wh_front_convergence.py` printed without content).

## Scheme corrections after the second review (2026-10-12)

Six points left open by the second review were changed, by decision of the authors, and
every production run of the paper was redone with the corrected code (the campaign of
Sec. 5 of README.md, 113 runs; the paper's figures and numbers were regenerated from it).

1. Moving mesh: the ICN corrector is evaluated on the mesh at t+dt with or without
   refinement (before, only with refinement; without it the corrector used the mesh at t,
   a first-order-in-dt inconsistency of the ALE scheme).
2. Lapse: (ln N)' = a_r is integrated outward and shifted so that N(r_max) = 1, the exact
   solution of the linear boundary problem, instead of a Newton shooting with a
   finite-difference derivative (round-off only); a non-finite lapse is a FATAL stop
   instead of being replaced by N = 1.
3. Horizon diagnostic: the bisection of the linear interpolant is replaced by its zero,
   printed together with the index of the containing cell and the time of the state
   (end of step) instead of the start-of-step time. stdout only.
4. The non-finite check runs on the state after every step (the last step included) and
   once on the initial data; `FIRST NaN` names the step that produced the value.
5. Zone A acquires the front only when the |P| peak exceeds 1e-6 (a run with `ref_ton`
   before the front's arrival used to lock onto a 1e-13 noise peak).
6. Scalar outer boundary: the condition d_t[P + ψ' + cψ] = 0 is applied in full; the term
   ψ d_t c (c = Θ'/Θ + Q/2 − K/3) was dropped before. d_t Θ' and d_t Θ come from the
   current right-hand sides, d_t K from the gated trK right-hand side and d_t Q from a
   backward difference of the recovered Q: over the previous step in the predictor
   stage, between the ICN stages in the corrector.

Effect on the six cases of the second review (same commands, previous binary vs
corrected binary, maximum relative difference of each field at the last dump): fixed
grid, refined grids and the probe case differ at round-off (1e-15 to 1e-12); the moving-mesh
case (`--ale=1`, t = 6) differs by 3e-7 in Θ and 6e-5 in ΘC (point 1); the horizon lines
agree in number and index and differ in the printed radius and time as intended; no FATAL
stop in either binary. Convergence orders of a coarse test ladder (A0 = 0.04, N = 1601 /
3201 / 6401, r_max = 150, fixed and moving mesh, t = 10 and 20) are the same to two digits
with the previous and the corrected scheme. The effect of point 6 appears only once the
scalar pulse reaches the outer boundary (t ≳ r_max), i.e. in the long energy and
luminosity runs; it is covered by the redone production runs and their convergence ladders.

## Effect of the initial-data stage fix on the ingoing production runs (2026-10-11)

Thirteen ingoing runs of the paper were repeated with the fixed initial data (P evaluated at the
RK4 stage radii, Q/a_r/K interpolated) and compared with the originals, which used
the grid-point values at the stages: the white-hole formation grids `box40_N12800` and
`box40pp_N12800` (t = 30 table: identical to three digits, e.g. |P| 0.253, width 13 points,
rho 0.0655 → 0.0656), the inner-shell ladders `fam_s3_A0.06_N19200` and `fam_r15_A0.020_N19200`
(peak position, height, width and constraint at t = 30 and 36 identical), the family runs
`fam_s3_A0.04_N19200`, `fam_r15_A0.028_N6400`, `fam_s3_A0.08_N6400` (every trapped/anti-trapped
transition time and extent identical), the frozen-centre runs `freeze_N1601/3201/6401`
(freeze times identical), `wh_wh_N25600` and `wh_bh_N25600` (identical timelines; central lapse at
t = 88 differs in the second digit of 1e-19), and the front run `whconv_N12800`. The finest grids
of the white-hole ladders, which were started before the fix, were at that point quoted as run;
they were redone in the full rerun below.

## Full rerun with the corrected scheme (2026-10-12/13)

Every recipe of README.md §5 was run again with the release source at the state of the
scheme corrections above (the `collapse.cpp` shipped here): 113 runs, started on
2026-10-12 at 12:55, up to 80 at once under `nice` on one 112-core machine, finished on
2026-10-13 at 10:14 (about 310 core-hours of measured wall time; the wall times now
quoted in README §5 are the ones measured in this campaign). The 21 analysis scripts were run on the new
output and every printed number compared with the paper (the scripts print the paper's
value next to their own); the 14 figures were regenerated and replaced in the paper.

Reproduced to the quoted digits:

* scattering ladder (9 runs): all 21 convergence orders of Table conv-scatter;
* amplitude scan at N = 6400 (11 runs): the same classification (0.05–0.06 disperse,
  0.07–0.11 and 0.13 central throat, 0.14–0.15 detached horizon; 0.12 disperses, as before);
* black-hole threshold (28 runs on four grids): early horizon from A0 = 0.135 / 0.130 /
  0.125 / 0.125, the same Θ_AH to three digits, `bh_threshold.png` pixel-identical;
* collapse ladder, mass function, moving mesh (6 runs): p = (2.62, 2.55, 2.34, 2.38, 2.62,
  2.42, 2.24), constraint ratios 3.69–4.69 / 3.95–6.26, Θ_AH(28) = 12.293 / 12.604 /
  12.677, M_ADM drift 3.70836 → 3.71171; Θ_AH equal to 8.0e-4 on the fixed and moving
  grids, Θ(r_max) 284 → 124; M_hor = 6.338, M_hor/M_ADM = 1.708. The one digit that moved:
  M_ADM(Θ = 250) of the moving-mesh t = 28 run is 3.7111 (3.7109 before the corrections;
  the moving-mesh corrector change, point 1 above), still 3.711 as the paper rounds it in
  the baseline recipe;
* spin-0 (10 runs): baseline pulse +0.0606 / −0.0718 and drift 10.7 / 11.3 in u;
  flux-identity residual 0.020 % → 0.002 % at dr = 0.023 and 0.098 % / 0.006 % at
  dr = 0.047; clock speed 1.2405 (1.241); Coulomb field 19.75 → 40.5, late 4.4–5.0, peak
  m − M_ADM 13.8, m_c 0.62; invariant waveform identical to three digits; probe medians
  1.00012 and 1.26706;
* white holes (49 runs): Θ_AH of the A0 = 0.12 black hole 4.846–5.068 on both grids; the
  A0 = 0.10 timeline (anti-trapped 2.7–17.6 at t = 36, 4.1–22.2 at t = 40, bulk gone by
  t = 44); the seven-grid inner-shell ladder at t = 30 (|P| = 0.253 / 0.355 / 0.509 on the
  uniform grids, |ΘC| = 28 / 49 / 85, exponents width ∝ dx^0.75, |P| ∝ dx^−0.5 per family;
the final rerun at the end of this file fits width ∝ dr^0.79 / dr^0.84);
  the front ladder (Θ_f, width, |P|, |ΘC| at t = 45–80 identical to three digits); the
  freezing times 55 / 63 / 77 / 101 / 135; every trapped / anti-trapped transition of the
  σ = 3 and r0 = 15 scans and ladders, including the A0 = 0.04 central collapse at
  t = 22 / 24 / 27 / 32.5.

Changed:

* the radiated-energy ladder of the c2 = 0 runs (r_max = 600, t ≤ 600): E_∞ = 0.6231 /
  0.5443 / 0.5323 on dr = 0.0469 / 0.0234 / 0.0117 (was 0.6219 / 0.5435 / 0.5317), observed
  order 2.71 (2.73), Richardson 0.5302 (0.530), E_∞ = 0.532 ± 0.023 = (14.4 ± 0.6) % of
  M_ADM = 3.708 (0.5323/3.708 = 14.36 %, was quoted as 14.3 %); E_ψ = 0.197 and
  E_æ = 0.335 unchanged. This is point 6 above (the ψ d_t c term of the scalar outer
  boundary condition), which acts once the pulse reaches r_max = 600, i.e. only in these
  long runs; the flux identity, which is local, is unaffected.
* one new run, σ = 3, A0 = 0.03, N = 76800 (dr = 0.0008, t ≤ 35, the finest grid of the
  A0 = 0.03 series, same command as the A0 = 0.04 run): the cell-sized central collapse is
  present again, at t ≈ 22 (it is present at dr = 0.0031 at t ≈ 17 and absent at
  dr = 0.0016 up to t = 35). The paper now states this non-monotonic behaviour.

Seven runs stop before their t_final on the documented guards, at the same step as in the
earlier campaigns and after the phases the paper uses: `bis_N25601_A0.135` at t = 32
(non-finite lapse; early horizon at t = 16), `box40_N12800` / `N25600` / `N51200` at
t = 55.3 / 55.6 / 55.9 (edge guard `dTheta/dr at rmax` once the front reaches r_max = 40;
profiles to t = 54 complete), `fam_r15_A0.028_N21334` at t = 92 (non-finite lapse; its
transitions are at t ≤ 60), `fam_s3_A0.03_N6400` at t = 15.6 (non-finite lapse; last dump
t = 14) and `fam_s3_A0.07_N42667` at t = 69.4 (fold of the moving mesh; anti-trapped phase
t = 16–34). Four of them used to end on a NaN one step later; with point 2 above the
non-finite lapse is reported first.

glibc malloc. The code allocates large temporaries every step, and with the default
glibc settings they are returned to the kernel each time (19–31 % of the CPU time in the
kernel on the authors' machine). With `MALLOC_MMAP_THRESHOLD_=4000000000
MALLOC_TRIM_THRESHOLD_=4000000000 MALLOC_TOP_PAD_=268435456` in the environment the
output is bit-identical and the run is about 18 % faster (393 s → 322 s for an
N = 25600 white-hole test to t = 1.5). The ten longest runs of the campaign
(`sl600_N51201`, `box40_N25600`, `box40_N51200`, `box40pp_N51200`, `box40pp_N51200_F16`,
`fam_s3_A0.06_N76800`, `fam_r15_A0.020_N76800`, `fam_s3_A0.03_N76800`,
`fam_s3_A0.04_N76800`, `fam_r15_A0.028_N42667`) were restarted with these settings on
2026-10-12 at 16:02 (at most 20 at once); README §5 marks them (m).

## Third review pass: predictor-stage d_t Q (2026-10-13)

A last read of the code before the release found that the predictor stage of point 6
above always used d_t Q = 0: the "previous" edge value of Q was recorded at the end of
the step, which is the value the next step starts from, so the backward difference
vanished and the ψ d_t c term entered the trapezoidal ICN average at about half weight.
The fix stores the start-of-step edge value instead (one line in `ICN_step`). The
output is bit-identical until the scalar pulse reaches r_max (three short cases: fixed
grid, zone B refinement, recorder + probe); afterwards ψ near the edge changes by about
1e-3 relative, so only the long energy runs (r_max = 600) can be affected. Those five
runs, the baseline spin-0 run and the fourteen dr = 0.0094 family scans were redone
with the fixed code: every quoted number and figure came out unchanged to the quoted
digits (the flux-identity figure differs in 11 pixels). In no other production run does
the scalar pulse reach r_max within the run's time span, so their output is bit-identical
to that of the fixed code.


## Final rerun with the shipped code (2026-11-04)

Code change since the full rerun above: the outer boundary condition for the incoming
spin-0 field. The incoming characteristic field w_in = (bN/v0) K + a_r is evolved as for
the outgoing linear spin-0 wave, d_t w_in = β̃ ∂_r w_in − N bN K/Θ (README §2), instead of
being held fixed, d_t w_in = 0; the outgoing field keeps its interior equation. Checked at
the commit: the flat-space spin-0 reflection drops from 3.5 % (r_max = 100) and 2.3 %
(r_max = 150) to ≤ 1e-6 for v0 = 1 and 1.267, and before the wave reaches r_max the output
changes by ≤ 2e-11 on the fixed grid and ≤ 2e-6 at the edge of the r_max = 40 white-hole
runs. The campaign binary (g++ 8.5.0, `-O2 -std=c++11`) is byte-identical to the one rebuilt
from that source (md5 e718004465beb8017bdf803dcbb2504e); the `collapse.cpp` shipped here
differs from it only by the readability pass below, which gives bit-identical output.

Runs: every recipe of README §5 again (113 runs, the command lines of the recipes; each
run directory keeps its `cmdline.txt`, checked against them) plus two scalar-reflection
runs (fixed grid, A0 = 1e-4, dr = 0.0117, r_max = 150 and a reference box r_max = 300, t ≤
175), 115 runs, started on 2026-10-29 at 19:48, up to 40 at once under `nice` on one
112-core machine, all with the malloc settings above, the last finished on 2026-10-31.
The disk of the run directories filled on 2026-10-30 while four of them were still
running, and these lost some profile dumps: `box40_N51200` (4 of 26), `fam_r15_A0.020_N76800`
(5 of 26), `fam_s3_A0.03_N76800` (8 of 55) and `whconv_N51200` (7 of 81). They were run
again from the start with the same command lines (`*_r2`, 2026-10-30/31, complete): every
dump that the two copies have in common, 164 in all, is byte-identical (`cmp`), and the
analyses use the `_r2` copies. The analysis scripts of README §5 were then run on the
output (the white-hole ones on 2026-11-04).

Result: every figure and every number of the paper was reproduced to the quoted digits, the
same as in the full rerun above unless stated:

* convergence: the 21 orders of Table conv-scatter; collapse p = (2.62, 2.55, 2.34, 2.38,
  2.62, 2.42, 2.24), constraint ratios 3.69–4.69 / 3.95–6.26, Θ_AH(28) = 12.293 / 12.604 /
  12.677 (Richardson 12.70), M_ADM 3.70836 → 3.71107 at Θ = 250 on the moving mesh
  (+0.07 %), 3.71155 / 3.71170 / 3.71171 at the fixed-grid edge; M_hor = 6.338, M_ADM =
  3.7111, M_hor/M_ADM = 1.708; moving mesh: Θ_AH equal to 8.0e-4, 12.60 at t = 28 → 11.8 at
  t = 150, Θ(r_max) 284 → 124 on the fixed grid;
* amplitude scan and threshold: the same classification of the eleven N = 6400 runs; early
  horizon from A0 = 0.135 / 0.130 / 0.125 / 0.125 on the four grids;
* spin-0: baseline pulse +0.0606 / −0.0718, drift 10.7 / 11.3 in u; flux-identity residual
  0.020 % → 0.002 % (dr = 0.023) and 0.098 % / 0.006 % (dr = 0.047); clock speed 1.2406;
  E_∞ = 0.6231 / 0.5443 / 0.5323, order 2.71, Richardson 0.5302, E_∞ = 0.532 ± 0.023 =
  (14.4 ± 0.6) % of M_ADM, E_ψ = 0.197; probe medians 1.00012 and 1.26706; scalar
  reflection 2.8e-4 of the incident amplitude (the paper's 0.03 %);
* white holes: Θ_AH of the A0 = 0.12 black hole 4.846–5.068; the A0 = 0.10 timeline
  (anti-trapped 2.7–17.6 at t = 36, 4.1–22.2 at t = 40, bulk gone by t = 44); freezing
  times 53 / 63 / 77 / 101 / 135 (N0 < 0.1 from t = 52, 64, 78, 102, 136; 55 instead of 53
  for dr = 0.0625 before); the front ladder;
  the family ladders (13 / 14 / 17 and 14 / 15 / 18 points at formation); every trapped /
  anti-trapped transition of the σ = 3 and r0 = 15 runs, including the A0 = 0.04 central
  collapse at t = 22 / 24 / 27 / 32.5 and the A0 = 0.03 one at dr = 0.0031 and 0.0008 but
  not 0.0016;
* inner shell at t = 30: |P| = 0.253 / 0.354 / 0.508 on the uniform grids (0.355 / 0.509
  before) and 0.694 / 0.925 / 1.249 / 1.194 on the refined ones, |ΘC| = 28 / 49 / 85;
  per-family power laws: height ∝ dr^−0.50 (uniform) / dr^−0.47 (refined), energy density
  ∝ dr^−0.99 / dr^−0.95, width ∝ dr^0.84 / dr^0.79.

The one change to a number of the paper: the width exponent of the inner shell, quoted as
dx^0.75 after the full rerun above, is now fitted as dr^0.84 (uniform) and dr^0.79
(refined), and the paper now says "roughly as dr^0.8". The other differences from the
full rerun (the last digit of |P|, the coarsest freezing time) are below the precision the
paper quotes.
The same seven runs as before stop on the documented guards, after the phases the paper
uses (README §5).

Boundary-reflection tests (README §5, `analysis/boundary_reflection.py`, 2026-11-04). The
flat-space spin-0 tests (pulse trK ∝ exp[−(r−50)²/16], dr = 0.05, fixed grid, boundary at
R = 100 and 150 against a reference box 2R) were run with the shipped binary: the reflected
w_in at r = R − 10 is 5.6e-7 / 4.1e-7 of the incident w_out for v0 = 1 and 1.5e-6 / 7.5e-7
for v0 = 1.267 (paper: below 1e-5); with the code before the boundary fix (w_in held fixed,
d_t w_in = 0) the same runs give 3.45 % / 2.33 % and 3.53 % / 2.33 % (paper:
3.5 % and 2.3 %). The scalar reflection is the 2.8e-4 above. A collapse-background
comparison ("the late signal agrees to about 5 % with a run whose boundary is twice as far",
from an earlier test with the code before the spin-0 fix) was repeated (`reflbg_R150` against `basef`): the difference of Θ∂_Θψ in Θ = 15–110 at
t = 236–300 is 7.3 % of the reference's maximum, but the two moving meshes do not share the
spacing at large Θ (about 2× coarser at Θ = 60 with r_max = 150) and the runs already differ
by 0.1–2.3 % before the pulse reaches r_max, so the number mixes resolution with reflection.
The paper does not quote it; the flat-space numbers carry the claim.

## Scripts and documents against the paper (2026-11-13)

The scripts and documents were checked against the paper: every figure has its script and
recipe (`analysis/README.md`, "Figures of the paper", with the paper's figure numbers;
`fig_invariant.py` makes a supplementary figure that is not in the paper), docstrings cite
the paper's section titles, and every simulation number of the paper is printed by a
script, from the final-rerun runs (no change to `collapse.cpp`):

* probe: largest |v − v0|/v0 for t ≥ 40, 2.96e-4 (v0 = 1) and 1.89e-4 (v0 = 1.267) (paper:
  below 3e-4 after t ≃ 40);
* luminosity clock: the burst (deepest minimum of L) at u = 9.7–14.7 on Θ = 120–400
  (paper: u ≃ 10–15);
* Coulomb field: peak Θ²a_r ∝ Θ^0.73 (paper: roughly Θ^0.7);
* energy: peak |∂_u m| / peak |∂_u m_c| on the same sphere 25 / 15 / 8.1 / 5.3 at Θ = 200 /
  150 / 100 / 75 (paper: "more than an order of magnitude smaller" than the luminosity of
  Fig. luminosity_clock(a), whose spheres Θ = 120–400 have peaks 2.3–7.5 against ≤ 0.27);
* mass function: M_ADM − M_hor = −2.639 with the Richardson Θ_AH (paper −2.64); ρ_æ > 0
  outside the horizon only at Θ = 18.3–20.3, inside the scalar shell (Θ = 12.7–20.6 at
  t = 28; paper 18 ≲ Θ ≲ 20);
* collapse ladder: N_min = 0.979 / 0.703 / 0.383 and, for A0 = 0.06 at t = 14, constraint
  ratios 1.79 / 2.46 (paper 0.98, 0.70, 0.4; 1.8, 2.5); Θ_AH differences shrinking by 4.3;
  Θ_AH = 8.47, r_AH = 9.2, Θ(r = 150) = 136.8 at t = 2 and central lapse 7.1e-7, r_AH = 64.8,
  Θ(r = 150) = 88.3 at t = 28 (paper ≃ 8.5, 9, 137; < 1e-6, 65, 90); at the fixed-grid edge
  Θ = 234.5 the fixed grids agree to 6.7e-5 and fixed grid and moving mesh to 6.9e-6 (paper
  1e-4, 2e-5; with the second-order end points of the next section); the change of m between
  t = 2 and 28 on Θ = 100 / 150 / 200 / 250 is 1.217 / 0.346 / 0.144 / 0.073 %, fitted by
  A + B/Θ³ with A = −0.0102 % (power free: p = 3.12, A = +0.0037 %) (paper: Θ^−3, 1.2 %,
  extrapolation about 0.01 %); Θ = 250 leaves the fixed grid between the dumps t = 18 and 20
  (Θ(r_max) = 253.9 / 249.8; paper t ≃ 20);
* white hole: the light ray from the shell is 3.16 / 0.30 / 0.09 behind the escaping pulse
  at t = 30 / 34 / 40 (paper 3.2, 0.3, 0.1); the inner front carries a jump of 0.28–0.32 in
  P with N = 0.125–0.154 for t = 20–30, moves at −β − N to within 0.005 up to t = 36 and
  stalls at Θ ≃ 3.2 (paper: about 0.3, N ≃ 0.13, Θ ≃ 3.2); the anti-trapped layer on the
  escaping pulse is 0.89 / 0.90 / 0.84 wide at t = 40 and 0.059 / 0.062 / 0.067 at t = 80 on
  the three zone-A grids, spread 7–18 % at every time t = 40–80 (paper: about 0.9 → 0.06,
  within 20 %); the pulse FWHM decays as e^(−kt) with k = 0.14 / 0.15 / 0.16 over t = 50–75
  and 0.12 / 0.13 / 0.14 over t = 50–80 (paper: roughly e^(−0.15t)); energy density of the
  inner shell ∝ dr^−0.99 (uniform) / dr^−0.95 (refined) (paper dr^−1); the A0 = 0.10
  anti-trapped region splits between t = 36 and 38 and its inner part is gone at t = 42 on
  wh_wh_N25600 (paper: splits at t ≃ 36, gone by t ≃ 44).

Every paper figure whose script changed was regenerated from the final-rerun runs with
the changed scripts and compared with the paper's PDF rendered at 150 dpi: all are
pixel-identical (white_hole as described under "Readability pass" below).

## Readability pass (2026-11-13)

Names, comments, docstrings, paper labels (`eq:…`, `fig:…`, `sec:…` instead of numbers),
shared helpers in `analysis/spin0_common.py` (`dump_path`, `load_dump`,
`requested_dumps`, `norm_l2`, `run_amplitude`, `steepest_slope`) in place of duplicated
code, consistent option and output names (`bisect_bh_threshold.py --t_early`,
`boundary_reflection.py --theta_window/--t_window`, `fig_moving_mesh.py --t`;
`fig_flux_identity.py` now writes `flux_identity_dumps.png`, `flux_identity_rec.py`
writes `flux_identity.png` by default; the unused `fig_bh_threshold.py --t_last` and
`fig_energy.py --title_note` removed), and in `collapse.cpp` the computational spacing
renamed `dx`, the coupling combinations set once in `main`, the dead coupling clamps
removed, `make smoke` and `--help-legacy` added.

* `collapse.cpp`: the binaries before and after give bit-identical output (every file:
  profiles, `phi.dat`, `psi_times.dat`, auxiliary dumps, `spheres.bin/.txt`, stdout,
  stderr) on four short runs: moving mesh + zone A + sphere recorder + probe (c2 = 0);
  ingoing white-hole box with zone B pre-positioned and auxiliary dumps; the legacy
  positional form on the fixed grid; ingoing shell with zones A and B acquiring and
  releasing.  `make warn` (and `-Wshadow`) is silent.
* Every analysis script was rerun on the final-rerun runs before and after: printed
  numbers, figures (PDF at 150 dpi and PNG, pixel by pixel) and the saved arrays are
  identical, except for changed labels in the printout and fig_white_hole (next item).
* `fig_white_hole.py` (a result change): the panel-(a) dumps are opened by their
  requested times (`--times`, default the README recipe) and placed at their exact times
  from `phi.dat`.  Before, the nominal one-decimal times of the file names were used
  (29.25 → 29.2, 29.75 → 29.8, …), so the escaping-pulse trace started at 29.5 instead of
  `--t_pulse` = 29.25 and the quarter-time points sat up to 0.05 off.  The recipe and the
  default panel-(a) run are now the complete `box40_N51200_r2`.  Printed numbers that
  change: only the inner-front lines at t = 29.25–31 (times now 29.25 / 29.75 / 30.25 /
  30.75, and dr/dt between dumps −0.019 … −0.016 instead of −0.023 … −0.013, now within
  0.002 of −β − N) and the end point of the light ray at t = 54 (Θ = 36.39 instead of 36.58,
  from the four dumps t = 46–52 that `box40_N51200` lacks).  None of the paper's quoted
  values changes (light ray 3.16 / 0.30 / 0.09 behind the pulse, jump 0.28–0.32, N, the
  per-family power laws 0.84 / 0.79 and −0.50 / −0.47).  Panel (a) gains the t = 46–52
  dumps and the exact quarter times; the paper's Fig. white_hole is this version.

## Analysis fixes after the readability review (2026-11-13)

Each was checked by rerunning the affected scripts on the final-rerun runs before and after
(printed numbers; figures as PDF at 150 dpi and PNG, pixel by pixel; saved arrays).

* Second-order end points for θ₊: `spin0_common.apparent_horizon` and `fig_white_hole.py`
  (θ₋) now take Θ' from `spin0_common.ddr` (`np.gradient(…, edge_order=2)`) instead of the
  default `edge_order=1`.  Only the two end points change and the horizon and θ₋ zeros are
  interior: the outputs of `fig_white_hole.py` (both recipes), `fig_mass_function.py`,
  `fig_moving_mesh.py`, `freeze_and_horizon.py --horizon`, `scan_bh_threshold.py`,
  `bisect_bh_threshold.py` and `fig_bh_threshold.py` are identical.  The same change in
  the grid-edge mass of `fig_conv_collapse.py` (`m_at`) moves the printed fixed-grid edge masses at t = 28 (3.71155 /
  3.71170 / 3.71171 → 3.71143 / 3.71164 / 3.71168, drift +0.09 % → +0.08–0.09 %) and the
  agreements the paper quotes as upper bounds 1e-4 and 2e-5 (printed 4.3e-5 / 1.5e-5 →
  6.7e-5 / 6.9e-6); every other output of the script is identical.
* One horizon criterion: every script now finds horizons with
  `spin0_common.apparent_horizon` at proper radius r > 0.5 (`bisect_bh_threshold.py` and
  `fig_bh_threshold.py` used r > 0; the black-hole verdict of `bisect_bh_threshold.py`
  keeps its Θ_AH > 0.5 cut), and `scan_bh_threshold.py --rmin` (minimum of θ₊) defaults to
  the same 0.5 instead of 0.3.  `bisect_bh_threshold.py` (with and without `--t_early 20`)
  and `fig_bh_threshold.py` (printout, figure, arrays) are identical: zeros at r < 0.5
  occur only at six intermediate dumps (t = 26–38, central throats with Θ ≈ 0.09–0.15,
  never Θ > 0.5), while the early (t ≤ 20) and last-dump values are unchanged.  `scan_bh_threshold.py`
  keeps its classification; only the minimum of θ₊ (A0 = 0.06: 0.677 → 0.742; 0.07:
  −1.167 → −0.998; 0.10: −1.069 → −1.058) and the first time with θ₊ < 0 (0.08: 14 → 15;
  0.13: 28 → 29) change, none of them in the paper.
* `bisect_bh_threshold.py` prints the central lapse at the last dump of the run as a new
  column N0(end), next to N0(t<=t_early) (the last dump searched).  All other columns and
  the brackets are unchanged (0.134/0.135, 0.125/0.130, 0.120/0.125, 0.120/0.125 with
  `--t_early 20`); N0(end) is between 4e-13 and 5e-203 for every run except
  `bis_N12801_A0.120` (0.957, dispersal) and `bis_N25601_A0.135` (0.257 at its last
  dump, t = 32).
* `boundary_reflection.py spin0` computes a_N, b_N, v0 from new `--c1/--c2` options
  (`spin0_common.couplings`; defaults 0.7 / −0.5119403, `--c2 0` for the v0 = 1.267 runs)
  instead of parsing `err.txt`, which is now optional: if present, its 6-decimal values
  are printed next to the computed ones, with a warning if they differ beyond its printed
  precision (tested with a wrong `--c2`).  On the eight flat-space pairs of the 2026-11-04
  test (v0 = 1 and 1.267, R = 100 and 150, code before and after the boundary fix) the
  printed results are identical: 5.59e-7 / 4.13e-7 (v0 = 1), 1.47e-6 / 7.49e-7
  (v0 = 1.267); before the fix 3.45e-2 / 2.33e-2 and 3.53e-2 / 2.33e-2.  The `scalar`
  (2.8e-4) and `late` (7.3 %) tests are unchanged.
