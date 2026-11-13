# collapse — spherical scalar collapse in Einstein–æther theory

Companion code of E. Albertini, M. Vaglio, E. Barausse, [arXiv number to be filled by the
authors] (2026); please cite that paper if you use it. Copyright (c) 2026 Emma Albertini,
Massimo Vaglio, Enrico Barausse; released under the MIT licence (file `LICENSE`).

This is the public version of the code used for the paper: the numerical scheme of the
development version plus the corrections recorded in `VALIDATION.md` (the ingoing initial
data, six points of the scheme and the exact outer boundary condition for the incoming
spin-0 field; every production run of the paper was redone with the code shipped here, §6),
with only the features used in the paper, a documented single-file
source and a named-option command line. The
positional command line of the development version is still accepted (unless it uses one
of the removed features, see "Command line" below). `analysis/` holds the scripts that
produce every figure and number of the paper from the code's output (Python ≥ 3.9 with
numpy ≥ 2.0 — `np.trapezoid` — and matplotlib ≥ 3.5; the spin-0 scripts fork worker
processes, so they need Linux or macOS).

```
make            # g++ -O2 -std=c++11 -o collapse collapse.cpp
make warn       # -Wall -Wextra, no warnings expected
make smoke      # a few-second sanity run (N = 1601, A0 = 0.35, t = 8) in smoke_run/
./collapse --help            # named options (--help-legacy: the positional form)
```

## 1. Physics

The formulation is that of Garfinkle, Eling & Jacobson (GEJ, gr-qc/0703093): a massless
scalar ψ coupled to Einstein–æther theory with couplings `c1`, `c2` (run time) and
`c3 = c4 = 0`. Spherical symmetry; the æther is hypersurface orthogonal and unit, so it is
the unit normal of the `t = const` slices. The metric is

    ds² = −N² dt² + (dr + β dt)² + Θ(t,r)² dΩ²,

with γ_rr = 1: `r` is proper radial distance along the slice and the **areal radius Θ is an
evolved field**. Signature (−,+,+,+), u·u = −1, geometric units with the GEJ normalisation
of the scalar.

**Freely evolved fields:** ψ (`Psi`), its momentum along the æther P = u^a∂_aψ (`P`), the
trace of the extrinsic curvature K (`trK`), the areal radius Θ (`Theta`) and the radial
æther acceleration a_r (`Acr`). Their equations, in words:

* ψ is advected by the shift and changes at the rate N P;
* P obeys the spherical wave equation on the slice: ψ'' + 2Θ'ψ'/Θ plus the coupling to the
  acceleration (a_r ψ') and to the expansion (K P), times N;
* K obeys the trace of the evolution equation (the Raychaudhuri equation of the æther
  congruence): N K²/3 plus the æther term (c14−2)(a_r' + 2Θ'a_r/Θ + a_r²), the scalar
  source 2P² and the shear term 3(1−c13)Q², all divided by 2+c13+3c2;
* a_r obeys the radial æther field equation combined with the momentum constraint (GEJ,
  Appendix); in unitary gauge the æther equation follows from the metric equations through
  the Bianchi identity (Jacobson, arXiv:1108.1496), i.e. it keeps the Hamiltonian constraint
  satisfied in time; its principal part couples it to K' with coefficient −c123/(c14(1−c13)),
  which together with the K equation makes up the spin-0 wave;
* Θ changes at the rate N Θ(Q/2 − K/3) (the areal expansion).

Every equation also contains the advection term β̃ ∂_r f, with β̃ the shift relative to the
moving mesh. The exact expressions are in the header of `collapse.cpp`.

**Recovered on every slice** (and at every ICN stage):

* the shear `Q` from the momentum constraint
  Q' = −3QΘ'/Θ + [−Pψ' + (2+c13+3c2)K'/3]/(1−c13), integrated outward from Q(0)=0 (the
  first cell from the regular series, the ODE being singular at r=0);
* the lapse `NN` from ∂_r ln N = a_r (GEJ eq. 28), integrated outward from the centre and
  shifted so that N(r_max) = 1 — the æther fixes the slicing, so the lapse is not a gauge
  choice. (The paper writes the same condition as an inward integration from N(r_max) = 1;
  the equation is linear in ln N, so both give the same function to round-off.)
* the shift `NR` from β' = N(Q + K/3), integrated outward from β(0)=0.

The slicing is singularity avoiding (a_r grows, N → 0 and the interior freezes), so there
is **no excision** and no inner boundary.

**Spin-0 mode.** The pair (K, a_r) has characteristic speeds β̃ ± N v0, with
aN = (c14−2)/(2+c13+3c2), bN = −c123/(c14(1−c13)), v0 = √(aN bN). The defaults c1 = 0.7,
c2 = −0.5119403 give v0 = 1 (GEJ's values); c2 = 0 gives v0 = 1.2669. The code prints v0 at
start-up.

**Initial data.** A variant of the GEJ shell (GEJ use ψ = a0 exp[−(r²−r0²)²/s⁴]; here it is
multiplied by r²/100, so the peak value is ≃ A0)

    ψ(0,r) = (A0/100) r² exp[−((r²−r0²)/σ²)²],  r0 = 10, σ = 4 by default,

time symmetric (P = 0; r0 = 10, σ = 4 are the GEJ values), or purely ingoing, P = ψ' + ψ/r (`--ingoing=1`, as in
arXiv:1608.06970). K = a_r = 0 (as in GEJ) and N = 1 at t = 0; Θ is obtained from the Hamiltonian
constraint by RK4 from Θ(0)=0, Θ'(0)=1, with ψ' (and P for the ingoing shell) evaluated
analytically at the RK4 stage radii. For time-symmetric data Q = β = 0 and the solve is
done once. For the ingoing shell P ≠ 0 sources Q through the momentum constraint and Q²
enters the Hamiltonian constraint, so the Θ and Q solves are iterated to a fixed point
(|ΔΘ| < 1e-13, 6 iterations; Q is interpolated linearly to the stage radii) and β ≠ 0
follows from Q; the t = 0 constraint residual is then second order in dx (VALIDATION.md).

## 2. Numerical method

* Uniform computational grid x_i = i dx, i = 0…N−1, dx = r_max/(N−1). Second-order centred
  differences, one-sided second-order stencils at the outer point, parity at r = 0
  (ψ, P, K, Q even; Θ, a_r, β odd).
* Iterated Crank–Nicolson (Euler predictor, up to 12 corrector iterations, tolerance
  1e-10 on ψ and Θ), dt = 0.1 dx (times min R_x with the refinement). On the moving mesh
  the corrector is evaluated on the mesh at t+dt (mesh velocities frozen for the step),
  with or without refinement. A non-finite field value after a step stops the run with a
  FATAL line (the state of every step, the last one included, is checked).
* Fourth-order Kreiss–Oliger dissipation after every step, f −= σ D⁴f with
  σ = ε (dt/dx)/16 = ε CFL/16 = 1/16 for ε = 10 (the strength per step, independent of
  N since dt/dx is fixed; 1/(16F) per step with the refinement, i.e. the same rate per
  unit time). At the outer edge the
  wave fields use one-sided interior fourth differences (non-reflecting); Θ keeps a linear
  ghost extrapolation that pins Θ_xx = 0 at the edge, which is its boundary datum.
* Outer boundary: characteristic conditions on the incoming spin-0 mode (with the
  v0-dependent eigenvectors, the incoming field w_in = (bN/v0) K + a_r is evolved as for
  the outgoing linear spin-0 wave, d_t w_in = β̃ ∂_r w_in − N bN K/Θ, while the outgoing
  field keeps its interior equation) and on the incoming scalar mode (the spherical
  condition d_t[P + ψ' + cψ] = 0 with c = Θ'/Θ + Q/2 − K/3, i.e. Θψ constant along
  outgoing rays, including the ψ d_t c term with d_t Q taken as a backward difference:
  over the previous step in the predictor stage, between the ICN stages in the
  corrector), each applied only while its mode actually enters the grid (the
  *gate*, re-evaluated every step and logged on stderr when it changes). Ψ and Θ keep the
  interior equations.
* Lapse: (ln N)' = a_r integrated outward by the trapezoid rule and shifted so that
  N(r_max) = 1 (the problem is linear, so this is its exact solution); a non-finite lapse
  stops the run. Shift: β' = N(Q + K/3) from β(0) = 0.
* Apparent horizon (diagnostic): outermost zero of θ₊ = Θ' + (3Q−2K)Θ/6 at r ≤ 45, printed
  on stdout as the index of the cell that contains it, the radius of the zero of the linear
  interpolant and the time of the state. The flag latches monotonically and is unreliable near the black-hole
  threshold (it reacts to transient focusing and to sign changes near the r = 45 search
  edge): horizons, the Misner–Sharp mass and θ₋ are evaluated by the analysis scripts from
  the profiles.

### Moving mesh (ALE), `--ale=1` (default)

In the γ_rr = 1 gauge the grid falls inward (β ≃ −1 outside the matter), so the areal radius
of the outer edge shrinks by about one unit per unit time. The mesh is moved with the map

    r = R(x,t) = x + g(x) s(t),

g a C² quintic taper, 0 for x ≤ 0.4 r_max, u³(10−15u+6u²) in between, 1 at x = r_max;
ds/dt is chosen every step so that Θ(r_max) stays fixed. Derivatives become
∂_r = ∂_x/R_x, ∂_r² = ∂_x²/R_x² − R_xx∂_x/R_x³, and every advection coefficient gets β → β + R_t.
The taper must be C² at r_max (checked at start-up). `--ale=0` is the plain fixed grid
(the same binary with R = x), used for the convergence ladders.

**The ALE grid is non-uniform in the outer 60%: analysis scripts must differentiate with
`np.gradient(f, r)` using the `r` column, never with a constant spacing.**

### Moving refinement (zones A and B)

Optional, requires `--ale=1`. Before the taper acts, the computational coordinate is
reparametrised, ξ = x − a w [B(u) − B(u/λ)], u = (x − x_c)/w, b(u) = (1−u²)⁴, B = ∫b, λ = 4:
a narrow bump (half-width w) compresses the grid by a factor F at the centre x_c and a wide
bump (half-width 4w) re-stretches it, so that ξ = x exactly for |x − x_c| > 4w. Strength
a = (1−1/F)/(1−1/λ), ramped in over 5 time units. The time step shrinks by F (not F²).
Up to two zones; their displacements add and their cores are kept disjoint.

* **Zone A** (`--ref_F`, `--ref_ton`, `--ref_w`) follows the outgoing front: from
  `ref_ton` on it acquires the peak of |P| with Θ > 15 (and Θ < Θ(r_max) − 10) once that
  peak exceeds 1e-6 (so that round-off noise before the front arrives is never acquired), then moves
  its centre with a PI controller (time constants 0.5 and 1, |ẋ_c| ≤ 2) to follow the peak.
* **Zone B** (`--refB_F`, `--refB_ton`, `--refB_w`) follows an inner feature: it acquires
  the |P| peak with 1 < Θ < 15 and then follows it wherever it goes. With `--refB_x0=x0 > 0`
  it is instead **pre-positioned**: switched on at `refB_ton` at the fixed centre x0, held
  there until `--refB_ttrack`, then tracks the peak. This is how the white-hole inner shell,
  which forms abruptly at t ≈ 29.5, r ≈ 9.96, was resolved from before its formation.
* A zone is released (ramped out, centre frozen) if its centre would leave the allowed range
  (support clear of the origin and of the last 20 points), if the peak is lost or falls
  below 5% of its value at acquisition, or on a clash with the other zone. Events are
  logged on stderr. The factor F is fixed for the whole run (F = 8 and 16 in the paper).
* Requirements (checked): w ≥ 16 dx; a pre-positioned centre must satisfy
  1 + 4w ≤ x0 ≤ r_max − 20 dx − 4w.

## 3. Command line

Named options, any order, all optional:

| option | default | meaning |
|---|---|---|
| `--N=` | 25600 | number of grid points (r = 0 … r_max inclusive) |
| `--A0=` | 0.35 | shell amplitude |
| `--t_final=` | none | stop time. Without it the run stops after 60000 steps |
| `--times=t1,t2,…` | none | times at which a full profile is written |
| `--rmax=` | 300 | outer edge of the grid |
| `--c1=` | 0.7 | coupling c1 |
| `--c2=` | −0.5119403 | coupling c2 (v0 = 1 with c1 = 0.7) |
| `--ingoing=` | 0 | ≠0: purely ingoing shell, P = ψ' + ψ/r |
| `--ale=` | 1 | 1: moving mesh; 0: fixed grid |
| `--ref_F=` | 0 | zone A refinement factor (≤1: off) |
| `--ref_ton=` | 50 | zone A: earliest acquisition time |
| `--ref_w=` | 5 | zone A: half-width of the fine core (x units) |
| `--refB_F=` | 0 | zone B refinement factor (≤1: off) |
| `--refB_ton=` | 20 | zone B: acquisition / pre-positioning time |
| `--refB_w=` | 2 | zone B: half-width |
| `--refB_x0=` | 0 | zone B: pre-positioned centre (0: search for the peak) |
| `--refB_ttrack=` | 0 | zone B: hold at x0 until this time, then track |
| `--r0=` | 10 | shell radius |
| `--sigma=` | 4 | shell width parameter (also the width of the probe Gaussian) |
| `--probe_amp=` | 0 | linear spin-0 probe: trK += amp exp(−((r−probe_r)/σ)²) at t=0 |
| `--probe_r=` | 150 | centre of the probe |
| `--dump_every=` | 0 | auxiliary per-field dumps every k steps (env `COLLAPSE_DUMP_EVERY` overrides, with a note on stderr) |
| `--rec_spheres=Θ1,Θ2,…` | none | sphere recorder (off by default): stencils around these areal radii, see §4 |
| `--rec_every=` | 1 | sphere recorder: record every k steps |

### Legacy positional form and mapping

`./collapse N A0 t_final times rmax c1 c2 ingoing g_ale ale_inner ref_F ref_ton ref_w refB_F
refB_ton refB_w refB_Fmax qform refB_x0 refB_ttrack r0 sigma` (`./collapse --help-legacy`)
is still accepted (trailing
slots may be omitted; pass an empty string `""` for `times` to reach later slots). The
slots of removed features must be absent or 0, otherwise the run stops with an error.

| slot | original name | release option |
|---|---|---|
| 1 | N | `--N` |
| 2 | A0 | `--A0` |
| 3 | t_final | `--t_final` |
| 4 | t1,t2,… | `--times` |
| 5 | rmax | `--rmax` |
| 6 | c1 | `--c1` |
| 7 | c2 | `--c2` |
| 8 | ingoing | `--ingoing` |
| 9 | g_ale | `--ale` |
| 10 | ale_inner | **removed** (inner taper; must be 0) |
| 11 | ref_F | `--ref_F` |
| 12 | ref_ton | `--ref_ton` |
| 13 | ref_w | `--ref_w` |
| 14 | refB_F | `--refB_F` |
| 15 | refB_ton | `--refB_ton` |
| 16 | refB_w | `--refB_w` |
| 17 | refB_Fmax | **removed** (adaptive zone-B factor; must be 0) |
| 18 | qform | **removed** (W-form momentum constraint; must be 0) |
| 19 | refB_x0 | `--refB_x0` |
| 20 | refB_ttrack | `--refB_ttrack` |
| 21 | r0 | `--r0` |
| 22 | sigma | `--sigma` |
| env | COLLAPSE_DUMP_EVERY | `--dump_every` (env still honoured) |
| — | (not in the development version) | `--probe_amp`, `--probe_r` |

Example: `./collapse 3200 0.10 35 "20,30" 40 0.8 -0.7111111 1 1 0 0 50 5 8 20 1.5 0 0 9.9 30.5`
is `./collapse --N=3200 --A0=0.10 --t_final=35 --times=20,30 --rmax=40 --c1=0.8
--c2=-0.7111111 --ingoing=1 --refB_F=8 --refB_ton=20 --refB_w=1.5 --refB_x0=9.9
--refB_ttrack=30.5`.

The development-tree command lines that **cannot** be run with this binary are those with
`ale_inner ≠ 0` (the t = 2000 "inner taper" runs, `... 1 1 3`), `refB_Fmax ≠ 0` or
`qform = 1`; none of them enters the paper.

## 4. Output

All files go to the working directory (run each job in its own directory).

* `prof_tTTT.T.csv` — one per requested time, written at the first step with t ≥ the
  requested time (the actual time is later by less than one step) and named after the
  requested time with format `%05.1f` (so two requested times that round to the same
  tenth overwrite each other; e.g. 29.25 is written as `prof_t029.2.csv`; this and the
  absence of a t = 0 profile — the first dump is the first step at or after the smallest
  requested time, e.g. `--times=0.01` — are known limitations kept for compatibility with
  the paper's data sets; the analysis scripts open such dumps by their requested time and
  read the exact time from `phi.dat`, `analysis/spin0_common.py`). One header line,
  then N rows, 17 significant digits:

  | column | content |
  |---|---|
  | `r` | physical radius R(x_i,t) (non-uniform with the mesh/refinement on) |
  | `Psi` | ψ |
  | `P` | P = u^a∂_aψ |
  | `trK` | K |
  | `Acr` | a_r |
  | `Theta` | areal radius Θ |
  | `NN` | lapse N |
  | `Q` | shear Q |
  | `ThetaC` | Θ × (Hamiltonian-constraint residual); ignore the last 3 points (Θ edge pin) |

  The shift is not written; it is β(r) = ∫₀^r N(Q+K/3) dr' (the analysis scripts rebuild it
  that way). Useful combinations: Misner–Sharp mass m = Θ(1−h)/2,
  h = Θ'² − (3Q−2K)²Θ²/36; expansions θ± = Θ' ± (3Q−2K)Θ/6; outgoing spin-0 characteristic
  w_out = −(bN/v0) K + a_r (= 0.8955 K + a_r for v0 = 1).
* `phi.dat` — `step t dt phi` every step, phi = ψ(r=0).
* `psi_times.dat` — index (`dump_id step t dt`) of the auxiliary dumps; header only when
  they are off.
* `psi_dump_NNNNNN.csv` (`i,r,Psi,P,trK,Acr,Theta,NN,NR,Q,DTheta`) and one
  `<field>_dump_NNNNNN.csv` per field — only with `--dump_every=k`. They are large; the
  paper's analyses use the profiles.
* `spheres.bin`, `spheres.txt` — only with `--rec_spheres=Θ1,Θ2,…` (sphere recorder, for
  time derivatives on spheres of fixed areal radius at the step cadence; used by
  `analysis/flux_identity_rec.py`). At t = 0 and after every `rec_every`-th completed step
  (fields, mesh and the Q/N/β solves all at the new time), for each Θ_s the code finds the
  OUTERMOST index i with Θ[i] ≤ Θ_s < Θ[i+1] and stores the 13 grid points i−6 … i+6 of
  r (the physical radius R(x_i,t) of the moving mesh), ψ, P, K, a_r, Θ, N, Q and β (the
  shift). `spheres.bin` is a sequence of little-endian float64 records

      step, t, dt,  then for each sphere:  i, r[13], Psi[13], P[13], trK[13], Acr[13],
                                           Theta[13], NN[13], Q[13], NR[13]

  (3 + 118·n_spheres doubles; dt = the step just completed, 0 in the t = 0 record; i = −1
  when Θ never reaches Θ_s, the stencil is then the outermost 13 points; the stencil is
  shifted to stay inside 0 … N−1). `spheres.txt` lists the layout, the spheres, rec_every,
  N, r_max, couplings and ID parameters. Size: 8·(3+118 n_s) bytes per record, e.g. 9
  spheres at N = 25601 to t = 600 (256 000 steps) is 2.2 GB with `--rec_every=1` and
  440 MB with `--rec_every=5`. Writes are buffered (4 MB); the recorder only reads the
  state, so all other output is bit-identical with it on or off (VALIDATION.md).
* stdout: `t= dt= hor=` every 500 steps (plus the mesh amplitude and zone centres with the
  refinement), `Horizon detected/updated: i r t` (diagnostic only, see §2), and
  `FIRST NaN` if a field goes non-finite, after which the run stops with a FATAL line on
  stderr (the profiles, `phi.dat` and `psi_times.dat` written so far are kept and
  complete). stderr: v0, gate changes, refinement
  events, the convergence of the Θ/Q iteration of the ingoing initial data, fatal errors.
* The code writes stdout and stderr to the terminal; the analysis scripts expect them
  saved in the run directory as `log.txt` and `err.txt` where they use them
  (`./collapse … > log.txt 2> err.txt`; `analysis/boundary_reflection.py spin0` checks
  its `--c1/--c2` against the couplings in `err.txt`), and
  optionally a copy of the command line as `cmd.txt` (`--A0=` is read from it by the
  threshold scripts, otherwise from the directory name). `analysis/README.md`, "Input".

## 5. Reproducing the runs of the paper

Every run in its own directory. `T(a,b,h)` below stands for the comma list a, a+h, …, b,
e.g. `$(python3 -c "print(','.join(f'{0.5*i:.1f}' for i in range(1,1201)))")` for 0.5…600
in steps of 0.5. The dump cadence matters: Δt = 0.5 is needed for the flux, luminosity
and spin-0 analyses (Δt = 2 smears the burst). The command lines below are the ones used
for the public-code reproduction of the paper (October 2026) and, unchanged, for its final
rerun with the code shipped here (29–31 October 2026, `VALIDATION.md`, "Final rerun"); each
recipe names the script in `analysis/` that produces the figure or the numbers from the
output (`analysis/README.md` lists them all). Every figure and number of the paper comes
from these commands.

Wall times are those measured in the reproduction campaign of October 2026 (113
single-core runs, g++ -O2, up to 80 of them at once under `nice` on one 112-core machine,
so they are loaded-machine times; the runs marked (m) were run with the malloc settings
below and at most 20 at once). In the final rerun (115 runs, up to 40 at once, all with the
malloc settings) the short runs took about as long or less and the longest ones up to about
twice as long, e.g. 38 h for the refined N = 51200 white-hole grid and 24 h for its F = 16
variant, while the four jobs repeated there with few neighbours took 7.3–12.7 h, close to
the times quoted here. The disk use is that of the profile dumps (17-digit CSV,
about 2.2 MB per 10 000 grid points per dump). The code itself needs memory for a few
dozen arrays of N doubles.

For long runs set, before starting the binary,

    export MALLOC_MMAP_THRESHOLD_=4000000000 MALLOC_TRIM_THRESHOLD_=4000000000 MALLOC_TOP_PAD_=268435456

The code allocates large temporaries every step and glibc otherwise returns them to the
kernel each time (19–31 % of the CPU time was spent in the kernel on the authors'
machine). With these settings the output is bit-identical and the runs are about 18 %
faster (393 s → 322 s for an N = 25600 white-hole test to t = 1.5).

**Scattering convergence ladder** (Sec. Convergence/Scattering; fixed grid, r_max = 150):

    for A in 0.01 0.04 0.06; do for N in 6401 12801 25601; do
      ./collapse --N=$N --A0=$A --t_final=14.5 --times=8,9,9.5,10,10.5,11,11.5,12,13,14 \
                 --rmax=150 --ale=0           # in a directory named A${A}_N${N}
    done; done                                # 3 min / 10-15 min / 45-60 min per run, 14 / 28 / 56 MB
    # then (Figs. conv_scatter_fields_A001/_A006 and conv_scatter_constraint, Table conv-scatter;
    # also N_min = 0.98 / 0.70 / 0.38 and the A0 = 0.06 constraint ratios 1.8 / 2.5 at t = 14)
    #   python3 analysis/fig_conv_scatter.py WORKDIR
    # (norms over r < 20, orders averaged over t = 12-14; equivalent positional form
    # "N A 14.5 TIMES 150 0.7 -0.5119403 0 0")

**Amplitude scan** (Sec. Convergence: on the dr = 0.023 grid A0 = 0.05–0.06 disperse,
A0 = 0.07–0.11 and 0.13 collapse the centre into a cell-sized frozen throat (0.12 disperses
again), and only A0 = 0.14–0.15 form a genuine detached horizon; `analysis/scan_bh_threshold.py astar_A*`; 7–13 min and
84 MB per run):

    for A in 0.05 0.06 0.07 0.08 0.09 0.10 0.11 0.12 0.13 0.14 0.15; do
      ./collapse --N=6400 --A0=$A --t_final=60 --times=T(1,60,1) --rmax=150 > log.txt
    done

**Black-hole threshold** (Fig. bh_threshold: bisection in A0 on four grids, r_max = 300,
t ≤ 60, dumps every 2; the figure plots the outermost trapped surface found at t ≤ 20 (the
early horizon) and at t = 60; the smallest amplitude with an early horizon is 0.135, 0.130,
0.125, 0.125 for N = 12801, 25601, 51201, 102401, so A* lies in 0.120–0.125 on the two
finest grids; `analysis/fig_bh_threshold.py bis_N*` for the figure,
`analysis/bisect_bh_threshold.py bis_N* --t_early 20` for the per-run table; 20–50 min /
1.6–1.9 h / 4.3–4.7 h / 8.5–8.9 h per run at N = 12801 / 25601 / 51201 / 102401, 84 / 168 /
335 / 670 MB):

    for N in 12801 25601 51201 102401; do for A in 0.105 0.110 0.115 0.120 0.125 0.130 0.135 0.140 0.145 0.150; do
      ./collapse --N=$N --A0=$A --rmax=300 --t_final=60 --times=T(2,60,2)   # in bis_N${N}_A${A}
    done; done
    # the paper's runs: N = 12801 at 0.120–0.150 (0.130–0.134 in steps of 0.001), N = 25601 at
    # 0.105–0.150, N = 51201 at 0.115–0.130, N = 102401 at 0.115–0.125; the N = 25601,
    # A0 = 0.135 run stops on the non-finite-lapse guard at t = 32 (after its early horizon)

**Collapse convergence ladder** (Sec. Convergence/Collapse and Appendix "Convergence
tests", Fig. conv_collapse, Tables I and II; fixed grid, norms outside the AH, orders
averaged over t = 4–28; paper: p = (2.62, 2.55, 2.34, 2.38, 2.62, 2.42, 2.24),
constraint-norm ratios 3.69–4.69 on the coarse pair and 3.95–6.26 on the fine pair,
Θ_AH(t = 28) = 12.293 / 12.604 / 12.677, differences shrinking by 4.3, Richardson 12.70;
M_ADM drift between t = 2 and 28 +0.07 %, 3.70836 → 3.71107 at Θ = 250 on the moving-mesh
run mass_hi below (pass it with --mesh_run); on the fixed grid the outer edge falls to
Θ ≈ 235 by t = 28, so the script reads the fixed-grid value at the grid edge, 3.70836 →
3.71143 / 3.71164 / 3.71168, +0.08–0.09 %; with --mesh_run it also prints the agreement there
(6.7e-5 between the fixed grids, 6.9e-6 with the moving mesh; paper 1e-4, 2e-5), the
dump at which Θ = 250 leaves the fixed grid (t = 20; paper t ≃ 20), the Θ^−3 falloff of
the change of m (1.2 % at Θ = 100) and its extrapolation (|A| ≲ 0.01 %), and the
t = 2 / 28 values of Θ_AH, r_AH, the central lapse and Θ(r = 150) quoted in the
appendix):

    for N in 12801 25601 51201; do            # 16 min / 45 min / 2.7 h; 39 / 78 / 156 MB
      ./collapse --N=$N --A0=0.35 --t_final=28 --times=T(2,28,2) --rmax=300 --ale=0
    done
    python3 analysis/fig_conv_collapse.py col_N12801 col_N25601 col_N51201 --mesh_run mass_hi   # Fig. conv_collapse + numbers

**Moving-mesh validation** (Fig. moving_mesh: Θ(r_max) and Θ_AH(t) of the two runs):

    ./collapse --N=25600 --A0=0.35 --t_final=150 --times=T(2,150,2) --rmax=300 --ale=0
    ./collapse --N=25600 --A0=0.35 --t_final=150 --times=T(2,150,2) --rmax=300 --ale=1
    python3 analysis/fig_moving_mesh.py mesh_ale0 mesh_ale1      # 2.7 h / 2.8 h, 420 MB per run

**Mass function and horizon** (Fig. mass_function, Sec. "The mass function and the
horizon", profile at t = 28, finest grid; paper: M_hor = 6.35 from the Richardson limit
12.70 of Θ_AH, M_ADM = 3.7111, M_ADM − M_hor = −2.64, M_hor/M_ADM = 1.71, the æther energy
density negative outside the horizon except across the outgoing scalar shell; the script
prints Θ_AH = 12.677 on this grid, M_hor = 6.338, M_ADM(Θ = 250) = 3.7111, 1.708, −2.639 with
the Richardson Θ_AH, and ρ_æ > 0 only at Θ = 18.3–20.3, inside the shell; the public-code
run before the moving-mesh corrector correction of VALIDATION.md gave 3.7109):

    ./collapse --N=51201 --A0=0.35 --t_final=28 --times=T(2,28,2) --rmax=300 --ale=1   # 2.7 h, 156 MB
    python3 analysis/fig_mass_function.py mass_hi          # Fig. mass_function + numbers

(the fixed-grid N = 51201 run of the collapse ladder gives the same t = 28 profile to
within the convergence error; the Θ_AH(t) history comes from the moving-mesh runs above).

**Spin-0 wave, v0 = 1 baseline** (Fig. wout_background; `analysis/fig_wout_background.py`;
4.5 h, 3.4 GB):

    ./collapse --N=25600 --A0=0.35 --t_final=300 --times=T(0.5,300,0.5) --rmax=300
    python3 analysis/fig_wout_background.py basef --madm 3.711

**Outer-boundary reflection** (Sec. "The spin-0 wave far from the source", "These
conditions reflect very little";
`analysis/boundary_reflection.py`; each test pairs a run with the boundary at R with a
reference at 2R and the same dr):

    # linear spin-0 pulse trK ∝ exp[-(r-50)^2/16] in flat space, v0 = 1 and 1.267 (minutes each)
    for c2 in -0.5119403 0; do for R in 100 150; do
      ./collapse --N=$((R*20+1))   --rmax=$R          --A0=0 --probe_amp=1e-6 --probe_r=50 --c2=$c2 --ale=0 --t_final=$((2*R)) --times=T(1,$((2*R)),1) 2> err.txt
      ./collapse --N=$((R*40+1))   --rmax=$((2*R))    --A0=0 --probe_amp=1e-6 --probe_r=50 --c2=$c2 --ale=0 --t_final=$((2*R)) --times=T(1,$((2*R)),1) 2> err.txt
      python3 analysis/boundary_reflection.py spin0 <run R> <run 2R> --R $R --c2=$c2   # checks <run R>/err.txt
    done; done
    # paper: below 1e-5 (public code: 0.4-1.5e-6). The 3.5 % (R = 100) / 2.3 % (R = 150)
    # with w_in set to zero come from the same runs with the code before the boundary
    # fix (w_in held fixed, d_t w_in = 0; VALIDATION.md, "Final rerun"): 3.45/3.53 % and 2.33/2.33 %
    # weak scalar shell in flat space (1.2 h / 3.7 h; 214 / 417 MB)
    ./collapse --N=12801 --rmax=150 --ale=0 --A0=1e-4 --t_final=175 --times=T(100,175,1)
    ./collapse --N=25601 --rmax=300 --ale=0 --A0=1e-4 --t_final=175 --times=T(100,175,1)
    python3 analysis/boundary_reflection.py scalar screfl_R150 screfl_R300   # paper 0.03 % (2.8e-4)
    # collapse background (not in the paper): basef above (r_max = 300) against the same run
    # with r_max = 150 (20 min, 1.6 GB). The moving mesh coarsens toward the boundary, so the
    # two runs differ in resolution at large Theta (spacing at Theta = 60 about 2x larger with
    # r_max = 150) and the difference, 0.1-2.3 % before the pulse reaches r_max and 7.3 % in
    # Theta = 15-110, t = 236-300, is an upper bound on the reflection, not a measurement of it
    ./collapse --N=12800 --A0=0.35 --t_final=300 --times=T(0.5,300,0.5) --rmax=150
    python3 analysis/boundary_reflection.py late reflbg_R150 basef --theta_window 15,110

**Superluminal spin-0, v0 = 1.267, energy and flux** (c1 = 0.7, **c2 = 0**; Figs.
luminosity_clock, flux_identity, aether_coulomb, spin0_flux_infinity, and the supplementary
invariant_wout_dm, which is not in the paper):

    # dr = 0.047 / 0.023 / 0.0117: luminosity clock, invariant waveform, Coulomb field,
    # E(Theta,u) and the E_inf resolution ladder (u in [200,250] at Theta <= 200 needs t >~ 430)
    for N in 12801 25601 51201; do            # 1.5 h / 4.5 h / 5.9 h (m); 3.4 / 6.7 / 13.4 GB of dumps
      ./collapse --N=$N --A0=0.35 --c2=0 --rmax=600 --t_final=600 --times=T(0.5,600,0.5)
    done
    # flux identity (Fig. flux_identity, analysis/flux_identity_rec.py): the same runs with
    # the sphere recorder instead of the profile dumps (dr = 0.023 is the paper's figure;
    # dr = 0.047 for the convergence of the residual; no --times needed)
    for N in 12801 25601; do                  # 1.7 h / 4.4 h; spheres.bin 218 / 436 MB
      ./collapse --N=$N --A0=0.35 --c2=0 --rmax=600 --t_final=600 \
                 --rec_spheres=40,60,80,100,120,140,160,180,200 --rec_every=5
    done
    python3 analysis/flux_identity_rec.py sl600_N25601_rec --c2 0 --plot 40,100,200   # Fig. flux_identity
    python3 analysis/flux_identity_rec.py sl600_N12801_rec --c2 0                     # dr = 0.047 residuals

Keep c1 = 0.7 and change only c2 (c2 = 0 with c1 = 0.8 gives a different v0). Paper
values from these runs: flux identity within 0.02 % of the peak flux on every sphere
Θ = 40–200 at dr = 0.023, |F_ψ| ≲ 0.02 against æther-flux peaks 0.9–4
(`flux_identity_rec.py`); luminosity burst at u ≃ 10–15 and clock speed 1.241
(`fig_luminosity_clock.py` on the dr = 0.023 run: u = 9.7–14.7, 1.2406); Θ²a_r peak growing
roughly as Θ^0.7 (`fig_coulomb.py`: Θ^0.73); M_ADM = 3.708 from the first dump; radiated
energy (eq. Einf, `fig_energy.py`) E_∞ = 0.532 ± 0.023 = (0.144 ± 0.006) M_ADM, of which the
scalar carries E_ψ = 0.197 (5.3 %) and the æther E_æ = 0.335 (0.090 M_ADM); the corrected
rate −∂_u m_c more than an order of magnitude below the luminosity (`fig_energy.py` prints
the ratio of the peaks on each sphere: 25 at Θ = 200, 15 at 150, 8.1 and 5.3 at 100 and 75).
Also printed, not quoted in the paper: the residual 0.020 % at Θ = 40 → 0.002 % at Θ = 200
(0.098 % / 0.006 % at dr = 0.047); the E_∞ ladder 0.6231 / 0.5443 / 0.5323 on the three
grids (observed order 2.71, Richardson 0.5302; before the scheme corrections of
VALIDATION.md 0.6219 / 0.5435 / 0.5317, order 2.73, 14.3 %: the full scalar outer boundary
condition acts once the pulse reaches r_max = 600).

**Probe runs** (Fig. probe_speed; run-differencing, the probe must be linear):

    for c2 in -0.5119403 0; do                # 2.5 h and 840 MB per run
      ./collapse --N=25600 --A0=0.35 --c2=$c2 --rmax=350 --t_final=150 --times=T(1,150,1)
      ./collapse --N=25600 --A0=0.35 --c2=$c2 --rmax=350 --t_final=150 --times=T(1,150,1) \
                 --probe_amp=1e-6 --probe_r=150
    done
    # analysis/fig_probe_speed.py differences the two runs dump by dump; speed
    # (dr_pk/dt + beta)/N with beta rebuilt from the background run.  Paper: deviation
    # below 3e-4 after t ~ 40, medians within 0.01 % (v0 = 1) and 0.02 % (v0 = 1.267) of
    # v0; printed: largest deviation for t >= 40 3.0e-4 / 1.9e-4, medians 1.00012 and
    # 1.26706 (v0 = 1.2669), MAD 4e-5 / 8e-5

**White hole vs black hole at c1 = 0.8** (Sec. "Large c1: transient white holes";
c2 = −0.7111111 keeps v0 = 1; ingoing shell):

    W="--c1=0.8 --c2=-0.7111111 --ingoing=1"
    E="0.5,1,1.5,2,2.5,3"      # early dumps: min theta_+ = 0.34 at t = 0.5, trapped by t = 2.5
    # black hole at A0 = 0.12 (Theta_AH = 4.85-5.07 for t > 30, still slowly varying at
    # t = 100: ~4% in time, ~2% between the two grids)
    for N in 25600 51200; do ./collapse --N=$N --A0=0.12 $W --rmax=300 --t_final=100 --times=$E,T(4,100,2); done
    # (directories wh_bh_N25600, wh_bh_N51200: 2.0 h / 5.7 h, 308 / 615 MB)
    # transient white hole at A0 = 0.10 (wh_wh_N25600; theta_+ < 0 at t = 6-16, anti-trapped
    # from t = 20, 2.7-17.6 at t = 36, split at t = 38, inner part gone by t = 42: paper
    # "dip at t ~ 6-17, appears at t ~ 20, splits at t ~ 36, gone by t ~ 44")
    ./collapse --N=25600 --A0=0.10 $W --rmax=300 --t_final=100 --times=$E,T(4,100,2)   # 2.5 h, 307 MB
    # Theta_AH(t) of the black-hole runs; trapped / anti-trapped history of both
    python3 analysis/freeze_and_horizon.py --horizon wh_bh_N25600 wh_bh_N51200
    python3 analysis/family_transitions.py wh_bh_N25600 wh_wh_N25600 --extent

**Escaping pulse with zone A locked on before the compression** (Fig. sharp_features(b) and
the escaping-pulse numbers: width ∝ e^(−0.15t), the anti-trapped layer 0.9 → 0.06 between
t = 40 and 80 agreeing to 20 % on the three grids; refinement 8× moving with the pulse; base
grids 12800/25600/51200 on r_max = 150):

    for N in 12800 25600 51200; do            # 2.6 h / 6.5 h / 13.4 h; 227 / 453 / 905 MB
      ./collapse --N=$N --A0=0.10 $W --rmax=150 --t_final=80 --times=T(40,80,0.5) \
                 --ref_F=8 --ref_ton=40 --ref_w=3
    done
    python3 analysis/wh_front_convergence.py whconv_N12800 whconv_N25600 whconv_N51200

**Inner shell: the seven-grid ladder at t = 30** (Fig. white_hole; r_max = 40; the shell
width, |P| peak and density at t = 30 against the grid spacing at the shell):

    Q="29,29.25,29.5,29.75,30.25,30.5,30.75,31"    # quarter-time dumps around the formation
    # three uniform grids (dx = 0.0031 / 0.0016 / 0.0008; 1.9 h / 2.6 h (m) / 7.3 h (m) up to the
    # stop at t = 55-56, 73 / 145 / 291 MB); N = 51200 also gives panel (a) (fig_white_hole.py
    # opens these dumps by the requested times: 29.25 is written as prof_t029.2.csv)
    for N in 12800 25600 51200; do
      ./collapse --N=$N --A0=0.10 $W --rmax=40 --t_final=56 --times=T(20,54,2),$Q
    done
    # three grids with zone B pre-positioned on the formation radius, F = 8
    # (4.7 h / 10.4 h / 18.2 h (m); 145 / 290 / 581 MB)
    for N in 12800 25600 51200; do
      ./collapse --N=$N --A0=0.10 $W --rmax=40 --t_final=50 \
                 --times=20,24,26,28,T(28.5,50,0.5),29.25,29.75,30.25,30.75 \
                 --refB_F=8 --refB_ton=20 --refB_w=1.5 --refB_x0=9.9 --refB_ttrack=30.5
    done
    # the finest, F = 16 (12.1 h (m), 101 MB)
    ./collapse --N=51200 --A0=0.10 $W --rmax=40 --t_final=31 --times=29,29.25,29.5,29.75,30,30.25,30.5,30.75,31 \
               --refB_F=16 --refB_ton=20 --refB_w=1.5 --refB_x0=9.9 --refB_ttrack=30.5

    # Fig. white_hole (panel (a) from the finest uniform run) and the scaling exponents; also
    # prints the light ray behind the escaping pulse (3.2 / 0.3 / 0.1 at t = 30 / 34 / 40)
    # and the inner-front track (jump ~0.3 in P, N ~ 0.13, stalling at Theta ~ 3.2)
    python3 analysis/fig_white_hole.py --uniform box40_N12800,box40_N25600,box40_N51200_r2 \
            --refined box40pp_N12800,box40pp_N25600,box40pp_N51200,box40pp_N51200_F16
    python3 analysis/wh_formation_ladder.py box40_N12800 box40_N25600 box40_N51200 \
            box40pp_N12800 box40pp_N25600 box40pp_N51200 box40pp_N51200_F16 --series box40pp_N51200_F16
    # Fig. sharp_features: the three sharp features (escaping pulse from whconv_N51200 above;
    # inner shell at t = 30; inner ingoing scalar front at Theta ~ 3.2-3.5, steepest slope
    # 4.9 / 6.2 / 6.9 at t = 24 and 7.1 / 14.8 / 30.7 at t = 40 on the three uniform grids)
    python3 analysis/fig_sharp_features.py --uniform box40_N12800,box40_N25600,box40_N51200 \
            --refined box40pp_N51200,box40pp_N51200_F16 --escaping whconv_N51200

(the uniform runs stop at t ≈ 55–56 on the FATAL edge guard `dTheta/dr at rmax`, after the
outgoing front has reached r_max = 40; the profiles up to t = 54 are unaffected). In the
final rerun `box40_N51200` and `whconv_N51200` lost dumps when the disk filled and were
repeated as `box40_N51200_r2` and `whconv_N51200_r2` (VALIDATION.md, "Final rerun"); the
paper's numbers and Fig. sharp_features come from the complete `_r2` copies, while the
paper's Fig. white_hole was drawn with the incomplete `box40_N51200` as the finest uniform
run (four dumps missing at t = 46–54 in panel (a)) and with the one-decimal nominal times
of its quarter-time dumps (e.g. t = 29.2 for the dump requested at 29.25, which dropped
that dump from the escaping-pulse trace). The `fig_white_hole.py` call above uses the
complete `box40_N51200_r2` and the exact dump times.

**Frozen centre vs resolution** (the freezing time of the central lapse moves later on
every finer grid, roughly as dr^−1/3, with no sign of converging to a finite time; the
script prints t ≃ 53, 63, 77, 101, 135 for dr = 0.0625, 0.031, 0.016, 0.0078, 0.0039):

    for N in 1601 3201 6401 12801; do         # 1 min / 7 min / 26 min / 2.4 h; 26 / 53 / 105 / 210 MB
      ./collapse --N=$N --A0=0.10 $W --rmax=100 --t_final=150 --times=T(2,150,2)
    done
    ./collapse --N=38401 --A0=0.10 $W --rmax=150 --t_final=150 --times=T(2,150,2)   # dr = 0.0039; 7.3 h, 629 MB
    python3 analysis/freeze_and_horizon.py --freeze freeze_N1601:0.0625 freeze_N3201:0.031 \
            freeze_N6401:0.016 freeze_N12801:0.0078 freeze_N38401_r150:0.0039

(at dr = 0.0039 the r_max = 100 box hits the outer-edge mesh guard at t ≈ 117 once the
front reaches the edge, hence the larger box for the finest grid).

**Other initial data** (the paragraph on `(r0,σ) = (10,3)` and `(15,4)` of Sec. "Large c1:
transient white holes"; r_max = 60 unless stated; directories
`fam_s3_A<A0>_N<N>` and `fam_r15_A<A0>_N<N>`):

    S3="--r0=10 --sigma=3"; R15="--r0=15 --sigma=4"
    # scans at dr = 0.0094 (13-27 min and 42 MB each): sigma = 3, A0 = 0.01 ... 0.08; r0 = 15, A0 = 0.020 ... 0.036
    for A in 0.01 0.015 0.02 0.03 0.04 0.05 0.06 0.07 0.08; do
      ./collapse --N=6400 --A0=$A $W $S3  --rmax=60 --t_final=60 --times=T(2,60,2); done
    for A in 0.020 0.024 0.028 0.032 0.036; do
      ./collapse --N=6400 --A0=$A $W $R15 --rmax=60 --t_final=60 --times=T(2,60,2); done
    # inner-shell ladders, dr = 0.0031 / 0.0016 / 0.0008 (2.4-2.5 h / 6.1-6.2 h / 8.4-8.7 h (m); 109 / 218 / 436 MB)
    for N in 19200 38400 76800; do
      ./collapse --N=$N --A0=0.06  $W $S3  --rmax=60 --t_final=45 --times=T(20,45,1)
      ./collapse --N=$N --A0=0.020 $W $R15 --rmax=60 --t_final=45 --times=T(20,45,1)
    done
    # sigma = 3, A0 = 0.03 and 0.04: central collapse vs resolution (dumps every 0.5 from t = 8;
    # 2.4 h / 5.8 h / 7.4-8.3 h (m); 231 / 461 / 922 MB). A0 = 0.04: the centre collapses at
    # t = 22 / 24 / 27 / 32.5 for dr = 0.0094 / 0.0031 / 0.0016 / 0.0008 (the first from the
    # scan above); A0 = 0.03: at t = 17 for dr = 0.0031, not up to t = 35 for dr = 0.0016, and
    # again at t = 22 for dr = 0.0008 (non-monotonic in the resolution; the dr = 0.0094 scan
    # run stops on the non-finite-lapse guard at t = 15.6, last dump t = 14)
    for A in 0.03 0.04; do for N in 19200 38400 76800; do
      ./collapse --N=$N --A0=$A $W $S3 --rmax=60 --t_final=35 --times=T(8,35,0.5); done; done
    # r0 = 15, A0 = 0.028 (trapped -> anti-trapped) on three grids
    ./collapse --N=21334 --A0=0.028 $W $R15 --rmax=200 --t_final=100 --times=T(2,100,2)   # dr = 0.0094; 2.4 h, 214 MB (stops at t = 92, see below)
    ./collapse --N=42667 --A0=0.028 $W $R15 --rmax=200 --t_final=100 --times=T(2,100,2)   # dr = 0.0047; 3.2 h (m), 466 MB
    ./collapse --N=25641 --A0=0.028 $W $R15 --rmax=60  --t_final=40  --times=T(6,40,1)    # dr = 0.0023; 3.1 h, 196 MB
    # sigma = 3, A0 = 0.07 and 0.08 at dr = 0.0047 in the larger box (4.9 h, 316 MB / 5.7 h, 466 MB)
    for A in 0.07 0.08; do
      ./collapse --N=42667 --A0=$A $W $S3 --rmax=200 --t_final=100 --times=T(2,100,2); done
    python3 analysis/family_transitions.py fam_*                       # classification of every run
    python3 analysis/wh_family_ladder.py fam_s3_A0.06_N19200 fam_s3_A0.06_N38400 fam_s3_A0.06_N76800
    python3 analysis/wh_family_ladder.py fam_r15_A0.020_N19200 fam_r15_A0.020_N38400 fam_r15_A0.020_N76800

(two of these runs stop before their t_final on the documented guards, after the phases
the paper uses: the r0 = 15, A0 = 0.028, dr = 0.0094 run on the non-finite-lapse guard at
t = 92 — its trapped → anti-trapped transition is at t ≈ 24 — and the σ = 3, A0 = 0.07,
dr = 0.0047 run on a fold of the moving mesh at t = 69.6, after its anti-trapped phase
t = 16–34.)

## 6. What was removed relative to the development version

* the inner (throat) ALE taper (`ale_inner`, argv[10], its switch logic, `INNER_T_HOLD`,
  `THROAT_THETA_MAX`): used only by exploratory t = 2000 runs whose frozen throat was later
  found to be a resolution artefact; not in the paper;
* the W-form momentum constraint (`qform`, argv[18], `QW_R0`): a test that changed nothing
  (~1%); not in the paper;
* the adaptive zone-B factor (`refB_Fmax`, argv[17], the "raising F" logic): exploratory;
  the paper's refined runs all use a fixed F (8 or 16);
* the unused `*_prev` copies of the fields, the per-step debugging print of Q and β for the
  first 20 steps, and the inner-taper columns of the progress line.

Added: the linear probe (`--probe_amp`, `--probe_r`), ported unchanged from the
run-differencing binary used for the probe analysis (it was never in the development
version); off by default, so the arithmetic of every other run is untouched. The sphere
recorder (`--rec_spheres`, `--rec_every`) is also new; it only reads the state
(`VALIDATION.md`).

Changed: the ingoing initial data (`--ingoing=1`). The development code solved the
Hamiltonian constraint for Θ with Q = 0 and with P, Q taken from the grid point at all
four RK4 stages, which left a constraint violation at t = 0 that did not shrink with
resolution. The public code iterates the Θ and Q solves to a fixed point and evaluates P
analytically (and Q, a_r, K by linear interpolation) at the stage radii; the t = 0 residual
is now second order (VALIDATION.md, "Ingoing initial data"). Ingoing runs therefore differ
from the development runs at the level of the old initial-data error (at t = 10 in the
r_max = 40 box at N = 1600: ≤ 4e-4 of the maximum of P, 2e-4 of K, 3e-5 of Θ); every
white-hole result of the paper was recomputed
with the corrected data. Time-symmetric data were unchanged and bit-identical at that point.
Changed afterwards, for every kind of data: the six scheme points of `VALIDATION.md`,
"Scheme corrections" (ICN corrector on the t+dt mesh also without refinement, lapse by
direct integration, horizon diagnostic, finiteness check after every step, zone-A
acquisition floor, full scalar outer boundary condition), and last the exact outer boundary
condition for the incoming spin-0 field (§2; `VALIDATION.md`, "Final rerun"); all
production runs of the paper were redone with the code shipped here. Also new:
non-numeric option values (including the entries of `--times` and `--rec_spheres`) and
`--ale` outside {0, 1} are rejected, a non-finite field stops the run instead of letting
NaNs propagate, and a FATAL stop no longer loses the buffered tail of `phi.dat` and
`psi_times.dat`.
