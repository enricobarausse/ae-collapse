# paper_checks: SymPy equation audit of the paper

`check_equations.py` checks every displayed and labelled equation, every inline relation that
carries physics, and every number that can be checked by arithmetic, in the paper
*Gravitational collapse in Lorentz-violating gravity and the spin-0 mode* (Albertini, Vaglio,
Barausse). Paper source version: `main.tex` of the paper's LaTeX source (arXiv v1); all line
numbers below refer to that file.
Sections follow the order of the paper. Each section is one function that quotes the LaTeX
verbatim with its `main.tex` line number ("..." marks omitted text; caption and long-paragraph
quotes are marked "abridged"), says what is checked and how, and prints one status line per
check. The script reads nothing from disk (neither the paper nor the simulation code need be
present) and writes no files.

(`check_equations.wl`, `check_equations.nb` and `README_mathematica.md` are the independent
Mathematica audit; `CHANGELOG.md` summarises what both audits found and what was corrected.)

## How to run

```
python3 check_equations.py            # full audit, about 1.5 min (the three slow derivations run in parallel worker processes)
python3 check_equations.py --quick    # skips S04, S10 and S39 (the slow derivations, reported as INFO), about 1 min
```

Requirements: Python 3 with SymPy (tested with 1.14) and NumPy (tested with 2.5); nothing
else. The output ends with a SUMMARY block that lists every FAIL and CAVEAT and gives the
counts by status; the exit code is 0 if no check FAILs and 1 otherwise.

Status codes:

* **PASS**: derived from first principles, or proven as an exact identity.
* **FAIL**: the statement as printed is wrong, or inconsistent with another statement of the paper.
* **CAVEAT**: correct only under an assumption the paper does not state, worded too strongly, or a
  minor misprint.
* **DATA**: a simulation number. Only its arithmetic and its consistency with other numbers in
  the paper are checked.
* **INFO**: not a check of a statement of the paper: a literature input or a definition (not checked here), a slow
  derivation skipped in `--quick` mode, or a passed check of the simulation code that the paper does not state (the
  "(code check)" lines of S35; a failed code check is FAIL). The Mathematica audit reports its passed code checks as PASS.

Some PASS identities are verified in exact rational arithmetic at random points (S02, S10); a few PASS lines are
statements reproduced numerically (S36, S51, S57).

## Refreshing the numbers

Every simulation number quoted in the paper is registered once, in the `DATA` table at the top
of the script (`_d('key', value, line, 'what')`), with the `main.tex` line it comes from; the checks read it
through `dval(key)` (the value) and `dline(key)` (the line). Measured values that a check needs but the paper does not print
(records of the figure runs, made by the `release/README.md` recipes and `release/analysis` scripts) are in the
`*_UNPRINTED` tables right after the DATA table, each with a comment saying where it comes from. The
DATA checks read from that table only. When a number changes in the paper, update the table
entry and the line number; section S00 prints the whole table at the start of every run.

Conventions are given in the script docstring: signature (−,+,+,+), G = c = 1, the normalisation
G_ab = T^æ_ab + T^ψ_ab with T^ψ = ∂ψ∂ψ − ½g(∂ψ)² (stated in the paper at L313-316), the aether
u_a = −N∂_a t, ∇_a u_b = −K_ab − u_a a_b, K^r_r = K/3 + Q, K^θ_θ = K/3 − Q/2. The code equations
used by the script were transcribed by hand from `release/collapse.cpp` (the file in this repository)
and checked against it line by line (the docstring lists the line ranges); the
script never reads the code. Lines of the code are always written `collapse.cpp:L…` (for example
`collapse.cpp:L875-886`, the trK and a_r right-hand sides in `compute_RHS`); a bare `L…` is always a line of
`main.tex`.

## What is derived, not restated

* **S08b, S08c: Eqs. (evolution) and (radial)**: every printed equation equals the code system (and S10 tests the printed
  forms directly against the field equations of the action, with negative controls); the Bonetti–Barausse 3+1
  equations, mapped with (α, β, λ) = (c14, c13, c2) and f = αa² and reduced to spherical symmetry (their K_ij has the
  opposite sign), give Eq. (hamiltonian) and the Q equation of Eq. (radial).
* **S10: the whole GEJ system in `collapse.cpp` follows from Eq. (action).**
  1. Reduced action, varied with respect to N, β, γ_rr, Θ, ψ and an aether tilt η.
  2. The six field equations are solved for the second-order jets at random rational points,
     with general c1…c4.
  3. The code's P, trK and a_r evolution equations, its Hamiltonian constraint (theta_C) and its
     momentum constraint (compute_Q) are reproduced with residual exactly 0.
  4. A negative control, with 3(1−c13)Q² changed to 3Q², fails as it should.
  S11 and S25 show that the displayed Eqs. (hamiltonian) and (evolution) (with S08b) are these equations.
* **S04: v2², v1², v0² from Eq. (speeds)**, from the second-order action about Minkowski: the
  plane-wave Euler–Lagrange equations of all fields, gauge-fixed, and the gcd of the maximal
  minors of the Fourier matrix.
* **S39: the Misner–Sharp identity** behind "the time derivative of the Misner–Sharp mass, evaluated at fixed areal radius,
  equals the flux of the stress-energy tensor through the sphere" (the identity is not printed), proven for a completely
  general 2-metric g_AB(t,r), with a negative control (the 4-d trace in place of G^(2)).
* **S12: constraint propagation at the centre**, regular Taylor series with general c_i and ᾱ.
* **S17: the null expansions θ±** as Θ/2 times the covariant expansions of u ± e, and the
  Misner–Sharp mass m = Θ(1 − θ₊θ₋)/2.
* **S23, S40: ∂_Θ m and the luminosity formula**, from the identity by frame algebra with a
  generic stress tensor (ρ, j, p); θ± = B ± A and the dominant-energy-condition argument. S40 also checks the wording
  "Using the field equations G_ab = T_ab, one obtains the luminosity L": sufficient, and the full tensor equation is needed.
* **S41: the precursor shear (not printed in the paper)**, Q = −2Mt/((1−c13)Θ³) + O(t²),
  derived from the code equations for any time-symmetric data, together with m − M ∝ t²/Θ³; S48
  uses it for "E slightly negative before the front".
* **S43: the Coulomb term**, from the Hamiltonian and momentum constraints with Θ′ kept general:
  with the Θ′ factor of Eq. (dm-coulomb) the linear statement is exact.
* **S25–S33: the spin-0 sector**: Eq. (evolution) as a first-order system, principal part, eigenvalues, left and right eigenvectors,
  the couplings that keep the characteristic variables from being exact (to Θ through a_r/Θ, to ψ through P, both
  falling off in the wave zone, where a_r ~ P ~ 1/r do not vanish), the static frame and the rigid translation r → r + c(t)
  (β = ċ, γ_rr = 1, c ≃ −t, ∂_t at fixed areal radius = ∂_t − β∂_r; exact in flat space, noted as INFO on the collapse
  background), linearisation, wave equation, the asymptotic solution with F(u) (the 1/r² part builds up and remains),
  w_out = O(1/r), w_in = O(1/r²) and its rate.
* **S35, S36: the outer boundary conditions, Eq. (bc-win)**: parities at r = 0 preserved; the entering fields of an
  outgoing wave, w_in = −b_N F/r² and P + ∂_rψ = −g/r² (one power of r below the radiation, → 0 at infinity); for the
  linear spin-0 system ∂_t w_in + b_N trK/r = (b_N/(v₀r))(∂_t + v₀∂_r)(r trK) identically (a_N b_N = v₀²), so the
  outgoing-characteristic condition is equivalent to ∂_t w_in = −b_N trK/r; Eq. (bc-win) exact for an outgoing linear
  wave on flat space with constant lapse and shift, and equal to the code's G; the scalar condition (u+e)(Θψ) = 0 and,
  in flat space, ∂_t(P + ∂_rψ) + P/r = ∂_t[r⁻¹(∂_t + ∂_r)(rψ)] identically (equivalent given a quiescent start);
  the reflection with w_in set to zero computed from the linear theory of the paper's test pulse (3.48 % / 2.33 %).
  The projector algebra of the code update, the scalar update through ∂_t P and the outflow gate are kept as INFO
  code checks (the paper does not state them).
* **S49, S59–S61: Sec. V.E (sec:pn_comparison) and App. B**: rows 1 and 3 of D₀ from the printed reduced
  Taherasghari–Will equations, det D₀, s₋₁ by elimination, Γ(x) exactly, the truncated
  redefinition T₁ = D₀⁻¹D₁ reproducing Eq. (N_transform), s₁ and κ₁, and the single-field example.
* **S51: the ingoing initial data**, a direct numerical solve of the t = 0 constraints showing
  the violation of a single Q = 0 solve and its removal by the iterated solve used in the code.

## History

The chronology of the audit runs, and of the changes to the paper and to the code that they
prompted, is in `CHANGELOG.md`.

## Sections

| Sec | main.tex line | What is checked | Statuses |
|---|---|---|---|
| S00 | — | DATA table: simulation numbers quoted in the paper | INFO |
| S01 | 128 | Introduction: c1 + c3 must vanish (GW170817); Hořava dispersion relation | PASS×2 |
| S02 | 236 | Eq. (action) and the tensor M^{ab}_{mn} | PASS |
| S03 | 250 | c_ij combinations; hypersurface-orthogonal aether; c13 = c14 = 0.7 | PASS×2 |
| S04 | 259 | Eq. (speeds): squared speeds of the spin-2, spin-1, spin-0 modes | PASS×3 (full run; INFO in --quick) |
| S05 | 267 | Sign of the spin-0 energy ~ c14(2 − c14); positive only if 0 < c14 < 2 | INFO, PASS |
| S06 | 277 | GW170817 bound | PASS |
| S07 | 309 | Eq. (lapse): kinematic identity a_i = D_i ln N | PASS |
| S08 | 317 | Eq. (metric); evolved variables; Q; d_r ln N = a_r; d_r β = N(Q + trK/3); β(0) = 0 | PASS×2 |
| S08b | 342 | Eqs. (evolution), (radial): printed GEJ system = code system (and, in S10, = field equations of the action); negative controls | PASS×3 |
| S08c | 339 | Bonetti–Barausse 3+1 equations with (α, β, λ) = (c14, c13, c2): action, Hamiltonian and momentum constraints in spherical symmetry | PASS×2 |
| S09 | 377 | a_r fixed on each slice by the lapse: D_a u_b = −K_ab − u_a a_b contains no d_t N (derived) | PASS |
| S10 | 378 | The GEJ system evolved in release/collapse.cpp follows from Eq. (action) | PASS×11, INFO (full run; INFO×2 in --quick) |
| S11 | 367 | Eq. (hamiltonian): Θ C = Θ × the code's Hamiltonian constraint | PASS |
| S12 | 384 | Eq. (constraint-prop): propagation of the Hamiltonian constraint at r = 0 | PASS |
| S13 | 395 | Singularity-avoiding slicing: ln N, monotonic N, sources ∝ N, the bound on \|β\| (from Eq. (radial) and β(0) = 0) | PASS×3, DATA |
| S14 | 438 | Discretization; grid falls inward (~1 areal unit per unit time on average); ALE map; Fig. 1 | PASS×3, DATA |
| S15 | 476 | Couplings: c2 = −0.5119403 gives v0 = 1; c2 = 0 gives v0 = 1.267 | PASS×3 |
| S16 | 482 | Eq. (id): time-symmetric shell (GEJ profile × r²/100, peak ≃ A0) | PASS×3 |
| S17 | 497 | Eqs. (expansions), (ms): null expansions; Misner–Sharp mass; M_hor = Θ_AH/2 | PASS×4, DATA |
| S18 | 531 | Convergence order; factor 4 for p = 2; dt = 0.1 h halved with dr keeps the error O(h²) (symbolic, dt factor from the DATA table) | PASS×2 |
| S19 | 548 | Θ C at the centre (parities) | PASS |
| S20 | 557 | Scattering: grid spacings, Table 1, black-hole threshold; App. A (incl. the ladder parameters of the conv-fields / conv-constraint captions vs the run record) | DATA×6 |
| S21 | 630 | Collapse: orders, Table 2, Θ_AH ladder and Richardson limit, m on Θ = 250 (+0.07 %, ∝ Θ⁻³, extrapolation ≈ 0.01 %: p = 3 fit −0.0101 %, free p +0.0037 %; Θ = 250 leaves the fixed grid at t ≃ 20); App. A (incl. the Table 2 / Fig. conv-collapse caption parameters and the N_r = 51201 moving mesh vs the run record) | DATA×10 |
| S22 | 679 | Masses of Fig. 4: M_hor = 6.35, M_ADM = 3.7111, −2.64 | DATA×2 |
| S23 | 689 | Eq. (dmdTheta); θ± = B ± A; m decreases outward only where the DEC fails; aether ≤ −2.64 | PASS×5, DATA |
| S24 | 711 | Penrose inequality from Eq. (dmdTheta) under the DEC; Eq. (mratio); area theorem; DEC vs WEC | PASS×2, DATA×2, INFO |
| S25 | 755 | Eq. (evolution) as a first-order system: trK, a_r lines = code; ψ, P, Θ equations free of ∂_r trK, ∂_r a_r (block-triangular symbol) | PASS×2 |
| S26 | 762 | Eq. (principal): principal part of the (trK, a_r) pair | PASS |
| S27 | 768 | Eq. (aNbN): a_N, b_N and a_N b_N = v0² | PASS×2 |
| S28 | 776 | Advection with the aether, eigenvalues and characteristic speeds; Franchini et al. speeds with (α, β, λ) = (c14, c13, c2) | PASS×3 |
| S29 | 791 | Eq. (wout); exact far out: the only couplings, to Θ (∝ a_r/Θ) and ψ (∝ P), fall off; scalar pair P ∓ ∂_rψ | PASS×3 |
| S30 | 809 | Static frame and the rigid translation r → r + c(t) (β = ċ, c ≃ −t, ∂_t at fixed areal radius = ∂_t − β∂_r, flat space); linearised (trK, a_r) system | PASS×3 |
| S31 | 823 | Spherical wave equation for trK | PASS |
| S32 | 826 | Eq. (asymptotic): outgoing solution | PASS×2 |
| S33 | 834 | Eq. (wout-asym): the 1/r² part builds up and remains; waveform; ∂_t w_in = −b_N trK/r | PASS×2 |
| S34 | 911 | Eq. (retarded) (Extraction paragraph): tortoise coordinate, slowing near the black hole | PASS |
| S35 | 849 | Boundary conditions: parities; entering fields −b_N F/r² and P + ∂_rψ = −g/r²; (∂_t + v0∂_r)(r trK) = 0 ⇔ ∂_t w_in = −b_N trK/r; Eq. (bc-win); scalar ∂_t(P + ∂_rψ) = −P/r; code checks (INFO) | PASS×7, INFO×4 |
| S36 | 877 | Boundary-condition tests: reflection with w_in set to zero from the linear theory (3.5 %, 2.3 %); reflection numbers vs the release VALIDATION record | PASS×2, DATA |
| S37 | 886 | Results: the two runs (v0 = 1 and 1.267; the two steps); Eq. (wbar): normalised Gaussian running mean (constant preserved, δ of a constant = 0, high-pass transfer), against release/analysis/spin0_common.py gaussian_running_mean (3σ truncation, renormalised ends) and the σ = 10 of the figure data; extraction of the pulse | PASS×3, DATA |
| S38 | 936 | Propagation speed: the probe ((ṙ + β)/N = v0 exactly for the outgoing characteristic, S28; medians 0.012 % / 0.015 % → 0.01 % / 0.02 %; max \|v − v0\|/v0 for t ≥ 40 = 2.96e-4 / 1.89e-4 < 3×10⁻⁴, larger early scatter, vs data/probe_speed.npz); the luminosity clock (Θ_* = Eq. (retarded), lag fit v = 1.2406 → 1.241, −2.07 % → 2 %, rms 0.31, scalar speed 1 excluded by 56; burst = narrow pulse of w_out within 0.16 in u; scalar speed 1; the burst's u = 9.7–14.7 = "u ≃ 10–15") | PASS×2, DATA×5 |
| S39 | 997 | Opening of "Mass conservation and radiated energy" (the two steps → S39–S48); Misner–Sharp flux identity (exact): time derivative at fixed Θ = flux of T through the sphere | PASS×2 (full run; INFO in --quick), INFO |
| S40 | 1004 | Eq. (lum): luminosity through fixed-Θ spheres; "using the field equations G_ab = T_ab" | PASS×4 |
| S41 | 1120 | The precursor shear ahead of the front, Q = −2Mt/((1−c13)Θ³), and m − M ∝ t²/Θ³ (derivation behind S48's "E slightly negative") | PASS×3 |
| S42 | 1016 | Flux identity: residuals and flux magnitudes | DATA |
| S43 | 1052 | Eq. (dm-coulomb): Coulomb term from the Hamiltonian constraint; −¾c13Q² | PASS×3 |
| S44 | 1037 | The "confusing" luminosity: swing of m grows with Θ (7.05 → 18.76), m settles above M_ADM, ∫L du = −1.24 … −1.09 < 0 on every sphere (aether_coulomb.npz); the aether Coulomb field a_r ∝ Θ⁻² and its name | PASS, DATA×2, INFO |
| S45 | 1067 | Eq. (mc): corrected mass and radiated energy | INFO |
| S46 | 1075 | Eq. (Mae): Einstein–aether energy and its spherical reduction | PASS, INFO |
| S47 | 1098 | Eq. (Einf): radiated energy and its split (M_ADM = 3.708 of the c2 = 0 run) | DATA×3, INFO |
| S48 | 1113 | Fig. fluxinf in time: corrected rate vs the uncorrected L of Fig. lum(a) (25× at Θ = 200, 11× at Θ = 120, 28× between the panel peaks; not printed: 5.3× / 8.3× on the inner spheres Θ = 75 / 100); E nearly the same on all spheres once the pulse has passed (spread ≤ 0.044 for u ≥ 110, 6.8 % for u = 200–250; up to 0.43 while the pulse crosses, excluded by the qualifier); E_inf(u) settles in 0.522–0.540 ⊂ Eq. (Einf); E < 0 before the front; overshoot | DATA×3, PASS, INFO |
| S49 | 1140 | Sec. V.E (sec:pn_comparison): v_L, Eqs. (pn_first_surviving), (G_exact), (inverse_propagator) | PASS×6, INFO |
| S50 | 1255 | White-hole couplings | PASS |
| S51 | 1257 | Eq. (wh-id): ingoing shell; initial-data constraint solve | PASS×2 |
| S52 | 1265 | θ₋ and the anti-trapped condition | PASS, INFO |
| S53 | 1272 | White hole or black hole: numbers; split of the anti-trapped region at t ≃ 36 (θ₋ = 0 zeros of box40_N51200); thin layer 0.9 → 0.06 on the refined grids, agreement within 20 % (whconv ladder: all-grid spread 8/18/11 % at t = 40/64/80), depth grows with resolution | DATA×3 |
| S54 | 1291 | Three sharp features: scalings of the unresolved inner shell; dr^0.8, dr^−0.5, dr^−1 vs the rerun fits; refinement factors 8 and 8 or 16, hold t = 20–30.5 vs the release recipe; anti-trapped layer < 1 unit wide (Fig. wh caption); outgoing ray 0.3 / 0.1 behind the pulse at t = 34 / 40 (gaps 0.302 / 0.091); inner ingoing front (t ≃ 7, jump ≃ 0.3, ṙ = −β − N, N ≃ 0.13, Θ ≃ 3.2 near the inner θ₋ = 0 edge, slope ratios at t = 24 / 40, t ≃ 44 / 54) and the Fig. sharp caption vs fig_sharp_features.py | DATA×12 |
| S55 | 1343 | No unique continuation; frozen-centre artefact | DATA |
| S56 | 1358 | Other initial data | DATA |
| S57 | 1377 | Comparison with Akhoury–Garfinkle–Gupta | DATA, PASS |
| S58 | 1436 | Conclusion (statements against the body of the paper; the far-zone flux sentence = the introduction; lapse collapse and freezing inside the horizon = Sec. method L410-412, the universal-horizon reading and Franchini et al. = L419-421 and the introduction L169-172) | INFO×2 |
| S59 | 1614 | App. B.1: single-field example (Eqs. simple_L, Tp, Lprime, Lprime_p1) | PASS×5 |
| S60 | 1670 | App. B.2: B̃^rr_LTF, reduction, D₀ rows, det D₀, v_L, s₋₁, the 1PN term | PASS×7, INFO×2 |
| S61 | 1787 | App. B.3: T₁, Eq. (M1), v_T, s₁, κ₁, Γ series, det identity; the summary (no v_L in s₁ through 1PN, v_L at 2PN) | PASS×10 |

## SymPy and Mathematica sections side by side

The two audits are independent and number their sections separately. The table pairs them by topic (a dash: the
topic is checked in one audit only). Status vocabularies also differ: the SymPy audit uses PASS / DATA / INFO /
CAVEAT / FAIL (its passed code checks and literature inputs are INFO, its simulation numbers DATA), the Mathematica
audit PASS / SKIP / CAVEAT / FAIL (its passed code checks and arithmetic on simulation numbers are PASS, literature
inputs and bare simulation numbers SKIP).

| Topic (main.tex label or section) | SymPy (`check_equations.py`) | Mathematica (`check_equations.wl`) |
|---|---|---|
| DATA table / snippet cross-check against main.tex | S00 | 0b, S81 |
| Introduction: GW170817 sentence, Hořava dispersion | S01, S06 | S01, S02, S09 |
| eq:action, the tensor M | S02 | S03 |
| c_ij combinations, production couplings | S03, S15 | S04, S21 |
| eq:speeds | S04 | S05, S06 |
| spin-0 energy sign, 0 < c14 < 2 | S05 | S07, S08 |
| eq:lapse | S07 | S10 |
| eq:metric, evolved variables, Q | S08 | S11, S12, S13 |
| eq:evolution, eq:radial (printed system = code system) | S08b | S13b |
| Bonetti–Barausse 3+1 equations | S08c | S13c |
| ∇_a u_b = −K_ab − u_a a_b (a_r fixed by the lapse) | S09 | S14 |
| code system from eq:action (derived) | S10 | — (transcribed, not derived) |
| eq:hamiltonian | S11 | S15 |
| eq:constraint-prop | S12 | S16 |
| singularity-avoiding slicing, \|β\| bound | S13 | S17, S18 |
| discretisation, moving mesh, Fig. mesh | S14 | S19, S20 |
| eq:id | S16 | S22 |
| eq:expansions, eq:ms, M_hor = Θ_AH/2 | S17 | S23 |
| convergence order, dt = 0.1 h, factor 4 | S18 | S24 |
| Θ C at the centre | S19 | S25 |
| scattering convergence, threshold, App. A ladders | S20 | S26–S28, S30c |
| collapse convergence, Θ_AH ladder, M_ADM drift | S21 | S29, S30, S30b, S30c |
| masses, eq:mratio | S22, S24 | S31 |
| eq:dmdTheta, θ± = B ± A, energy conditions | S23 | S32–S34 |
| Penrose inequality, area theorem | S24 | S35 |
| first-order system, block-triangular symbol | S25 | S36 |
| eq:principal | S26 | S37 |
| eq:aNbN | S27 | S38 |
| eigenvalues, advection with the aether | S28 | S39 |
| eq:wout, characteristic variables | S29 | S40, S41 |
| linearised system far from the source | S30 | S42 |
| wave equation | S31 | S43 |
| eq:asymptotic | S32 | S44 |
| eq:wout-asym | S33 | S45 |
| eq:retarded | S34 | S46 |
| outer boundary, eq:bc-win, scalar condition, code update and gate | S35 | S47–S52 |
| boundary-condition tests (reflection) | S36 | S53 |
| Results opening, eq:wbar, extraction | S37 | S54, S54a |
| probe and luminosity clock (propagation speed) | S38 | S54b–S54d |
| Misner–Sharp flux identity | S39 | S55 |
| eq:lum | S40 | S56, S56b |
| precursor shear ahead of the front (not printed) | S41 | — |
| flux-identity numbers | S42 | S58 |
| aether stress tensor → constraints | — | S59 |
| eq:dm-coulomb, −¾ c13 Q² | S43 | S61, S67 |
| confusing luminosity, aether Coulomb field, Θ² a_r | S44 | S60, S60b, S62 |
| eq:mc | S45 | S63 |
| eq:Mae | S46 | S64 |
| eq:Einf | S47 | S65, S66 |
| Fig. fluxinf in time, E < 0 ahead of the front | S48 | S67b |
| sec:pn_comparison (Sec. V.E) | S49 | S68 |
| white-hole couplings | S50 | S69 |
| eq:wh-id, ingoing initial data | S51 | S70 |
| θ₋ < 0, anti-trapped | S52 | S71 |
| white-hole timeline, split of the anti-trapped region | S53 | S72, S72b |
| sharp features (escaping pulse, inner shell, inner front) | S54 | S73, S73b–S73e |
| frozen centre, other initial data | S55, S56 | S74 |
| Akhoury–Garfinkle–Gupta comparison | S57 | S75 |
| Conclusion | S58 | S76 |
| App. B.1 single-field example | S59 | S77 |
| App. B.2 spherical reduction, D₀, s₋₁ | S60 | S78, S79 |
| App. B.3 first-order redefinition, s₁, κ₁ | S61 | S80 |

## Remaining CAVEATs

None against the arXiv v1 source. Earlier findings and how the paper resolved them are listed in `CHANGELOG.md`.

## Last run

Full run of 2026-11-13 against the arXiv v1 `main.tex`:

```
============================================================================================
SUMMARY
============================================================================================
  counts: DATA=58, INFO=20, PASS=140
  paper main.tex of the LaTeX source of arXiv v1; run time 86 s (FULL mode)
```
