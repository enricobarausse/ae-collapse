# Mathematica equation checks for the paper

An independent Wolfram Language audit of the equations and numbers of "Gravitational collapse in Lorentz-violating
gravity and the spin-0 mode" (`main.tex`, the arXiv source). It complements the SymPy audit in this directory
(`check_equations.py`, `README.md`): the two are independent and share no code.

Files (this directory):

* `check_equations.wl` - Wolfram Language script, one section per equation / inline relation / piece of arithmetic
  of `main.tex` in paper order (section 0a = shared setup, 0b = the DATA association, S01-S81 (with S13b, S13c, S30b, S30c, S54a, S54b, S54c, S54d, S56b, S60b, S67b, S72b, S73b, S73c, S73d, S73e; there is no S57) = the checks, last =
  summary). Each section has the LaTeX label and line numbers, the verbatim LaTeX, a comment on what is checked and
  how, the check, and exactly one `PASS:` / `FAIL:` / `CAVEAT:` / `SKIP:` line with the residual.
* `check_equations.nb` - the same code as a notebook (one Input cell per section), generated from the `.wl` by
  `make_notebook.wl` (run `wolframscript -file make_notebook.wl` and commit the regenerated `.nb` after every edit of the `.wl`).
* `make_notebook.wl` - the notebook generator.
* `CHANGELOG.md` - what the two audits (SymPy and Mathematica) found and what was corrected.
* `check_equations.py` and `README.md` (the SymPy audit) are independent of these files.

## How to run

Requirements: Mathematica 14.2 with `wolframscript` (no packages beyond the standard distribution).

```
wolframscript -file check_equations.wl [path/to/main.tex]
```
(or `math -script check_equations.wl`, or `Get["check_equations.wl"]` in a kernel). Run time about 40 s, most of it
the geometry setup of section 0a and the plane-wave derivation S06. The path of `main.tex` is used only by the snippet cross-check S81: first
command-line argument, else the environment variable `PAPER_TEX`, else `main.tex` in the current directory; if the
file is absent S81 is SKIP and everything else runs unchanged (exit code 0). Run as a script the exit code is 1 if
any section FAILs, else 0 (CAVEATs do not change it) (a stale snippet gives `FAIL: S81` and exit 1). In the front end open `check_equations.nb`
and use Evaluation > Evaluate Notebook (the `Exit` at the end is skipped inside a front end; set `PAPER_TEX` or put
`main.tex` in the kernel's current directory for S81).

The `.nb` was tested by importing all 99 Input cells and evaluating them in order in a kernel (`ToExpression`,
with `PAPER_TEX` set): status lines identical to the script run below, 93 PASS / 0 FAIL / 0 CAVEAT / 3 SKIP (only
the printed run time differs).

## What is transcribed and what is derived

**Transcribed from the code, not derived.** The evolution equations of P, trK, a_r and Theta (`codeRHS`), the
Hamiltonian-constraint diagnostic (`thetaCcode`), the momentum-constraint ODE (`momConCode`), the coefficients
a_N, b_N, the outer-boundary update, the initial-data profile and the ALE taper are transcribed from the simulation
code `release/collapse.cpp` of this repository (code lines are always written `collapse.cpp:L…`, bare
`L…` and the "main.tex lines" column are lines of `main.tex`): header equations collapse.cpp:L44-76; `compute_RHS`
collapse.cpp:L817-930 (interior loop L844-890: P L872-873, trK L875-881, a_r L883-886, Theta L888-889; outer boundary
L892-929: spin-0 pair L901-915 with G at L911, scalar pair L917-928); outer-boundary gate in `ICN_step`
collapse.cpp:L1001-1030; `theta_C` collapse.cpp:L1305-1315; `compute_Q` collapse.cpp:L699-735; couplings
`c13`, `c14`, `c123`, `den_K`, `den_a` and `spin0_aN/bN` in `main` collapse.cpp:L1607-1616; `shell_psi`/`shell_dpsi`
collapse.cpp:L163-171; `set_initial_data` collapse.cpp:L1208-1298 (ingoing P L1252-1261); `quintic_taper`
collapse.cpp:L325-332. The transcriptions were checked against that file line by line. This script does not derive those
equations from the action; the SymPy audit (`check_equations.py`, its S04/S10) does.

**Names of the fields in the script.** The script's symbols are not the paper's letters: `Nf` = lapse N, `betaf` =
shift beta, `Thf` = areal radius Theta, `psi` = psi, `Pf` = P, `trk` = trK, `Qf` = Q, `arf` = a_r (collapse.cpp:
`NN`, `NR`, `Theta`, `Psi`, `P`, `trK`, `Q`, `Acr`); `ThfDot`, `psiDot`, `PfDot`, `trkDot`, `arfDot` are their code time
derivatives at fixed r. The legend is repeated at their definition in section 0a.

**Differences from the SymPy audit.** The precursor shear ahead of the front, `Q = -2Mt/((1-c13)Theta^3)` (not printed
in the paper; behind "E slightly negative"), is derived only in the SymPy audit (its S41); this script has no
counterpart. Checks of the simulation code that are not statements of the paper (the code update, the scalar pair and
the outflow gate, S48-S52) are PASS here, while the SymPy audit reports the same kind of check as INFO (its S35). The
section numbers of the two audits differ (each follows the paper with its own sub-numbering); the table "SymPy and
Mathematica sections side by side" below (also in `README.md`) pairs them.

**Derived here from first principles** (symbolically, for the metric `ds^2 = -N^2dt^2 + (dr+beta dt)^2 +
Theta^2 dOmega^2` with `u_a = -N d_a t`):
* Christoffels, Einstein tensor, `nabla u`, `a_a`, covariant `K_ab` (S10-S14): eq:lapse, eq:metric,
  `beta' = N(Q+trK/3)`, `nabla_a u_b = -K_ab - u_a a_b`, `trK = -div u`, the Theta evolution equation.
* eq:action: `M` as a 4-index array, `J^a_m`, `M nabla u nabla u = c13 K.K + c2 K^2 - c14 a^2` for a
  hypersurface-orthogonal aether (S03); the three mode speeds of eq:speeds from the second-order action (S06).
* The aether stress tensor from arXiv:1501.07274 eq. (2.34), projected on `uu` and `ue` (S59): with kinematic
  definitions only, `G_uu = rho_ae + rho_psi` IS the transcribed `theta_C` = the printed eq:hamiltonian (S15), and
  `G_ue = T_ue` IS the transcribed `compute_Q` (this one needs the a_r evolution equation).
* The normalisation `G_ab = kappa T_ab` actually used, kappa = 1 (S33).
* eq:constraint-prop (with `\bar\alpha`) by Taylor expansion at the centre (S16); the parities at r = 0 preserved by
  the equations (S47); `Theta C -> (Theta'(0)^2-1)/2` (S25); the `|beta|` bound (S18).
* eq:expansions as the rescaled expansions of the radial null congruences, eq:ms, `M_hor = Theta_AH/2` (S23);
  eq:dmdTheta, `theta_pm = B pm A` and the dominant-energy-condition argument (S32-S34); the Misner-Sharp identity and
  eq:lum (S55, S56), and the wording "using the field equations G_ab = T_ab" for eq:lum (S56b); eq:dm-coulomb with its `Theta'` factor and its `O(M/Theta)` remainder (S61, S67).
* The principal symbol of the full first-order system: eq:evolution (trK, a_r lines) = code, block-triangular (S36); the only couplings that keep
  the characteristic variables from being exact (Theta' entry `~ a_r/Theta`, psi' entry `~ P`, both falling off for
  `a_r ~ P ~ 1/r`), scalar combinations `P -+ psi'` (S41); the linearised trK, a_r system far out (inline, main.tex L818-823, no label), the wave equation,
  eq:asymptotic, eq:wout-asym with the right eigenvector and `d_t w_in` (S42-S45).
* The outer boundary (S48-S53): projector algebra of the entering-field condition (main.tex L851-853) and its equality with the code's update for the
  spin-0 pair; the boundary datum `(d_t - beta d_r) w_in = -N b_N trK/Theta` exact for the outgoing wave in flat space
  with constant lapse and shift; the scalar datum `(u+e)(Theta psi) = 0`; the outflow conditions; the reflection of
  the frozen entering field from the linear solution of the paper's test pulse.
* The post-Newtonian argument (S68) and App. B (S77-S80): `v_L`, `v_T`, Gamma exact, the inverse propagator, the
  toy model, rows 1 and 3 of `D0` from the printed reduced Eqs. (4.5a), (4.6), `det D0`, `s_(-1)`, `T1 = D0^-1 D1` to
  O(x) = eq:N_transform, eq:M1, `s_1`, `kappa_1`.

**Literature inputs, not re-derived:** the sign of the spin-0 mode energy (S07, SKIP; its consequence
`0 < c14 < 2` is S08), the coefficient `c14/(8 pi G)` of eq:Mae (S64 checks only its spherical reduction), the
observational value `1e-15` of the GW170817 bound (S09), the reduced equations of Taherasghari-Will (only their
printed spherical reduction is used; row 2 of `D0` is not checked).

**Arithmetic on the paper's own numbers** (DATA, below): S17, S19, S20, S24, S26-S31, S35, S53, S54 (S54b-S54d read the probe and the luminosity clock against the
figure-run records probe.log, probe_speed.npz, lum_clock.log and invariant_wout_dm.npz), S60b and S67b (the swing of m, m above
M_ADM and int L du < 0, the corrected rate against the uncorrected luminosity, E on all spheres and E_inf(u), against the
figure-run records aether_coulomb.npz, spin0_flux_infinity.npz, luminosity_clock.npz and energy.log), S65, S66,
S72-S76 (S72b reads the split of the anti-trapped region and the width and depth of the thin layer on the escaping
pulse against the theta_- = 0 zeros of the runs; S73b also fits the rerun table of the white-hole inner shell; S73c
compares the refinement factors and hold window with the release recipe and the width of the outermost anti-trapped
layer with its measurement on the dr = 0.0008 run; S73d reads the inner ingoing front and the Fig. sharp caption
against the record of fig_sharp_features.py and the front tracked on box40_N51200; S73e reads the convergence of
the outgoing ray on the escaping pulse against the fig_white_hole.py record).

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

## DATA: simulation numbers in one place

Every simulation number quoted in the paper sits in the association `DATA` at the top of the script (section
0b), as `key -> <|"v" -> value, "L" -> "main.tex line(s)", "tex" -> "snippet as printed"|>` (152 entries). The
arithmetic-consistency sections read from it (`dat["key"]`), so when a run is redone the number is changed in one
place. Section S81 cross-checks every snippet against `main.tex` (whitespace, backslashes, braces, `$` and `~`
ignored): FAIL with the stale keys listed if any snippet is not found at its quoted line(s), SKIP if `main.tex` is
not found. Claims of the text whose wording matters are DATA entries too: the numbers `dmClaim` = coefficient 1/2
(S33), `coulombClaim` = remainder order 1 and `q2Claim` = -3/4 (S67), and the quoted relations and sentences (`"v" ->
Null`), each read or quoted by one section: `dispersion` S01 (quoted in its comment), `gwCombination` S02,
`evolutionEqs` S13b, `bonettiMap` S13c, `decClaim` S34, `penroseClaim` S35, `charExact` S41, `rigidShift` S42,
`outgoingUnchanged` S48 and S51, `noIncoming`, `outgoingChar`, `binRate` S49, `scalarRate` S50, `resultsSteps` S54,
`wbarEq` S54a, `probeExact` S54b, `lumDef`, `lagLine`, `burstSpin0` S54c, `fluxOpening`, `fluxHayward` S55, `lumRoute`
S56b, `lumConfusing`, `coulombName` S60b, `Esame` S67b, `whTracking` S73c (read against the release recipe),
`conclFreeze` S76 (quoted only: literature and Sec. method), `ltf` S78; `lumCapCoulomb` (caption of Fig. lum) is
guarded only by S81. S81 guards that the wording of every one of them is the quoted one.

Measured values that the paper does NOT print but that the arithmetic checks compare with the printed numbers (the
fixed-grid edge 124 at t = 150 read off Fig. mesh, S20; the m drift on four spheres and its fit, S30b; the luminosity
clock Theta_*, arrival times and fit, S54c and S54d; the narrow w_out pulse, S54c) sit in a second association,
`RECORD`, right after DATA; S81 does not cover them.

## Execution log (2026-11-13, `$Version` = 14.2.0 for Linux x86 (64-bit) (December 26, 2024), main.tex of the arXiv v1 source)

The script was run as `wolframscript -file check_equations.wl main.tex` in a directory holding the arXiv v1 `main.tex`.

```
setup done, 7.0 s
DATA: 152 entries; main.tex for S81: main.tex
PASS: S01 UV dispersion [main.tex L140-143] omega^2 = k^2 + alpha k^4 + beta k^6: d omega/dk -> Infinity (~ 3 sqrt(beta) k^2)   | residual: {Infinity, 3*Sqrt[be], Infinity}
PASS: S02 L130 c1+c3 must vanish [main.tex L130] a part in 1e15 on v2 <=> c1 + c3 = c13 ~ 1e-15 (v2 - 1 = c13/2 + O(c13^2))   | residual: 0
PASS: S03 eq:action/M [main.tex L236-248] J from M; M grad u grad u = c13 K.K + c2 K^2 - c14 a^2 (HO aether, twist-free)   | residual: 0
PASS: S04 c13=c14=0.7 [main.tex L250-251,477] c13 = c14 = 7/10 for the production couplings   | residual: {7/10, 7/10}
PASS: S05 eq:speeds, v0 values [main.tex L259-266,478-479] v0^2 > 0; v0 = 1 (base couplings), 1.266862 printed as 1.267 (c2 = 0): arithmetic only   | residual: {0.999999997954823, 1.2668615834434864}
PASS: S06 eq:speeds derived [main.tex L259-266] v2^2, v1^2, v0^2 from the second-order action (plane waves; spin-2, spin-1, spin-0 sectors)   | residual: <large expression>
SKIP: S07 spin-0 energy sign [main.tex L267-271,1110-1111] literature input (Eling:2005zq, Jacobson:2008aj), not re-derived; c14(2-c14) = 0.91 > 0 at the couplings used   | residual: 91/100
PASS: S08 0 < c14 < 2 [main.tex L271] c14(2-c14) > 0 <=> Inequality[0, Less, cc, Less, 2]   | residual: Inequality[0, Less, cc, Less, 2]
PASS: S09 |c13|<1e-15 [main.tex L128-131,277-278] v2 - 1 = c13/2 + O(c13^2) (eq:speeds); the 1e-15 is an observational input   | residual: 0
PASS: S10 eq:lapse [main.tex L304,309-312] u.u=-1, a_r = d_r ln N, a.u=0 for u_a=-N grad t   | residual: {0, 0, 0, 0, 0, 0, 0}
PASS: S11 eq:metric [main.tex L317-328] line element = setup metric; gamma_rr = 1; sphere area 4 pi Theta^2; g^tt = -1/N^2   | residual: {0, 0, 0, 0}
PASS: S12 evolved variables, Q [main.tex L331-339] five evolved fields, no d_t of N, beta, Q; K^th_th = trK/3 - Q/2 from the trace   | residual: {5, {}, 0}
PASS: S13 beta'=N(Q+K/3) [main.tex L359,363-364] K_rr = beta'/N from the covariant K_ab, hence beta' = N(Q+trK/3)   | residual: {0, 0}
PASS: S13b eq:evolution, eq:radial [main.tex L351-359] printed psi, Theta, P, trK, a_r equations == code system; printed Q equation == momentum constraint (S59: G_ue = T_ue); 3 negative controls differ   | residual: {{{0, 0, 0, 0, 0}, 0}, {(3*(c1 + c3)*Nf[t, r]*Qf[t, r]^2)/(2 + c1 + 3*c2 + c3), -2*Pf[t, r]*Derivative[0, 1][psi][t, r], (-4*Nf[t, r]*Derivative[0, 1][psi][t, r]*Derivative[0, 1][Thf][t, r])/Thf[t, r]}}
PASS: S13c Bonetti-Barausse map [main.tex L339-342] mapped action = ours; -(Theta^2/4) x their Hamiltonian constraint = eq:hamiltonian (theta_C); (1-beta) x their momentum constraint = code momentum constraint = -(Q eq. of eq:radial); flipped matter sign fails   | residual: {{0, 0, 0}, 2*Pf[t, r]*Derivative[0, 1][psi][t, r]}
PASS: S14 grad u = -K - u a [main.tex L377-378] decomposition, K symmetric, trK=-div u, d_t Theta equation, a free of d_t N   | residual: 0
PASS: S15 eq:hamiltonian [main.tex L365-375] printed Theta C == code theta_C identically (field equation: S59)   | residual: 0
PASS: S16 eq:constraint-prop [main.tex L384-393] dC1/dt|_{C=0} = (1-alphabar) a_r'(0) trK(0) N(0) (Taylor expansion at r=0)   | residual: 0
PASS: S17 ln N = -Int a_r [main.tex L395-402] d_r ln N = a_r with N->1 at infinity; d_r N = a_r N >= 0 for a_r >= 0; Theta 8 < Theta_AH 12 (= 11.8 rounded)   | residual: {0, 0, 0}
PASS: S18 sources ~ N, |beta| bound [main.tex L403-410] RHS at N = 0 is pure advection; |beta| <= N(r) Int|Q+K/3| from beta(0)=0, beta'=N(Q+K/3), N monotone (example: max(|beta|-bound) = -0.0233)   | residual: {-0.023280140142374478, -0.1473198205360167, -0.4739852374526406, -0.9368217179590936, -1.7496150326301287, -2.6229413360154386}
PASS: S19 ALE C^2 quintic [main.tex L446-456,461-472] taper is C^2 (not C^3), outer 60%; Theta_AH maximum at t ~ 28 as in L679, L731   | residual: {0, 0, 0, 0, 0, 0, 0, 0}
PASS: S20 boundary fall rate [main.tex L447-448,662-668] "roughly one areal unit per unit time" = run average of Fig. 1(a), (284-124)/150 = 1.07; early rates (284-235)/28 = 1.75, (284-250)/20 = 1.7   | residual: {1.0666666666666667, 1.75, 1.7}
PASS: S21 v0(c2) values [main.tex L477-479] c2 = -0.5119403 = -343/670, the exact v0 = 1 root at c1 = 0.7 (v0 values: S05)   | residual: {-343/670, -0.5119402985074627}
PASS: S22 eq:id [main.tex L479-490] profile = code = (r^2/100) x GEJ profile; peak at r = 10.063 with value 1.006 A0; Q = C/Theta^3 => Q = 0 regular   | residual: {0, 0, 10.063002656856709, 1.0063397940042738, C[1]/th0[r]^3}
PASS: S23 eq:expansions, eq:ms, M_hor [main.tex L492-512] theta_pm = (Theta/2) x expansion of u pm e (ingoing sign reversed) = printed; theta+ theta- = (grad Theta)^2; m(theta+=0) = Theta/2; flat theta_pm = 1   | residual: {0, 0, 0, 0, 0, 0, 0, 0}
PASS: S24 order formula [main.tex L530-541,1551-1552] Log2 ratio = p; factor 4 <-> p = 2; N_r^2 factors 1, 4, 16   | residual: {{-4.440892098500626*^-16, 0., 0.}, 2, {1., 3.9993751220474283, 15.996250805503875}}
PASS: S25 Theta C at r=0 [main.tex L548-555] (Theta C)(0) = (Theta'(0)^2-1)/2; the remaining terms are O(r^2)   | residual: {0, 0, 0}
PASS: S26 grid spacings [main.tex L440,559,607-608,631-632] dr = r_max/(N_r-1) = {0.0234, 0.0117, 0.00586} (scattering = collapse); threshold ladder 0.023 ... 0.003 (factor 8)   | residual: {{0.0234375, 0.01171875, 0.005859375}, {0.0234375, 0.01171875, 0.005859375}, {0.0234375, 0.0029296875}, 8}
PASS: S27 scattering table/text [main.tex L557-604,586-588] weak/strong rows >= 1.9; A0=0.06: psi, Theta >= 2, aether + constraint in [1.4,1.5]; N_min decreasing; t = 14 factors 1.8 < 2.5 < 4   | residual: {1.95, 1.9, {2.81, 2.42}, {1.45, 1.47, 1.42, 1.49}, {0.84799690655495, 1.3219280948873624}}
PASS: S28 threshold bracket [main.tex L606-617,622-626] switch 0.135 -> 0.125 with resolution; A* in [0.120, 0.125]; 0.06 = 0.5 x lower end   | residual: {True, True, True, 0.5}
PASS: S29 collapse orders / Table 2 [main.tex L590,642-654,1583-1586,1598-1599] 2.24-2.62 = min/max of the row; coarse ratios ~4; fine ratios ~4 to t=8 then ~6 (t~10 in between); mean log2 = 2.22   | residual: {{2.24, 2.62}, {3.69, 4.69}, {3.95, 3.97}, {5.59, 6.26}, 2.2232708054285224}
PASS: S30 Theta_AH ladder / ADM [main.tex L655-670,1571-1573] factor 4.26, Richardson 12.701, m(Theta=250) drift +0.073%, 3.71107 -> 3.7111   | residual: {4.2602739726027155, 12.701333333333332, 0.07307812617976654, 3.7111}
PASS: S30b m drift ~ Theta^-3 [main.tex L671-674,666] printed: p = 3.1, two-point extrapolation -0.0073%; measured p = 3.07, p=3 fit A = -0.0101% (= "about 0.01%" to two digits), free-p fit (p = 3.12) +0.0037%   | residual: {3.1011790198015454, -0.007264957264957264, 3.0718707774276295, -0.010115914612962407}
PASS: S30c App. A ladder parameters [main.tex L1529-1598] scattering captions A0 = 0.01/0.04/0.06, N_r = 6401/12801/25601, r_max = 150, 0 < r < 20, factor 4, t >~ 12; collapse Table conv-ratio / Fig. conv-collapse A0 = 0.35, r_max = 300, fixed grid, r_AH + 2 <= r <= 150, N_r = 12801/25601/51201; moving mesh N_r = 51201: = body DATA and the run record   | residual: {True, True, True, True, True, True}
PASS: S31 masses, eq:mratio [main.tex L679-687,714-715,724-727,743] 12.70/2 = 6.35; 3.7111 - 6.35 = -2.64; 6.35/3.7111 = 1.711 (1.71)   | residual: {6.35, -2.6388999999999996, 1.7110829673142733}
PASS: S32 d_r m identity [main.tex L689-697] d_r m = Theta^2/2 (rho_G B + j_G A), exact (geometry only)   | residual: 0
PASS: S33 eq:dmdTheta as written [main.tex L692,693-701] paper: coefficient 1/2 of Theta^2 [rho + (A/B) j]; derived kappa = 1, coefficient kappa/2; A = Theta(Q/2 - trK/3)   | residual: {1/2, 0}
PASS: S34 theta = B pm A, DEC, deficit [main.tex L706-712,702-718] theta_pm = B pm A; untrapped <=> B > |A|; -T u = rho u + j e, causal+future <=> rho >= |j|; => bracket >= 0; scalar obeys DEC; aether <= total = -2.64   | residual: {0, 0, {0, 0, 0, 0}, 0, 0}
PASS: S35 Penrose / area theorem [main.tex L719-722,725-735] 1.71 > 1; Theta_AH 12.6 -> 11.8 (area -12.%) after t ~ 28; 12.6 = 12.604 rounded   | residual: {True, True, True, 0.12295288485764677, True}
PASS: S36 eq:evolution block-triangular [main.tex L755-761,348-353] trK, a_r lines of eq:evolution == code RHS; Theta, psi, P equations free of d_r trK, d_r a_r   | residual: {0, 0, {{0, 0}, {0, 0}, {0, 0}}}
PASS: S37 eq:principal [main.tex L762-767] principal matrix of the code's (trK,a_r) equations   | residual: {{0, 0}, {0, 0}}
PASS: S38 eq:aNbN [main.tex L768-779] aN bN = v0^2 (eq:speeds); aN, bN = code spin0_aN, spin0_bN   | residual: {0, 0, 0}
PASS: S39 advection, eigenvalues [main.tex L776-789] aether rdot = -beta; lambda = beta +- N v0; rdot = -lambda; speeds -+ v0 relative to the aether; Franchini et al. c0^2 (alpha, beta, lambda = c14, c13, c2) = v0^2   | residual: {0, 0, 0, 0, 0, 0, 0, 0}
PASS: S40 eq:wout [main.tex L791-799] left eigenvectors and their speeds   | residual: {{0, 0}, {0, 0}}
PASS: S41 exact char. variables [main.tex L799-804] padded w_out, w_in: only a Theta' entry N a_r/Theta x -2. and a psi' entry N P x 3.333; ~ r^-2, r^-1 for a_r ~ P ~ 1/r (neither zero): exact far out; P -+ psi' <-> beta -+ N   | residual: <large expression>
PASS: S42 linearised system [main.tex L812-818,809-823] static frame and the rigid translation r -> r + c(t) (beta = c', gamma_rr = 1; beta = -1 => c = -t + const, dR/dt = -1; d_t|_R = d_t - beta d_r, flat space); linearised code equations   | residual: {{0, 0}, {0, 0, 0, 0}}
PASS: S43 wave equation [main.tex L823-824] d_t^2 trK = v0^2 Laplacian trK   | residual: 0
PASS: S44 eq:asymptotic [main.tex L824-834] outgoing solution solves the linearised system and the wave equation exactly   | residual: {0, 0, 0}
PASS: S45 eq:wout-asym, waveform, d_t w_in [main.tex L835-847] r w_out -> -2 bN f/v0; w_in has no 1/r part, d_t w_in = -bN trK/r (changes while the pulse goes by); the 1/r^2 part F grows during the pulse and stays after it   | residual: {{0, 0, {0, 0}, {0, 0}, 0, 0}, {0.054596713366184324, 0.2219969081188155, 0.38939710319330495, 0.443993816237631, 0.443993816237631, 0.443993816237631}}
PASS: S46 eq:retarded [main.tex L911-921,931] Schwarzschild tortoise derivative; dTheta/dt = v0(1 - 2M/Theta) along u = const (slowing near the black hole)   | residual: {0, 0}
PASS: S47 parities at r=0 [main.tex L849-850] RHS of the odd fields (Theta, a_r) odd and of the even fields (psi, P, trK) even; beta odd   | residual: {{0, 0, 0}, {0, 0, 0}, {0, 0}, {0, 0}, {0, 0, 0}, {0, 0, 0}}
PASS: S48 entering fields; code update [main.tex L851-853] w_in = theta+.u is the field with lambda_+ (enters iff lambda_+ > 0); code check: eigenpairs, theta^J.e_I = delta (theta+.e+ = 2 in the code), P^2 = P, P e_in = 0, theta^in P = 0, theta^out P = theta^out; d_t(theta^in.u) = B; code update = P R + e_in G   | residual: {{0, 0}, {0, 0}, {0, 0}, {0, 0}, 0, 0, 0, {{0, 0}, {0, 0}}, {0, 0}, {0, 0}, {0, 0}, 0, 0, {0, 0}, {{0, 0}, {0, 0}}}
PASS: S49 eq:bc-win [main.tex L853-860,860-864,868,collapse.cpp:L911] w_in = -bN F/r^2 -> 0 at infinity, |w_in|/|w_out| = v0 F/(2 f r) (one power of r below); linear system: d_t w_in + bN K/r = (bN/(v0 r))(d_t + v0 d_r)(r K) identically (equivalence); N^-1(d_t - beta d_r) = u.grad; flat space, constant N, beta: eq:bc-win exact for the outgoing wave; code G = beta d_r w_in - N bN K/Theta = d_t w_in   | residual: {{0, 0, 0, 0, 0, 0, 0, 0}, 0}
PASS: S50 scalar boundary datum [main.tex L803-805,852,857,861-862,871-875,collapse.cpp:L917-928] (u+e)(Theta psi) = Theta (P + psi' + c psi), c = Theta'/Theta + Q/2 - trK/3 = the code's c; exact for g(t-r)/r; entering field P + psi' = -g/r^2; flat: d_t X = d_t(P + psi') + P/r with X = r^-1 (d_t + d_r)(r psi), so d_t(P + psi') = -P/r <=> (d_t + d_r)(r psi) = 0 given X = 0 initially; g(t-r)/r obeys it   | residual: {0, 0, 0, {0, 0}, 0, 0, 0, 0}
PASS: S51 scalar pair in the code [main.tex L851-853,collapse.cpp:L917-928] code imposes d_t(P+psi') = B through d_t P alone: same entering condition; d_t(P-psi') differs from the projector form by -(theta_in.R - B), zero for outgoing waves (B = R_P + R_psi'): the outgoing scalar equation is unchanged for outgoing waves (code check)   | residual: {0, 0, 0, 0}
PASS: S52 outflow gate (code) [collapse.cpp:L1001-1030] both eigenvalues negative <=> beta < -v0 N (spin-0), beta < -N (scalar)   | residual: {True, True}
PASS: S53 w_in set to zero: reflection [main.tex L857-860,877-882] linear theory of the paper's pulse gives 3.48% (r_max=100, printed 3.5%) and 2.33% (r_max=150, printed 2.3%); ~ width/r_max as stated; release record: spin-0 <= 1.5e-6 < 1e-5, pre-fix 3.45-3.53% / 2.33%, scalar 2.8e-4 -> 0.03%   | residual: {3.477369351736214, 2.331889137442018}
PASS: S54 extraction / speed numbers [main.tex L887-977,989] the two runs: v0 = 1 and 1.267 (c2 = 0, > 1: the two signals separate); sigma = 10 in [5,15]; 1.241 is 2.04% below v0 (2%); probe medians 0.01-0.02%   | residual: {True, 2.0413898235978833, 0.02}
PASS: S54a eq:wbar running mean [main.tex L897-909] w_bar[c] = c on any interval (delta of a constant = 0), a + b t kept, t^2 -> t^2 + sigma^2, transfer Exp[-sigma^2 omega^2/2]; code (3-sigma truncation, renormalised ends, sigma = 10): constant to 2.e-15, test series |code - eq| = 5.5e-4 < 0.27% max|w|; figure record <= 0.51% of the pulse   | residual: {{0, 0, 0, 0}, 1.7763568394002505*^-15, 0.0005478391709112795, 0.002699796063260207, 0.00509915014164306}
PASS: S54b probe speed [main.tex L940-947,976-977] (rdot + beta)/N = v0 exactly for the outgoing spin-0 characteristic (no background approximation); record: medians 0.0122% / 0.0154% -> 0.01% / 0.02%; max |v - v0|/v0 for t >= 40 = 2.96e-4 / 1.89e-4 < 3e-4, larger before (3.84e-4 / 3.22e-4)   | residual: {0, 0, {0.01, 0.02}, 0.000296, {0.00008800000000000002, 0.000133}}
PASS: S54c luminosity clock [main.tex L952-967,989] Theta_* = eq:retarded to 0.0049; lag = Delta Theta_*/v, fit v = 1.2406 -> 1.241, -2.07% (2%); rms 0.31; scalar (speed 1) misses by 56; burst = narrow w_out pulse within 0.16 in u   | residual: {0.0048901377655283795, 1.2406102880266892, -2.0721518246249904, 0.31352508341390434, 56.10249946273463, 0, 0.1594813242834956}
PASS: S54d burst at u ~ 10-15 [main.tex L955-956] u of the burst (min of L) = {9.7, 11., 12., 12.7, 13.6, 14.2, 14.7} on Theta = 120 ... 400: 10 to 15 to the nearest unit, the drift = the 2% lag (4.8)   | residual: {{9.700000000000001, 11., 12., 12.700000000000001, 13.600000000000001, 14.200000000000001, 14.700000000000001}, 4.833518224318718}
PASS: S55 MS identity [main.tex L1004] d_a m = Th^2/2 (G_a^b - delta G2) d_b Th (identity)   | residual: {0, 0}
PASS: S56 eq:lum [main.tex L1004-1015] flux formula (identity, kappa=1), A, scalar rho/j/p, F>0 outward   | residual: {{0, 0}, 0, {0, 0, 0}}
PASS: S56b eq:lum from G_ab = T_ab [main.tex L1004] "Using the field equations G_ab = T_ab, one obtains L": S55 identity + G = T gives eq:lum exactly (S56); the full tensor equation is needed (dF/dp = N Theta^2 A/2, the rr equation; dm/dQ = Theta^2 A/2 needs d_t Q, absent from eq:evolution)   | residual: {((-2*Kq + 3*Qq)*Thq^3)/12, (aa1x56*nnx56*Thq^2)/2}
SKIP: S58 flux-test numbers [main.tex L1016-1021,1027] 0.02% residual, |F_psi| <~ 0.02 vs 0.9-4: simulation data (caption spheres 40, 100, 200 within 40-200)   | residual: "n/a"
PASS: S59 T^ae -> constraints [main.tex L365-375,1044-1061] rho_ae, j_ae from 1501.07274 (2.34); Ham. constraint = code theta_C = eq:hamiltonian (kinematics only); momentum constraint = compute_Q (needs the a_r equation)   | residual: {{0, 0, 0, 0, 0, 0}, {}, {arf}}
PASS: S60 Coulomb field a ~ M/Theta^2 [main.tex L1044-1049] static aether: Theta^2 a_r -> M (a_r ~ Theta^-2), the F_inf/r^2 term   | residual: MM
PASS: S60b luminosity swing, m above M_ADM [main.tex L1037-1042,1049] record: swing of m = 7.05, 9.38, 14.05, 18.76 (grows with radius); late m - M_ADM = 1.05-1.17 > 0 (> 0 from u = 107); int L du = -Delta m = {-1.24, -1.19, -1.12, -1.09} < 0 on every sphere   | residual: {0, {-1.239, -1.1869999999999998, -1.119, -1.085}, {-1.26, -1.263, -1.286, -1.32}}
PASS: S61 eq:dm-coulomb [main.tex L1051-1061] exact: d_Theta m = (c14/2) Theta' d_Theta(Theta^2 a) + quadratic = eq:dm-coulomb (Theta' factor included)   | residual: {0, 0}
SKIP: S62 Theta^2 a_r ~ Theta^0.7 [main.tex L1090-1091] simulation data (exponent 0.7)   | residual: "n/a"
PASS: S63 eq:mc [main.tex L1067-1074] vacuum a=K=Q=0 data: theta_C = 0 <=> m = const, so m_c = M_ADM, E = 0 initially   | residual: 0
PASS: S64 eq:Mae [main.tex L1075-1082] surface term c14/(8 pi) oint a dA = c14/2 Theta^2 a_r: the coefficient of eq:mc   | residual: 0
PASS: S65 eq:Einf arithmetic [main.tex L1098-1111,1135,1452-1453] 0.532 - 0.197 = 0.335 exactly; fractions share M_ADM = 0.532/0.144 = 3.694; 9% + 5% = 14%; caption = eq:Einf rounded   | residual: {0.335, 0.09067669172932331, 0.006225563909774435, 0.05332330827067669, 3.6944444444444446}
PASS: S66 eq:Einf M_ADM [main.tex L1098,1104-1105] with the c2=0 run's M_ADM = 3.708 (L1098): 0.532/3.708 = 0.1435 -> 0.144   | residual: {0.14347357065803668, {0.14331940137521912, 0.14362778152393796}}
PASS: S67 Coulomb cancellation, -3/4 c13 Q^2 [main.tex L1121-1122,1124] paper: cancellation up to O(1/Theta^1); derived remainder (c14/2)(Theta'-1) d(Theta^2 a) = O(1/Theta^1); Q^2 coefficient (-3*c1)/4 - (3*c3)/4   | residual: {-1/4*((c1 + c4)*MM^2)/Th^2 - ((c1 + c4)*MM)/(2*Th), (-3*c1)/4 - (3*c3)/4}
PASS: S67b fluxinf: rate, E, E_inf(u) [main.tex L1114-1119] record: uncorrected/corrected max |L| = 25 (Theta = 200), 11.3 (Theta = 120), panel peaks 28 (> 10); E spread <= 0.044 = 7.7% of the mean for u >= 110 (once the pulse has passed), 6.8% for u = 200-250; E_inf(u) in 0.522-0.540 for u >= 200, inside 0.532 +- 0.023   | residual: {25.364238410596027, 11.279296875, 28.33962264150943, 0.06771378708551494, {0.522, 0.54}, 0.07678883071553229}
PASS: S68 PN: v_L, Gamma, inverse propagator [main.tex L1140-1243] eq:pn_first_surviving, eq:G_exact, eq:inverse_propagator: v0^2(c13=0) = v_L^2; v1^2(c13=0) = c1/c14; Gamma = -(1-c14/2)(1-x)/s_{-1} = (1-x/x*)/(1-x/vL^2) = 1 + (c14/2) y/(1-y) exactly; inverse propagator = y/(1-y); 2PN lapse operator = k^2 s_1; pole at x = vL^2 not cancelled   | residual: {0, 0, 0, 0, 0, 0, 0, 0}
PASS: S69 c2(c1) for v0=1 [main.tex L1255-1256] c2 = -c1^3/(2-4c1+3c1^2); = -0.711111 at c1=0.8, v0 = 1.   | residual: {0, -0.7111111111111111, 1.0000000000000002}
PASS: S70 eq:wh-id [main.tex L1257-1265] ingoing data P = psi' + psi/r   | residual: 0
PASS: S71 theta_- anti-trapped [main.tex L1265-1270] theta_- = -(ingoing null expansion); theta_- < 0 <=> ingoing rays diverge; 1 in flat space   | residual: {0, 0}
PASS: S72 white-hole timeline [main.tex L1272-1287] BH at 0.12 > WH at 0.10; 2.5 < 6 < 17 < 20 < 36 < 40 < 44; AGG t ~ 17 vs 20   | residual: {True, True, True}
PASS: S72b anti-trapped split and layer [main.tex L1277-1284] one region at t = 34, two at t = 36; inner remnant 0.18 wide at t = 44; layer width 0.84-0.91 (t = 40) -> 0.061-0.068 (t = 80); all-grid spread {8, 18, 11}% (<= 20%; successive grids {8, 11, 8}%); min theta_- grows with resolution   | residual: {True, True, True, True, True, True, {8, 11, 8}, {8, 18, 11}}
PASS: S73 white-hole scalings [main.tex L1291-1330] three features; rho ~ |P|^2: 2(-0.5) = -1; width x height ~ dr^0.3 -> 0; 7 = 3 + 4 grids; t = 30 just after 29.5   | residual: {0., True, True, True, True, 6.666666666666667}
PASS: S73b white-hole exponents vs rerun [main.tex L1324-1326] width dr^0.84 (uniform) / dr^0.8 (refined), |P| dr^-0.5 / dr^-0.48, rho dr^-0.99 / dr^-0.95: each within 0.06 of the printed dr^0.8, dr^-0.5, dr^-1   | residual: {{0.84, -0.5, -0.99}, {0.8, -0.48, -0.9500000000000001}}
PASS: S73c white-hole refinements and layer [main.tex L1297-1323,1423-1426] factors 8 (pulse) and 8/16 (shell), hold t = 20-30.5 = release recipe, brackets 29.5 and 30; x0 = 9.9 vs r_form = 9.96; layer widths {0.67, 0.43} < 1 with the |P| peak inside   | residual: {True, True, True, True, True, {0.67, 0.43}, True, True}
PASS: S73d inner front and Fig. sharp [main.tex L1332-1341,1400-1408] jump 0.26-0.35 ~ 0.3; dr/dt = -beta - N within 0.0015; N 0.125-0.137; Theta 3.15-3.22, just inside the inner theta_- = 0 edge; slope ratios {1.27, 1.11} (t = 24) and {2.08, 2.07} (t = 40); caption peaks 21.17 -> 21, 61.51 -> 62, ladder and factors   | residual: {True, True, True, True, True, True, True, True, True, True, True, {1.27, 1.11}, {2.08, 2.07}}
PASS: S73e rays converge on the pulse [main.tex L1309-1312,1423-1424] pulse - ray gap 3.158 (t = 30) -> 0.302 (t = 34) -> 0.091 (t = 40), printed 3.2, 0.3, 0.1; joins by t ~ 34   | residual: {True, True, True, True, True, True}
PASS: S74 frozen centre / other data [main.tex L1343-1375] factor-2 ladders (0.009/0.0047, 0.0031/0.0016/0.0008); dispersing amplitudes < 0.08; trapped (t~20) before anti-trapped (t~24)   | residual: {True, True, True, True, True, True, True}
PASS: S75 AGG comparison [main.tex L1377-1395,1415-1416] width ratio 2/0.8 = 2.5 (FWHM ratio 2.55335); 40,000 -> 80,000 halving; Theta ~ 7 = 6.9 rounded   | residual: {5/2, 2.5533536958926075, 2, 7}
PASS: S76 conclusion [main.tex L1450-1453] 14% = 9% + 5% = eq:Einf   | residual: {True}
PASS: S77 App. B.1 toy [main.tex L1614-1668,1875-1884] eq:simple_scalar, eq:simple_L, eq:Tp, eq:Lprime, eq:Lprime_p1: L = -k^2(1-y), y <-> (eta/v^2) d_t^2 nabla^-2; L T_p = -k^2(1-y^{p+1}); p=1 form; term eta^{p+1} d_t^{2p+2} (nabla^-2)^p; y = 1 kept, p spurious roots   | residual: {0, 0, 0, 0, 0}
PASS: S78 B^rr_LTF [main.tex L1682] (d_j d_k gamma2) trace-free rr part = (2/3)(gamma2'' - gamma2'/r), as printed   | residual: {0, (Derivative[1][gam2][rr0] - rr0*Derivative[2][gam2][rr0])/(3*rr0)}
PASS: S79 App. B.2 reduction [main.tex L1670-1785] eq:x_def, eq:det_M0, eq:vL_def, eq:sminus1, eq:pn_1_term: vector-Laplacian identity; sin(kr)/r; gauge => K_hat = N; Eqs. (4.5a), (4.6) on the mode = rows 1, 3 of D0; Box -> -k^2(1-x); det D0, s_{-1}, g1 and the 1PN term as printed   | residual: {0, 0, Nn, 0, 0, 0, 0, 0, 0}
PASS: S80 App. B.3 redefinition [main.tex L1787-1907] T1 = D0^-1 D1 to O(x): row 1 = eq:N_transform, column 1 at x^0 as printed; D0 T1 = eq:M1 + O(x^2); s_1 = eq:s1, kappa_1 (both forms); s_1 = s_{-1}(1 + c14 y/2) + O(x^3); Gamma series; det identity; summary: no v_L in s_1 through 1PN, v_L at 2PN   | residual: {{{0, 0, 0}, {0, 0, 0}, {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}}, 0, 0, 0, 0, 0, 0}, {0, 0, True}}
PASS: S81 DATA snippets vs main.tex [main.tex Lall] 152/152 snippets found at the quoted lines of main.tex   | residual: {}

==================== SUMMARY ====================
status  id                                          main.tex lines (collapse.cpp:L... = code lines)
PASS    S01 UV dispersion                           140-143
PASS    S02 L130 c1+c3 must vanish                  130
PASS    S03 eq:action/M                             236-248
PASS    S04 c13=c14=0.7                             250-251,477
PASS    S05 eq:speeds, v0 values                    259-266,478-479
PASS    S06 eq:speeds derived                       259-266
SKIP    S07 spin-0 energy sign                      267-271,1110-1111
PASS    S08 0 < c14 < 2                             271
PASS    S09 |c13|<1e-15                             128-131,277-278
PASS    S10 eq:lapse                                304,309-312
PASS    S11 eq:metric                               317-328
PASS    S12 evolved variables, Q                    331-339
PASS    S13 beta'=N(Q+K/3)                          359,363-364
PASS    S13b eq:evolution, eq:radial                351-359
PASS    S13c Bonetti-Barausse map                   339-342
PASS    S14 grad u = -K - u a                       377-378
PASS    S15 eq:hamiltonian                          365-375
PASS    S16 eq:constraint-prop                      384-393
PASS    S17 ln N = -Int a_r                         395-402
PASS    S18 sources ~ N, |beta| bound               403-410
PASS    S19 ALE C^2 quintic                         446-456,461-472
PASS    S20 boundary fall rate                      447-448,662-668
PASS    S21 v0(c2) values                           477-479
PASS    S22 eq:id                                   479-490
PASS    S23 eq:expansions, eq:ms, M_hor             492-512
PASS    S24 order formula                           530-541,1551-1552
PASS    S25 Theta C at r=0                          548-555
PASS    S26 grid spacings                           440,559,607-608,631-632
PASS    S27 scattering table/text                   557-604,586-588
PASS    S28 threshold bracket                       606-617,622-626
PASS    S29 collapse orders / Table 2               590,642-654,1583-1586,1598-1599
PASS    S30 Theta_AH ladder / ADM                   655-670,1571-1573
PASS    S30b m drift ~ Theta^-3                     671-674,666
PASS    S30c App. A ladder parameters               1529-1598
PASS    S31 masses, eq:mratio                       679-687,714-715,724-727,743
PASS    S32 d_r m identity                          689-697
PASS    S33 eq:dmdTheta as written                  692,693-701
PASS    S34 theta = B pm A, DEC, deficit            706-712,702-718
PASS    S35 Penrose / area theorem                  719-722,725-735
PASS    S36 eq:evolution block-triangular           755-761,348-353
PASS    S37 eq:principal                            762-767
PASS    S38 eq:aNbN                                 768-779
PASS    S39 advection, eigenvalues                  776-789
PASS    S40 eq:wout                                 791-799
PASS    S41 exact char. variables                   799-804
PASS    S42 linearised system                       812-818,809-823
PASS    S43 wave equation                           823-824
PASS    S44 eq:asymptotic                           824-834
PASS    S45 eq:wout-asym, waveform, d_t w_in        835-847
PASS    S46 eq:retarded                             911-921,931
PASS    S47 parities at r=0                         849-850
PASS    S48 entering fields; code update            851-853
PASS    S49 eq:bc-win                               853-860,860-864,868,collapse.cpp:L911
PASS    S50 scalar boundary datum                   803-805,852,857,861-862,871-875,collapse.cpp:L917-928
PASS    S51 scalar pair in the code                 851-853,collapse.cpp:L917-928
PASS    S52 outflow gate (code)                     collapse.cpp:L1001-1030
PASS    S53 w_in set to zero: reflection            857-860,877-882
PASS    S54 extraction / speed numbers              887-977,989
PASS    S54a eq:wbar running mean                   897-909
PASS    S54b probe speed                            940-947,976-977
PASS    S54c luminosity clock                       952-967,989
PASS    S54d burst at u ~ 10-15                     955-956
PASS    S55 MS identity                             1004
PASS    S56 eq:lum                                  1004-1015
PASS    S56b eq:lum from G_ab = T_ab                1004
SKIP    S58 flux-test numbers                       1016-1021,1027
PASS    S59 T^ae -> constraints                     365-375,1044-1061
PASS    S60 Coulomb field a ~ M/Theta^2             1044-1049
PASS    S60b luminosity swing, m above M_ADM        1037-1042,1049
PASS    S61 eq:dm-coulomb                           1051-1061
SKIP    S62 Theta^2 a_r ~ Theta^0.7                 1090-1091
PASS    S63 eq:mc                                   1067-1074
PASS    S64 eq:Mae                                  1075-1082
PASS    S65 eq:Einf arithmetic                      1098-1111,1135,1452-1453
PASS    S66 eq:Einf M_ADM                           1098,1104-1105
PASS    S67 Coulomb cancellation, -3/4 c13 Q^2      1121-1122,1124
PASS    S67b fluxinf: rate, E, E_inf(u)             1114-1119
PASS    S68 PN: v_L, Gamma, inverse propagator      1140-1243
PASS    S69 c2(c1) for v0=1                         1255-1256
PASS    S70 eq:wh-id                                1257-1265
PASS    S71 theta_- anti-trapped                    1265-1270
PASS    S72 white-hole timeline                     1272-1287
PASS    S72b anti-trapped split and layer           1277-1284
PASS    S73 white-hole scalings                     1291-1330
PASS    S73b white-hole exponents vs rerun          1324-1326
PASS    S73c white-hole refinements and layer       1297-1323,1423-1426
PASS    S73d inner front and Fig. sharp             1332-1341,1400-1408
PASS    S73e rays converge on the pulse             1309-1312,1423-1424
PASS    S74 frozen centre / other data              1343-1375
PASS    S75 AGG comparison                          1377-1395,1415-1416
PASS    S76 conclusion                              1450-1453
PASS    S77 App. B.1 toy                            1614-1668,1875-1884
PASS    S78 B^rr_LTF                                1682
PASS    S79 App. B.2 reduction                      1670-1785
PASS    S80 App. B.3 redefinition                   1787-1907
PASS    S81 DATA snippets vs main.tex               all
CAVEATs (findings for the authors):
PASS: 93   FAIL: 0   CAVEAT: 0   SKIP: 3   total: 96   time: 33.3 s
```

STATUS meaning: PASS = derived / verified; FAIL = the statement as written in main.tex is wrong; CAVEAT = correct
only with a qualification, or a printed detail does not match the paper's own equations or numbers (a finding for
the authors; exit code unchanged); SKIP = simulation data or literature input, not checkable from equations (or
main.tex absent, S81).

## Equation list and status (main.tex of the arXiv v1 source)

| id | main.tex lines (collapse.cpp:L… = code lines) | item | status |
|---|---|---|---|
| S01 | 140-143 | Introduction: omega^2 = k^2 + alpha k^4 + beta k^6, group velocity without bound | PASS |
| S02 | 130 | Introduction: "the combination of the coupling parameters $c_1+c_3$" must vanish (v2 - 1 = c13/2) | PASS |
| S03 | 236-248 | eq:action and the tensor M | PASS |
| S04 | 250-251, 477 | Inline: c_ij combinations and the production couplings | PASS |
| S05 | 259-266, 478-479 | eq:speeds: v0 arithmetic | PASS |
| S06 | 259-266 | eq:speeds: v2, v1, v0 derived from the action | PASS |
| S07 | 267-271, 1110-1111 | Inline: the spin-0 energy has the sign of c14(2-c14) (literature) | SKIP |
| S08 | 271 | Inline: positive spin-0 energy "only if $0<c_{14}<2$" | PASS |
| S09 | 128-131, 277-278 | Inline: GW170817 bound /c13/ < 1e-15 | PASS |
| S10 | 304, 309-312 | eq:lapse and inline u_a=-N nabla_a t | PASS |
| S11 | 317-328 | eq:metric: gamma_rr = 1, areal radius | PASS |
| S12 | 331-339 | Displayed and inline: evolved variables; Q defined by K^r_r = trK/3 + Q | PASS |
| S13 | 359, 363-364 | eq:radial: d_r beta = N(Q + trK/3) from gamma_rr = 1 | PASS |
| S13b | 351-359 | eq:evolution and eq:radial: printed equations == code system; momentum constraint == G_ue = T_ue; negative controls | PASS |
| S13c | 339-342 | Bonetti-Barausse 3+1 equations: map (alpha, beta, lambda) = (c14, c13, c2), f = alpha a^2; action, Hamiltonian, momentum constraint | PASS |
| S14 | 377-378 | Inline: a_r fixed on each slice by the lapse; nabla_a u_b = -K_ab - u_a a_b contains no d_t N (derived) | PASS |
| S15 | 365-375 | eq:hamiltonian: the printed Theta C is the code's theta_C | PASS |
| S16 | 384-393 | eq:constraint-prop, with d_r ln N = alphabar a_r | PASS |
| S17 | 395-402 | Inline: ln N(r) = - Int_r^infinity a_r dr'; N monotone where a_r >= 0 | PASS |
| S18 | 403-410 | Displayed: sources proportional to N; /beta(r)/ <= N(r) Int_0^r /Q + trK/3/ dr' | PASS |
| S19 | 446-456, 461-472 | ALE map and the moving-mesh caption | PASS |
| S20 | 447-448, 662-668 | Inline: "roughly one areal unit per unit time" = run average of Fig. 1(a), (284-124)/150 = 1.07 | PASS |
| S21 | 477-479 | Inline: c2=-0.5119403 is the exact v0=1 root -343/670 (the v0 values themselves: S05) | PASS |
| S22 | 479-490 | eq:id | PASS |
| S23 | 492-512 | eq:expansions, eq:ms, M_hor = Theta_AH/2 | PASS |
| S24 | 530-541, 1551-1552 | Convergence order, dt = 0.1 h, the factor 4, and the rescaled-constraint factors | PASS |
| S25 | 548-555 | Inline: Theta C at r -> 0 is the regularity condition Theta'(0) = 1 | PASS |
| S26 | 440, 559, 607-608, 631-632 | Grid spacings | PASS |
| S27 | 557-604, 586-588 | Scattering (L557-604, Table 1 L586-588, appendix L1520-1556) | PASS |
| S28 | 606-617, 622-626 | Black-hole threshold (L606-617, Fig. 2 caption L622-626) | PASS |
| S29 | 590, 642-654, 1583-1586, 1598-1599 | Collapse orders (Table 1 last row, L590; L642-654) and Table 2 | PASS |
| S30 | 655-670, 1571-1573 | Theta_AH ladder, Richardson limit, ADM conservation | PASS |
| S30b | 671-674, 666 | m on Theta = 250 grows by 0.07% but falls off as Theta^-3 (1.2% at Theta = 100), extrapolation about 0.01% (p = 3 fit -0.0101%, record -0.0102%; free-p fit +0.0037%) | PASS |
| S30c | 1529-1598 | App. A ladder parameters: captions of Figs. conv-fields-A001/A006 and conv-constraint (A0 = 0.01/0.04/0.06, N_r = 6401/12801/25601, r_max = 150, 0 < r < 20, factor 4, t >~ 12), Table conv-ratio and Fig. conv-collapse (A0 = 0.35, r_max = 300, fixed grid, r_AH + 2 <= r <= 150, N_r = 12801/25601/51201), the N_r = 51201 moving mesh, against the body DATA and the run record   [DATA["appScatFig"], ["appScatA006"], ["appConstraintCap"], ["appAHMesh"], ["appMeshN"], ["tabConvRatioCap"], ["convCollapseCap"]] | PASS |
| S31 | 679-687, 714-715, 724-727, 743 | Masses and eq:mratio | PASS |
| S32 | 689-697 | The exact relation behind eq:dmdTheta: d_r m = (Theta^2/2) (rho_G B + j_G A) | PASS |
| S33 | 692, 693-701 | eq:dmdTheta as written, with A = u(Theta), B = Theta' and G_ab = T_ab   [claim in DATA["dmClaim"]] | PASS |
| S34 | 706-712, 702-718 | theta_pm = B pm A, untrapped spheres, m decreases outward only where the DEC fails, aether deficit | PASS |
| S35 | 719-722, 725-735 | Penrose inequality from eq:dmdTheta under the DEC; eq:mratio; area theorem | PASS |
| S36 | 755-761, 348-353 | eq:evolution as a first-order system: trK, a_r lines == code RHS; Theta, psi, P equations free of d_r trK, d_r a_r | PASS |
| S37 | 762-767 | eq:principal | PASS |
| S38 | 768-779 | eq:aNbN | PASS |
| S39 | 776-789 | Inline: advection rdot = -beta; eigenvalues lambda_pm = beta pm N v0; speeds pm v0 relative to the aether; Franchini et al.: c0^2 with (alpha, beta, lambda) = (c14, c13, c2) equals v0^2 | PASS |
| S40 | 791-799 | eq:wout: left eigenvectors and their speeds | PASS |
| S41 | 799-804 | Inline: exact characteristic variables far from the source (couplings to Theta ~ a_r/Theta and to psi ~ P fall off); the scalar combinations P -+ d_r psi | PASS |
| S42 | 812-818, 809-823 | Static frame; rigid translation r -> r + c(t) (beta = c-dot, gamma_rr = 1, beta = -1 => c = -t, d_t at fixed areal radius = d_t - beta d_r in flat space); linearised (trK, a_r) system | PASS |
| S43 | 823-824 | Wave equation | PASS |
| S44 | 824-834 | eq:asymptotic | PASS |
| S45 | 835-847 | eq:wout-asym: the 1/r^2 part builds up as the pulse goes by and remains; waveform; w_in changes, d_t w_in = -b_N trK/r | PASS |
| S46 | 911-921, 931 | eq:retarded (Extraction paragraph): tortoise coordinate, slowing near the black hole | PASS |
| S47 | 849-850 | Inline: parities at r = 0 are preserved by the evolution equations | PASS |
| S48 | 851-853 | Inline: the entering fields w_in and P + d_r psi; the code update (projector algebra, code check) | PASS |
| S49 | 853-860, 860-864, 868, collapse.cpp:L911 | eq:bc-win: entering fields of an outgoing wave (w_in = -b_N F/r^2 -> 0); (d_t + v0 d_r)(r trK) = 0 <=> d_t w_in = -b_N trK/r; exact on flat space; code G | PASS |
| S50 | 803-805, 852, 857, 861-862, 871-875, collapse.cpp:L917-928 | Inline: P + d_r psi = -g/r^2; Theta psi constant along the outgoing light rays; flat space: d_t(P + d_r psi) = -P/r | PASS |
| S51 | 851-853, collapse.cpp:L917-928 | The scalar condition in the code: through d_t P alone, equivalent up to terms vanishing for outgoing waves (code check) | PASS |
| S52 | collapse.cpp:L1001-1030 | Outflow gate of the code: beta < -v0 N (spin-0), beta < -N (scalar) (not described in the paper; code check) | PASS |
| S53 | 857-860, 877-882 | Inline: reflection with w_in set to zero, 3.5% and 2.3% from the linear theory; the reflection numbers against the release VALIDATION record | PASS |
| S54 | 887-977, 989 | Results opening (v0 = 1 and 1.267; the two steps), extraction and propagation speed   [DATA["resultsIntro"], ["resultsSteps"]] | PASS |
| S54a | 897-909 | eq:wbar: the Gaussian running mean is normalised (a constant is returned on any interval, delta of a constant = 0), keeps a + b t, t^2 -> t^2 + sigma^2, transfer Exp[-sigma^2 omega^2/2] (delta w_out a high-pass filter); against release/analysis/spin0_common.py gaussian_running_mean (3-sigma truncation, renormalised ends) on a test series and the sigma = 10 figure record   [DATA["wbarEq"]] | PASS |
| S54b | 940-947, 976-977 | The probe: (rdot + beta)/N = v0 exactly for the outgoing spin-0 characteristic (no background approximation); record (probe.log, probe_speed.npz): medians 0.0122 % / 0.0154 % -> 0.01 % / 0.02 % (body and caption), max \|v - v0\|/v0 for t >= 40 = 2.96e-4 / 1.89e-4 < 3e-4, larger before t = 40   [DATA["probeExact"], ["probeScatter"], ["probeMedText"]] | PASS |
| S54c | 952-967, 989 | The luminosity clock: Theta_* = eq:retarded (M_ADM = 3.7084), lag = Delta Theta_*/v, least-squares v = 1.2406 -> 1.241, 2.07 % -> 2 %, rms 0.31, the scalar speed 1 misses by 56; scalar speed 1 from its eigenvalues; burst = narrow w_out pulse within 0.16 in u   [DATA["lumDef"], ["lumSpheres"], ["lagLine"], ["lagSpeed"], ["burstSpin0"]] | PASS |
| S54d | 955-956 | "on every sphere a sharp burst at u ~ 10-15" (an earlier wording, "at u ~ 10", was a CAVEAT; `CHANGELOG.md`): the burst's u = t_arr - Theta_*/v0 is 9.7 ... 14.7 on Theta = 120 ... 400, 10 to 15 to the nearest unit (the drift 4.8 of the 2 % lag)   [DATA["lumBurst"]] | PASS |
| S55 | 997-1004 | Opening of "Mass conservation and radiated energy" (the two steps: S55-S58, S59-S66); Misner-Sharp flux identity: time derivative at fixed areal radius = flux through the sphere   [DATA["fluxOpening"], ["fluxHayward"]] | PASS |
| S56 | 1004-1015 | eq:lum and its definitions | PASS |
| S56b | 1004 | Inline: "Using the field equations G_ab = T_ab, one obtains the luminosity L" (sufficient; the full tensor equation is needed) | PASS |
| S58 | 1016-1021, 1027 | Flux-test numbers | SKIP |
| S59 | 365-375, 1044-1061 | Aether stress tensor -> Hamiltonian and momentum constraints (behind eq:hamiltonian and eq:dm-coulomb) | PASS |
| S60 | 1044-1049 | Inline: a_r prop Theta^-2 around a black hole; Theta^2 a_r relaxes to a constant | PASS |
| S60b | 1037-1042, 1049 | "the luminosity gives a confusing answer": record (aether_coulomb.npz): swing of m 7.05 / 9.38 / 14.05 / 18.76 on Theta = 75 ... 200 (grows with radius), m - M_ADM 1.05-1.17 > 0 behind the pulse (> 0 from u = 107), int L du = -Delta m = -1.24 ... -1.09 < 0 on every sphere; "We call such a contribution the aether Coulomb field" (a name)   [DATA["lumConfusing"], ["coulombName"]] | PASS |
| S61 | 1051-1061 | eq:dm-coulomb | PASS |
| S62 | 1090-1091 | Caption: the amplitude of Theta^2 a_r grows roughly as Theta^0.7 | SKIP |
| S63 | 1067-1074 | eq:mc | PASS |
| S64 | 1075-1082 | eq:Mae | PASS |
| S65 | 1098-1111, 1135, 1452-1453 | eq:Einf (L1098-1111), caption L1135, conclusion L1450-1453 | PASS |
| S66 | 1098, 1104-1105 | eq:Einf with the c2 = 0 run's M_ADM = 3.708 | PASS |
| S67 | 1121-1122, 1124 | Inline: the linear Coulomb term cancels up to O(M/Theta); -3/4 c13 Q^2   [DATA["coulombClaim"], DATA["q2Claim"]] | PASS |
| S67b | 1114-1119 | Fig. fluxinf in time, record: corrected max \|L\| 0.265 ... 0.151 (Theta = 75 ... 200) vs the uncorrected 2.31 ... 7.51 of Fig. lum(a) (Theta = 120 ... 400): ratio 25 at Theta = 200, 11.3 at Theta = 120, 28 between the panel peaks ("more than an order of magnitude"; not printed: 5.3 / 8.3 on the inner spheres 75 / 100); E "nearly the same on all spheres once the pulse has passed": spread <= 0.044 for u >= 110, 0.5945 ... 0.5557 (6.8 %) for u = 200-250; up to 0.43 during the crossing, u = 10-95, excluded by the qualifier; E_inf(u) in 0.522-0.540 for u = 200-255, inside eq:Einf   [DATA["corrRate"], ["Esame"]] | PASS |
| S68 | 1140-1243 | Sec. V.E (sec:pn_comparison): the spin-0 mode in the post-Newtonian expansion; eq:pn_first_surviving, eq:G_exact, eq:inverse_propagator | PASS |
| S69 | 1255-1256 | Inline: c2 = -c1^3/(2-4c1+3c1^2) = -0.7111 keeps v0 = 1 | PASS |
| S70 | 1257-1265 | eq:wh-id | PASS |
| S71 | 1265-1270 | Inline: theta_- < 0 is anti-trapped | PASS |
| S72 | 1272-1287 | Inline: the A0 = 0.12 black hole and the A0 = 0.10 white-hole timeline (split at t = 36) | PASS |
| S72b | 1277-1284 | Inline: the anti-trapped region splits at t = 36; the thin layer narrows from 0.9 (t = 40) to 0.06 (t = 80), grids agree to within 20% (all-grid spread 8/18/11%), depth grows with resolution (theta_- = 0 zeros of box40_N51200, whconv ladder) | PASS |
| S73 | 1291-1330 | Inline: the three sharp features; scalings of the escaping pulse and the inner shell | PASS |
| S73b | 1324-1326 | Inline: dr^0.8, dr^-0.5, dr^-1 against the per-family fits of the rerun (0.84/0.80, -0.50/-0.48, -0.99/-0.95) | PASS |
| S73c | 1297-1323, 1423-1426 | Inline and Fig. wh caption: refinement factors 8 and 8 or 16, hold t = 20-30.5 (release recipe); anti-trapped layer less than one unit wide | PASS |
| S73d | 1332-1341, 1400-1408 | Inline and Fig. sharp caption: the inner ingoing front (t = 7, jump 0.3, rdot = -beta - N, N = 0.13, Theta = 3.2 near the inner theta_- = 0 edge, slope ratios at t = 24 / 40, t = 44 / 54) and the caption numbers (fig_sharp_features.py record) | PASS |
| S73e | 1309-1312, 1423-1424 | Inline and Fig. wh caption: the outgoing ray leaving the inner shell at t = 30 converges on the escaping pulse (3.2 behind at t = 30, 0.3 at t = 34, 0.1 at t = 40; gaps 3.158, 0.302, 0.091) | PASS |
| S74 | 1343-1375 | Inline: frozen centre (dr^-1/3) and the other initial data | PASS |
| S75 | 1377-1395, 1415-1416 | Inline: comparison with Akhoury-Garfinkle-Gupta | PASS |
| S76 | 1450-1453 | Conclusion: 14% = 9% + 5% (the new freezing / universal-horizon sentence, L1474-1479, is literature and Sec. method, guarded by S81)   [DATA["conclSplit"], ["conclFreeze"]] | PASS |
| S77 | 1614-1668, 1875-1884 | App. B.1: the single-field example; eq:simple_scalar, eq:simple_L, eq:Tp, eq:Lprime, eq:Lprime_p1 | PASS |
| S78 | 1682 | App. B.2: B^rr_LTF = (2/3)(gamma_2'' - gamma_2'/r) | PASS |
| S79 | 1670-1785 | App. B.2: spherical reduction, D0, det D0, s_{-1}, the 1PN v_L term; eq:x_def, eq:det_M0, eq:vL_def, eq:sminus1, eq:pn_1_term | PASS |
| S80 | 1787-1907 | App. B.3: the first-order redefinition, D1, s_1, kappa_1; the summary (L1903-1905: no v_L in s_1 through 1PN, v_L at 2PN) | PASS |
| S81 | all | DATA snippets vs main.tex (is the text audited here the text of the paper?) | PASS |

## Open findings

None against the arXiv v1 source. The findings raised on the drafts of the paper (FAILs and CAVEATs), how each was
resolved, and the sections that now check the corrected statements are listed in `CHANGELOG.md`.

## Not checked

Simulation outputs (scattering and collapse orders, threshold ladder, residual percentages, probe medians, Theta_AH
history 12.6 -> 11.8, WH times and widths), literature inputs (the c14(2-c14) energy-sign criterion, the
observational value of the GW bound, Garfinkle-Jacobson positivity, the field equations of Taherasghari-Will and
row 2 of D0).
