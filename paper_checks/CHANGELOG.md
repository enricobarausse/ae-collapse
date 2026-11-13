# paper_checks: history

The two audits in this directory check the paper *Gravitational collapse in Lorentz-violating
gravity and the spin-0 mode* against first principles, against its own numbers and against the
simulation code `release/collapse.cpp`:

* `check_equations.py` (SymPy, `README.md`) derives every displayed and labelled equation and every
  inline relation that carries physics from the action or from the geometry of the code's metric,
  checks the arithmetic of every simulation number the paper quotes, and compares the measured
  values behind the paper's statements with the records of the production runs.
* `check_equations.wl` / `check_equations.nb` (Mathematica, `README_mathematica.md`) does the same
  independently (no shared code), and in addition checks that every quoted number and sentence
  occurs at its line of `main.tex` (S81).

Both audits were run on the drafts of the paper while it was written, and on the final text. Every
FAIL and CAVEAT they raised was resolved before the arXiv v1 source, either by a correction of the
paper or of the code (below) or, in one case, by the authors' reading of the statement, which the
audits now check. Section ids are those of the current scripts (SymPy / Mathematica).

## Findings that corrected the paper

* **Mass function and energy conditions** (S23 / S32-S34). The draft wrote ∂_Θ m = 4πΘ²ρ; with
  the normalisation of the paper and of the code (G_ab = T_ab) the exact relation is
  ∂_Θ m|_t = ½Θ²[ρ + (A/B) j], and a negative ρ + (A/B) j violates the *dominant*, not the weak,
  energy condition (a WEC-respecting counterexample is constructed). Now Eq. (dmdTheta) and the
  surrounding text.
* **Aether Coulomb term** (S43 / S61, S67). "The linear Coulomb term cancels exactly" holds only
  up to (c14/2)(Θ′ − 1)∂_Θ(Θ²a_r) = O(M/Θ), 2-5 % of the Coulomb swing at Θ = 75-200. Now
  Eq. (dm-coulomb) keeps the Θ′ factor and the text says "up to O(M/Θ)".
* **Speed claim** (S38 / S54b, S54c). The precision of the spin-0 speed was claimed for both the
  probe and the luminosity lag; the lag is 2 % low. The paper now quotes the probe's precision and
  the lag's 2 % separately. A later wording, "afterwards the deviation stays below ~10⁻⁴", was
  exceeded by most v0 = 1 samples after t = 40 (largest 2.96 × 10⁻⁴); the paper prints 3 × 10⁻⁴.
* **Luminosity burst** (S38 / S54d). "On every sphere a sharp burst at u ≃ 10": on the record the
  burst is at u = 9.7 … 14.7 on Θ = 120 … 400 (the drift is the 2 % lag itself). Now "u ≃ 10-15".
* **Derivation of the luminosity** (S39-S40 / S55, S56b). "Using the evolution equations and the
  Hamiltonian constraint with Eq. (ms) gives L" is not sufficient: m depends on Q, so dm/dt needs
  ∂_t Q (absent from the evolution equations), and Eq. (lum) contains the radial pressure. Now
  "Using the field equations G_ab = T_ab".
* **Wording and misprints** (one finding each): a stray "=" in "the combination of the coupling
  parameters c1 + c3" (S01 / S02); the entering condition for the scalar pair, which the code imposes
  through ∂_t P alone, equivalent up to terms that vanish for outgoing waves (S35 / S51); a sentence
  attributing the eigenvectors of the spin-0 matrix to the scalar pair as well, removed;
  the reflection with the entering field set to zero, 2.4 % → 2.3 % (linear theory 2.33 %; S36 / S53);
  the fraction 0.144 M_ADM of Eq. (Einf), which needs the M_ADM = 3.708 of the c2 = 0 run, now printed
  (S47 / S66); a convergence claim of the conclusion restricted to where the solution is resolved (the
  claims are now made and checked in the body, S20-S21 / S27-S30); the factor 2/3 of the
  trace-free part B̃^rr_LTF (S60 / S78); the white-hole width sentence, which concluded from "levels
  off near 0.003" that a narrower structure could not be excluded (a non sequitur), rewritten;
  C₁ defined; "grows roughly linearly" for an amplitude growing as Θ^0.7 (S44 / S62);
  "radius-independent radiated energy" qualified.
* **Numbers aligned with the run records** (DATA checks in both audits): the moving-mesh and fixed-grid
  M_ADM (3.7111 vs 3.7117) explained in the text; the frozen-centre ladder given with its five points
  (dr^(−1/3)); the radiated energy split into spin-0 and scalar parts with the subtraction exact in
  the printed digits (E_∞ = 0.532 ± 0.023, E_ψ = 0.197, E_æ = 0.335); the A0 = 0.03 central artefact
  given per grid; the Θ = 250 sphere leaving the fixed grid at t ≃ 20, ρ_æ > 0 only for
  18 ≲ Θ ≲ 20 at t = 28 and the extrapolated mass drift "about 0.01 %" (S21 / S30b); the
  anti-trapped layer widths "within 20 %" and the light ray "0.3 behind at t = 34" (S53-S54 / S72b,
  S73e); the cumulative energy "nearly the same on all spheres once the pulse has passed"
  (S48 / S67b).

One statement was kept by the authors and is checked as they read it: "roughly one areal unit per
unit time" for the fall of the fixed-grid outer edge is the average over the run of Fig. mesh(a),
(284 − 124)/150 = 1.07 (S14 / S20; the early rate, 1.70 up to t = 20, is reported as information).

## Finding that corrected the code

* **Ingoing initial data** (SymPy S51). `set_initial_data` solved the Hamiltonian constraint for Θ
  with Q = 0 and computed Q afterwards; for the ingoing shell P ≠ 0 sources Q, so the data violated
  the constraint at t = 0 by an amount that did not shrink with resolution. The Θ and Q solves are
  now iterated to a fixed point (and P, Q, a_r, K are evaluated at the RK4 stage radii), which makes
  the t = 0 residual second order; every white-hole result of the paper was recomputed with the
  corrected data (`release/VALIDATION.md`, "Initial-data fix" and "Ingoing initial data").

## Text against code

The code equations used by the audits are transcribed by hand from `release/collapse.cpp` and
re-checked against it line by line whenever the code changed (code lines are written
`collapse.cpp:L…`). Two text-versus-code inconsistencies inside the audits themselves were fixed: a
check of dt = 0.1 h that was a tautology is now a symbolic check that the error C₁dr² + C₂dt² is
O(h²) (SymPy S18), and the convergence bounds stated in the SymPy S20 text are now the ones the code
tests.

## Status

Against the `main.tex` of the arXiv v1 source: SymPy full run PASS 140, DATA 58, INFO 20,
CAVEAT 0, FAIL 0; Mathematica 93 PASS, 0 FAIL, 0 CAVEAT, 3 SKIP over 96 sections, with all 152
DATA snippets found at their quoted lines (S81). The notebook `check_equations.nb` gives the same
status lines as the script.
