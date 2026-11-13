#!/usr/bin/env python3
r"""
check_equations.py -- equation-by-equation audit of the paper
    "Gravitational collapse in Lorentz-violating gravity and the spin-0 mode"
    (Albertini, Vaglio, Barausse).

    Paper source version: main.tex of the paper's LaTeX source (arXiv v1).  Every line
    number quoted below (L...) refers to that file.

PURPOSE
-------
Every displayed equation, every labelled equation (eq:action, eq:speeds, eq:lapse,
eq:metric, eq:evolution, eq:radial, eq:hamiltonian, eq:constraint-prop, eq:id, eq:expansions, eq:ms, eq:dmdTheta,
eq:mratio, eq:principal, eq:aNbN, eq:wout, eq:asymptotic,
eq:wout-asym, eq:bc-win, eq:wbar, eq:retarded, eq:lum, eq:dm-coulomb,
eq:mc, eq:Mae, eq:Einf, eq:pn_first_surviving, eq:G_exact, eq:inverse_propagator,
eq:wh-id, eq:simple_scalar, eq:simple_L, eq:Tp, eq:Lprime, eq:Lprime_p1, eq:x_def,
eq:det_M0, eq:vL_def, eq:sminus1, eq:pn_1_term, eq:N_transform, eq:M1, eq:s1,
eq:kappa1), every inline relation that carries physics and every number that can be
checked by arithmetic gets one section below (S00-S61, with S08b, S08c), in the order of the paper.  Each
section is one function.  It quotes the LaTeX verbatim with its main.tex line number
("..." marks omitted text; quotes of figure captions and of long paragraphs are marked
"abridged"), says what is checked and how (in the docstring and, for most sections, in a
printed HOW paragraph), and prints one status line per check:

  PASS    derived from first principles with SymPy or proven as an exact identity (some
          identities in exact rational arithmetic at random points: S02, S10); also used
          for a statement reproduced numerically (the constraint solve of S51, the
          reflection estimate of S36, the shell widths of S57)
  FAIL    the statement as printed is wrong, or inconsistent with another statement
          of the paper; the correct form is printed
  CAVEAT  correct only under an assumption the paper does not state, worded too
          strongly, or a minor misprint; the precise statement is printed
  DATA    a simulation number: only its arithmetic and its consistency with the other
          numbers of the paper are checked (all such numbers live in the DATA table); a
          DATA line without a comparison just records the number and the arithmetic
          derived from it
  INFO    not a check of a statement of the paper: a literature input or a definition
          (not checked here), a slow derivation skipped in --quick mode, or a check of
          the simulation code against the paper's equations that passed (S35, marked
          "(code check)"; such a check prints FAIL if it fails)

What is NOT checked: the simulation numbers themselves (they are registered in the DATA
table and cross-checked against each other, never re-computed), and literature inputs
(the mode energies, the GW170817 bound, the results of Akhoury-Garfinkle-Gupta and of
Taherasghari-Will), which are recorded as INFO.

THE DATA TABLE
--------------
Every simulation number quoted in the paper is registered once, in the DATA table right
below the helpers, as  _d(key, value, main.tex line, description).  The DATA checks read
from that table only (dval(key) = the value, dline(key) = its main.tex line), so a number that changes in the paper is updated in one place;
section S00 prints the whole table at the start of every run.

HOW TO RUN
----------
    python3 check_equations.py            full audit
    python3 check_equations.py --quick    skips the three slow derivations (S04, S10,
                                          S39), which then report INFO

Requirements: Python 3 with SymPy (tested with 1.14) and NumPy (tested with 2.5); no
other dependency.  The script is self-contained: it reads nothing from disk (neither the
paper nor the simulation code) and writes no files.  The three slow derivations (S04
plane-wave speeds, S10 field equations from the action, S39 Misner-Sharp identity) run in
parallel worker processes while the fast sections print; the output of each is shown in
its own section.  Run time: about 1.5 min for the full audit (three cores), about 1 min
with --quick.  The output ends with a SUMMARY block that lists every FAIL and
CAVEAT and gives the counts by status.  Exit code: 0 if no check FAILs, 1 otherwise.

CONVENTIONS (identical to release/collapse.cpp and to Garfinkle-Eling-Jacobson gr-qc/0703093)
----------------------------------------------------------------------------------------
* signature (-,+,+,+), units G = c = 1.
* The aether is unit and timelike, u.u = -1, and hypersurface orthogonal: u_a = -N d_a t,
  i.e. it is the future unit normal n of the t = const slices.  The action is
      S = 1/(16 pi G) int sqrt(-g) [ R - M^{ab}_{mn} D_a u^m D_b u^n + lambda (u.u + 1) ]
          + S_matter,
      M^{ab}_{mn} = c1 g^{ab} g_{mn} + c2 d^a_m d^b_n + c3 d^a_n d^b_m - c4 u^a u^b g_{mn}.
* Normalisation of the field equations used by the code, by Eq. (lum) and by every
  derivation below (stated in the paper at L313-316):
      G_ab = T^ae_ab + T^psi_ab ,   T^psi_ab = d_a psi d_b psi - (1/2) g_ab (d psi)^2,
  i.e. S_matter = -(1/16 pi G) int sqrt(-g) (d psi)^2, so the code's psi is sqrt(8 pi G)
  times a canonically normalised scalar.  With this normalisation
  rho_psi = (P^2 + psi'^2)/2 and the Hamiltonian constraint, Eq. (hamiltonian), has
  (Theta^2/4)(P^2 + psi'^2).
* Metric in the gamma_rr = 1 gauge, Eq. (metric):  ds^2 = -N^2 dt^2 + (dr + beta dt)^2 + Theta^2 dOmega^2,
  so  d_t = N u + beta d_r,  e = d_r is the unit radial vector, and for any field
  d_t f = beta f' + N u(f)   (this is how every RHS of release/collapse.cpp is written).
* Extrinsic curvature: D_a u_b = -K_ab - u_a a_b,  K_ij = -(1/2N)(d_t g_ij - D_i b_j - D_j b_i),
  trK = K = -D_a u^a (K > 0 = converging), and the shear Q is
      K^r_r = K/3 + Q ,   K^theta_theta = K/3 - Q/2 ,
  so  A == u(Theta) = Theta (Q/2 - K/3)  (the Theta evolution equation) and  beta' = N(Q + K/3)
  (the shift equation).  P = u(psi);  a_r = d_r ln N.
* Couplings of the runs: c1 = 0.7, c2 = -0.5119403, c3 = c4 = 0 (v0 = 1); c2 = 0 (v0 = 1.266862, printed 1.267);
  white holes c1 = 0.8, c2 = -0.7111.  All symbolic derivations keep c1..c4 GENERAL, except
  Sec. V.E (sec:pn_comparison) and App. B (S49, S59-S61), which the paper writes for c13 = 0 (c3 = -c1).

CODE TRANSCRIPTION
------------------
The "code equations" used below (R_P, R_K, R_a, C_H, C_M in the block after the spherical
machinery; the spin-0 coefficients in S27; the outer-boundary formula in S35; the quintic
taper in S14; the initial-data profile in S16 and S51) were transcribed by hand from
release/collapse.cpp (the file in this repository) and checked against it
line by line.  Code lines are written collapse.cpp:L..., never bare L... (which always
means main.tex): header comment collapse.cpp:L44-76 (radial relations L44-49, evolution
equations L51-60, aN/bN L62-63, Hamiltonian constraint L71-76); compute_Q
collapse.cpp:L699-735 (momentum constraint); compute_RHS collapse.cpp:L817-930 (couplings
den_K = 2+c13+3c2, den_a = c14(1-c13), set once in main at L1607-1609; P: L872-873, trK:
L875-880, a_r: L883-886, Theta: L888-889; outer boundary L892-929, spin-0 pair L901-915,
scalar pair L917-928); the outer-boundary gate in ICN_step collapse.cpp:L1001-1030;
theta_C collapse.cpp:L1305-1315 (Hamiltonian constraint); spin0_aN/spin0_bN
collapse.cpp:L1615-1616; quintic_taper collapse.cpp:L325-332; shell_psi/shell_dpsi
collapse.cpp:L163-171; set_initial_data collapse.cpp:L1208-1298 (ingoing P L1252-1261,
the Theta/Q iteration L1265-1280).  S10 then shows that these equations are
the field equations of Eq. (action); the script never reads the code.
"""
import contextlib
import io
import math
import multiprocessing
import random
import sys
import time
from collections import Counter, OrderedDict, namedtuple
from itertools import combinations

import mpmath
import numpy as np
import sympy as sp
from sympy.calculus.euler import euler_equations

QUICK = '--quick' in sys.argv
T_START = time.time()
RESULTS = []          # (section, line, status, one-line description)
PAPER_SOURCE = 'the LaTeX source of arXiv v1'   # version of main.tex whose line numbers are quoted


def report(sec, line, status, desc, detail=''):
    """Record and print one status line of section `sec` for main.tex line `line`."""
    RESULTS.append((sec, line, status, desc))
    tag = f'[{status}]'
    print(f'  {tag:9s} {desc}' + (f'\n            {detail}' if detail else ''))


def header(sec, line, title, latex):
    """Print the section banner and the LaTeX of main.tex quoted verbatim."""
    print('\n' + '=' * 92)
    print(f'{sec}  main.tex L{line}:  {title}' if line else f'{sec}  (all main.tex lines):  {title}')
    for ln in latex.strip('\n').split('\n'):
        print('    | ' + ln)
    print('-' * 92)


def how(text):
    """Print the HOW paragraph of a section."""
    print('  HOW: ' + text.strip('\n').replace('\n', '\n  '))


def is_zero(e):
    """True if the expression simplifies to zero."""
    return sp.simplify(sp.expand(e)) == 0


def elapsed():
    """Wall-clock time since the start of the run, as a short string."""
    return f'(time so far {time.time() - T_START:.0f} s)'


# =======================================================================================
# DATA TABLE: every simulation number quoted in the paper, with its main.tex line.
# The DATA checks read from this table only.  Refresh here when the paper changes.
# =======================================================================================
Entry = namedtuple('Entry', 'value line what')
DATA = OrderedDict()


def _d(key, value, line, what):
    DATA[key] = Entry(value, line, what)


# --- Sec. III.A, formulation, Fig. 1 (moving mesh) ----------------------------------------
_d('lapse_collapse',   (150, 8, 12), 401, 'by t ~ 150 the lapse collapses at Theta ~ 8, with Theta_AH ~ 12 (L401-402)')
_d('taper_start',      0.4, 452, 'ALE taper g supported for x > 0.4 r_max (the outer 60%)')
_d('mesh_run',         (300, 25600, 0.35), 461, 'Fig. 1 caption: r_max, N_r, A0 of the moving-mesh validation')
_d('t_ThAH_max',       28, 470, 'Fig. 1 caption: Theta_AH reaches its maximum at t ~ 28 (L470-471; also L731)')
_d('Theta_ADM_read',   250, 508, 'M_ADM is read at Theta >~ 250 (also L662, L686)')
# --- Sec. III.B, convergence: scattering ------------------------------------------------------
_d('dt_factor',        0.1, 535, 'dt = 0.1 h')
_d('scatter_A',        (0.01, 0.04, 0.06), 558, 'scattering amplitudes')
_d('scatter_N',        (6401, 12801, 25601), 559, 'scattering ladder (also L575, L1530)')
_d('scatter_rmax',     150, 559, 'scattering box')
_d('scatter_dr',       (0.023, 0.012, 0.006), 559, 'scattering spacings')
_d('N_min',            (0.98, 0.70, 0.4), 563, 'minimum of the central lapse during the bounce (also Table 1, L587-588)')
_d('scatter_window',   20, 567, 'norms over 0 < r < 20 (also L576, L1533, L1550)')
_d('scatter_t',        (12, 14), 568, 'orders averaged over t = 12-14 (also L576)')
_d('tab1', {0.01: (2.35, 2.52, 2.17, 2.16, 2.92, 2.16, 1.95),
            0.04: (2.62, 2.46, 2.20, 2.14, 3.36, 2.12, 1.90),
            0.06: (2.81, 1.99, 1.45, 1.47, 2.42, 1.42, 1.49)}, 586,
   'Table 1: p(psi, P, trK, a_r, Theta, Q, ||Theta C||) per amplitude (L586-588)')
_d('p_low',            (1.4, 1.5), 598, 'A0 = 0.06: aether variables and constraint converge with p ~ 1.4-1.5')
_d('ratio_t14',        (1.8, 2.5, 14), 602, 'A0 = 0.06, t = 14: constraint-norm factor on the coarser / finer pair (L601-603)')
_d('thr_N',            (12801, 102401), 608, 'threshold bisection: coarsest and finest of the four grids')
_d('thr_dr',           (0.023, 0.003), 608, 'their spacings')
_d('A_switch',         (0.135, 0.125), 614, 'smallest amplitude with an early horizon: coarsest grid / two finest (L614-615)')
_d('A_star_bracket',   (0.120, 0.125), 615, 'A_* lies between these (L615-616)')
_d('thr_t_form',       20, 624, 'Fig. 2 caption: horizon that forms when the shell reaches the centre (t <= 20)')
_d('sep_t',            12, 1543, 'App. A: aether curves of A0 = 0.06 separate after the bounce, t >~ 12 (also L1556)')
_d('conv_rescale',     (1, 4, 16), 1552, 'Fig. conv-constraint caption: rescaling factors N_r^2 normalised to the coarsest grid (L1551-1552)')
# --- Sec. III.B, convergence: collapse; App. A ------------------------------------------------
_d('collapse_N',       (12801, 25601, 51201), 631, 'collapse ladder (also L578, L1578, L1597)')
_d('collapse_rmax',    300, 631, 'collapse box')
_d('collapse_dr',      (0.023, 0.012, 0.006), 632, 'collapse spacings ("as in the scattering test")')
_d('collapse_tmax',    28, 633, 'the three runs are compared every two time units up to t = 28')
_d('AH_t2',            2, 634, 'an apparent horizon is present at t = 2 (also L1566)')
_d('collapse_window',  (2, 150), 639, 'norms over r_AH + 2 <= r <= 150 (also L579, L1577, L1596)')
_d('collapse_t_avg',   (4, 28), 642, 'orders averaged from t = 4 to 28 (also L580)')
_d('p_collapse',       (2.62, 2.55, 2.34, 2.38, 2.62, 2.42, 2.24), 590,
   'Table 1, last row: collapse orders p(psi, P, trK, a_r, Theta, Q, ||Theta C||)')
_d('p_collapse_range', (2.24, 2.62), 644, 'text: the collapse orders lie between 2.24 and 2.62')
_d('ratio_fine_text',  (8, 6), 650, 'finer pair: ratios close to 4 up to t ~ 8, then ~ 6')
_d('ThAH_factor',      4.3, 658, 'differences of Theta_AH between successive grids shrink by 4.3 at t = 28')
_d('ThAH_richardson',  12.70, 659, 'Richardson limit of Theta_AH (also L681)')
_d('ADM_t_leave',      20, 664, 'the sphere Theta = 250 leaves the fixed grid at t ~ 20 (record: Theta(r_max) = 253.9 at t = 18, 249.8 at t = 20)')
_d('mesh_N_fine',      51201, 665, 'moving-mesh run used for the ADM mass (also L1572)')
_d('Theta_edge_mesh',  284, 666, 'areal radius of the outer edge held by the moving mesh')
_d('ADM_drift_pct',    0.07, 666, 'growth of m on the sphere Theta = 250 between t = 2 and t = 28, per cent')
_d('drift_falloff',    (3, 100, 250), 672, 'the change of m between t = 2 and 28 falls off as Theta^-3 on Theta = 100-250 (L671-672)')
_d('drift_100_pct',    1.2, 672, 'the change of m on Theta = 100, per cent')
_d('drift_inf_pct',    0.01, 674, 'its extrapolation to Theta -> infinity is about this, per cent (L673-674)')
_d('Theta_edge_fixed28', 235, 668, 'areal radius of the outer edge of the fixed grid at t = 28')
_d('ADM_agree',        (1e-4, 2e-5), 669, 'on that sphere: three fixed grids agree to 1e-4, fixed vs moving mesh to 2e-5 (L669-670)')
_d('ThAH_t2',          8.5, 1566, 'App. A: Theta_AH at t = 2')
_d('lapse_t28',        1e-6, 1567, 'App. A: central lapse below 1e-6 by t = 28')
_d('rAH_range',        (9, 65), 1568, 'App. A: r_AH grows from ~ 9 (t = 2) to ~ 65 (t = 28) (L1568-1569)')
_d('Theta_r150',       (137, 90), 1570, 'App. A: areal radius at r = 150 from ~ 137 to ~ 90 (L1569-1570)')
_d('ThAH_ladder',      (12.293, 12.604, 12.677), 1571, 'App. A: Theta_AH at t = 28 on the three grids')
_d('ADM_t2_t28',       (3.70836, 3.71107), 1572, 'App. A: ADM mass at Theta = 250, moving mesh, t = 2 and t = 28 (L1572-1573)')
_d('ratio_t',          (4, 8, 12, 16, 20, 24, 28), 1583, 'Table 2: output times')
_d('ratio_row_coarse', (3.69, 3.80, 4.69, 4.53, 4.26, 4.05, 3.92), 1585, 'Table 2: ||ThC||_12801/||ThC||_25601')
_d('ratio_row_fine',   (3.95, 3.97, 5.59, 6.26, 6.12, 5.95, 5.82), 1586, 'Table 2: ||ThC||_25601/||ThC||_51201')
_d('finest_below_t',   10, 1598, 'Fig. conv-collapse caption: from t ~ 10 the finest-grid curve lies below the others (L1598-1599)')
# --- Sec. IV, mass function --------------------------------------------------------------------
_d('t_mass',           28, 679, 'Fig. 4 time (largest apparent horizon; also L717, L741)')
_d('Mhor',             6.35, 681, 'M_hor = Theta_AH/2 with the Richardson limit of Theta_AH')
_d('MADM',             3.7111, 683, 'M_ADM of the baseline run (also Fig. 4 caption, L743)')
_d('M_deficit',        -2.64, 687, 'net exterior energy M_ADM - M_hor (also L715)')
_d('ae_positive_band', (18, 20), 717, 'aether energy density positive only across the outgoing shell, 18 <~ Theta <~ 20 at t = 28 (record: rho_ae > 0 only at Theta = 18.3-20.3)')
_d('Mratio',           1.71, 725, 'Eq. (mratio)')
_d('ThAH_shrink',      (12.6, 11.8, 150), 732, 'Theta_AH decreases from 12.6 to 11.8 by t = 150 (grid of Fig. 1(b))')
_d('spike_Theta',      5, 743, 'Fig. 4 caption: discretisation spike inside the horizon at Theta ~ 5')
# --- Sec. V.B, boundary-condition tests --------------------------------------------------------
_d('refl_scalar_pct',  0.03, 882, 'linear scalar pulse reflected with 0.03% of its amplitude (flat space)')
_d('probe_bc',         (50, 16, 100, 150), 878, 'spin-0 test pulse trK ~ exp[-(r-50)^2/16], r_max = 100 and 150 (L878-879)')
_d('refl_new',         1e-5, 880, 'reflection into w_in with the outgoing-wave boundary condition: below 1e-5 of the incident w_out')
_d('refl_B0_pct',      (3.5, 2.3), 880, 'reflection with w_in set to zero (held at 0) at r_max = 100 and 150, per cent (L880-881)')
# --- Sec. V.C, results ---------------------------------------------------------------------------
_d('v0_super',         1.267, 479, 'spin-0 speed for c2 = 0 as printed (exact value 1.266862; also L879, L888, L947, L955, L963, L974, L985, L989)')
_d('sigma_filter',     10, 908, 'Gaussian running-mean width (also Fig. 3 caption, L929)')
_d('sigma_range',      (5, 15), 909, 'filter widths for which the pulse amplitude is the same on every sphere')
_d('probe_amp',        (1e-6, 150), 936, 'linear probe pulse of trK: amplitude, centre (L936-937)')
_d('probe_scatter_t',  40, 944, 'probe: scatter at early times, t <~ 40, from the splitting of the initial trK pulse')
_d('probe_bound',      3e-4, 947, 'probe: afterwards the deviation |v - v0|/v0 stays below 3e-4')
_d('probe_dev_text',   (0.01, 0.02), 947, 'body: the median agrees with v0 to 0.01% and 0.02% for v0 = 1 and 1.267')
_d('lum_burst_u',      (10, 15), 956, 'luminosity burst on every sphere at u ~ 10-15 (L955-956)')
_d('lum_spheres',      (120, 400), 958, 'clock spheres Theta = 120-400')
_d('lag_speed',        1.241, 963, 'luminosity-lag speed in tortoise distance (also Fig. lum caption, L989)')
_d('lag_dev_pct',      2, 963, 'quoted deviation of the lag speed from v0, per cent (also L989)')
_d('probe_dev',        (0.01, 0.02), 977, 'Fig. probe caption: medians agree with v0 to 0.01% (v0 = 1) and 0.02% (v0 = 1.267)')
# --- Sec. V.D, energy flux ---------------------------------------------------------------------
_d('flux_resid_pct',   0.02, 1017, 'flux identity: agreement within 0.02% of the peak flux')
_d('flux_spheres',     (40, 200), 1018, 'on every sphere from Theta = 40 to 200')
_d('Fpsi_max',         0.02, 1020, 'bound on |F_psi|')
_d('Fae_peak',         (0.9, 4), 1021, 'aether-flux peaks')
_d('fluxid_spheres',   (40, 100, 200), 1027, 'Fig. fluxid caption: spheres shown')
_d('corr_rate_factor', 10, 1115, '"more than an order of magnitude": corrected rate L = -d_u m_c of Fig. fluxinf(a) against the uncorrected\n'
   '    luminosity of Fig. lum(a) (L1114-1116)')
_d('Theta2ar_exp',     0.7, 1091, 'Fig. coulomb caption: amplitude of Theta^2 a_r grows roughly as Theta^0.7')
_d('Theta_fit',        (75, 200), 1100, 'spheres of the E_inf + b/Theta + O(1/Theta^2) fit')
_d('E_psi',            0.197, 1100, 'scalar part of E')
_d('Einf',             (0.532, 0.023), 1104, 'Eq. (Einf): E_inf and its uncertainty')
_d('Einf_frac',        (0.144, 0.006), 1104, 'Eq. (Einf): E_inf / M_ADM and its uncertainty')
_d('E_ae',             (0.335, 0.090), 1105, 'Eq. (Einf): E_ae = E_inf - E_psi and E_ae / M_ADM')
_d('pct_ae',           9, 1109, 'the spin-0 wave carries off 9% of the ADM mass (also L1453)')
_d('MADM_super',       3.708, 1098, 'M_ADM of the c2 = 0 radiated-energy run')
_d('Einf_caption',     (0.53, 0.02), 1135, 'Fig. fluxinf caption: E_inf = 0.53 +/- 0.02')
# --- Sec. VI, white holes ------------------------------------------------------------------------
_d('wh_c',             (0.8, -0.7111), 1255, 'white-hole couplings c1, c2 (L1255-1256)')
_d('wh_A_bh',          0.12, 1273, 'black-hole amplitude')
_d('wh_t_cross',       2.5, 1274, 'theta_+ crosses zero by t ~ 2.5 at A0 = 0.12')
_d('wh_A_wh',          0.10, 1275, 'white-hole amplitude')
_d('wh_timeline',      (6, 17, 20, 36, 40, 44), 1276, 'A0 = 0.10: dip of theta_+ t ~ 6-17, anti-trapped from t ~ 20, splits at t ~ 36, inner part shrinks after t ~ 40, gone by t ~ 44 (L1275-1279)')
_d('wh_layer_narrow',  (0.9, 40, 0.06, 80, 20), 1282, 'thin anti-trapped layer on the escaping pulse: width about 0.9 at t = 40 to 0.06 at t = 80 on the refined grids, which agree on its width to within 20% (L1282-1283)')
_d('agg_at_t',         17, 1287, 'AGG: anti-trapped surface first appears at t ~ 17')
_d('wh_compress_rate', 0.15, 1297, 'escaping pulse width shrinks as exp(-0.15 t)')
_d('wh_refine_factor', 8, 1301, 'moving refinement around the escaping pulse: spacing eight times finer, fixed factor, centre on the |P| peak (L1301-1305)')
_d('wh_Theta_shell',   6.9, 1316, 'areal radius of the inner shell')
_d('wh_t_steep',       29.5, 1316, 'time at which the shell steepens (also L1327, L1348, L1422)')
_d('wh_fig_t',         30, 1318, 'Fig. wh(b,c) time (also L1421, L1428)')
_d('wh_ngrids',        (7, 3, 4), 1318, 'seven grids: three uniform, four with a refined region (L1318-1319)')
_d('wh_shell_factor',  (8, 16), 1320, 'refined region on the inner shell: eight or sixteen times finer')
_d('wh_shell_hold',    (20, 30.5), 1322, 'that region is held at the shell formation radius from t = 20 until t = 30.5 (L1322-1323)')
_d('wh_width_exp',     0.8, 1325, 'peak width ~ dr^0.8 in each group of runs (uniform and refined)')
_d('wh_height_exp',    -0.5, 1326, 'peak height ~ dr^-0.5')
_d('wh_rho_exp',       -1, 1326, 'energy density ~ dr^-1')
_d('wh_nfeatures',     3, 1291, 'three sharp features (L1291-1292; also the Fig. sharp caption, L1400)')
_d('wh_ray',           (30, 3.2, 34, 0.3, 40, 0.1), 1310, 'the outgoing ray leaving the inner shell at t = 30, 3.2 behind the escaping-pulse peak, is 0.3 behind it at t = 34 and 0.1 at t = 40 (L1309-1311)')
_d('wh_ray_join',      34, 1424, 'Fig. wh caption: the ray joins the escaping pulse by t ~ 34 (L1423-1424)')
_d('wh_front',         (7, 0.3, 0.13, 3.2), 1332, 'third feature, inner ingoing scalar front: from t ~ 7, a jump of about 0.3 in P, along rdot = -beta - N, lapse N ~ 0.13 where it travels, stalls at Theta ~ 3.2, close to the inner edge of the anti-trapped region (L1332-1336)')
_d('wh_front_slope_t', (24, 40), 1338, 'its steepest slope changes little with resolution at t = 24 but doubles with each halving of dr by t = 40 (L1336-1338)')
_d('wh_front_fall',    (44, 54), 1341, 'the front falls to the centre at t ~ 44 and has dispersed by t ~ 54 (L1340-1341)')
_d('wh_freeze_exp',    -1 / 3, 1356, 'freezing time of the frozen-centre artefact ~ dr^(-1/3)')
_d('wh_other_data',    ((10, 3), (15, 4)), 1359, 'other (r0, sigma)')
_d('wh_disperse_s3',   (0.01, 0.015, 0.03), 1360, 'dispersing amplitudes for sigma = 3')
_d('wh_s3_A08',        (0.08, 0.009, 0.0047), 1365, 'sigma = 3, A0 = 0.08: trapped region that disappears at dr ~ 0.009, black hole at 0.0047 (L1365-1366)')
_d('wh_r15',           (0.028, 20, 24, 3), 1367, 'r0 = 15, A0 = 0.028: trapped gone by t ~ 20, anti-trapped from t ~ 24, on three grids (L1367-1370)')
_d('wh_s3_A03',        (0.03, 0.0031, 0.0016, 0.0008), 1372, 'sigma = 3, A0 = 0.03 central artefact at dr = 0.0031, not at 0.0016, again at 0.0008 (L1372-1373)')
_d('wh_s3_A04',        0.04, 1373, 'sigma = 3, A0 = 0.04: the artefact happens later on each finer grid')
_d('agg_t',            32, 1385, 'AGG: sharp features near t ~ 32')
_d('agg_N',            (40000, 80000), 1386, 'AGG grid points')
_d('agg_sigma',        2, 1391, 'AGG shell sigma')
_d('agg_width_ratio',  2.5, 1391, 'AGG shell is about 2.5 times wider than ours (L1391-1392)')
_d('agg_Theta_shell',  7, 1393, 'our feature sits at Theta ~ 7')
_d('wh_fig_dr',        0.0008, 1416, 'Fig. wh(a): uniform grid dr = 0.0008')
_d('wh_layer_width',   1, 1424, 'Fig. wh(a): the outermost theta_- = 0 branch is an anti-trapped layer less than one unit wide (L1423-1426)')
_d('sharp_fig',        (40, 0.0008), 1401, 'Fig. sharp(a): |P| at t = 40 on the uniform grid dr = 0.0008 (L1401-1402)')
_d('sharp_esc',        (40, 80, 21, 62, 8), 1402, 'Fig. sharp(b): escaping pulse at t = 40-80, peak from Theta ~ 21 to 62, region eight times finer (L1402-1404)')
_d('sharp_shell',      (30, (0.0031, 0.0016, 0.0008), (8, 16)), 1405, 'Fig. sharp(c): inner shell at t = 30 on dr = 0.0031, 0.0016, 0.0008 and F = 8, 16 (L1405-1407)')
_d('sharp_front_t',    40, 1408, 'Fig. sharp(d): inner front at t = 40 on the three uniform grids (L1407-1408)')
# --- Conclusion ------------------------------------------------------------------------------------
_d('concl_split',      (14, 9, 5), 1452, 'conclusion: 14% of the ADM mass, 9% in the spin-0 wave, 5% in the escaping scalar (L1452-1453)')

# Not printed in the paper (used by S14): the areal radius of the fixed-grid outer edge at t = 150 in Fig. 1(a),
# 124 (read off the figure; the value is printed neither in the text nor in the caption).
FIG1_EDGE_T150_UNPRINTED = (124, 150)
# Not printed (used by S21): the measurement behind L666-674 (analysis of the profiles of the N_r = 51201
# moving-mesh run): m(t=28)/m(t=2) - 1 in per cent on Theta = 100, 150, 200, 250, and the extrapolations of the
# fits A + B/Theta^3 and A + B/Theta^p (p free).
DRIFT_MEASURED_UNPRINTED = {'pct': {100: 1.2166, 150: 0.3460, 200: 0.1436, 250: 0.0729}, 'A_p3': -0.0101, 'A_pfree': 0.0037}
# Not printed (used by S54): the refinement recipe of the white-hole runs in release/README.md (zone A on the escaping
# pulse --ref_F=8; zone B on the inner shell --refB_F=8 / 16 --refB_ton=20 --refB_w=1.5 --refB_x0=9.9 --refB_ttrack=30.5;
# the shell forms at t ~ 29.5, r ~ 9.96), and the theta_- = 0 zeros of the outermost anti-trapped layer and the |P| peak
# of the escaping pulse on run box40_N51200 (the uniform dr = 0.0008 grid of Fig. wh(a); a figure run) at t = 40 and 44.
WH_REFINE_RECIPE_UNPRINTED = {'zoneA_F': 8, 'zoneB_F': (8, 16), 'refB_ton': 20, 'refB_ttrack': 30.5, 'refB_x0': 9.9,
                              'refB_w': 1.5, 'r_form': 9.96}
WH_LAYER_MEASURED_UNPRINTED = {40: (20.74, 21.41, 21.10), 44: (25.47, 25.90, 25.69)}   # t: (Theta_in, Theta_out, Theta of max |P|)
# Not printed (used by S53): the theta_- = 0 zeros (r > 1) on box40_N51200 at t = 34-44 (the region splits between t = 34
# and 36), and on the refined escaping-pulse runs whconv_N12800 / 25600 / 51200 (zone A, --ref_F=8, r_max = 150) the width
# of the layer (outermost theta_- = 0 pair at Theta > 15) and min theta_- there, at t = 40, 64, 80 (coarse to fine).
# Not printed (used by S54): release/analysis/fig_white_hole.py on box40_N51200, data/white_hole.npz: Theta of the
# escaping-pulse peak (largest local max of |P| beyond the shell) minus Theta of the outgoing ray leaving the inner shell at
# t = 30, at t = 30, 34, 40; and the inner edge of the anti-trapped region (innermost theta_- = 0, r > 1) on box40_N51200.
RAY_GAP_UNPRINTED = {30: 3.158, 34: 0.302, 36: 0.171, 40: 0.091}
# App. A ladders (S20, S21), from the command lines of the figure runs (the release/README.md recipes) and the
# defaults of release/analysis/fig_conv_scatter.py / fig_conv_collapse.py:
#   A0.01|0.04|0.06_N6401|12801|25601: --N=<N> --A0=<A0> --t_final=14.5 --rmax=150 --ale=0;  fig_conv_scatter.py: norm window
#   --rlo 0 --rhi 20, fields ||f_h - f_h/2|| against 4 ||f_h/2 - f_h/4||, Theta C rescaled by stride^2 (1, 4, 16);
#   col_N12801|25601|51201: --N=<N> --A0=0.35 --t_final=28 --rmax=300 --ale=0;  fig_conv_collapse.py: window [r_AH + 2,
#   --rmax_norm 150], --ratio_times 4,8,...,28;  mass_hi (ADM drift, --mesh_run): --N=51201 --A0=0.35 --rmax=300 --ale=1.
CONV_RUNS_UNPRINTED = {'scatter': {'A0': (0.01, 0.04, 0.06), 'N': (6401, 12801, 25601), 'rmax': 150, 'ale': 0, 'window': (0, 20), 'field_factor': 4,
                                   'rescale': (1, 4, 16)},
                       'collapse': {'A0': 0.35, 'N': (12801, 25601, 51201), 'rmax': 300, 'ale': 0, 'window_lo_plus': 2, 'window_hi': 150,
                                    'ratio_times': (4, 8, 12, 16, 20, 24, 28)},
                       'mesh': {'N': 51201, 'A0': 0.35, 'rmax': 300, 'ale': 1}}
# Eq. (wbar) against the code (S37): release/analysis/spin0_common.py gaussian_running_mean is
#   w_bar(t_i) = sum_s exp(-(t_s - t_i)^2/2 sigma^2) dt_s w(t_s) / sum_s exp(-(t_s - t_i)^2/2 sigma^2) dt_s,
#   s over the dumps with |t_s - t_i| <= truncate*sigma (truncate = 3), dt = np.gradient(t);
# fig_wout_background.py calls it per sphere (fixed Theta) with --sigma (default 10; the release README recipe passes none;
# --sigma_scan default 5,7.5,10,12.5,15) and plots only where t <= t_end - 3 sigma.  Record (data/wout_background.npz of
# the figure runs, sigma = 10 stored, dumps t = 0.50-300.00 every 0.50, background plotted for t <= 269.5; recomputed here off
# line): max over t of |Theta (w_bar_code - w_bar_eq)|, w_bar_eq = Eq. (wbar) with untruncated trapezoidal integrals over
# the whole data range, on Theta = 50, 70, 90, 110, 130, 150, 180, against the pulse amplitude max |Theta delta w_out|.
WBAR_CODE_UNPRINTED = {'sigma': 10.0, 'truncate': 3.0, 'dt': 0.5,
                       'maxdiff': (3.6e-4, 3.3e-4, 3.2e-4, 2.9e-4, 2.9e-4, 3.4e-4, 2.6e-4),
                       'pulse': (0.0727, 0.0706, 0.0706, 0.0709, 0.0716, 0.0723, 0.0734)}
# Not printed (used by S38): the probe record (release/analysis/fig_probe_speed.py on the figure runs: logs/probe.log and
# data/probe_speed.npz): (median - v0)/v0 in per cent, max |v - v0|/v0 for t >= 40 and for t < 40, first sample time and
# number of samples, for v0 = 1 and v0 = 1.266862.
PROBE_MEASURED_UNPRINTED = {'med_pct': (0.0122, 0.0154), 'late_max': (2.96e-4, 1.89e-4), 'early_max': (3.84e-4, 3.22e-4),
                            't_first': 16.0, 'n': (134, 127)}
# Not printed (used by S38): the luminosity clock (release/analysis/fig_luminosity_clock.py on the figure run sl600_N25601,
# logs/lum_clock.log): spheres, Theta_* (M_ADM = 3.7084 from the first dump), burst arrival t_arr = time of the deepest
# minimum of L, the fitted speed and the rms residual of the lag fit; and, from data/invariant_wout_dm.npz (the data of the
# supplementary Fig. invariant of fig_invariant.py, not in the paper; sl600_N51201), the retarded time of the steepest point of the narrow pulse of Theta w_out.
LUM_CLOCK_UNPRINTED = {'Theta': (120, 160, 200, 240, 300, 350, 400), 'MADM': 3.7084,
                       'Theta_star': (140.17, 182.43, 224.15, 265.55, 327.26, 378.43, 429.44),
                       't_arr': (120.385, 155.026, 188.910, 222.333, 271.941, 312.941, 353.711), 'v_fit': 1.2406, 'rms': 0.313,
                       'wout_pulse_u': {120: 9.85, 180: 11.49, 240: 12.88, 300: 13.68, 360: 14.23}}
# Not printed (used by S44, S48): the energy record of the c2 = 0 figure run.  data/aether_coulomb.npz (fig_coulomb.py,
# sl600_N51201, Fig. coulomb(b) thin lines m - M_ADM): the swing (max - min over u = -60 ... 260) of m, m - M_ADM averaged over
# u = 200-250, its minimum for u >= 120, the last u at which m < M_ADM, and the integral of L = -d_u m over the record,
# -(m - M_ADM)|_end + (m - M_ADM)|_start; the uncorrected max |L| = max |d_u m| on the same spheres (np.gradient over u).
# data/spin0_flux_infinity.npz (fig_energy.py, Fig. fluxinf): max |L| of the corrected rate L = -d_u m_c (key L); logs/energy.log:
# E averaged over u = 200-250 on each sphere and the b of E_inf + b/Theta; the spread max - min of E over the spheres for
# u = 10-95 (during the crossing, largest at u = 10.5) and for u >= 110; E_inf(u) (key Einf_u), plotted for u >= u_fit = 95:
# its maximum (the overshoot) and its band for u = 200-255 (u_max of the record).  data/luminosity_clock.npz (Fig. lum(a),
# sl600_N25601): max |L| of the uncorrected luminosity on Theta = 120 ... 400.
ENERGY_TIME_UNPRINTED = {
    'Theta': (75, 100, 150, 200),
    'm_swing': (7.05, 9.38, 14.05, 18.76), 'm_late': (1.171, 1.140, 1.094, 1.052), 'm_min_u120': (1.045, 0.962, 0.813, 0.673),
    'm_last_below_u': (97.4, 101.3, 105.3, 107.1), 'intL': (-1.239, -1.186, -1.119, -1.085),
    'L_uncorr_same': (1.409, 1.834, 2.732, 3.620), 'L_corr': (0.265, 0.222, 0.179, 0.151),
    'E_avg': (0.5945, 0.5786, 0.5632, 0.5557), 'E_b': 4.650,
    'E_spread_cross': (0.430, 10.5), 'E_spread_late': 0.044,
    'Einf_u_fit': 95, 'Einf_overshoot': (0.585, 133.0), 'Einf_band_late': (0.522, 0.540),
    'L_clock': {120: 2.31, 160: 3.09, 200: 3.83, 240: 4.54, 300: 5.77, 350: 6.72, 400: 7.51},
    'intL_from_t0': (-1.26, -1.32),          # -(m - M_ADM) at the end of the record, integrated from t = 0 (range over the spheres)
}
# Not printed (used by S36): release/VALIDATION.md (boundary-reflection tests, release/analysis/boundary_reflection.py):
# reflected w_in / incident w_out for (v0 = 1, R = 100, 150) and (v0 = 1.267, R = 100, 150); the same runs with the
# pre-fix code (w_in held fixed), in per cent, ((v0 = 1: R = 100, 150), (v0 = 1.267: R = 100, 150)); the flat-space
# scalar reflection.
BC_VALIDATION_UNPRINTED = {'spin0': (5.6e-7, 4.1e-7, 1.5e-6, 7.5e-7), 'prefix': ((3.45, 2.33), (3.53, 2.33)), 'scalar': 2.8e-4}
# Not printed (used by S54): per-family least-squares exponents of the rerun that made Fig. wh(b,c) (seven grids at
# t = 30: release/analysis/fig_white_hole.py and wh_formation_ladder.py, logs white_hole.log and formation.log):
# (width (FWHM), max |P|, max rho) ~ dr^(...) for the uniform and the refined family.
WH_EXPONENT_FITS_UNPRINTED = {'uniform': (0.84, -0.50, -0.99), 'refined': (0.79, -0.47, -0.95)}
FRONT_EDGE_UNPRINTED = {30: (3.22, 3.62), 34: (3.15, 3.29)}   # t: (Theta of the front, Theta of the inner theta_- = 0 edge)
WH_SPLIT_ZEROS_UNPRINTED = {34: (3.29, 14.49), 36: (2.87, 14.82, 15.17, 16.82), 40: (2.34, 14.39, 20.74, 21.41),
                            44: (13.28, 13.46, 25.47, 25.90)}
WH_LAYER_LADDER_UNPRINTED = {'width': {40: (0.897, 0.908, 0.839), 64: (0.112, 0.100, 0.095), 80: (0.061, 0.063, 0.068)},
                             'min_theta_minus': {40: (-0.21, -0.25, -0.28), 64: (-1.19, -1.50, -1.69), 80: (-0.85, -1.12, -1.36)}}
# Not printed (used by S54): release/analysis/fig_sharp_features.py on the figure runs box40_N12800 / 25600 / 51200
# (uniform, r_max = 40: dr = 0.0031 / 0.0016 / 0.0008), box40pp_N51200(_F16) (zone B, F = 8 / 16) and whconv_N51200 (zone A,
# F = 8, --ref_ton=40, t_final = 80): steepest |dP/dTheta| of the inner front in 1.5 < Theta < 5 at t = 24 and 40 (coarse to
# fine) and the |P| peak of the escaping pulse at t = 40 and 80.  The inner front tracked on box40_N51200 (steepest |dP/dTheta|
# in 1.5 < Theta < 5; spin0_common.load_profile, beta_shift): the jump of P across it (P extremum within 0.3 outside minus
# extremum within 0.3 inside), its dr/dt (centred differences over t +- 2) against -beta - N at the front, the lapse at the
# front, its Theta, and the steepest slope at t = 44 / 54; and on wh_wh_N25600 (r_max = 300, dr = 0.0117, dumps from t = 0.5)
# the early front.
SHARP_SCRIPT_UNPRINTED = {'slope': {24: (4.9, 6.2, 6.9), 40: (7.1, 14.8, 30.7)}, 'esc_peak': {40: 21.17, 80: 61.51}}
FRONT_MEASURED_UNPRINTED = {
    'jump': {20: -0.286, 22: -0.294, 24: -0.304, 26: -0.316, 28: -0.331, 30: -0.352},
    'jump_t32_36': (-0.386, -0.438, -0.513),
    'rdot': {22: -0.0268, 24: -0.0230, 26: -0.0203, 28: -0.0185, 30: -0.0180},
    'minus_beta_minus_N': {22: -0.0264, 24: -0.0226, 26: -0.0201, 28: -0.0183, 30: -0.0174},
    'N_front': {24: 0.125, 26: 0.128, 28: 0.137},
    'Theta_stall': {30: 3.220, 32: 3.180, 34: 3.153, 36: 3.164},
    'Theta_t44': 1.870, 'slope_t54': 0.002,
    'early_wh_wh': {4: (6.64, -0.64, -0.258), 6: (5.61, -0.88, -0.267), 8: (4.92, -1.17, -0.262)},   # t: (Theta, dP/dTheta, jump)
    'early_rdot': {6: (-0.375, -0.363), 8: (-0.240, -0.230), 10: (-0.155, -0.148)},                  # t: (dr/dt, -beta - N)
}


def dval(key):
    """Value of a DATA entry."""
    return DATA[key].value


def dline(key):
    """main.tex line of a DATA entry."""
    return DATA[key].line


# =======================================================================================
# Common symbols
# =======================================================================================
c1, c2, c3, c4 = sp.symbols('c1 c2 c3 c4', real=True)
c13, c14, c123 = c1 + c3, c1 + c4, c1 + c2 + c3
DEN = 2 + c13 + 3 * c2                      # 2 + c13 + 3 c2
aN_sym = (c14 - 2) / DEN                    # Eq. (aNbN)
bN_sym = -c123 / (c14 * (1 - c13))
v0sq_sym = c123 * (2 - c14) / (c14 * (1 - c13) * (2 + c13 + 3 * c2))   # Eq. (speeds)
v2sq_sym = 1 / (1 - c13)
v1sq_sym = (2 * c1 - c1 ** 2 + c3 ** 2) / (2 * c14 * (1 - c13))
BASE = {c1: sp.Rational(7, 10), c2: sp.Float('-0.5119403', 30), c3: 0, c4: 0}
SUPER = {c1: sp.Rational(7, 10), c2: 0, c3: 0, c4: 0}
WHITE = {c1: sp.Rational(4, 5), c2: sp.Rational(-32, 45), c3: 0, c4: 0}   # c2 = -c1^3/(2-4c1+3c1^2) at c1 = 4/5 (v0 = 1), printed -0.7111 (S50)
D = sp.diff


def aether_lagrangian_M(gl, gi, Du, u_up):
    """sum_{abmn} M^{ab}_{mn} D_a u^m D_b u^n with M built LITERALLY from Eq. (action)."""
    tot = 0
    for a in range(4):
        for b in range(4):
            for m in range(4):
                for n in range(4):
                    M = (c1 * gi[a, b] * gl[m, n]
                         + (c2 if (a == m and b == n) else 0)
                         + (c3 if (a == n and b == m) else 0)
                         - c4 * u_up[a] * u_up[b] * gl[m, n])
                    if M != 0:
                        tot += M * Du[a, m] * Du[b, n]
    return tot


def aether_lagrangian(gl, gi, Du, u_up):
    """Expanded form c1 Du.Du + c2 (div u)^2 + c3 D_a u^m D_m u^a - c4 a.a (faster)."""
    t1 = sum(gi[a, b] * gl[m, n] * Du[a, m] * Du[b, n]
             for a in range(4) for b in range(4) for m in range(4) for n in range(4)
             if gi[a, b] != 0 and gl[m, n] != 0)
    dv = sum(Du[a, a] for a in range(4))
    t3 = sum(Du[a, m] * Du[m, a] for a in range(4) for m in range(4))
    acc = [sum(u_up[a] * Du[a, m] for a in range(4)) for m in range(4)]
    t4 = sum(gl[m, n] * acc[m] * acc[n] for m in range(4) for n in range(4))
    return c1 * t1 + c2 * dv ** 2 + c3 * t3 - c4 * t4


def christoffel(g, gi, X, simp=True):
    """Christoffel symbols Gam[a][b][c] of the metric g (inverse gi) in coordinates X."""
    f = sp.simplify if simp else (lambda e: e)
    return [[[f(sum(gi[a, d] * (sp.diff(g[d, b], X[c]) + sp.diff(g[d, c], X[b])
                                 - sp.diff(g[b, c], X[d])) for d in range(4)) / 2)
              for c in range(4)] for b in range(4)] for a in range(4)]


def ricci(Gam, X, a, b):
    """Component R_ab of the Ricci tensor from the Christoffel symbols."""
    s = 0
    for m in range(4):
        s += sp.diff(Gam[m][a][b], X[m]) - sp.diff(Gam[m][a][m], X[b])
        for l in range(4):
            s += Gam[m][m][l] * Gam[l][a][b] - Gam[m][b][l] * Gam[l][a][m]
    return s


# =======================================================================================
# Spherical machinery (S03, S07-S11, S17, S23, S28, S35, S52): general spherically symmetric
# metric with lapse N, shift beta, radial scale L (gamma_rr = L^2; L -> 1 is the code
# gauge, kept general so that its field equation can be varied), areal radius Theta, and
# an aether tilted by eta relative to the slice normal: u = n + eta e (eta -> 0 after
# variation; only O(eta) is needed for the radial aether equation).
# =======================================================================================
t, r, th, ph = sp.symbols('t r theta phi', real=True)
XS = [t, r, th, ph]
Nf = sp.Function('N')(t, r)
bf = sp.Function('beta')(t, r)
Lf = sp.Function('L')(t, r)
Tf = sp.Function('Theta')(t, r)
pf = sp.Function('psi')(t, r)
ef = sp.Function('eta')(t, r)
EPS = sp.Symbol('EPS')
gS = sp.Matrix([[-Nf ** 2 + Lf ** 2 * bf ** 2, Lf ** 2 * bf, 0, 0],
                [Lf ** 2 * bf, Lf ** 2, 0, 0],
                [0, 0, Tf ** 2, 0],
                [0, 0, 0, Tf ** 2 * sp.sin(th) ** 2]])
giS = sp.simplify(gS.inv())
GamS = christoffel(gS, giS, XS)
n_up = sp.Matrix([1 / Nf, -bf / Nf, 0, 0])
e_up = sp.Matrix([0, 1 / Lf, 0, 0])


def cov_up(v):            # Du[a, b] = D_a v^b
    return sp.Matrix(4, 4, lambda a, b: sp.diff(v[b], XS[a])
                     + sum(GamS[b][a][c] * v[c] for c in range(4)))


u_tilt = n_up + EPS * ef * e_up
Du_tilt = cov_up(u_tilt)
Du0 = Du_tilt.subs(EPS, 0)

# code variables written through the metric (gamma_rr = 1 => L = 1)
Krr_c = D(bf, r) / Nf                                        # K^r_r = beta'/N
Ktt_c = -(D(Tf, t) - bf * D(Tf, r)) / (Nf * Tf)              # K^th_th
K_c = Krr_c + 2 * Ktt_c
Q_c = Krr_c - K_c / 3
P_c = (D(pf, t) - bf * D(pf, r)) / Nf
a_c = D(Nf, r) / Nf
A_c = Tf * (Q_c / 2 - K_c / 3)                               # u(Theta)
B_c = D(Tf, r)                                               # Theta'

# =======================================================================================
# The GEJ system as coded in release/collapse.cpp (compute_RHS collapse.cpp:L817-930, theta_C
# collapse.cpp:L1305-1315, compute_Q collapse.cpp:L699-735; transcribed by hand, see the docstring),
# written as residuals in the metric functions: used by S08b, S10 and S11.
# =======================================================================================
R_P = D(P_c, t) - (bf * D(P_c, r) + Nf * (K_c * P_c + a_c * D(pf, r) + D(pf, r, 2) + 2 * D(Tf, r) * D(pf, r) / Tf))
R_K = D(K_c, t) - (bf * D(K_c, r) + Nf * K_c ** 2 / 3 + Nf * ((c14 - 2) * (D(a_c, r) + 2 * D(Tf, r) * a_c / Tf + a_c ** 2)
                                                         + 2 * P_c ** 2 + 3 * (1 - c13) * Q_c ** 2) / DEN)
R_K_negctl = D(K_c, t) - (bf * D(K_c, r) + Nf * K_c ** 2 / 3 + Nf * ((c14 - 2) * (D(a_c, r) + 2 * D(Tf, r) * a_c / Tf + a_c ** 2)
                                                            + 2 * P_c ** 2 + 3 * Q_c ** 2) / DEN)
R_a = D(a_c, t) - (bf * D(a_c, r) + Nf * ((2 * K_c / 3 - Q_c) * a_c + c13 / (c14 * (1 - c13)) * P_c * D(pf, r)
                                          - c123 / (c14 * (1 - c13)) * D(K_c, r)))
C_H = (D(Tf, r, 2) + (D(Tf, r) ** 2 - 1) / (2 * Tf) + c14 * a_c * D(Tf, r)
       + Tf / 4 * (c14 * (2 * D(a_c, r) + a_c ** 2) + P_c ** 2 + D(pf, r) ** 2
                   + sp.Rational(3, 2) * (1 - c13) * Q_c ** 2 - DEN / 3 * K_c ** 2))
C_M = (1 - c13) * (D(Q_c, r) + 3 * Q_c * D(Tf, r) / Tf) - DEN / 3 * D(K_c, r) + P_c * D(pf, r)

# Eqs. (evolution) and (radial) of the paper (L344-362), typed in independently, as residuals (S08b, S10)
dtb = lambda f: D(f, t) - bf * D(f, r)                           # (d_t - beta d_r)
PR_psi = dtb(pf) - Nf * P_c
PR_Th = dtb(Tf) - Nf * Tf * (Q_c / 2 - K_c / 3)
PR_P = dtb(P_c) - Nf * (K_c * P_c + D(pf, r, 2) + (a_c + 2 * D(Tf, r) / Tf) * D(pf, r))
PR_P_negctl = dtb(P_c) - Nf * (K_c * P_c + D(pf, r, 2) + (a_c - 2 * D(Tf, r) / Tf) * D(pf, r))
PR_K = dtb(K_c) - Nf * (K_c ** 2 / 3 + ((c14 - 2) * (D(a_c, r) + a_c ** 2 + 2 * a_c * D(Tf, r) / Tf)
                                        + 2 * P_c ** 2 + 3 * (1 - c13) * Q_c ** 2) / (2 + c13 + 3 * c2))
PR_K_negctl = dtb(K_c) - Nf * (K_c ** 2 / 3 + ((c14 - 2) * (D(a_c, r) + a_c ** 2 + 2 * a_c * D(Tf, r) / Tf)
                                            + 2 * P_c ** 2 + 3 * Q_c ** 2) / (2 + c13 + 3 * c2))
PR_a = dtb(a_c) - Nf * ((2 * K_c / 3 - Q_c) * a_c + (c13 * P_c * D(pf, r) - c123 * D(K_c, r)) / (c14 * (1 - c13)))
PR_Q = (1 - c13) * (D(Q_c, r) + 3 * D(Tf, r) / Tf * Q_c) - ((2 + c13 + 3 * c2) / 3 * D(K_c, r) - P_c * D(pf, r))
PR_Q_negctl = (1 - c13) * (D(Q_c, r) + 3 * D(Tf, r) / Tf * Q_c) - ((2 + c13 + 3 * c2) / 3 * D(K_c, r) + P_c * D(pf, r))
PR_beta = D(bf, r) - Nf * (Q_c + K_c / 3)

# The same RHS in plain symbols (field values and first r-derivatives at a point)
Nsym, bsym = sp.symbols('N beta')
F = {nm: sp.Symbol(nm) for nm in ('Th', 'Thr', 'psr', 'psrr', 'P', 'K', 'Kr', 'a_r', 'dr_a_r', 'Q')}
RHS = {
    'psi': Nsym * F['P'],
    'P': Nsym * (F['K'] * F['P'] + F['a_r'] * F['psr'] + F['psrr'] + 2 * F['Thr'] * F['psr'] / F['Th']),
    'trK': Nsym * F['K'] ** 2 / 3 + Nsym * ((c14 - 2) * (F['dr_a_r'] + 2 * F['Thr'] * F['a_r'] / F['Th'] + F['a_r'] ** 2)
                                           + 2 * F['P'] ** 2 + 3 * (1 - c13) * F['Q'] ** 2) / DEN,
    'a_r': Nsym * ((2 * F['K'] / 3 - F['Q']) * F['a_r'] + c13 / (c14 * (1 - c13)) * F['P'] * F['psr']
                   - c123 / (c14 * (1 - c13)) * F['Kr']),
    'Theta': Nsym * F['Th'] * (F['Q'] / 2 - F['K'] / 3)}
ADV = {'psi': bsym * F['psr'], 'P': 0, 'trK': bsym * F['Kr'], 'a_r': bsym * F['dr_a_r'], 'Theta': bsym * F['Thr']}
# full RHS of trK and a_r including advection (S25, S26, S30)
rhsK = ADV['trK'] + RHS['trK']
rhsA = ADV['a_r'] + RHS['a_r']

# Generic principal-part matrix (S28, S29, S32, S33, S35, S38)
aN_, bN_ = sp.symbols('a_N b_N')
Mab = sp.Matrix([[bsym, Nsym * aN_], [Nsym * bN_, bsym]])
v0_ab = sp.sqrt(aN_ * bN_)

# Frame stress tensor symbols (S23, S40)
rho, jj, pr_, AA, BB = sp.symbols('rho j p A B', real=True)
Th_s = sp.symbols('Theta', positive=True)
Nn = sp.symbols('N', positive=True)
Pp, psp = sp.symbols('P psi_r', real=True)
rho_psi = (Pp ** 2 + psp ** 2) / 2
j_psi = -Pp * psp
rr = sp.symbols('r', positive=True)

# =======================================================================================
# The three slow derivations run in forked worker processes; their printed output is
# captured and shown, in order, in their own sections (S04, S10, S39).
# =======================================================================================


def _captured(fn):
    """Run fn, returning its printed output and the status lines it recorded."""
    buf = io.StringIO()
    n0 = len(RESULTS)
    with contextlib.redirect_stdout(buf):
        fn()
    return buf.getvalue(), RESULTS[n0:]


def _job_ms():
    """S39 worker: prove d_a m = (Theta^2/2)(G_a^b - delta_a^b G^(2)) d_b Theta for a general 2-metric."""
    ga, gb, gc = [sp.Function(s)(t, r) for s in ('g_tt', 'g_tr', 'g_rr')]
    g2 = sp.Matrix([[ga, gb, 0, 0], [gb, gc, 0, 0], [0, 0, Tf ** 2, 0], [0, 0, 0, Tf ** 2 * sp.sin(th) ** 2]])
    gi2 = sp.simplify(g2.inv())
    Gam2 = christoffel(g2, gi2, XS)
    Rl = sp.Matrix(4, 4, lambda p, q: ricci(Gam2, XS, p, q) if (p < 2 and q < 2) or p == q else 0)
    Rsc = sum(gi2[p, q] * Rl[p, q] for p in range(4) for q in range(4))
    Gm = gi2 * (Rl - Rsc * g2 / 2)                          # G^a_b
    G2 = Gm[0, 0] + Gm[1, 1]
    G4 = G2 + Gm[2, 2] + Gm[3, 3]
    dT = [D(Tf, xv_) for xv_ in XS]
    m2 = Tf / 2 * (1 - sum(gi2[p, q] * dT[p] * dT[q] for p in range(2) for q in range(2)))
    oks = []
    for A_ in (0, 1):
        res = sp.simplify(D(m2, XS[A_]) - Tf ** 2 / 2 * sum((Gm[B_, A_] - (G2 if A_ == B_ else 0)) * dT[B_] for B_ in range(2)))
        oks.append(res == 0)
    report('S39', 1004, 'PASS' if all(oks) else 'FAIL', 'd_a m = (Theta^2/2)(G_a^b - delta_a^b G^(2)) d_b Theta, a = t, r, any 2-metric\n'
           '            (at fixed Theta: dm/dt|_Theta = (Theta^2/2) G^R_T, the flux of T_ab through the sphere once G_ab = T_ab)')
    test = {ga: -(1 + t * r / 7 + r ** 2 / 30), gb: t / 5 + r / 11, gc: 1 + r * t / 13, Tf: r + t * r ** 2 / 17 + r ** 3 / 50}
    bad = (D(m2, r) - Tf ** 2 / 2 * sum((Gm[B_, 1] - (G4 if B_ == 1 else 0)) * dT[B_] for B_ in range(2)))
    val = bad.subs(test).doit().subs({t: sp.Rational(1, 3), r: sp.Rational(7, 5), th: sp.Rational(1, 2)})
    report('S39', 1004, 'PASS' if abs(complex(sp.N(val))) > 1e-8 else 'FAIL',
           f'negative control: 4-d trace instead of G^(2) gives residual {float(sp.re(sp.N(val))):.3g} != 0')


def _job_gej():
    """S10 worker: the coded GEJ evolution equations and constraints follow from Eq. (action)."""
    Lae0 = sp.simplify(aether_lagrangian(gS, giS, Du0, n_up))
    Lae1 = sp.simplify(sp.diff(aether_lagrangian(gS, giS, Du_tilt, u_tilt), EPS).subs(EPS, 0))
    kin = sp.simplify(sum(giS[a, b] * D(pf, XS[a]) * D(pf, XS[b]) for a in range(4) for b in range(4)))
    Rs = sp.simplify(sum(giS[a, b] * ricci(GamS, XS, a, b) for a in range(4) for b in range(4)))
    sg = Nf * Lf * Tf ** 2                    # sqrt(-g)/sin(theta)
    Lag = sg * (Rs - Lae0 - kin) - sg * Lae1
    EL = euler_equations(Lag, [Nf, bf, Lf, Tf, ef, pf], [t, r])
    EL = [sp.together(e.lhs.subs(ef, 0).doit().subs(Lf, 1).doit()) for e in EL]
    funcs = {'N': Nf, 'b': bf, 'T': Tf, 'p': pf}
    jet = {(nm, i, j): sp.Symbol(f'{nm}_{i}{j}') for nm in funcs for i in range(3) for j in range(3 - i)}

    def tojet(e):
        reps = []
        for (nm, i, j), s in jet.items():
            if i + j == 0:
                continue
            f = funcs[nm]
            var = [(t, i)] if j == 0 else ([(r, j)] if i == 0 else [(t, i), (r, j)])
            reps.append((i + j, sp.Derivative(f, *var), s))
        for _, dd, s in sorted(reps, key=lambda z: -z[0]):
            e = e.subs(dd, s)
        for nm, f in funcs.items():
            e = e.subs(f, jet[(nm, 0, 0)])
        return e
    ELj = [tojet(e) for e in EL]
    codej = {nm: tojet(ex.doit()) for nm, ex in
             (('P eq.', R_P), ('trK eq.', R_K), ('a_r eq.', R_a), ('Ham. constraint', C_H),
              ('mom. constraint', C_M), ('printed P (evolution)', PR_P), ('printed trK (evolution)', PR_K),
              ('printed a_r (evolution)', PR_a), ('printed Q (radial)', PR_Q),
              ('NEG. CONTROL trK eq. with 3Q^2', R_K_negctl), ('NEG. CONTROL printed radial +P psi\'', PR_Q_negctl))}
    unk = [jet[('T', 0, 2)], jet[('b', 0, 2)], jet[('T', 2, 0)], jet[('b', 1, 1)], jet[('N', 1, 1)], jet[('p', 2, 0)]]
    random.seed(7)
    worst = {k_: 0 for k_ in codej}
    ntrial = 0
    while ntrial < 3:
        vals = {s: sp.Rational(random.choice([-1, 1]) * random.randint(1, 9), random.randint(1, 7))
                for s in jet.values() if s not in unk}
        vals[jet[('N', 0, 0)]] = sp.Rational(random.randint(1, 9), 5)
        vals[jet[('T', 0, 0)]] = sp.Rational(random.randint(1, 9), 2)
        cv = {c1: sp.Rational(random.randint(1, 8), 10), c2: sp.Rational(random.randint(-9, 9), 10),
              c3: sp.Rational(random.randint(-5, 5), 10), c4: sp.Rational(random.randint(-5, 5), 10)}
        if any(sp.simplify(x_.subs(cv)) == 0 for x_ in (c14, 1 - c13, DEN, c123)):
            continue                                   # degenerate couplings: resample
        sol = sp.solve([sp.numer(sp.together(e.subs(vals).subs(cv))) for e in ELj], unk, dict=True)
        if len(sol) != 1:
            continue                                   # singular jet point: resample
        ntrial += 1
        for k_, v in codej.items():
            res = sp.nsimplify(sp.simplify(v.subs(vals).subs(cv).subs(sol[0])))
            worst[k_] = max(worst[k_], abs(float(res)))
    for k_, v in worst.items():
        if k_.startswith('NEG'):
            report('S10', 378, 'PASS' if v > 0 else 'FAIL', f'{k_}: residual {v:.3g} != 0 (the test has teeth)')
        else:
            report('S10', 378, 'PASS' if v == 0 else 'FAIL',
                   f'{"" if k_.startswith("printed") else "code "}{k_:24s} implied by the field equations of Eq. (action): max residual {v:.1e}')
    ratioH = sp.factor(sp.cancel(ELj[0] / codej['Ham. constraint']))
    print(f'    E_N / C_H = {ratioH}   (so C_H is exactly the Hamiltonian constraint, G = T normalisation)')
    print(f'    {elapsed()}')


def _job_speeds():
    """S04 worker: v2^2, v1^2, v0^2 from the quadratic action about Minkowski."""
    tt, xx, yy, zz = sp.symbols('t x y z', real=True)
    XP = [tt, xx, yy, zz]
    eps, w, k = sp.symbols('epsilon omega k')

    def lagr2(g, usp):
        eta0 = sp.diag(-1, 1, 1, 1)
        hh = g - eta0
        gi = sp.expand(eta0 - eta0 * hh * eta0 + eta0 * hh * eta0 * hh * eta0)   # exact to O(eps^2)
        Gam = christoffel(g, gi, XP, simp=False)
        R = 0
        for a in range(4):
            for b in range(4):
                if gi[a, b] != 0:
                    R += gi[a, b] * ricci(Gam, XP, a, b)
        ut, s1_, s2_ = sp.symbols('ut s1 s2')
        uu = sp.Matrix([ut] + list(usp))
        nrm = (uu.T * g * uu)[0] + 1
        ns = sp.series(nrm.subs(ut, 1 + eps * s1_ + eps ** 2 * s2_), eps, 0, 3).removeO()
        S1 = sp.solve(ns.coeff(eps, 1), s1_)[0]
        S2 = sp.solve(ns.coeff(eps, 2).subs(s1_, S1), s2_)[0]
        uu = uu.subs(ut, 1 + eps * S1 + eps ** 2 * S2)
        Du = sp.Matrix(4, 4, lambda a, b: sp.diff(uu[b], XP[a]) + sum(Gam[b][a][c] * uu[c] for c in range(4)))
        Lag = sp.series(sp.sqrt(-g.det()), eps, 0, 3).removeO() * (R - aether_lagrangian(g, gi, Du, uu))
        return sp.expand(sp.diff(Lag, eps, 2).subs(eps, 0) / 2)

    def dispersion(L2, fields, gauge0):
        E = [sp.expand(e.lhs) for e in euler_equations(L2, fields, [tt, zz])]
        live = [f for f in fields if f not in gauge0]
        amps = sp.symbols('A0:%d' % len(live))
        phs = sp.exp(sp.I * (k * zz - w * tt))
        rep = {f: Aa * phs for f, Aa in zip(live, amps)}
        rep.update({f: 0 for f in gauge0})
        rows = []
        for e in E:
            ee = sp.expand(sp.powsimp(sp.expand(e.subs(rep).doit() / phs)))
            rows.append([ee.coeff(Aa) for Aa in amps])
        Mx = sp.Matrix(rows)
        n = len(live)
        mins = [sp.factor(Mx.extract(list(rr_), list(range(n))).det()) for rr_ in combinations(range(Mx.rows), n)]
        mins = [m for m in mins if m != 0]
        gg = mins[0]
        for m in mins[1:]:
            gg = sp.gcd(gg, m)
        return sp.factor(gg)

    def obeys(poly, vsq):
        return is_zero(sp.numer(sp.together(poly.subs(w, k * sp.sqrt(vsq)))))

    hxy = sp.Function('h')(tt, zz)
    g = sp.diag(-1, 1, 1, 1); g[1, 2] = g[2, 1] = eps * hxy
    d2 = dispersion(lagr2(g, [0, 0, 0]), [hxy], [])
    report('S04', 261, 'PASS' if obeys(d2, v2sq_sym) else 'FAIL', f'spin-2: dispersion {d2} = 0  <=>  w^2 = k^2/(1-c13)')
    Af, Bf, vf = [sp.Function(s)(tt, zz) for s in ('A', 'B', 'v')]
    g = sp.diag(-1, 1, 1, 1); g[0, 1] = g[1, 0] = eps * Af; g[1, 3] = g[3, 1] = eps * Bf
    d1 = dispersion(lagr2(g, [eps * vf, 0, 0]), [Af, Bf, vf], [Bf])
    report('S04', 262, 'PASS' if obeys(d1, v1sq_sym) else 'FAIL', 'spin-1: dispersion relation  <=>  w^2/k^2 = (2c1 - c1^2 + c3^2)/(2 c14 (1-c13))')
    h00, h0z, hzz, hs = [sp.Function(s)(tt, zz) for s in ('h00', 'h0z', 'hzz', 'hs')]
    g = sp.diag(-1 + eps * h00, 1 + eps * hs, 1 + eps * hs, 1 + eps * hzz); g[0, 3] = g[3, 0] = eps * h0z
    d0 = dispersion(lagr2(g, [0, 0, eps * vf]), [h00, h0z, hzz, hs, vf], [h0z, hzz])
    report('S04', 263, 'PASS' if obeys(d0, v0sq_sym) else 'FAIL',
           'spin-0: dispersion relation  <=>  w^2/k^2 = c123(2-c14)/(c14(1-c13)(2+c13+3c2))', elapsed())


_FUT = {}


def start_slow_jobs():
    """Fork the three slow derivations so that they run while the fast sections print."""
    if QUICK:
        return
    pool = multiprocessing.get_context('fork').Pool(3)
    for fn in (_job_speeds, _job_gej, _job_ms):
        _FUT[fn.__name__] = pool.apply_async(_captured, (fn,))


def show_slow(name, sec, line, skipped_desc):
    """Print the captured output of a slow job (or an INFO line in --quick mode)."""
    if QUICK:
        report(sec, line, 'INFO', f'skipped in --quick mode ({skipped_desc})')
        return
    out, res = _FUT[name].get()
    print(out, end='')
    RESULTS.extend(res)


# =======================================================================================
# Sections, in the order of the paper
# =======================================================================================


def S00_data_table():
    """List the DATA table (every simulation number of the paper with its main.tex line)."""
    header('S00', 0, 'DATA table: simulation numbers quoted in the paper', 'see the DATA block at the top of the script')
    for k, e in DATA.items():
        print(f'    L{e.line:<5d} {k:22s} = {e.value!s:30s} {e.what}')
    report('S00', 0, 'INFO', f'{len(DATA)} simulation numbers registered; DATA checks below read from this table')


def S01_introduction():
    """Introduction, L128-131 and L140-143: the GW170817 sentence (the combination c1 + c3 = c13
    must vanish) and the anisotropic Horava dispersion relation, whose group velocity grows
    without bound."""
    header('S01', 128, 'Introduction: GW170817 sentence; Horava dispersion relation', r"""
GRB~170817A~\cite{Abbott:2017gw} fixed the speed of the spin-2 mode to
that of light to within a part in $10^{15}$. This requires the combination of the coupling parameters $c_1+c_3$ to
vanish~\cite{Oost:2018tcv,EmirGumrukcuoglu:2017cfa,Ramos:2018oku,Barausse:2019yuk}.            (L129-131)
gravity the higher spatial derivatives make the dispersion relations
anisotropic, $\omega^2=k^2+\alpha k^4+\beta k^6$ schematically, so the
group velocity grows without bound in the ultraviolet and no finite speed
limits the propagation of signals                                               (L140-143)
""")
    xx = sp.symbols('x')
    report('S01', 130, 'PASS' if sp.expand(sp.series(1 / sp.sqrt(1 - xx), xx, 0, 2).removeO() - 1 - xx / 2) == 0 else 'FAIL',
           'a part in 1e15 on v2 requires c1 + c3 = c13 ~ 1e-15: v2 = (1 - c13)^(-1/2) = 1 + c13/2 + O(c13^2) (see S06)')
    kk, al, be = sp.symbols('k alpha beta', positive=True)
    om = sp.sqrt(kk ** 2 + al * kk ** 4 + be * kk ** 6)
    vg = sp.diff(om, kk)
    lim = sp.limit(vg, kk, sp.oo)
    lead = sp.limit(vg / kk ** 2, kk, sp.oo)
    report('S01', 141, 'PASS' if lim == sp.oo and is_zero(lead - 3 * sp.sqrt(be)) else 'FAIL',
           f'omega^2 = k^2 + alpha k^4 + beta k^6 (beta > 0): v_g = d omega/dk ~ {lead} k^2 -> {lim}: "grows without bound in the ultraviolet"')


def S02_action():
    """Eq. (action), L236-248: the literal contraction M^{ab}_{mn} D_a u^m D_b u^n equals the
    expanded Lagrangian c1 (Du)^2 + c2 (D.u)^2 + c3 D_a u^m D_m u^a - c4 a^2, at random
    rational metrics, gradients and vectors (exact arithmetic).  Everything else is derived
    from this action, so this fixes the reading of M."""
    header('S02', 236, 'Eq. (action) and the tensor M^{ab}_{mn}', r"""
S={}&\frac{1}{16\pi G}\int d^4x\,\sqrt{-g}\,
\Bigl(R - M^{ab}{}_{mn}\nabla_a u^m\nabla_b u^n\\
&+ \lambda\,(u^a u_a+1)\Bigr) + S_{\rm matter},                                  (L238-240, split)
M^{ab}{}_{mn}=c_1 g^{ab}g_{mn}+c_2\delta^a_m\delta^b_n
+c_3\delta^a_n\delta^b_m-c_4u^au^bg_{mn},                                         (L246-247)
""")
    how('the literal contraction M^{ab}_{mn} D_a u^m D_b u^n is compared with the expanded form\n'
        'c1 D_a u_m D^a u^m + c2 (D.u)^2 + c3 D_a u^m D_m u^a - c4 a_m a^m (used by every derivation\n'
        'below) at 3 random rational metrics / gradients / vectors, exactly.')
    random.seed(12345)
    ok = True
    for _ in range(3):
        while True:
            A = sp.Matrix(4, 4, lambda i, j: sp.Rational(random.randint(-4, 4), random.randint(1, 3)))
            gl = sp.diag(-3, 2, 2, 2) + (A + A.T) / 10
            if gl.det() != 0:
                break
        gi = gl.inv()
        Du = sp.Matrix(4, 4, lambda i, j: sp.Rational(random.randint(-9, 9), random.randint(1, 5)))
        u_up = sp.Matrix([sp.Rational(random.randint(1, 9), 3)] + [sp.Rational(random.randint(-3, 3), 4) for _ in range(3)])
        ok &= (sp.expand(aether_lagrangian_M(gl, gi, Du, u_up) - aether_lagrangian(gl, gi, Du, u_up)) == 0)
    report('S02', 239, 'PASS' if ok else 'FAIL',
           'M-contraction == c1(Du)^2 + c2(D.u)^2 + c3 D_a u^m D_m u^a - c4 a^2 (exact, 3 random points)')


def S03_combinations():
    """L250-254, L477: with a hypersurface-orthogonal aether (u = n) the Lagrangian reduces to
    c13 K_ab K^ab + c2 K^2 - c14 a^2, so only c13, c14, c2 (hence c123) enter; the baseline
    couplings give c13 = c14 = 0.7."""
    header('S03', 250, 'c_ij combinations; hypersurface-orthogonal aether; c13 = c14 = 0.7', r"""
The four dimensionless couplings $c_i$ are also indicated through the combinations $c_{ij}\equiv c_i+c_j$ and
$c_{123}\equiv c_1+c_2+c_3$. When the \ae ther is hypersurface-orthogonal,
as it necessarily is in spherical symmetry, the twist term can be
removed ...                                                                       (L250-253)
$c_1=0.7$, $c_2=-0.5119403$, $c_3=c_4=0$ (so $c_{13}=c_{14}=0.7$)                 (L477)
""")
    how('build the aether Lagrangian from M for u = n in the general spherical metric and show\n'
        'it equals c13 K_ab K^ab + c2 K^2 - c14 a^2 (c3 is absorbed: every equation of the paper\n'
        'depends on c13, c14, c123 only).')
    Lae0 = sp.simplify(aether_lagrangian(gS, giS, Du0, n_up))
    Lae0_L1 = sp.simplify(Lae0.subs(Lf, 1).doit())
    KK = Krr_c ** 2 + 2 * Ktt_c ** 2
    target = c13 * KK + c2 * K_c ** 2 - c14 * a_c ** 2
    report('S03', 252, 'PASS' if is_zero(Lae0_L1 - target) else 'FAIL',
           'L_ae(u = n) = c13 K_ab K^ab + c2 K^2 - c14 a_r^2  (gamma_rr = 1, general c_i)')
    okc = (BASE[c1] + BASE[c3] == sp.Rational(7, 10) and BASE[c1] + BASE[c4] == sp.Rational(7, 10))
    report('S03', 477, 'PASS' if okc else 'FAIL', 'c1 = 0.7, c3 = c4 = 0  =>  c13 = c14 = 0.7')


def S04_speeds():
    """Eq. (speeds), L259-266: v2^2, v1^2, v0^2 derived from the second-order action about
    Minkowski (plane waves along z; Euler-Lagrange equations of all fields, gauge-fixed; the
    dispersion relation is the gcd of the maximal minors of the Fourier matrix).  Slow;
    runs in a worker."""
    header('S04', 259, 'Eq. (speeds): squared speeds of the spin-2, spin-1, spin-0 modes', r"""
v_2^2=\frac{1}{1-c_{13}},\qquad
v_1^2=\frac{2c_1-c_1^2+c_3^2}{2c_{14}(1-c_{13})},\\
v_0^2=\frac{c_{123}(2-c_{14})}{c_{14}(1-c_{13})(2+c_{13}+3c_2)} .                  (L261-263, gathered)
""")
    how('(derived) expand sqrt(-g)(R - M Du Du) to 2nd order about Minkowski with u = d_t,\n'
        'plane waves along z.  spin-2: h_xy.  spin-1: (h_tx, h_xz, v^x), gauge h_xz = 0.  spin-0:\n'
        '(h_tt, h_tz, h_zz, h_xx = h_yy, v^z), gauge h_tz = h_zz = 0; u^t from the unit norm.\n'
        'Euler-Lagrange equations of ALL fields (constraints included) -> Fourier matrix; the\n'
        'dispersion relation is the gcd of its maximal minors.')
    show_slow('_job_speeds', 'S04', 277, 'plane-wave derivation of v2, v1, v0')
    for lab, sub in (('c2 = -0.5119403', BASE), ('c2 = 0', SUPER), ('white hole c1 = 0.8, c2 = -32/45', WHITE)):
        print(f'    numerical values ({lab}): v2 = {float(sp.sqrt(v2sq_sym.subs(sub))):.4f}, '
              f'v1 = {float(sp.sqrt(v1sq_sym.subs(sub))):.4f}, v0 = {float(sp.sqrt(v0sq_sym.subs(sub))):.5f}')


def S05_spin0_energy_sign():
    """L267-271: the sign of the spin-0 energy is that of c14(2-c14) (literature); hence positive
    energy only if 0 < c14 < 2 (derived); value for the runs."""
    header('S05', 267, 'Sign of the spin-0 energy ~ c14 (2 - c14); positive only if 0 < c14 < 2', r"""
The energy densities of the linear modes were computed by
Eling~\cite{Eling:2005zq} with the Weinberg and Einstein energy--momentum
pseudotensors; that of the spin-0 mode has the sign of
$c_{14}(2-c_{14})$~\cite{Eling:2005zq,Jacobson:2008aj}, so the spin-0
mode carries positive energy only if $0<c_{14}<2$.                              (L267-271)
""")
    val = float((c14 * (2 - c14)).subs(BASE))
    val_wh = float((c14 * (2 - c14)).subs(WHITE))
    report('S05', 270, 'INFO', f'literature input (Eling:2005zq, Jacobson:2008aj), not re-derived; c14(2-c14) = {val:.2f} (c14 = 0.7), '
           f'{val_wh:.2f} (c14 = 0.8, white holes) > 0')
    xx = sp.symbols('x', real=True)
    sol = sp.solve_univariate_inequality(xx * (2 - xx) > 0, xx, relational=False)
    report('S05', 271, 'PASS' if sol == sp.Interval.open(0, 2) else 'FAIL', f'c14(2-c14) > 0  <=>  c14 in {sol}: "positive energy only if 0 < c14 < 2"')


def S06_gw170817():
    """L277-278 (and L128-131): GW170817 bounds v2 and hence |c13| <~ 1e-15, because v2 - 1 = c13/2 + O(c13^2)."""
    header('S06', 277, 'GW170817 bound', r"""
bounds the spin-2 speed and hence $|c_{13}|\lesssim10^{-15}$~\cite{Oost:2018tcv,EmirGumrukcuoglu:2017cfa,Gupta:2021vdj}.  (L278)
""")
    x = sp.symbols('x')
    ser_v2 = sp.series(1 / sp.sqrt(1 - x), x, 0, 2).removeO()
    report('S06', 278, 'PASS' if sp.expand(ser_v2 - (1 + x / 2)) == 0 else 'FAIL',
           'v2 - 1 = c13/2 + O(c13^2): |v2 - 1| = O(1e-15) <=> |c13| = O(1e-15) (order of magnitude; also "a part in 1e15", L130)')


def S07_lapse():
    """Eq. (lapse), L302-312: a_b = u^a D_a u_b computed covariantly for u = n in the general
    spherical metric (any gamma_rr) equals d_r ln N, and a.u = 0."""
    header('S07', 309, 'Eq. (lapse): kinematic identity a_i = D_i ln N', r"""
GEJ~\cite{Garfinkle:2007bk} write it as
$u_a=-N\nabla_a t$, i.e.,                                                        (L303-304)
a_i=D_i\ln N ,\qquad a^a\equiv u^b\nabla_b u^a .                                 (L310)
""")
    how('compute a_b = u^a D_a u_b covariantly for u = n (u_a = -N d_a t) in the general spherical\n'
        'metric (L arbitrary) and compare with d_r ln N; check a.u = 0.')
    u_lo0 = sp.simplify(gS * n_up)
    acc_lo = [sp.simplify(sum(gS[b, c] * sum(n_up[a] * Du0[a, c] for a in range(4)) for c in range(4))) for b in range(4)]
    okr = is_zero(acc_lo[1] - D(Nf, r) / Nf)
    oku = is_zero(sum(acc_lo[b] * n_up[b] for b in range(4)))
    report('S07', 310, 'PASS' if okr and oku else 'FAIL',
           'a_r = d_r ln N exactly (any gamma_rr), a.u = 0; u_a = ' + str(list(u_lo0)))
    return acc_lo


def S08_metric_and_variables():
    """Eq. (metric), L317-328, and L330-339, L363-365: gamma_rr = 1 and sphere area 4 pi Theta^2 for the
    metric; the five evolved fields; d_r beta = N(Q + trK/3) follows from K^r_r = beta'/N with
    K^r_r = K/3 + Q; beta(0) = 0 by parity."""
    header('S08', 317, 'Eq. (metric); evolved variables; Q; d_r ln N = a_r; d_r beta = N(Q + trK/3); beta(0) = 0', r"""
ds^2=-N^2(t,r)\,dt^2+(dr+\beta(t,r)\,dt)^2
+\Theta^2(t,r)\,d\Omega^2 ,                                                      (L318-319)
... the symmetry spheres have area $4\pi\Theta^2$ ... $\gamma_{rr}=1$.           (L323-328, abridged)
\bigl(\psi,\,P,\,\trK,\,\Theta,\,a_r\bigr),                                       (L332)
The shear $Q$ of
the extrinsic curvature, defined by $K^r{}_r=\trK/3+Q$ and
$K^\theta{}_\theta=K^\phi{}_\phi=\trK/3-Q/2$, the lapse $N$ and the
shift $\beta^r\equiv\beta$ are recovered on each slice.                          (L335-339)
\beta'=N\Bigl(Q+\frac{\trK}{3}\Bigr),                                              (L359)
and Eq.~\eqref{eq:lapse}, integrated outward from $Q(0)=0$ and from the
regularity condition $\beta(0)=0$, and inward from $N\to1$ at large
radius.                                                                           (L363-365)
""")
    how('read gamma_ij and the sphere metric off Eq. (metric) (the machinery metric at L = 1);\n'
        'with gamma_rr = 1 the rr extrinsic curvature is K_rr = -(1/2N)(0 - 2 D_r beta_r) = beta\'/N;\n'
        'the paper defines K^r_r = K/3 + Q, hence beta\' = N(Q + K/3) (checked from the covariant\n'
        'D_a u_b in S09).  beta is odd in r (radial component of a vector) -> beta(0) = 0.')
    g1 = gS.subs(Lf, 1).doit()
    metric_paper = sp.Matrix([[-Nf ** 2 + bf ** 2, bf, 0, 0], [bf, 1, 0, 0], [0, 0, Tf ** 2, 0], [0, 0, 0, Tf ** 2 * sp.sin(th) ** 2]])
    area = sp.integrate(sp.integrate(g1[2, 2] * sp.sin(th), (ph, 0, 2 * sp.pi)), (th, 0, sp.pi))   # sqrt(g_thth g_phph) = Theta^2 sin(theta)
    # spatial metric: gamma_ij = g_ij (shift form), gamma_rr = g_rr = 1
    report('S08', 318, 'PASS' if (sp.simplify(g1 - metric_paper) == sp.zeros(4, 4) and g1[1, 1] == 1 and sp.simplify(g1[3, 3] - g1[2, 2] * sp.sin(th) ** 2) == 0
                                  and sp.simplify(area - 4 * sp.pi * Tf ** 2) == 0) else 'FAIL',
           'Eq. (metric) = the machinery metric at L = 1: gamma_rr = g_rr = 1 (r is proper radial distance), sphere area 4 pi Theta^2')
    report('S08', 359, 'PASS' if is_zero((K_c / 3 + Q_c) - D(bf, r) / Nf) else 'FAIL',
           'K^r_r = K/3 + Q = beta\'/N  <=>  d_r beta = N(Q + trK/3)')


def S08b_evolution_radial():
    """Eqs. (evolution) and (radial), L342-365: every printed equation is compared, as a residual in the
    metric functions N, beta, Theta, psi (gamma_rr = 1), with the transcribed code system (R_P, R_K, R_a,
    C_M), which S10 derives from Eq. (action); the psi, Theta and beta' equations are identities of the
    definitions P = u(psi), A = u(Theta), K^r_r = beta'/N.  The worker of S10 also tests the printed
    forms directly against the field equations.  Negative controls: 3Q^2 for 3(1-c13)Q^2 in the trK
    equation, +P psi' for -P psi' in the momentum constraint, a_r - 2Theta'/Theta in the P equation."""
    header('S08b', 342, 'Eqs. (evolution), (radial): the printed GEJ system', r"""
symmetry, with a prime for $\partial_r$, the evolution equations are             (L343)
(\partial_t-\beta\partial_r)\psi&=NP,\nonumber\\
(\partial_t-\beta\partial_r)\Theta&=N\Theta\Bigl(\frac Q2-\frac{\trK}{3}\Bigr),\nonumber\\
(\partial_t-\beta\partial_r)P&=N\Bigl[\trK P+\psi''+\Bigl(a_r+\frac{2\Theta'}{\Theta}\Bigr)\psi'\Bigr],\nonumber\\
(\partial_t-\beta\partial_r)\trK&=\frac N3\,\trK^2\nonumber\\
&\quad+\frac{N}{2+c_{13}+3c_2}\Bigl[2P^2+3(1-c_{13})Q^2\nonumber\\
&\quad+(c_{14}-2)\Bigl(a_r'+a_r^2+\frac{2\Theta'}{\Theta}a_r\Bigr)\Bigr],\label{eq:evolution}\\
(\partial_t-\beta\partial_r)a_r&=N\Bigl(\frac{2\trK}{3}-Q\Bigr)a_r\nonumber\\
&\quad+N\,\frac{c_{13}P\psi'-c_{123}\trK'}{c_{14}(1-c_{13})}.\nonumber              (L345-352)
(1-c_{13})\Bigl(Q'+\frac{3\Theta'}{\Theta}Q\Bigr)
=\frac{2+c_{13}+3c_2}{3}\,\trK'-P\psi',\\
\beta'=N\Bigl(Q+\frac{\trK}{3}\Bigr),                                              (L357-359)
""")
    how('(derived) the printed equations, typed in independently, are written as residuals in the metric functions with\n'
        'P = (psi_t - beta psi\')/N, trK and Q from K^r_r = beta\'/N, K^th_th = -(Theta_t - beta Theta\')/(N Theta), a_r = N\'/N,\n'
        'and subtracted from the code residuals of the machinery (S10: those follow from Eq. action).')
    ok_id = all(is_zero(e) for e in (PR_psi, PR_Th, PR_beta))
    report('S08b', 345, 'PASS' if ok_id else 'FAIL',
           'printed psi, Theta and beta\' equations are identities of P = u(psi), u(Theta) = Theta(Q/2 - trK/3), K^r_r = trK/3 + Q = beta\'/N')
    dif = {k: sp.simplify(sp.expand(a - b)) for k, (a, b) in
           {'P': (PR_P, R_P), 'trK': (PR_K, R_K), 'a_r': (PR_a, R_a), 'Q (radial)': (PR_Q, C_M)}.items()}
    report('S08b', 348, 'PASS' if all(v == 0 for v in dif.values()) else 'FAIL',
           'printed P, trK and a_r equations (eq:evolution) and the momentum constraint (eq:radial) == the code system exactly; S10 shows\n'
           '            the same printed forms to be implied by the field equations of Eq. (action) (residual 0 at random jets)',
           'differences: ' + ', '.join(f'{k}: {v}' for k, v in dif.items()))
    bad = {'trK with 3Q^2': PR_K_negctl - R_K, 'radial with +P psi\'': PR_Q_negctl - C_M, 'P with a_r - 2Theta\'/Theta': PR_P_negctl - R_P}
    nz = {k: not is_zero(v) for k, v in bad.items()}
    report('S08b', 357, 'PASS' if all(nz.values()) else 'FAIL',
           'negative controls: ' + ', '.join(f'{k} differs from the code system ({"yes" if v else "NO"})' for k, v in nz.items()))


def S08c_bonetti_barausse():
    """L339-342: the general 3+1 unitary-gauge equations of Bonetti and Barausse (arXiv:1502.05554, Sec.
    "Modified field equations") with alpha, beta, lambda = c14, c13, c2 and f(a) = alpha a^2.  Their adapted
    action (1-beta)[K_ij K^ij - (1+lambda)/(1-beta) K^2 + 3R/(1-beta) + f(a)/(1-beta)] must equal ours,
    3R + K_ij K^ij - K^2 - (c13 K_ij K^ij + c2 K^2 - c14 a^2); their Hamiltonian constraint (variation
    with respect to N, chi = f'(a)/(2a)) and momentum constraint, reduced to spherical symmetry, must give
    Eq. (hamiltonian) and the Q equation of Eq. (radial).  Their K_ij = +(1/2N)(d_t gamma_ij - ...) is
    minus ours; their matter terms 16 pi G E and 8 pi G J^i are 2 rho and j in our normalisation G_ab = T_ab
    (E = N^2 T^00 = rho, J^i = N(T^0i + N^i T^00) = -T(n, e^i) = j)."""
    header('S08c', 339, 'Bonetti-Barausse 3+1 equations with alpha, beta, lambda = c14, c13, c2', r"""
The 3+1 field
equations in unitary gauge, for general matter and without symmetry, are
given in Ref.~\cite{Bonetti:2015oda} (their $\alpha$, $\beta$, $\lambda$
are our $c_{14}$, $c_{13}$, $c_2$, with $f(a)=\alpha a^2$).                    (L339-342)
""")
    al, be, la = c14, c13, c2
    KK_, Ks, aa_ = sp.symbols('KijKij K a2', real=True)
    L_BB = (1 - be) * (KK_ - (1 + la) / (1 - be) * Ks ** 2) + al * aa_          # + 3R in both
    L_us = KK_ - Ks ** 2 - (c13 * KK_ + c2 * Ks ** 2 - c14 * aa_)
    report('S08c', 341, 'PASS' if is_zero(L_BB - L_us) else 'FAIL',
           'their adapted action (1-beta)[K.K - (1+lambda)/(1-beta) K^2] + f(a), f = alpha a^2, with (alpha, beta, lambda) = (c14, c13, c2)\n'
           '            == R - (c13 K.K + c2 K^2 - c14 a^2) of Eq. (action) (S03), up to the total derivative in R = 3R + K.K - K^2')
    # spherical reduction
    Th = sp.Function('Theta3')(r)
    g3 = sp.diag(-1, 1, Th ** 2, Th ** 2 * sp.sin(th) ** 2)                     # ultrastatic: 4-d R = 3R
    gi3 = sp.simplify(g3.inv())
    Gm3 = christoffel(g3, gi3, XS)
    R3 = sp.simplify(sum(gi3[a, b] * ricci(Gm3, XS, a, b) for a in range(4) for b in range(4)))
    Kq, Qq, ar_, Pq = sp.Function('Kq')(r), sp.Function('Qq')(r), sp.Function('aq')(r), sp.symbols('P')
    psq = sp.Function('psq')(r)
    Krr_t, Ktt_t = -(Kq / 3 + Qq), -(Kq / 3 - Qq / 2)                             # theirs = minus ours
    KK_t = Krr_t ** 2 + 2 * Ktt_t ** 2
    K_t = Krr_t + 2 * Ktt_t
    chi = al                                                                       # f'(a)/(2a) for f = alpha a^2
    div_a = D(ar_, r) + 2 * D(Th, r) / Th * ar_
    H_BB = R3 - (1 - be) * KK_t + (1 + la) * K_t ** 2 + al * ar_ ** 2 - 2 * chi * ar_ ** 2 - 2 * chi * div_a - (Pq ** 2 + D(psq, r) ** 2)
    ThC_pr = (Th * D(Th, r, 2) + (D(Th, r) ** 2 - 1) / 2 + c14 * ar_ * Th * D(Th, r)
              + Th ** 2 / 4 * (c14 * (2 * D(ar_, r) + ar_ ** 2) + Pq ** 2 + D(psq, r) ** 2
                               + sp.Rational(3, 2) * (1 - c13) * Qq ** 2 - sp.Rational(1, 3) * (2 + c13 + 3 * c2) * Kq ** 2))
    okH = is_zero(-Th ** 2 / 4 * H_BB - ThC_pr)
    divK_r = D(Krr_t, r) + 2 * D(Th, r) / Th * (Krr_t - Ktt_t)                   # D_j K^{rj} for diagonal K, gamma_rr = 1
    M_BB = (1 - be) * (divK_r - (1 + la) / (1 - be) * D(K_t, r)) + (-Pq * D(psq, r))    # + 8 pi G J^r, J^r -> j_psi = -P psi'
    PRQ = (1 - c13) * (D(Qq, r) + 3 * D(Th, r) / Th * Qq) - ((2 + c13 + 3 * c2) / 3 * D(Kq, r) - Pq * D(psq, r))
    okM = is_zero(M_BB + PRQ)
    M_bad = (1 - be) * (divK_r - (1 + la) / (1 - be) * D(K_t, r)) + (Pq * D(psq, r))
    report('S08c', 342, 'PASS' if okH and okM and not is_zero(M_bad + PRQ) else 'FAIL',
           'spherical reduction (gamma = dr^2 + Theta^2 dOmega^2, K^r_r = -(trK/3 + Q), K^th_th = -(trK/3 - Q/2) in their sign, chi = alpha,\n'
           '            16 pi G E -> P^2 + psi\'^2, 8 pi G J^r -> -P psi\'): -(Theta^2/4) x their Hamiltonian constraint == Eq. (hamiltonian) and\n'
           '            (1-beta) x their momentum constraint == -(the Q equation of Eq. radial); with the matter sign flipped it fails'
           f' (Hamiltonian {okH}, momentum {okM})')


def S09_nabla_u(acc_lo):
    """L377-378: "a_r is fixed on each slice by the lapse through Eq. (lapse)".  The reason (not printed
    in the paper): D_a u_b = -K_ab - u_a a_b holds covariantly
    (gamma_rr = 1) and contains no d_t N; the ADM K_ij maps onto K^r_r = K/3 + Q,
    K^th_th = K/3 - Q/2 (L337-338)."""
    header('S09', 377, 'a_r fixed on each slice by the lapse: D_a u_b = -K_ab - u_a a_b contains no d_t N', r"""
\emph{Free evolution of $a_r$.} Although $a_r$ is fixed on each slice
by the lapse through Eq.~\eqref{eq:lapse}                                       (L377-378)
""")
    how('compute D_a u_b covariantly (gamma_rr = 1); build K_ab from the ADM formula\n'
        'K_ij = -(1/2N)(d_t g_ij - D_i b_j - D_j b_i) extended by K_ab n^b = 0, and a_b from S07.')
    g1 = gS.subs(Lf, 1).doit()
    Du1_lo = sp.simplify(sp.Matrix(4, 4, lambda a, b: sum(g1[b, c] * Du0[a, c] for c in range(4))).subs(Lf, 1).doit())
    gam = sp.diag(1, Tf ** 2, Tf ** 2 * sp.sin(th) ** 2)
    Y = [r, th, ph]
    gami = gam.inv()
    G3 = [[[sum(gami[a, d] * (D(gam[d, b], Y[c]) + D(gam[d, c], Y[b]) - D(gam[b, c], Y[d])) for d in range(3)) / 2
            for c in range(3)] for b in range(3)] for a in range(3)]
    blo = [bf, 0, 0]
    Db = sp.Matrix(3, 3, lambda i, j: D(blo[j], Y[i]) - sum(G3[l][i][j] * blo[l] for l in range(3)))
    Kij = sp.Matrix(3, 3, lambda i, j: -(D(gam[i, j], t) - Db[i, j] - Db[j, i]) / (2 * Nf))
    bup = [bf, 0, 0]
    K4 = sp.zeros(4, 4)
    for i in range(3):
        for j in range(3):
            K4[i + 1, j + 1] = Kij[i, j]
    for i in range(3):
        K4[0, i + 1] = K4[i + 1, 0] = sum(bup[j] * Kij[j, i] for j in range(3))
    K4[0, 0] = sum(bup[i] * bup[j] * Kij[i, j] for i in range(3) for j in range(3))
    ulo1 = sp.Matrix([-Nf, 0, 0, 0])
    alo1 = [sp.simplify(a.subs(Lf, 1).doit()) for a in acc_lo]
    resid = sp.simplify(Du1_lo + K4 + sp.Matrix(4, 4, lambda a, b: ulo1[a] * alo1[b]))
    okD = resid == sp.zeros(4, 4)
    noNt = all(sp.simplify(sp.diff(Du1_lo[a, b], D(Nf, t))) == 0 for a in range(4) for b in range(4))
    okmap = is_zero(Kij[0, 0] - (K_c / 3 + Q_c)) and is_zero(Kij[1, 1] / Tf ** 2 - (K_c / 3 - Q_c / 2))
    report('S09', 377, 'PASS' if okD and noNt and okmap else 'FAIL',
           'D_a u_b + K_ab + u_a a_b = 0 for all a, b; d(D_a u_b)/d(N_t) = 0; K^r_r = K/3 + Q, K^th_th = K/3 - Q/2 (L337-338)')


def S10_gej_from_action():
    """L377-384 (and Eqs. evolution, hamiltonian): the evolution equations of P, trK, a_r, the
    Hamiltonian constraint (theta_C) and the momentum constraint (compute_Q) of
    release/collapse.cpp (transcription checked by hand, see the docstring) follow from
    Eq. (action) by a reduced-action variation, with general c_i.  S11, S08b and S25 show that the
    displayed Eqs. (hamiltonian) and (evolution) are these equations.  Negative control included.
    Slow; runs in a worker."""
    header('S10', 399, 'The GEJ system evolved in release/collapse.cpp follows from Eq. (action)', r"""
we evolve it with the radial
\ae ther equation combined with the momentum
constraint~\cite{Garfinkle:2007bk}. Because in unitary gauge the
\ae ther equation follows from the metric
equations~\cite{Jacobson:2011cc}, this keeps the Hamiltonian constraint
satisfied as long as $N$ is rebuilt from $a_r$ with
Eq.~\eqref{eq:lapse}.                                                             (L378-384)
(checked are the RHS functions of release/collapse.cpp: compute_RHS, compute_Q (momentum
constraint), theta_C (Hamiltonian constraint), with P = u(psi), a_r = N'/N, beta' = N(Q+K/3),
d_t Theta = beta Theta' + N Theta (Q/2 - K/3); Eqs. (hamiltonian) and (evolution) are compared with
them in S11 and S08b)
""")
    how('(derived) reduced action.  Symmetric criticality holds for SO(3), so varying\n'
        'int N L Theta^2 [R - L_ae - (d psi)^2] w.r.t. N, beta, L, Theta, psi and the aether tilt\n'
        'eta (u = n + eta e; the unit norm holds to O(eta) so lambda drops out) gives the full\n'
        'field equations: E_N (Hamiltonian), E_beta (momentum), E_L (rr), E_Theta (angular),\n'
        'E_eta (radial aether eq.), E_psi.  At L = 1, eta = 0, the six equations are solved for\n'
        'the second-order jets (Theta_rr, beta_rr, Theta_tt, beta_tr, N_tr, psi_tt) at random\n'
        'rational points (all other jets and c1..c4 random), and the code\'s 3 evolution\n'
        'equations (P, trK, a_r) and 2 constraints are evaluated: they must vanish EXACTLY.\n'
        'Negative control: the same test with 3(1-c13)Q^2 -> 3Q^2 in the trK equation must fail.')
    show_slow('_job_gej', 'S10', 399, 'field equations from the action')
    report('S10', 381, 'INFO', 'the Bianchi-identity statement (aether equation implied by the metric equations in unitary gauge) is a\n'
           '            literature input; the derivation above uses all six field equations')


def S11_hamiltonian():
    """Eq. (hamiltonian), L364-375: the printed Theta C is Theta times the Hamiltonian constraint
    of the code (theta_C), which S10 derives from the action (E_N / C_H)."""
    header('S11', 367, 'Eq. (hamiltonian): the Hamiltonian constraint as Theta C', r"""
The Hamiltonian constraint is not imposed; multiplied by
$\Theta$, it reads                                                                 (L365-366)
\Theta C\equiv{}&\Theta\Theta''+\tfrac12(\Theta'^2-1)+c_{14}\,a_r\Theta\Theta'\\
&+\tfrac{\Theta^2}{4}\Bigl[c_{14}(2a_r'+a_r^2)+P^2+\psi'^2\\
&+\tfrac32(1-c_{13})Q^2-\tfrac13(2+c_{13}+3c_2)\trK^2\Bigr]=0 ,                   (L369-371, split)
and we monitor it during the evolution.                                           (L375)
""")
    ThC_paper = (Tf * D(Tf, r, 2) + (D(Tf, r) ** 2 - 1) / 2 + c14 * a_c * Tf * D(Tf, r)
                 + Tf ** 2 / 4 * (c14 * (2 * D(a_c, r) + a_c ** 2) + P_c ** 2 + D(pf, r) ** 2
                                  + sp.Rational(3, 2) * (1 - c13) * Q_c ** 2 - sp.Rational(1, 3) * (2 + c13 + 3 * c2) * K_c ** 2))
    report('S11', 369, 'PASS' if is_zero(ThC_paper - Tf * C_H) else 'FAIL',
           'printed Theta C == Theta x C_H (release/collapse.cpp theta_C), which S10 shows to be the Hamiltonian constraint G_uu = T_uu')


def S12_constraint_propagation():
    """Eq. (constraint-prop), L384-393: regular Taylor series at r = 0 with general c_i and
    general bar-alpha; d/dt of the O(r) coefficient C_1 of the Hamiltonian constraint on C = 0
    equals (1 - bar-alpha) a_r'(0) trK(0) N(0)."""
    header('S12', 384, 'Eq. (constraint-prop): propagation of the Hamiltonian constraint at r = 0', r"""
We check this by relaxing Eq.~\eqref{eq:lapse} to
$\partial_r\ln N=\bar\alpha\,a_r$. Near the centre, the constraint of
Eq.~\eqref{eq:hamiltonian} behaves as $C=C_1\,r+O(r^3)$, and for data
with $C=0$                                                                         (L384-387)
\left.\frac{dC_1}{dt}\right|_{C=0}
=(1-\bar\alpha)\,a_r'(0)\,\trK(0)\,N(0),                                          (L389-390)
""")
    how('(derived) C_1 is the O(r) Taylor coefficient of the Hamiltonian constraint C = (Theta C)/Theta\n'
        '(C is odd in r, so C_1 r is its leading term).  Regular Taylor series at r = 0 (Theta = r +\n'
        'th3 r^3, a_r odd, psi, P, K even), Q from the momentum constraint, N = n0 exp(alpha int a_r),\n'
        'beta from beta\' = N(Q+K/3), d_t of each Taylor coefficient from the code RHS (S10), chain\n'
        'rule, then impose C_1 = 0.  c1..c4 general.')
    lam = sp.symbols('alphabar', real=True)
    th3, a1, a3, k0, k2, p0, p2, ps2, ps4, n0, q2, q4 = sp.symbols('th3 a1 a3 k0 k2 p0 p2 ps2 ps4 n0 q2 q4', real=True)
    ser = lambda e: sp.series(e, rr, 0, 5).removeO().expand()
    Th_ = rr + th3 * rr ** 3
    Ac_ = a1 * rr + a3 * rr ** 3
    K_ = k0 + k2 * rr ** 2
    P_ = p0 + p2 * rr ** 2
    ps_ = ps2 * rr ** 2 + ps4 * rr ** 4
    Qa_ = q2 * rr ** 2 + q4 * rr ** 4
    mom = ser(D(Qa_, rr) * Th_ + 3 * D(Th_, rr) * Qa_ - Th_ * (-P_ * D(ps_, rr) + DEN / 3 * D(K_, rr)) / (1 - c13))
    Q_ = Qa_.subs(sp.solve([mom.coeff(rr, kk) for kk in (2, 4)], [q2, q4], dict=True)[0])
    N_ = ser(n0 * sp.exp(lam * sp.integrate(Ac_, rr)))
    b_ = ser(sp.integrate(ser(N_ * (Q_ + K_ / 3)), rr))
    dTh, dps, dK, dA = D(Th_, rr), D(ps_, rr), D(K_, rr), D(Ac_, rr)
    Th_t = ser(b_ * dTh + N_ * Th_ * (Q_ / 2 - K_ / 3))
    P_t = ser(b_ * D(P_, rr) + N_ * (K_ * P_ + Ac_ * dps + D(ps_, rr, 2) + 2 * sp.cancel(dTh * dps / Th_)))
    K_t = ser(b_ * dK + N_ * K_ ** 2 / 3 + N_ * ((c14 - 2) * (dA + 2 * sp.cancel(dTh * Ac_ / Th_) + Ac_ ** 2)
                                             + 2 * P_ ** 2 + 3 * (1 - c13) * Q_ ** 2) / DEN)
    A_t = ser(b_ * dA + N_ * ((2 * K_ / 3 - Q_) * Ac_ + c13 / (c14 * (1 - c13)) * P_ * dps - c123 / (c14 * (1 - c13)) * dK))
    psi_t = ser(b_ * dps + N_ * P_)
    Cs = ser(sp.cancel(Th_ * (D(Th_, rr, 2) + c14 * Ac_ * dTh + Th_ / 4 * (c14 * (2 * dA + Ac_ ** 2) + P_ ** 2 + dps ** 2
                                                                       + sp.Rational(3, 2) * (1 - c13) * Q_ ** 2 - DEN / 3 * K_ ** 2))
                       + (dTh ** 2 - 1) / 2) / Th_)
    C1 = Cs.coeff(rr, 1)
    dts = {th3: Th_t.coeff(rr, 3), a1: A_t.coeff(rr, 1), k0: K_t.coeff(rr, 0), p0: P_t.coeff(rr, 0), ps2: psi_t.coeff(rr, 2)}
    dC1 = sum(D(C1, v) * e for v, e in dts.items())
    src = sp.factor(sp.expand(dC1.subs(th3, sp.solve(C1, th3)[0])))
    report('S12', 389, 'PASS' if is_zero(src - (1 - lam) * a1 * k0 * n0) else 'FAIL',
           f'dC1/dt|_(C=0) = {src}  ==  (1-alphabar) a_r\'(0) trK(0) N(0)  (general c_i)',
           'with C_1 = the O(r) Taylor coefficient of the Hamiltonian constraint C about r = 0, on data with C = 0 (L385-387);\n'
           '            it vanishes only for alphabar = 1 (L393)')


def S13_singularity_avoidance():
    """L395-412: ln N = -int_r^inf a_r; with a_r >= 0 inside the horizon N decreases inward;
    every non-advective source is proportional to N; the displayed bound
    |beta(r)| <= N(r) int_0^r |Q + trK/3| follows from beta(0) = 0, beta' = N(Q + K/3) and the
    monotonicity of N; the quoted numbers (Theta ~ 8 inside Theta_AH ~ 12 at t ~ 150)."""
    header('S13', 395, 'Singularity-avoiding slicing: ln N, monotonic N, sources ~ N, |beta| bound', r"""
$\ln N(r)=-\int_r^\infty a_r\,dr'$. Once the horizon has formed, $a_r$ is
non-negative inside it, so $N$ decreases monotonically inward there. ... (by $t\simeq150$ this happens
at $\Theta\simeq8$, with $\ThAH\simeq12$) ... Every source term in the
evolution equations is proportional to $N$. The advection terms are
proportional to $\beta$ instead, but Eq.~\eqref{eq:radial}
and $\beta(0)=0$ together with the monotonicity of $N$ give                                   (L397-406, abridged)
|\beta(r)|\le N(r)\int_0^r\Bigl|Q+\frac{\trK}{3}\Bigr|\,dr' ,                     (L408)
""")
    aa = sp.Function('a')
    s_ = sp.symbols('s')
    lnN = -sp.Integral(aa(s_), (s_, rr, sp.oo))
    report('S13', 397, 'PASS' if sp.simplify(D(lnN, rr).doit() - aa(rr)) == 0 else 'FAIL',
           'd_r[-int_r^inf a_r] = a_r and -> 0 as r -> inf  (N -> 1); a_r = (ln N)\' >= 0 => N non-decreasing outward')
    report('S13', 403, 'PASS' if all(sp.simplify(v.subs(Nsym, 0)) == 0 for v in RHS.values()) else 'FAIL',
           'all non-advective RHS terms of the 5 evolved fields (code RHS, S10) vanish at N = 0')
    how('(derived) beta(r) = int_0^r N(r\') g(r\') dr\' with g = Q + K/3 (beta(0) = 0).  If a_r >= 0 on [0, r], N(r\') <= N(r)\n'
        'there, so |beta(r)| <= int_0^r N |g| <= N(r) int_0^r |g|.  Illustrated on random monotone N and signed g; a\n'
        'non-monotone N (a_r < 0 somewhere) can violate the bound, so the monotonicity is needed.')
    rng = np.random.default_rng(3)
    xg = np.linspace(0, 1, 2001)
    ok_b, viol = True, False
    for _ in range(200):
        ar_ = np.abs(rng.normal(size=8)) @ np.array([np.sin((k_ + 1) * xg) ** 2 for k_ in range(8)]) * 3   # a_r >= 0
        lnN_ = np.concatenate([[0], np.cumsum(0.5 * (ar_[1:] + ar_[:-1]) * np.diff(xg))]) - 10
        Nn_ = np.exp(lnN_)
        g_ = rng.normal(size=6) @ np.array([np.cos((k_ + 1) * 3 * xg) for k_ in range(6)])
        beta_ = np.concatenate([[0], np.cumsum(0.5 * ((Nn_ * g_)[1:] + (Nn_ * g_)[:-1]) * np.diff(xg))])
        bound = Nn_ * np.concatenate([[0], np.cumsum(0.5 * (np.abs(g_)[1:] + np.abs(g_)[:-1]) * np.diff(xg))])
        ok_b &= bool(np.all(np.abs(beta_) <= bound * (1 + 1e-9) + 1e-15))
        Nbad = np.exp(-10 + 3 * np.sin(6 * xg))             # non-monotone lapse
        bb_ = np.abs(np.concatenate([[0], np.cumsum(0.5 * ((Nbad * np.abs(g_))[1:] + (Nbad * np.abs(g_))[:-1]) * np.diff(xg))]))
        bnd_ = Nbad * np.concatenate([[0], np.cumsum(0.5 * (np.abs(g_)[1:] + np.abs(g_)[:-1]) * np.diff(xg))])
        viol |= bool(np.any(bb_ > bnd_ * 1.01))
    report('S13', 408, 'PASS' if ok_b and viol else 'FAIL',
           '|beta(r)| <= N(r) int_0^r |Q + trK/3| dr\' follows from beta(0) = 0, beta\' = N(Q+K/3) and a_r >= 0 on [0, r]\n'
           '            (200 random monotone examples obey it; a non-monotone N violates it, so "the monotonicity of N" is the needed\n'
           '            hypothesis, valid inside the horizon where a_r >= 0, L397-398)')
    tl, Th8, Th12 = dval('lapse_collapse')
    report('S13', dline('lapse_collapse'), 'DATA' if Th8 < Th12 and round(dval('ThAH_shrink')[1]) == Th12 and tl == dval('ThAH_shrink')[2] else 'FAIL',
           f'lapse collapse at Theta ~ {Th8} lies inside Theta_AH ~ {Th12} at t ~ {tl} ("at a smaller radius inside it"); '
           f'Theta_AH ~ 12 = the 11.8 of L732 at t = 150 rounded')


def S14_discretization_mesh():
    """L438-456 and Fig. 1 caption (L461-472): dr = r_max/(N_r - 1); the grid falls at ~ beta in
    areal radius because d_t Theta|_r -> beta at large radius; the quintic taper is C^2 and
    supported on the outer 60%.  "Roughly one areal unit per unit time" is the average over the run of
    Fig. 1(a): 284 -> 124 by t = 150 (read off the figure, not printed), 1.07 per unit
    time; the early rate, from the numbers of L666-668 (284 -> 235 by t = 28), is 1.75."""
    header('S14', 438, 'Discretization; grid falls inward; ALE map r = x + g(x) s(t); Fig. 1', r"""
Spatial derivatives are second-order centered finite
differences on a uniform grid of $N_r$ points,
$dr=r_{\max}/(N_r-1)$ ...                                                        (L438-440)
Because, in the $\gamma_{rr}=1$ gauge, the grid falls inward
($\beta$ is nearly uniform across the exterior and of order $-1$), the outer boundary shrinks in areal
radius at roughly one areal unit per unit time
(Fig.~\ref{fig:mesh}(a)). ... map $r=R(x,t)=x+g(x)s(t)$, with $g$ a $C^2$ quintic taper supported for
$x>0.4\,r_{\max}$ (the outer $60\%$ of the grid) and $s(t)$ chosen to
hold $\Theta(r_{\max})$ fixed.                                                   (L446-453)
Fig. 1 caption (abridged): Moving-mesh validation ($r_{\max}=300$, $N_r=25600$, $A_0=0.35$). ... Note that
$\ThAH$ reaches its maximum at $t\simeq28$                                       (L461-471)
""")
    Th_t_far = (ADV['Theta'] + RHS['Theta']).subs({F['Q']: 0, F['K']: 0, F['Thr']: 1})
    report('S14', 447, 'PASS' if Th_t_far == bsym else 'FAIL',
           f'd_t Theta|_r = beta Theta\' + N Theta(Q/2-K/3) -> {Th_t_far} when K, Q -> 0, Theta\' -> 1: the boundary falls at |beta| ~ 1 in areal radius')
    e0, e28, t28 = dval('Theta_edge_mesh'), dval('Theta_edge_fixed28'), dval('collapse_tmax')
    rate, rate22 = (e0 - e28) / t28, (e0 - dval('Theta_ADM_read')) / dval('ADM_t_leave')
    e150, t150 = FIG1_EDGE_T150_UNPRINTED
    rate150 = (e0 - e150) / t150
    report('S14', 448, 'PASS' if 0.9 < rate150 < 1.2 else 'FAIL',
           f'"roughly one areal unit per unit time" (L448) is the run average of Fig. 1(a): ({e0}-{e150})/{t150} = {rate150:.2f} (284 -> 124 by t = 150,\n'
           f'            read off the figure, not printed); the early rate is higher, ({e0}-{e28})/{t28} = {rate:.2f} over t = 0-28\n'
           f'            (L666-668; Theta = 250 leaves the grid at t ~ {dval("ADM_t_leave")}: {rate22:.2f})')
    u_ = sp.symbols('u')
    gq = u_ ** 3 * (10 - 15 * u_ + 6 * u_ ** 2)                       # release/collapse.cpp quintic_taper, collapse.cpp:L325-332
    conds = [gq.subs(u_, 0), D(gq, u_).subs(u_, 0), D(gq, u_, 2).subs(u_, 0),
             gq.subs(u_, 1) - 1, D(gq, u_).subs(u_, 1), D(gq, u_, 2).subs(u_, 1)]
    report('S14', 451, 'PASS' if all(c == 0 for c in conds) and D(gq, u_, 3).subs(u_, 1) != 0 and abs(1 - dval('taper_start') - 0.6) < 1e-12 else 'FAIL',
           'g = 10u^3-15u^4+6u^5 (release/collapse.cpp quintic_taper): g, g\', g\'\' continuous at u = 0, 1 (C^2, not C^3); support x > 0.4 rmax = outer 60%')
    rmx, Nr, A0 = dval('mesh_run')
    report('S14', dline('mesh_run'), 'DATA' if abs(A0 - 0.35) < 1e-12 and dval('t_ThAH_max') == dval('t_mass') else 'FAIL',
           f'Fig. 1 run: r_max = {rmx}, N_r = {Nr} (dr = {rmx / (Nr - 1):.5f}; not one of the nested 2^k 100 + 1 grids, which the moving\n'
           f'            mesh does not need), A0 = {A0} = the baseline shell; Theta_AH maximum at t ~ {dval("t_ThAH_max")} = the t of L679 and L731')


def S15_couplings():
    """L476-479: v0 = 1 for c3 = c4 = 0 iff c2 = -c1^3/(2 - 4c1 + 3c1^2) (used again at L1256);
    c1 = 0.7 gives c2 = -0.5119403; c2 = 0 gives v0 = 1.266862, printed 1.267."""
    header('S15', 476, 'Couplings: c2 = -0.5119403 gives v0 = 1; c2 = 0 gives v0 = 1.267', r"""
$c_1=0.7$, $c_2=-0.5119403$, $c_3=c_4=0$ (so $c_{13}=c_{14}=0.7$), the
values of GEJ~\cite{Garfinkle:2007bk} for which $v_0=1$. ... we also study the case
$c_2=0$ ($v_0=1.267$).                                                           (L477-479, abridged)
""")
    c2sol = sp.solve(sp.Eq(v0sq_sym.subs({c3: 0, c4: 0}), 1), c2)
    report('S15', 478, 'PASS' if len(c2sol) == 1 and is_zero(c2sol[0] + c1 ** 3 / (2 - 4 * c1 + 3 * c1 ** 2)) else 'FAIL',
           f'v0^2 = 1 with c3 = c4 = 0  <=>  c2 = {c2sol[0]}  (= -c1^3/(2-4c1+3c1^2), the formula of L1256)')
    c2_07 = float((-c1 ** 3 / (2 - 4 * c1 + 3 * c1 ** 2)).subs(c1, sp.Rational(7, 10)))
    v0b = math.sqrt(float(v0sq_sym.subs(BASE)))
    v0s = math.sqrt(float(v0sq_sym.subs(SUPER)))
    report('S15', 477, 'PASS' if abs(c2_07 + 0.5119403) < 1e-7 and abs(v0b - 1) < 1e-7 else 'FAIL',
           f'c1 = 0.7: c2 = {c2_07:.8f}; with c2 = -0.5119403, v0 = {v0b:.9f}')
    report('S15', 479, 'PASS' if round(v0s, 3) == dval('v0_super') else 'FAIL',
           f'c2 = 0: v0 = {v0s:.6f} exactly, printed as 1.267 throughout the paper (the exact value is used in every check below)')
    return v0s


def S16_initial_data():
    """Eq. (id), L479-490: the profile is the GEJ profile times r^2/100 and equals the code's
    (analytic psi_r compared); P = K = a_r = 0 makes the momentum constraint homogeneous so Q = 0
    (time symmetry); the shell peaks at r ~ 10 with peak value ~ A0."""
    header('S16', 482, 'Eq. (id): time-symmetric shell', r"""
shell, a variant of that of GEJ~\cite{Garfinkle:2007bk} (whose profile is
$a_0\exp[-(r^2-r_0^2)^2/s^4]$) multiplied by $r^2/100$,                         (L480-481)
\psi(0,r)=\frac{A_0}{100}\,r^2\exp\!\Bigl[-\Bigl(\frac{r^2-r_0^2}{\sigma^2}\Bigr)^{2}\Bigr],
\qquad P(0,r)=0,                                                                  (L483-484)
with $\trK=a_r=0$ as in GEJ~\cite{Garfinkle:2007bk} (hence $Q=0$), and
$r_0=10$ and $\sigma=4$ (the GEJ values), a shell centred at $r\simeq10$
with peak value $\simeq A_0$.                                                    (L487-489)
""")
    how('(i) same profile as release/collapse.cpp set_initial_data / shell_dpsi (analytic derivative\n'
        'compared) and = (r^2/100) x the GEJ profile with a0 = A0, s = sigma; (ii) P = 0 and trK = 0 make\n'
        'the momentum constraint homogeneous, so Q = 0 is the regular solution (time symmetry);\n'
        '(iii) location and value of the maximum of psi.')
    A0, r0, sig, a0s, ss = sp.symbols('A0 r0 sigma a0 s', positive=True)
    psi_id = A0 / 100 * rr ** 2 * sp.exp(-((rr ** 2 - r0 ** 2) / sig ** 2) ** 2)
    gej = a0s * sp.exp(-(rr ** 2 - r0 ** 2) ** 2 / ss ** 4)
    gg_ = (rr ** 2 - r0 ** 2) / sig ** 2
    psir_code = A0 / 100 * sp.exp(-gg_ ** 2) * 2 * rr * (1 - 2 * rr ** 2 * gg_ / sig ** 2)
    report('S16', 483, 'PASS' if is_zero(D(psi_id, rr) - psir_code) and is_zero(psi_id - rr ** 2 / 100 * gej.subs({a0s: A0, ss: sig})) else 'FAIL',
           'profile = (r^2/100) x GEJ profile (a0 = A0, s = sigma); the code\'s analytic psi_r (shell_dpsi) is its derivative')
    Qh = sp.Function('Qh')(rr)
    ode = (1 - c13) * (D(Qh, rr) + 3 * Qh / rr)
    solh = sp.dsolve(sp.Eq(ode, 0), Qh)
    consts = [s_ for s_ in solh.rhs.free_symbols if s_ != rr]
    ok_q = (len(consts) == 1 and sp.simplify(ode.subs(Qh, solh.rhs).doit()) == 0
            and sp.simplify(D(solh.rhs * rr ** 3, rr)) == 0 and sp.limit(solh.rhs.subs(consts[0], 1), rr, 0) == sp.oo)
    report('S16', 487, 'PASS' if ok_q else 'FAIL',
           f'P = K = 0: momentum constraint -> Q = {solh.rhs}, which solves it, is ~ 1/r^3 and diverges at r = 0 unless C1 = 0 => Q = 0 ("hence Q = 0")')
    f1 = psi_id.subs({A0: 1, r0: 10, sig: 4})
    rmax_ = sp.nsolve(D(f1, rr), rr, 10)
    pk = float(f1.subs(rr, rmax_))
    report('S16', 488, 'PASS' if abs(float(rmax_) - 10) < 0.5 and abs(pk - 1) < 0.02 else 'FAIL', f'psi peaks at r = {float(rmax_):.3f} (~10) with psi_max/A0 = {pk:.4f} ("peak value ~ A0")')


def S17_expansions_misner_sharp():
    """Eqs. (expansions), (ms), L492-512: theta_pm = (Theta/2) x the expansions of the null
    normals u +/- e (the ingoing one sign-reversed) = Theta' +/- (3Q-2K) Theta/6; m = Theta(1 -
    theta_+ theta_-)/2 is the Misner-Sharp mass Theta(1 - (grad Theta)^2)/2; theta_+ = 0 gives
    M_hor = Theta_AH/2, so M_hor and the horizon area grow or shrink together."""
    header('S17', 497, 'Eqs. (expansions), (ms): null expansions; Misner-Sharp mass; M_hor = Theta_AH/2', r"""
We use the expansions of
the outgoing and ingoing radial null congruences orthogonal to the
symmetry spheres, with null normals $u^a\pm e^a$ ($e^a$ the unit radial
vector), rescaled by $\Theta/2$ and with the sign of the ingoing one
reversed, so that both are positive on an ordinary sphere:                       (L492-496)
\theta_\pm\equiv\partial_r\Theta\pm\tfrac16(3Q-2\trK)\Theta .                    (L498)
m=\frac{\Theta}{2}\,\bigl(1-\theta_+\theta_-\bigr),                               (L503)
the apparent horizon is the outermost surface where $\theta_+=0$ ... On the apparent horizon $\theta_+=0$,
so Eq.~\eqref{eq:ms} gives $m=\ThAH/2$ ... the horizon mass and the horizon area
$4\pi\ThAH^2$ grow or shrink together.                                            (L506-512, abridged)
""")
    how('the expansion of a congruence k orthogonal to the spheres is q^{ab} D_a k_b = 2 k(Theta)/Theta (q the sphere metric),\n'
        'computed covariantly here; rescaled by Theta/2 it is k(Theta).  Misner-Sharp: 1 - 2m/Theta = g^{ab} d_a Theta d_b Theta.')
    g1 = gS.subs(Lf, 1).doit()
    gi1 = sp.simplify(giS.subs(Lf, 1).doit())
    Gam1 = christoffel(g1, gi1, XS)

    def expansion(kup):
        klo = g1 * kup
        Dk = sp.Matrix(4, 4, lambda a, b: sp.diff(klo[b], XS[a]) - sum(Gam1[c][a][b] * klo[c] for c in range(4)))
        return sp.simplify(sum(gi1[a, b] * Dk[a, b] for a in (2, 3) for b in (2, 3)))
    nu1 = n_up.subs(Lf, 1)
    eu1 = e_up.subs(Lf, 1)
    exp_out, exp_in = expansion(nu1 + eu1), expansion(nu1 - eu1)
    thp = B_c + (3 * Q_c - 2 * K_c) * Tf / 6
    thm = B_c - (3 * Q_c - 2 * K_c) * Tf / 6
    ok_e = is_zero(Tf / 2 * exp_out - thp.subs(Lf, 1)) and is_zero(-Tf / 2 * exp_in - thm.subs(Lf, 1))
    report('S17', 498, 'PASS' if ok_e else 'FAIL',
           'theta_+ = (Theta/2) x expansion of u + e, theta_- = -(Theta/2) x expansion of u - e (covariant, code gauge):\n'
           '            = Theta\' +/- (3Q - 2trK) Theta/6 as printed')
    gradTh2 = sp.simplify(sum(giS[a, b] * D(Tf, XS[a]) * D(Tf, XS[b]) for a in range(2) for b in range(2)).subs(Lf, 1).doit())
    report('S17', 503, 'PASS' if is_zero(gradTh2 - thp * thm) else 'FAIL',
           'theta_+ theta_- = g^{ab} d_a Theta d_b Theta  =>  m = (Theta/2)(1 - theta_+ theta_-) is the Misner-Sharp mass')
    flat = {Tf: r, Nf: 1, bf: 0}
    report('S17', 496, 'PASS' if sp.simplify(thp.subs(flat).doit()) == 1 and sp.simplify(thm.subs(flat).doit()) == 1 else 'FAIL',
           'flat space (Theta = r, N = 1, beta = 0): theta_+ = theta_- = 1 > 0 ("both are positive on an ordinary sphere")')
    ThAH = sp.symbols('Theta_AH', positive=True)
    report('S17', 510, 'PASS' if sp.simplify(ThAH / 2 - sp.sqrt(4 * sp.pi * ThAH ** 2 / (16 * sp.pi))) == 0 else 'FAIL',
           'theta_+ = 0 => m = Theta_AH/2 = sqrt(Area/16 pi): M_hor and the area 4 pi Theta_AH^2 grow or shrink together')
    report('S17', 508, 'DATA', f'M_ADM read at Theta >~ {dval("Theta_ADM_read")}: stated consistently at L508, L662 and L686')


def S18_convergence_order():
    """L525-541: with f_h = f + C h^p the ratio of successive differences is 2^p; dt = 0.1 h is
    halved with dr; for p = 2 the finer difference times 4 equals the coarser one."""
    header('S18', 531, 'Convergence order', r"""
p=\log_2\frac{\|f_{h}-f_{h/2}\|_2}{\|f_{h/2}-f_{h/4}\|_2},                       (L532)
where $f_h$ is $f$ computed with grid spacing $dr=h$ and time step
$dt=0.1\,h$, so that both are halved together ... Equivalently, the difference between the two finer grids, multiplied by
$4$, should coincide with the difference between the two coarser ones.            (L534-540, abridged)
""")
    hh, Cc, qq_, pp = sp.symbols('h C q p', positive=True)       # q = 2^(-p)
    fh = {1: Cc * hh ** pp, 2: Cc * hh ** pp * qq_, 4: Cc * hh ** pp * qq_ ** 2}   # f_h - f = C h^p
    ratio = sp.simplify((fh[1] - fh[2]) / (fh[2] - fh[4]))
    report('S18', 532, 'PASS' if is_zero(ratio - 1 / qq_) else 'FAIL',
           f'f_h = f + C h^p, q = 2^-p  =>  ratio = {ratio} = 2^p  =>  log2(ratio) = p')
    e2 = sp.simplify((4 * (fh[2] - fh[4]) - (fh[1] - fh[2])).subs(qq_, sp.Rational(1, 4)))
    # time step tied to the spacing, dt = dt_factor h (DATA, L535): the leading error C1 dr^2 + C2 dt^2 of the
    # second-order scheme is then (C1 + dt_factor^2 C2) h^2, homogeneous of degree 2 in h, so it obeys the same
    # factor-4 identity; the coefficient 1/100 quoted in the status line is checked against dt_factor.
    C1_, C2_ = sp.symbols('C1 C2', positive=True)
    cdt = sp.nsimplify(dval('dt_factor'))
    err_dt = lambda h: C1_ * h ** 2 + C2_ * (cdt * h) ** 2                   # dt halved with dr
    ok_dt = (is_zero(4 * (err_dt(hh / 2) - err_dt(hh / 4)) - (err_dt(hh) - err_dt(hh / 2)))
             and is_zero(err_dt(hh / 2) - err_dt(hh) / 4)
             and is_zero(err_dt(hh) - (C1_ + C2_ / 100) * hh ** 2))
    report('S18', 539, 'PASS' if e2 == 0 and ok_dt else 'FAIL',
           f'p = 2 (q = 1/4): 4 (f_{{h/2}} - f_{{h/4}}) = f_h - f_{{h/2}}; with dt = {dval("dt_factor")} h, halved with dr, the error '
           f'C1 dr^2 + C2 dt^2 = (C1 + 0.01 C2) h^2\n'
           '            is O(h^2), drops by 4 per halving and obeys the same identity (checked symbolically with dt_factor from the DATA table)')


def S19_regularised_constraint():
    """L543-555: (Theta'^2 - 1)/(2 Theta) is 0/0 at r = 0 and diverges unless Theta'(0) = 1;
    with Theta, a_r odd and the other fields even every term of Theta*C vanishes at r = 0 except
    (Theta'^2 - 1)/2, so Theta*C -> (Theta'(0)^2 - 1)/2."""
    header('S19', 548, 'Regularised constraint Theta*C at the centre', r"""
The constraint $C$ itself contains the term $(\Theta'^2-1)/(2\Theta)$,
which is $0/0$ at a regular centre and cannot be evaluated reliably there
on the grid. We therefore monitor $\Theta C$, Eq.~\eqref{eq:hamiltonian},
which is finite at $r=0$: on a regular slice $\Theta$ and $a_r$ are odd
in $r$ and the other fields even, so near the centre every term of
$\Theta C$ vanishes except $\tfrac12(\Theta'^2-1)$, and
$\Theta C\to\tfrac12[\Theta'(0)^2-1]$.                                            (L548-554)
""")
    th1, th3, a1, a3, p0, p2, s0, s2, k0, k2, q0, q2 = sp.symbols('theta1 th3 a1 a3 p0 p2 s0 s2 k0 k2 q0 q2', real=True)
    Th_reg = th1 * rr + th3 * rr ** 3
    ar_, P_, ps_, K_, Q_ = a1 * rr + a3 * rr ** 3, p0 + p2 * rr ** 2, s0 + s2 * rr ** 2, k0 + k2 * rr ** 2, q0 + q2 * rr ** 2
    term = (D(Th_reg, rr) ** 2 - 1) / (2 * Th_reg)
    lim_rC = sp.limit(rr * term, rr, 0)
    ThC = (Th_reg * D(Th_reg, rr, 2) + (D(Th_reg, rr) ** 2 - 1) / 2 + c14 * ar_ * Th_reg * D(Th_reg, rr)
           + Th_reg ** 2 / 4 * (c14 * (2 * D(ar_, rr) + ar_ ** 2) + P_ ** 2 + D(ps_, rr) ** 2
                                + sp.Rational(3, 2) * (1 - c13) * Q_ ** 2 - DEN / 3 * K_ ** 2))
    others = ThC - (D(Th_reg, rr) ** 2 - 1) / 2
    report('S19', 553, 'PASS' if is_zero(ThC.subs(rr, 0) - (th1 ** 2 - 1) / 2) and others.subs(rr, 0) == 0
           and is_zero(lim_rC - (th1 ** 2 - 1) / (2 * th1)) else 'FAIL',
           f'Theta = theta1 r + ...: (Theta\'^2-1)/(2Theta) ~ {lim_rC}/r (divergent unless theta1 = 1); with the stated parities every other\n'
           '            term of Eq. (hamiltonian) vanishes at r = 0 and Theta*C -> (theta1^2 - 1)/2: the constraint at r = 0 is Theta\'(0) = 1')


def S20_scattering():
    """L557-617, Table 1 (L586-588), Fig. 2 caption (L622-626) and App. A (L1517-1557): dr =
    rmax/(N-1); Table 1 agrees with the text; the t = 14 factors 1.8 / 2.5; the threshold ladder
    (r_max = 300 implied, switch 0.135 -> 0.125, A_* in 0.120-0.125, scattering amplitudes well
    below); the N_r^2 rescaling factors 1, 4, 16 of the constraint figure."""
    header('S20', 557, 'Scattering: grid spacings, Table 1, black-hole threshold; App. A', r"""
three amplitudes, $A_0=0.01$, $0.04$ and $0.06$, on the grids
$N_r=6401/12801/25601$ ($r_{\max}=150$, $dr=0.023/0.012/0.006$). ... which is
$N_{\min}\simeq0.98$, $0.70$ and $0.4$ for the three amplitudes. ... We compute
the norms over $0<r<20$ ... We average the orders over $t=12$--$14$          (L558-568, abridged)
[Table 1, rows L586-588: p(psi, P, trK, a_r, Theta, Q, ||Theta C||) for A_0 = 0.01, 0.04, 0.06 (DATA tab1)]
For $A_0=0.01$ and $0.04$, all fields and the constraint converge at
second order. For $A_0=0.06$, the scalar field and the areal radius
still converge at second order. The \ae ther variables $\trK$, $a_r$ and
$Q$ and the constraint, instead, converge only with $p\approx1.4$--$1.5$. ... At $t=14$, the
constraint norm decreases by a factor $1.8$ between the two coarser grids
and by $2.5$ between the two finer ones, approaching the factor $4$ of
second order.                                                                    (L595-604, abridged)
All three scattering amplitudes lie well below the black-hole threshold
$A_*$. We locate it by bisection in $A_0$ on four grids, from
$N_r=12801$ ($dr=0.023$) to $N_r=102401$ ($dr=0.003$). ... The smallest
amplitude that forms it decreases with resolution, from $0.135$ on the
coarsest grid to $0.125$ on the two finest. Therefore $A_*$ lies between
$0.120$ and $0.125$.                                                             (L606-616, abridged)
Fig. conv-constraint caption (abridged): rescaled by $N_r^2$ and normalised to
the coarsest grid (factors $1$, $4$ and $16$)                                     (L1551-1552)
Fig. conv-fields-A001 caption (abridged):
$A_0=0.01$, on the nested grids $N_r=6401/12801/25601$
($r_{\max}=150$). ... the $L_2$ norms over
$0<r<20$ of $f_h-f_{h/2}$
and $4(f_{h/2}-f_{h/4})$. The curves overlap for
second-order convergence.                                                        (L1530-1535)
Fig. conv-fields-A006 caption (abridged):
for $A_0=0.06$, the
strongest shell of the scattering test. The curves of $\trK$, $a_r$ and
$Q$ separate after the bounce ($t\gtrsim12$).                                    (L1541-1543)
Fig. conv-constraint caption (abridged):
$\|\Theta C\|_2$ over $0<r<20$ against
$t$ for $N_r=6401/12801/25601$ ... (a) Weak field, $A_0=0.01$; (b) strong field, $A_0=0.04$;
(c) $A_0=0.06$, the strongest shell. ... for $A_0=0.06$ they
separate after the bounce ($t\gtrsim12$).                                        (L1550-1556)
""")
    drs = [dval('scatter_rmax') / (n - 1) for n in dval('scatter_N')]
    report('S20', 559, 'DATA' if tuple(round(d, 3) for d in drs) == dval('scatter_dr') else 'FAIL',
           'dr = rmax/(N_r-1) = ' + ', '.join(f'{d:.5f}' for d in drs))
    tab = dval('tab1')
    lo, hi = dval('p_low')
    ok = (min(tab[0.01]) >= 1.9 and min(tab[0.04]) >= 1.9 and all(lo - 0.05 <= tab[0.06][i] <= hi for i in (2, 3, 5, 6))
          and tab[0.06][0] >= 1.95 and tab[0.06][1] >= 1.95 and tab[0.06][4] >= 1.95)
    report('S20', 586, 'DATA' if ok else 'FAIL', 'Table 1 consistent with the text: A0 = 0.01, 0.04 all p >= 1.9; '
           f'A0 = 0.06: trK, a_r, Q, Theta C in {lo - 0.05:.2f}-{hi} ("{lo}-{hi}"), psi, P, Theta >= 1.95 (second order; '
           f'values {tab[0.06][0]}, {tab[0.06][1]}, {tab[0.06][4]})')
    f1, f2, tq = dval('ratio_t14')
    report('S20', 602, 'DATA' if f1 < f2 < 4 and tq in range(dval('scatter_t')[0], dval('scatter_t')[1] + 1) else 'FAIL',
           f'A0 = 0.06 at t = {tq}: factors {f1} < {f2} < 4 (p = {math.log2(f1):.2f}, {math.log2(f2):.2f}): "improves with resolution, approaching 4"; '
           f'N_min {dval("N_min")} in text and Table 1')
    (n1, n4), (d1, d4) = dval('thr_N'), dval('thr_dr')
    rm_impl = (d1 * (n1 - 1), d4 * (n4 - 1))
    sw, br = dval('A_switch'), dval('A_star_bracket')
    okt = (all(abs(x - 300) / 300 < 0.03 for x in rm_impl) and (n4 - 1) == 8 * (n1 - 1) and sw[0] > sw[1] and br[1] == sw[1] and br[0] < br[1]
           and max(dval('scatter_A')) <= 0.5 * br[0])
    report('S20', 608, 'DATA' if okt else 'FAIL',
           f'threshold ladder: N_r = {n1} ... {n4} with dr = {d1} ... {d4} implies r_max = {rm_impl[0]:.0f} / {rm_impl[1]:.0f} (= 300, the collapse box); '
           f'switch {sw[0]} -> {sw[1]} with resolution; bracket {br} ends at the switch of the two finest grids; the strongest '
           f'scattering shell A0 = {max(dval("scatter_A"))} is {100 * max(dval("scatter_A")) / br[0]:.0f}% of the lower end ("well below")')
    fac = [(n / dval('scatter_N')[0]) ** 2 for n in dval('scatter_N')]
    report('S20', dline('conv_rescale'), 'DATA' if tuple(round(f) for f in fac) == dval('conv_rescale') else 'FAIL',
           f'(N_r/6401)^2 = {", ".join(f"{f:.3f}" for f in fac)} -> factors {dval("conv_rescale")}; aether curves separate for t >~ {dval("sep_t")} '
           f'(the bounce, within the averaging window t = 12-14)')
    rs = CONV_RUNS_UNPRINTED['scatter']
    ok_cap = (rs['A0'] == dval('scatter_A') and rs['N'] == dval('scatter_N') and rs['rmax'] == dval('scatter_rmax') and rs['window'] == (0, dval('scatter_window'))
              and rs['field_factor'] == 2 ** 2 and rs['rescale'] == dval('conv_rescale') and rs['ale'] == 0 and dval('sep_t') in range(dval('scatter_t')[0], dval('scatter_t')[1] + 1))
    report('S20', 1530, 'DATA' if ok_cap else 'FAIL',
           f'App. A captions (Figs. conv-fields-A001/A006, conv-constraint): A0 = {rs["A0"]}, N_r = {"/".join(map(str, rs["N"]))}, r_max = {rs["rmax"]}, '
           f'0 < r < {rs["window"][1]}\n            = the figure runs (--A0, --N, --rmax={rs["rmax"]} --ale={rs["ale"]}, fixed grid) and the fig_conv_scatter.py window '
           f'--rlo 0 --rhi 20; the factor 4 of 4(f_h/2 - f_h/4) = 2^p\n            for p = 2 (S18); t >~ {dval("sep_t")} = the bounce, where the A0 = 0.06 orders drop (Table 1)')


def S21_collapse_numbers():
    """L629-673, Table 1 (last row), Table 2 (L1581-1588), App. A (L1560-1602): collapse dr equal
    to the scattering ones; the orders' range; the Table 2 description; Theta_AH ladder factor
    4.3 and Richardson limit 12.70; m on Theta = 250 grows by +0.07% from 3.70836 -> 3.71107 (= the 3.7111 of
    Sec. IV); outer-edge numbers consistent."""
    header('S21', 630, 'Collapse: orders, constraint-norm ratios (Table 2), Theta_AH ladder, ADM drift', r"""
$A_0=0.35$, on the grids $N_r=12801/25601/51201$ with $r_{\max}=300$, so
that $dr=0.023/0.012/0.006$ as in the scattering test, and compare the
three runs every two time units up to $t=28$.                                   (L631-633)
We average the orders over the outputs from $t=4$ to $t=28$ ... They lie between
$2.24$ and $2.62$ ...                                                            (L642-644, abridged)
On the coarser pair of grids they stay close to $4$. On the finer pair
they are close to $4$ up to $t\simeq8$ and then rise to $\simeq6$         (L649-650)
At $t=28$, the differences between its values on successive grids
shrink by a factor $4.3$. Assuming second order, its Richardson limit is
$\ThAH\simeq12.70$.                                                              (L657-659)
We read it at
$\Theta=250$ ... this sphere leaves the
grid at $t\simeq20$. We therefore read the mass on the moving mesh at
the finest resolution ($N_r=51201$), where the outer edge stays at
$\Theta\simeq284$. Between $t=2$ and $t=28$, $m$ on this sphere grows by $0.07\%$. This is
neither a truncation error nor an effect of the outer boundary: at $t=28$
the outer edge of the fixed grid is at $\Theta\simeq235$, and on that
sphere the three fixed grids agree to $10^{-4}$ and the fixed grid and
the moving mesh agree to $2\times10^{-5}$. Nor is it a change of
$\MADM$, the limit of $m$ at large $\Theta$: on the spheres
$\Theta=100$--$250$ the change falls off as $\Theta^{-3}$ (it is $1.2\%$
at $\Theta=100$), and its extrapolation to $\Theta\to\infty$ is
about $0.01\%$.                                                                          (L661-674, abridged)
The apparent horizon
is present from $t=2$ ... On the moving
mesh at $N_r=51201$, the ADM mass read at $\Theta=250$ is                       (L1565-1572, abridged)
Table 2 caption:
Constraint-norm ratios outside the apparent horizon
($r_{\rm AH}+2\le r\le150$) on the nested fixed grid
$N_r=12801/25601/51201$ ($r_{\max}=300$, $A_0=0.35$); a ratio of $4$ is
exact second order.                                                              (L1576-1579)
Fig. conv-collapse caption (abridged):
($A_0=0.35$,
$r_{\max}=300$, fixed grid): $\|\Theta C\|_2$ outside the apparent
horizon, over $r_{\rm AH}+2\le r\le150$, against $t$ on the nested grids
$N_r=12801/25601/51201$, rescaled by $N_r^2$.                                   (L1594-1597)
App. A (L1566-1573, abridged): $\ThAH\simeq8.5$ ... central lapse is below $10^{-6}$ ... $r_{\rm AH}$ grows
from $\simeq9$ at $t=2$ to $\simeq65$ at $t=28$, while the areal radius at $r=150$ decreases from $\simeq137$
to $\simeq90$. ... $\ThAH=12.293$, $12.604$ and $12.677$ ... $3.70836$ at $t=2$ and $3.71107$ at $t=28$.
[Table 2, rows L1585-1586: 3.69 3.80 4.69 4.53 4.26 4.05 3.92 / 3.95 3.97 5.59 6.26 6.12 5.95 5.82]
""")
    drc = [dval('collapse_rmax') / (n - 1) for n in dval('collapse_N')]
    okd = all(abs(a - b) < 1e-12 for a, b in zip(drc, [150 / (n - 1) for n in dval('scatter_N')])) and tuple(round(d, 3) for d in drc) == dval('collapse_dr')
    report('S21', 632, 'DATA' if okd else 'FAIL',
           f'collapse ladder dr = {drc[0]:.5f}/{drc[1]:.5f}/{drc[2]:.5f} = the scattering spacings (nested, halving)')
    pc = dval('p_collapse')
    report('S21', 644, 'DATA' if (min(pc), max(pc)) == dval('p_collapse_range') and min(pc) >= 2 else 'FAIL',
           f'Table 1 last row: min {min(pc)}, max {max(pc)} = "between 2.24 and 2.62" (all >= 2: "second order")')
    rc, rf, tt_ = dval('ratio_row_coarse'), dval('ratio_row_fine'), dval('ratio_t')
    t8, six = dval('ratio_fine_text')
    i8 = tt_.index(t8)
    okr = (all(3.6 <= x <= 4.7 for x in rc) and all(abs(x - 4) < 0.1 for x in rf[:i8 + 1]) and all(abs(x - six) < 0.5 for x in rf[i8 + 1:])
           and tt_[i8] < dval('finest_below_t') < tt_[i8 + 1])
    report('S21', 649, 'DATA' if okr else 'FAIL',
           f'Table 2: coarse pair {min(rc)}-{max(rc)} ("close to 4"); fine pair {rf[:i8 + 1]} up to t = {t8}, then {min(rf[i8 + 1:])}-{max(rf[i8 + 1:])} '
           f'("rise to ~6"); the caption\'s "from t ~ {dval("finest_below_t")}" lies between t = {tt_[i8]} and {tt_[i8 + 1]}')
    r_p = 2 ** pc[-1]
    report('S21', 590, 'DATA' if min(rc + rf) <= r_p <= max(rc + rf) else 'FAIL',
           f'2^p(||Theta C||) = 2^{pc[-1]} = {r_p:.2f} lies within the Table 2 ratios [{min(rc + rf)}, {max(rc + rf)}]')
    a_, b_r, c_ = dval('ThAH_ladder')
    fac = (b_r - a_) / (c_ - b_r)
    rich = c_ + (c_ - b_r) / 3
    report('S21', 658, 'DATA' if round(fac, 1) == dval('ThAH_factor') and round(rich, 2) == dval('ThAH_richardson') else 'FAIL',
           f'Theta_AH ladder (L1571): ({b_r}-{a_})/({c_}-{b_r}) = {fac:.2f} -> {dval("ThAH_factor")} (p = {math.log2(fac):.2f}); '
           f'Richardson with p = 2: {rich:.4f} -> {dval("ThAH_richardson")}')
    m2, m28 = dval('ADM_t2_t28')
    cons = (m28 - m2) / m2 * 100
    report('S21', 666, 'DATA' if round(cons, 2) == dval('ADM_drift_pct') and round(m28, 4) == dval('MADM') else 'FAIL',
           f'm on Theta = 250: ({m28} - {m2})/{m2} = +{cons:.4f}% -> +{dval("ADM_drift_pct")}%; {m28} (t = 28, moving mesh, N_r = {dval("mesh_N_fine")}) '
           f'= the M_ADM = {dval("MADM")} of Sec. IV (t = 28) rounded')
    e0, e28, tl = dval('Theta_edge_mesh'), dval('Theta_edge_fixed28'), dval('ADM_t_leave')
    okE = e0 > dval('Theta_ADM_read') > e28 and tl < dval('collapse_tmax') and dval('ADM_agree')[1] < dval('ADM_agree')[0] < dval('ADM_drift_pct') / 100 * m2
    report('S21', 666, 'DATA' if okE else 'FAIL',
           f'outer edge {e0} (held by the mesh) > read sphere {dval("Theta_ADM_read")} > fixed-grid edge {e28} at t = 28, so Theta = 250 leaves '
           f'before t = 28 (t ~ {tl}); agreements {dval("ADM_agree")} are below the drift {dval("ADM_drift_pct") / 100 * m2:.1e}')
    pf3, t100, t250 = dval('drift_falloff')
    d100, d250, dinf = dval('drift_100_pct'), dval('ADM_drift_pct'), dval('drift_inf_pct')
    p_pr = math.log(d100 / d250) / math.log(t250 / t100)
    A_pr = (d250 * t250 ** 3 - d100 * t100 ** 3) / (t250 ** 3 - t100 ** 3)            # A + B/Theta^3 through the two printed points
    meas = DRIFT_MEASURED_UNPRINTED
    p_me = math.log(meas['pct'][100] / meas['pct'][250]) / math.log(2.5)
    S_ = np.array(sorted(meas['pct'])); y_ = np.array([meas['pct'][k] for k in S_])
    A_fit = np.linalg.lstsq(np.vstack([np.ones_like(S_, dtype=float), S_ ** -3.0]).T, y_, rcond=None)[0][0]
    ok_d = (abs(p_pr - pf3) < 0.2 and abs(A_pr) < dinf and round(abs(meas['A_p3']), 2) == dinf and round(meas['pct'][100], 1) == d100 and round(meas['pct'][250], 2) == d250
            and abs(p_me - pf3) < 0.1 and abs(A_fit - meas['A_p3']) < 2e-4)
    report('S21', 672, 'DATA' if ok_d else 'FAIL',
           f'printed: ln({d100}/{d250})/ln({t250}/{t100}) = p {p_pr:.2f} ~ {pf3} ("falls off as Theta^-3"); A + B/Theta^3 through (100, {d100}%) and\n'
           f'            (250, {d250}%) gives A = {A_pr:+.4f}%.  Measurement (not printed): {meas["pct"][100]}% / {meas["pct"][250]}% at\n'
           f'            Theta = 100 / 250 (ratio {meas["pct"][100] / meas["pct"][250]:.1f} vs 2.5^3 = 15.6, p = {p_me:.2f}); the four-sphere fit A + B/Theta^3\n'
           f'            gives A = {A_fit:+.4f}% (fig_conv_collapse.py record: -0.0102%), i.e. {dinf}% to two digits = "about {dinf}%"; the free-p fit\n'
           f'            (p = 3.12) gives {meas["A_pfree"]:+.4f}% and the two printed points {A_pr:+.4f}%, all of order {dinf}% or below')

    (r9, r65), (t137, t90) = dval('rAH_range'), dval('Theta_r150')
    okA = (r9 < r65 < 150 and t137 > t90 and dval('ThAH_t2') < a_ and dval('AH_t2') == 2 and dval('collapse_window')[1] == 150)
    report('S21', 1568, 'DATA' if okA else 'FAIL',
           f'App. A: window r_AH + 2 <= r <= 150 stays non-empty (r_AH {r9} -> {r65}); Theta(r = 150) {t137} -> {t90}; Theta_AH {dval("ThAH_t2")} at t = 2 '
           f'< {a_} at t = 28; central lapse < {dval("lapse_t28")} (frozen interior)')
    rc_, rm_ = CONV_RUNS_UNPRINTED['collapse'], CONV_RUNS_UNPRINTED['mesh']
    ok_cc = (rc_['N'] == dval('collapse_N') and rc_['rmax'] == dval('collapse_rmax') and rc_['A0'] == 0.35 and rc_['ale'] == 0
             and (rc_['window_lo_plus'], rc_['window_hi']) == dval('collapse_window') and rc_['ratio_times'] == dval('ratio_t')
             and rm_['N'] == dval('collapse_N')[-1] and rm_['ale'] == 1 and (rm_['rmax'], rm_['A0']) == (rc_['rmax'], rc_['A0'])
             and all(round((n / rc_['N'][0]) ** 2) == f for n, f in zip(rc_['N'], (1, 4, 16))))
    report('S21', 1577, 'DATA' if ok_cc else 'FAIL',
           f'App. A, Table 2 and Fig. conv-collapse captions: N_r = {"/".join(map(str, rc_["N"]))}, r_max = {rc_["rmax"]}, A0 = {rc_["A0"]}, fixed grid,\n'
           f'            r_AH + {rc_["window_lo_plus"]} <= r <= {rc_["window_hi"]}, t = {rc_["ratio_times"][0]}...{rc_["ratio_times"][-1]} = the runs col_N* (--ale=0) and the fig_conv_collapse.py window\n'
           f'            [r_AH + 2, --rmax_norm 150] and --ratio_times; "moving mesh at N_r = {rm_["N"]}" = the --mesh_run mass_hi (--N={rm_["N"]} --ale=1);\n'
           f'            rescaling by N_r^2 gives 1, 4, 16 (ratio 4 = exact second order, S18)')
    return cons


def S22_mass_numbers():
    """L679-687: M_hor = Theta_AH/2 = 12.70/2 = 6.35 (Richardson limit), M_ADM = 3.7111, net
    exterior energy -2.64."""
    header('S22', 679, 'Masses of Fig. 4: M_hor, M_ADM, exterior energy', r"""
Figure~\ref{fig:mass} shows the Misner--Sharp mass at $t=28$, when the
apparent horizon is largest. Outside the horizon, $m$ decreases outward
from $\Mhor=\ThAH/2=6.35$ (Richardson limit of $\ThAH$,
Sec.~\ref{sec:conv}) ... to $\MADM=3.7111$ at large radius, ... The exterior therefore
carries a net energy $\MADM-\Mhor=-2.64$.                                      (L679-687, abridged)
""")
    report('S22', 681, 'DATA' if round(dval('ThAH_richardson') / 2, 2) == dval('Mhor') and dval('t_mass') == dval('t_ThAH_max') else 'FAIL',
           f'M_hor = Theta_AH/2 = {dval("ThAH_richardson")}/2 = {dval("ThAH_richardson") / 2:.3f} -> {dval("Mhor")}; t = {dval("t_mass")} = time of the largest horizon (Fig. 1)')
    report('S22', 687, 'DATA' if round(dval('MADM') - dval('Mhor'), 2) == dval('M_deficit') else 'FAIL',
           f'M_ADM - M_hor = {dval("MADM")} - {dval("Mhor")} = {dval("MADM") - dval("Mhor"):.4f} -> {dval("M_deficit")}')


def S23_dm_dTheta():
    """Eq. (dmdTheta), L689-718: d_Theta m|_t = (Theta^2/2)[rho + (A/B) j] derived from the
    Misner-Sharp identity (S39) with G_ab = T_ab; A = u(Theta), B = Theta'; theta_pm = B pm A;
    untrapped <=> B > |A|; the dominant energy condition gives rho >= |j| and then a
    non-negative bracket; the scalar obeys it, so the aether carries at most -2.64."""
    header('S23', 689, 'Eq. (dmdTheta): d_Theta m; theta_pm = B pm A; dominant energy condition', r"""
\partial_\Theta m|_t=\frac{\Theta^2}{2}\Bigl[\rho+\frac{A}{B}\,j\Bigr],\\
A\equiv u^a\partial_a\Theta=\Theta\Bigl(\frac{Q}{2}-\frac{\trK}{3}\Bigr),
\qquad B\equiv\partial_r\Theta ,                                                  (L692-694, gathered)
where $\rho=T_{ab}u^au^b$ and $j=-T_{ab}u^ae^b$ ... $T_{ab}$ is the total (scalar plus \ae ther) stress
tensor. ... The expansions of Eq.~\eqref{eq:expansions} are $\theta_\pm=B\pm A$, so
a sphere is untrapped ($\theta_\pm>0$) when $B>|A|$. The dominant energy
condition requires the energy-flux vector
$-T^a{}_bu^b=\rho\,u^a+j\,e^a$ to be causal and future-directed, i.e.\
$\rho\ge|j|$. The bracket in Eq.~\eqref{eq:dmdTheta} is then
non-negative on untrapped spheres: outside the apparent horizon, $m$ can
decrease outward only where the dominant energy condition is
violated~\cite{Hayward:1994bu}.                                                  (L698-709, abridged)
... the scalar obeys the condition, so the
violation comes from the \ae ther, whose contribution to the exterior
energy is at most $-2.64$.                                                       (L713-715)
""")
    how('(derived) contract the exact identity of S39 with e^a (fixed t):\n'
        '   e(m) = (Theta^2/2)[T(e, grad Theta) - Theta\' T^(2)],  T^(2) = -rho + p,  grad Theta = -A u + B e\n'
        'with T_ab = rho u_a u_b + j(u_a e_b + e_a u_b) + p e_a e_b  (rho = T(u,u), j = -T(u,e), p = T(e,e)),\n'
        'in the normalisation G_ab = T_ab stated at L315 (code and Eq. lum).')
    T_ue = -jj
    T2 = -rho + pr_
    em = Th_s ** 2 / 2 * ((-AA) * T_ue + BB * pr_ - BB * T2)     # T(e, -A u + B e) - B T2
    dmdTheta = sp.simplify(em / BB)
    print(f'    derived:  d_Theta m|_t = {dmdTheta}')
    report('S23', 692, 'PASS' if is_zero(dmdTheta - Th_s ** 2 / 2 * (rho + AA / BB * jj)) and is_zero(A_c - Tf * (Q_c / 2 - K_c / 3))
           and is_zero(A_c - sum(n_up[a] * D(Tf, XS[a]) for a in range(4)).subs(Lf, 1)) else 'FAIL',
           'd_Theta m|_t = (Theta^2/2)[rho + (A/B) j] (G_ab = T_ab) from the identity of S39; A = u^a d_a Theta = Theta(Q/2 - trK/3)')
    thp = B_c + (3 * Q_c - 2 * K_c) * Tf / 6
    thm = B_c - (3 * Q_c - 2 * K_c) * Tf / 6
    report('S23', 702, 'PASS' if is_zero(thp - (B_c + A_c)) and is_zero(thm - (B_c - A_c)) else 'FAIL',
           'theta_pm = B pm A; theta_+ > 0 and theta_- > 0  <=>  B > |A| (untrapped)')
    # T^a_b u^b with u.u = -1, u.e = 0, e.e = 1: build in an orthonormal frame
    eta2 = sp.diag(-1, 1)                                                  # frame metric on (u, e)
    Tfr = sp.Matrix([[rho, -jj], [-jj, pr_]])                              # T(u,u) = rho, T(u,e) = -j, T(e,e) = p
    uvec = sp.Matrix([1, 0])
    flux = -(eta2.inv() * Tfr * uvec)                                       # -T^a_b u^b, components along (u, e)
    norm2 = (flux.T * eta2 * flux)[0]
    ok_f = sp.simplify(flux - sp.Matrix([rho, jj])) == sp.zeros(2, 1) and sp.simplify(norm2 - (jj ** 2 - rho ** 2)) == 0
    report('S23', 705, 'PASS' if ok_f else 'FAIL',
           '-T^a_b u^b = rho u^a + j e^a; its norm is j^2 - rho^2 and its time component rho: causal and future-directed <=> rho >= |j|')
    xx2 = sp.symbols('x', real=True)
    rv, jv = sp.symbols('rho_v j_v', real=True)
    convex = sp.expand((1 + xx2) / 2 * (rv + jv) + (1 - xx2) / 2 * (rv - jv) - (rv + xx2 * jv))
    report('S23', 706, 'PASS' if convex == 0 else 'FAIL',
           'rho + x j = ((1+x)/2)(rho + j) + ((1-x)/2)(rho - j) (identity): for |x| = |A|/B < 1 (untrapped) both weights are > 0 and,\n'
           '            with rho >= |j|, both brackets >= 0: outside the apparent horizon (B > 0) m can decrease outward only where the DEC fails')
    f_x = rho_psi + xx2 * j_psi
    ok_min = all(is_zero(f_x.subs(xx2, s) - (Pp - s * psp) ** 2 / 2) for s in (1, -1))
    report('S23', 713, 'PASS' if ok_min else 'FAIL',
           'scalar: rho_psi -/+ j_psi = (P -/+ psi\')^2/2 >= 0, i.e. rho_psi >= |j_psi| (DEC); since Eq. (dmdTheta) is linear in T_ab, the\n'
           '            scalar part of the exterior energy is >= 0 and the aether part is <= the total')
    report('S23', 715, 'DATA' if dval('M_deficit') < 0 else 'FAIL',
           f'E_ae = E_total - E_psi <= E_total = {dval("M_deficit")}: "at most -2.64"; aether density positive only across the shell '
           f'{dval("ae_positive_band")} at t = {dval("t_mass")}')


def S24_mass_ratio():
    """Eq. (mratio), L711-735: the Penrose inequality from integrating Eq. (dmdTheta) under the DEC;
    6.35/3.7111 = 1.711 -> 1.71; the area-theorem numbers; a negative rho + (A/B) j requires a
    violation of the dominant energy condition, not necessarily of the WEC (explicit
    counterexample), so "dominant" (L711-712, L1474) is the right condition."""
    header('S24', 711, 'Eq. (mratio), Penrose inequality, area theorem, energy condition', r"""
The outward decrease of $m$ in Fig.~\ref{fig:mass} therefore requires a
violation of the dominant energy condition outside the horizon.                 (L711-712)
In GR with matter obeying the dominant energy condition, integrating
Eq.~\eqref{eq:dmdTheta} from the horizon to infinity gives the
spherical Penrose inequality~\cite{Penrose:1973um,Malec:1994sy,Hayward:1994bu},
$\Mhor\le\MADM$. Here instead                                                     (L720-723)
\frac{\Mhor}{\MADM}=1.71 ,                                                        (L725)
... whereas in our run $\ThAH$ peaks at $t\simeq28$
and then decreases slowly, from $12.6$ to $11.8$ by $t=150$ on the grid
of Fig.~\ref{fig:mesh}(b). The area theorem relies only on the null
energy condition, which the \ae ther stress tensor does not satisfy in
general~\cite{Akhoury:2016mrw}.                                                   (L731-735)
""")
    Th_ah, Th_p, M_ = sp.symbols('Theta_AH Theta_p M', positive=True)
    gfun = sp.Function('g')                                  # d_Theta m = g >= 0 (DEC on untrapped spheres, S23)
    m_inf = Th_ah / 2 + sp.Integral(gfun(Th_p), (Th_p, Th_ah, sp.oo))
    ok_pen = sp.simplify(m_inf - Th_ah / 2 - sp.Integral(gfun(Th_p), (Th_p, Th_ah, sp.oo))) == 0
    report('S24', 720, 'PASS' if ok_pen else 'FAIL',
           'M_ADM = m(inf) = m(Theta_AH) + int_{Theta_AH}^inf d_Theta m dTheta = Theta_AH/2 + (>= 0 when the DEC holds on the untrapped\n'
           '            exterior, S23) >= M_hor: integrating Eq. (dmdTheta) gives the spherical Penrose inequality (needs the exterior untrapped)')
    report('S24', 725, 'DATA' if round(dval('Mhor') / dval('MADM'), 2) == dval('Mratio') and dval('Mratio') > 1 else 'FAIL',
           f'{dval("Mhor")}/{dval("MADM")} = {dval("Mhor") / dval("MADM"):.4f} -> {dval("Mratio")} > 1 (violates M_hor <= M_ADM)')
    a12, a118, t150 = dval('ThAH_shrink')
    report('S24', 732, 'DATA' if a118 < a12 and round(dval('ThAH_ladder')[1], 1) == a12 else 'FAIL',
           f'area theorem violated: Theta_AH {a12} -> {a118} by t = {t150} ({100 * (1 - a118 / a12):.0f}% decrease; 12.6 = the N_r ~ 25600\n'
           f'            value 12.604 of L1571 rounded, the resolution of Fig. 1(b))')
    rv, jv, pv, xv = 1.0, -5.0, 100.0, 0.5
    wec_ok = all(rv - 2 * tv * jv + pv * tv * tv >= 0 for tv in np.linspace(-1, 1, 2001))
    dec_val = rv + xv * jv
    report('S24', 711, 'PASS' if wec_ok and dec_val < 0 else 'FAIL',
           'd_Theta m < 0 on an untrapped sphere means T(u, v_Theta) < 0 for two future-directed timelike vectors: it "requires a violation\n'
           '            of the dominant energy condition" (L711-712, L1474); "weak" would not follow:',
           f'counterexample: rho = 1, j = -5, p_r = 100, p_tan = 0 obeys the WEC for every observer, yet rho + (A/B) j = {dec_val} < 0 at A/B = 0.5')
    report('S24', 734, 'INFO', 'the area theorem for apparent horizons rests on the null energy condition (literature, L729-735)')


def S25_evolution_block():
    """L755-762: Eqs. (evolution) are first order with d_r psi as an extra variable; their trK and a_r
    lines, written as d_t f = beta f' + N[...], are the code's (S08b, S10); the psi, P and Theta
    equations contain no derivatives of trK or a_r, so the principal symbol is block-triangular."""
    header('S25', 755, 'Eq. (evolution) as a first-order system; block-triangular symbol', r"""
With $\partial_r\psi$ as an extra variable, Eqs.~\eqref{eq:evolution}
are first order in time and space, and the eigenvalues and eigenvectors
of their principal part give the characteristic speeds and variables.
The equations for $\psi$, $P$ and $\Theta$ contain no derivatives of
$\trK$ or $a_r$, so the principal symbol is block-triangular and the
spin-0 speeds are those of the $(\trK,a_r)$ block                              (L755-761)
""")
    Kp = (bsym * F['Kr'] + Nsym / 3 * F['K'] ** 2 + Nsym / (2 + c13 + 3 * c2) * (2 * F['P'] ** 2 + 3 * (1 - c13) * F['Q'] ** 2
                                                                            + (c14 - 2) * (F['dr_a_r'] + 2 * F['Thr'] / F['Th'] * F['a_r'] + F['a_r'] ** 2)))
    ap = (bsym * F['dr_a_r'] + Nsym * (sp.Rational(2, 3) * F['K'] - F['Q']) * F['a_r']
          + Nsym / (c14 * (1 - c13)) * (c13 * F['P'] * F['psr'] - c123 * F['Kr']))
    report('S25', 755, 'PASS' if is_zero(Kp - rhsK) and is_zero(ap - rhsA) else 'FAIL',
           'the trK and a_r lines of Eq. (evolution) (L348-352), written as d_t f = beta f\' + N[...], == the code RHS (compute_RHS,\n'
           '            collapse.cpp:L875-886) used for the principal part below (S08b, S10: these follow from Eq. action)')
    free = all(sp.diff(ADV[k] + RHS[k], v) == 0 for k in ('Theta', 'psi', 'P') for v in (F['Kr'], F['dr_a_r']))
    report('S25', 758, 'PASS' if free else 'FAIL',
           'the psi, P and Theta equations contain no d_r trK, d_r a_r (and d_t psi\' = d_r(beta psi\' + N P) only through beta\' = N(Q+K/3),\n'
           '            N\' = N a_r): the (trK, a_r) block of the principal symbol has zeros below it, i.e. the symbol is block-triangular')

def S26_principal_part():
    """Eq. (principal), L762-766: the Jacobian of the code RHS of (trK, a_r) with respect to
    (trK', a_r') equals M with a_N, b_N of Eq. (aNbN)."""
    header('S26', 762, 'Eq. (principal): principal part of the (trK, a_r) pair', r"""
\partial_t\mathbf u = \mathsf M\,\partial_r\mathbf u+\dots,\qquad
\mathsf M=\begin{pmatrix}\beta & N a_N\\ N b_N & \beta\end{pmatrix},             (L763-764)
""")
    Mprin = sp.Matrix([[D(rhsK, F['Kr']), D(rhsK, F['dr_a_r'])], [D(rhsA, F['Kr']), D(rhsA, F['dr_a_r'])]])
    Mpaper = sp.Matrix([[bsym, Nsym * aN_sym], [Nsym * bN_sym, bsym]])
    report('S26', 764, 'PASS' if sp.simplify(Mprin - Mpaper) == sp.zeros(2, 2) else 'FAIL',
           'Jacobian d(RHS)/d(trK\', a_r\') == M with a_N, b_N of Eq. (aNbN)',
           'the trK and a_r equations also contain (2 a_r/Theta) Theta\' and P psi\' (first derivatives of Theta, psi with\n'
           '            field-dependent coefficients): the off-diagonal blocks of the block-triangular symbol (S25, S29)')


def S27_aNbN():
    """Eq. (aNbN), L768-775: a_N b_N = v0^2 identically; the code's spin0_aN, spin0_bN
    (release/collapse.cpp main, collapse.cpp:L1615-1616) are the same expressions."""
    header('S27', 768, 'Eq. (aNbN): a_N, b_N and a_N b_N = v0^2', r"""
a_N=\frac{c_{14}-2}{2+c_{13}+3c_2},\qquad
b_N=-\frac{c_{123}}{c_{14}(1-c_{13})},\\
a_Nb_N=v_0^2 .                                                                    (L770-772, gathered)
""")
    report('S27', 772, 'PASS' if is_zero(aN_sym * bN_sym - v0sq_sym) else 'FAIL', 'a_N b_N - v0^2(Eq. speeds) = 0 identically')
    # collapse.cpp:L1607-1609: c13 = c1+c3, c14 = c1+c4, c123 = c1+c2+c3, den_K = 2+c13+3c2, den_a = c14(1-c13)
    code_aN = (c1 + c4 - 2) / (2 + c1 + c3 + 3 * c2)            # collapse.cpp:L1615: spin0_aN = (c14 - 2)/den_K
    code_bN = -(c1 + c2 + c3) / ((c1 + c4) * (1 - c1 - c3))     # collapse.cpp:L1616: spin0_bN = -c123/den_a
    report('S27', 770, 'PASS' if is_zero(code_aN - aN_sym) and is_zero(code_bN - bN_sym) else 'FAIL',
           f'release/collapse.cpp spin0_aN, spin0_bN (collapse.cpp:L1615-1616) == a_N, b_N of Eq. (aNbN); '
           f'baseline a_N = {float(aN_sym.subs(BASE)):.6f}, b_N = {float(bN_sym.subs(BASE)):.6f}')


def S28_eigenvalues():
    """L776-788: the diagonal is advection with the aether (aether worldlines have rdot = -beta);
    eigenvalues beta +/- N v0; d_t w = lambda d_r w moves at rdot = -lambda; the speed relative to
    the aether, (rdot + beta)/N, is -/+ v0 whatever N and beta.  Franchini, Herrero-Valea and Barausse
    (arXiv:2103.00929, App. A): their c0^2 = (lambda+beta)(2-alpha)/(alpha(1-beta)(2+3lambda+beta)) is v0^2 with
    (alpha, beta, lambda) = (c14, c13, c2), and rdot = +/- v0 N/sqrt(gamma_rr) at zero shift is a speed v0 relative to
    the aether in proper distance and proper time."""
    header('S28', 776, 'Advection, eigenvalues and characteristic speeds', r"""
The diagonal entries describe advection with the \ae ther observers,
$\dot r=-\beta$; ... The eigenvalues
$\lambda_\pm=\beta\pm Nv_0$ give characteristics $\dot r=-\beta\mp Nv_0$,
i.e.\ speeds $\pm v_0$ relative to the \ae ther in its proper time and
proper radial distance ($\gamma_{rr}=1$), whatever the lapse and shift.          (L776-782)
Franchini, Herrero-Valea and Barausse~\cite{Franchini:2021bpt} found the
same speeds for the infrared limit of Ho\v{r}ava gravity in spherical
symmetry, with the same slicing but zero shift and $\gamma_{rr}\neq1$,
using derivatives of the metric functions as first-order variables:
their spin-0 characteristics, $\dot r=\pm v_0N/\sqrt{\gamma_{rr}}$,
again correspond to a speed $v_0$ relative to the
\ae ther~\cite{Jacobson:2010mx}.                                               (L783-789)
""")
    rdot_ae = sp.simplify(n_up[1] / n_up[0])
    report('S28', 777, 'PASS' if sp.simplify(rdot_ae + bf) == 0 else 'FAIL',
           f'aether worldlines u^a = (1/N, -beta/N, 0, 0): dr/dt = {rdot_ae} (the diagonal beta d_r is advection along them)')
    ev = list(Mab.eigenvals())
    ok = all(any(is_zero(e - (bsym + s * Nsym * v0_ab)) for e in ev) for s in (1, -1))
    rd = [-(bsym + s * Nsym * v0_ab) for s in (1, -1)]
    phys = [sp.simplify((x_ + bsym) / Nsym) for x_ in rd]
    report('S28', 780, 'PASS' if ok and all(is_zero(a - b) for a, b in zip(phys, [-v0_ab, v0_ab])) else 'FAIL',
           'eigenvalues beta +/- N sqrt(a_N b_N) = beta +/- N v0; w(r + lambda t) moves at rdot = -lambda; (rdot + beta)/N = -/+ v0')
    alF, beF, laF = c14, c13, c2
    c0sq = (laF + beF) * (2 - alF) / (alF * (1 - beF) * (2 + 3 * laF + beF))
    grr, Nz = sp.symbols('gamma_rr N_z', positive=True)
    vv_ = sp.Symbol('v0', positive=True)
    rdotF = vv_ * Nz / sp.sqrt(grr)
    speedF = sp.simplify(sp.sqrt(grr) * rdotF / Nz)                               # proper distance / proper time, zero shift
    report('S28', 787, 'PASS' if is_zero(c0sq - v0sq_sym) and speedF == vv_ else 'FAIL',
           'Franchini et al. (arXiv:2103.00929, App. A, Eq. c0): c0^2 = (lambda+beta)(2-alpha)/(alpha(1-beta)(2+3lambda+beta)) with\n'
           '            (alpha, beta, lambda) = (c14, c13, c2) equals v0^2 of Eq. (speeds); their rdot = +/- v0 N/sqrt(gamma_rr) at zero shift is\n'
           '            sqrt(gamma_rr) rdot/N = +/- v0 in proper distance per proper time')


def S29_characteristic_variables():
    """Eq. (wout), L791-804: the rows (-/+ b_N/v0, 1) are left
    eigenvectors of M for beta -/+ N v0; padded with zeros, the only obstruction to being exact left
    eigenvectors of the full principal symbol of the printed trK and a_r equations is the coupling to
    Theta (coefficient ~ a_r/Theta, from 2 Theta' a_r/Theta in the trK equation) and to psi (coefficient
    ~ P, from c13 P psi' in the a_r equation); with the wave-zone fall-offs a_r ~ P ~ 1/r (neither
    vanishes) both coefficients fall off; the scalar pair P -/+ d_r psi moves along rdot = -beta +/- N;
    the numerical coefficient -b_N/v0 is not printed (Fig. 3 caption)."""
    header('S29', 791, 'Eq. (wout): characteristic variables; exactness far out; scalar pair', r"""
\wout=-\frac{b_N}{v_0}\,\trK+a_r,\qquad
\win=+\frac{b_N}{v_0}\,\trK+a_r ,                                                 (L793-794)
which obey $\partial_t w-\lambda\,\partial_r w=(\text{lower order})$,
with $\lambda=\lambda_-$ for $\wout$, outgoing relative to the \ae ther,
and $\lambda=\lambda_+$ for $\win$. They become exact characteristic
variables of the full system far from the source, where the terms
coupling them to $\Theta$ and $\psi$, proportional to $a_r/\Theta$ and
$P$, fall off.                                                                   (L797-802)
... The outgoing and ingoing scalar combinations are
$P\mp\partial_r\psi$, which move along $\dot r=-\beta\pm N$.                   (L803-804)
""")
    lout = sp.Matrix([[-bN_ / v0_ab, 1]])
    lin = sp.Matrix([[bN_ / v0_ab, 1]])
    ok_out = sp.simplify(lout * Mab - (bsym - Nsym * v0_ab) * lout) == sp.zeros(1, 2)
    ok_in = sp.simplify(lin * Mab - (bsym + Nsym * v0_ab) * lin) == sp.zeros(1, 2)
    report('S29', 793, 'PASS' if ok_out and ok_in else 'FAIL',
           'w_out = l_out.u with l_out M = (beta - N v0) l_out (lambda_-, rdot = -beta + N v0, outgoing); w_in <-> lambda_+')
    how('full principal symbol in U = (trK, a_r, Theta, phi = d_r psi, P): rows d(RHS)/d(U\') of the code equations, with\n'
        'd_t phi = d_r(beta phi + N P) (principal part beta phi\' + N P\').  Off-diagonal entries: trK row, Theta\' column\n'
        '2N(c14-2) a_r/((2+c13+3c2) Theta); P row, Theta\' column 2 N phi/Theta.  With psi (not phi) as variable the a_r row\n'
        'has a psi\' entry N c13 P/(c14(1-c13)).')
    a_s, Th_s2, ph_s, P_s = sp.symbols('a Theta phi P')
    vv = sp.symbols('v0', positive=True)
    S = sp.zeros(5, 5)
    S[0, 0], S[0, 1], S[0, 2] = bsym, Nsym * aN_, 2 * Nsym * (c14 - 2) * a_s / (DEN * Th_s2)
    S[1, 0], S[1, 1] = Nsym * bN_, bsym
    S[2, 2] = bsym
    S[3, 3], S[3, 4] = bsym, Nsym
    S[4, 2], S[4, 3], S[4, 4] = 2 * Nsym * ph_s / Th_s2, Nsym, bsym
    sub = {aN_: vv ** 2 / bN_}
    Lout5 = sp.Matrix([[-bN_ / vv, 1, 0, 0, 0]])
    Lin5 = sp.Matrix([[bN_ / vv, 1, 0, 0, 0]])
    rout = sp.simplify((Lout5 * S - (bsym - Nsym * vv) * Lout5).subs(sub))
    rin = sp.simplify((Lin5 * S - (bsym + Nsym * vv) * Lin5).subs(sub))
    ok_ex = all(sp.simplify(x_.subs(a_s, 0)) == 0 for x_ in list(rout) + list(rin)) and any(x_ != 0 for x_ in rout)
    # the printed trK and a_r equations (eq:evolution), typed in again, in the variables (trK, a_r, Theta, psi, P):
    # rows d(RHS)/d(U') of the two lines; the padded rows of w_out, w_in only see these two rows.
    Kx, ax, Thx, Px, Qx, dK, da, dTh, dps, dPx = sp.symbols('K_ a_ Theta_ P_ Q_ dK da dTheta dpsi dP')
    rhsK_printed = (bsym * dK + Nsym / 3 * Kx ** 2
            + Nsym / DEN * (2 * Px ** 2 + 3 * (1 - c13) * Qx ** 2 + (c14 - 2) * (da + ax ** 2 + 2 * dTh / Thx * ax)))
    rhsA_printed = bsym * da + Nsym * (2 * Kx / 3 - Qx) * ax + Nsym * (c13 * Px * dps - c123 * dK) / (c14 * (1 - c13))
    dU = (dK, da, dTh, dps, dPx)

    def pad_residual(rk, ra, sign):
        rows = sp.Matrix([[sp.diff(rk, d_) for d_ in dU], [sp.diff(ra, d_) for d_ in dU]])
        l2 = sp.Matrix([[sign * bN_sym / sp.sqrt(v0sq_sym), 1]])
        lam = bsym + sign * Nsym * sp.sqrt(v0sq_sym)
        return sp.simplify(l2 * rows - lam * sp.Matrix([[l2[0], l2[1], 0, 0, 0]]))
    res_o, res_i = pad_residual(rhsK_printed, rhsA_printed, -1), pad_residual(rhsK_printed, rhsA_printed, +1)
    cth, cps = sp.simplify(res_o[2] * Thx / (Nsym * ax)), sp.simplify(res_o[3] / (Nsym * Px))
    ok_struct = all(sp.simplify(r_[k_]) == 0 for r_ in (res_o, res_i) for k_ in (0, 1, 4)) and \
        cth.free_symbols <= {c1, c2, c3, c4} and cps.free_symbols <= {c1, c2, c3, c4} and \
        sp.simplify(cth) != 0 and sp.simplify(cps) != 0 and \
        sp.simplify(res_i[2] * Thx / (Nsym * ax)).free_symbols <= {c1, c2, c3, c4} and \
        sp.simplify(res_i[3] / (Nsym * Px)).free_symbols <= {c1, c2, c3, c4}
    # each coupling comes from the term named in the paper's note: drop it and the column vanishes
    rk0 = rhsK_printed.subs(dTh, 0)
    ra0 = rhsA_printed.subs(dps, 0)
    ok_orig = sp.simplify(pad_residual(rk0, rhsA_printed, -1)[2]) == 0 and sp.simplify(pad_residual(rhsK_printed, ra0, -1)[3]) == 0
    # wave zone: Theta = r, a_r = alpha/r, P = pi/r (neither vanishes); couplings ~ r^-2 and r^-1 -> 0, principal entries O(1)
    al_, pi_, rw = sp.symbols('alpha pi_ r', positive=True)
    wz = {Thx: rw, ax: al_ / rw, Px: pi_ / rw, Nsym: 1}
    co_th, co_ps = sp.simplify(res_o[2].subs(wz)), sp.simplify(res_o[3].subs(wz))
    ok_wz = (sp.limit(co_th, rw, sp.oo) == 0 and sp.limit(co_ps, rw, sp.oo) == 0
             and sp.simplify(co_th * rw ** 2).free_symbols.isdisjoint({rw}) and sp.simplify(co_ps * rw).free_symbols.isdisjoint({rw})
             and sp.limit(co_th * rw ** 2, rw, sp.oo) != 0 and sp.limit(co_ps * rw, rw, sp.oo) != 0)
    ok_psi0 = sp.simplify(res_o[2].subs({Px: 0, dps: 0})) != 0      # psi = 0 (P = 0) leaves the Theta coupling
    report('S29', 799, 'PASS' if ok_ex and ok_struct and ok_orig and ok_wz and ok_psi0 else 'FAIL',
           f'(w_out, w_in) padded with zeros vs the full principal symbol of the printed trK, a_r equations: the residual has only a\n'
           f'            Theta\' column N a_r/Theta x {float(cth.subs(BASE)):.4f} (from 2 Theta\' a_r/Theta in the trK equation) and a psi\' column\n'
           f'            N P x {float(cps.subs(BASE)):.4f} (from c13 P psi\' in the a_r equation; c13 = c14 = 0.7, c2 = -0.5119403), i.e. "the terms coupling\n'
           f'            them to Theta and psi, proportional to a_r/Theta and P"; with a_r ~ P ~ 1/r (neither zero in the wave zone) the\n'
           f'            couplings fall off as r^-2 and r^-1: "exact far from the source"; psi = P = 0 alone would leave the Theta coupling')
    Msc = sp.Matrix([[bsym, Nsym], [Nsym, bsym]])               # (phi, P) block
    lm = sp.Matrix([[-1, 1]])                                   # P - phi
    lp = sp.Matrix([[1, 1]])                                    # P + phi
    ok_sc = (sp.simplify(lm * Msc - (bsym - Nsym) * lm) == sp.zeros(1, 2) and sp.simplify(lp * Msc - (bsym + Nsym) * lp) == sp.zeros(1, 2))
    report('S29', 803, 'PASS' if ok_sc else 'FAIL',
           'scalar block [[beta, N], [N, beta]] on (d_r psi, P): P - d_r psi <-> beta - N (rdot = -beta + N, outgoing), P + d_r psi <-> beta + N\n'
           '            (rdot = -beta - N, ingoing): "P -/+ d_r psi move along rdot = -beta +/- N"')


def S30_linear_system():
    """L809-823 (inline): the (c14 - 2) bracket of the trK equation is the divergence
    of Theta^2 N a_r; the O(eps) part of the code RHS on N = 1, beta = 0, Theta = r is
    d_t trK = a_N r^-2 d_r(r^2 a_r), d_t a_r = b_N d_r trK; the scalar enters only quadratically."""
    header('S30', 809, 'Linearised (trK, a_r) system far from the source (inline, L821-822)', r"""
We solve them in
the static frame of flat space, $N=1$ and zero shift, with the areal
radius $\Theta=r$ as radial coordinate. Our gauge differs from this frame
only by a rigid translation of the radial labels on each slice,
$r\to r+c(t)$, which keeps $\gamma_{rr}=1$ and gives the shift
$\beta=\dot c$. Since $\beta\simeq-1$ far out, $c(t)\simeq-t$: a point at
fixed $r$ falls inward in areal radius at about unit speed
(Sec.~\ref{sec:method}), and a time derivative at fixed areal radius
becomes $\partial_t-\beta\partial_r$ in our coordinates.                         (L810-818)
Keeping only the
terms linear in the fields, the $\trK$ and $a_r$ equations of
Eq.~\eqref{eq:evolution} become
$\partial_t\trK=a_N\,r^{-2}\partial_r(r^2a_r)$ and
$\partial_t a_r=b_N\,\partial_r\trK$; the scalar drops out, since it
enters them only quadratically.                                                  (L818-823)
""")
    Nfr = sp.Function('N')(rr)
    Thr = sp.Function('Theta')(rr)
    ar_exact = D(Nfr, rr) / Nfr
    lhs_div = D(Thr ** 2 * D(Nfr, rr), rr) / (Nfr * Thr ** 2)
    cfun = sp.Function('c')(t)
    Rr_ = r + cfun
    dR = sp.Matrix([D(Rr_, t), 1])                                       # dR = c' dt + dr
    g_flat_2d = sp.Matrix([[-1 + dR[0] ** 2, dR[0]], [dR[0], 1]])          # -dt^2 + dR^2 in (t, r)
    g_gauge = sp.Matrix([[-1 + D(cfun, t) ** 2, D(cfun, t)], [D(cfun, t), 1]])   # -N^2 dt^2 + (dr + beta dt)^2, N = 1, beta = c'
    Hf, s2_ = sp.Function('H'), sp.Symbol('s_2', real=True)
    w_ = Hf(Rr_, s2_)                     # field w(R, s2), s2 = t: d_t|_r w = D(w, t) + D(w, s2), d_t|_R w = D(w, s2)
    ok_fixR = sp.simplify((D(w_, t) + D(w_, s2_)) - D(cfun, t) * D(w_, r) - D(w_, s2_)) == 0
    c_unit = sp.dsolve(sp.Eq(D(cfun, t), -1), cfun).rhs                      # beta = c' = -1
    rfix_R = sp.simplify(D(r + c_unit, t))                                    # areal radius R = r + c of a point at fixed r
    report('S30', 814, 'PASS' if sp.simplify(g_flat_2d - g_gauge) == sp.zeros(2, 2) and g_gauge[1, 1] == 1 and ok_fixR and rfix_R == -1 else 'FAIL',
           'flat space -dt^2 + dR^2 + R^2 dOmega^2 with R = r + c(t): -dt^2 + (dr + c\' dt)^2 + (r + c)^2 dOmega^2, i.e. N = 1, gamma_rr = 1, beta = c\'\n'
           f'            ("a rigid translation of the radial labels ... keeps gamma_rr = 1 and gives the shift beta = c\'", L812-815); beta = -1 gives\n'
           f'            c = {c_unit}, so a point at fixed r has dR/dt = {rfix_R} ("falls inward in areal radius at about unit speed", L815-816); a field\n'
           '            w(t, R) has d_t|_R = d_t|_r - c\' d_r = d_t - beta d_r ("a time derivative at fixed areal radius becomes d_t - beta d_r", L817-818)')
    report('S30', 821, 'PASS' if is_zero(lhs_div - (D(ar_exact, rr) + 2 * D(Thr, rr) * ar_exact / Thr + ar_exact ** 2)) else 'FAIL',
           '(c14-2) bracket = (1/(N Theta^2)) d_r(Theta^2 N a_r) exactly (with a_r = N\'/N)')
    eps_ = sp.symbols('epsilon')
    Kf, af = sp.Function('K')(t, rr), sp.Function('a')(t, rr)
    sub_lin = {F['K']: eps_ * Kf, F['Kr']: eps_ * D(Kf, rr), F['a_r']: eps_ * af, F['dr_a_r']: eps_ * D(af, rr), F['Q']: 0, F['P']: 0,
               F['psr']: 0, F['Th']: rr, F['Thr']: 1, Nsym: 1, bsym: 0}
    linK = sp.expand(D(rhsK.subs(sub_lin), eps_).subs(eps_, 0))
    linA = sp.expand(D(rhsA.subs(sub_lin), eps_).subs(eps_, 0))
    okK = is_zero(linK - aN_sym * D(rr ** 2 * af, rr) / rr ** 2)
    okA = is_zero(linA - bN_sym * D(Kf, rr))
    report('S30', 821, 'PASS' if okK and okA else 'FAIL',
           'O(eps) part of the code RHS on N = 1, beta = 0, Theta = r == the linearised system of L821-822',
           'Q is O(eps) but enters trK, a_r only through Q^2 and Q a_r (second order); the scalar enters only quadratically (P psi\', P^2): it drops out')


def S31_wave_equation():
    """L823-824: eliminating a_r gives the spherical wave equation with speed^2 = a_N b_N = v0^2."""
    header('S31', 823, 'Spherical wave equation for trK', r"""
Together they give the wave equation
$\partial_t^2\trK=v_0^2\,r^{-2}\partial_r(r^2\partial_r\trK)$                  (L823-824)
""")
    Kf = sp.Function('K')(t, rr)
    lhs = aN_sym / rr ** 2 * D(rr ** 2 * (bN_sym * D(Kf, rr)), rr)
    report('S31', 824, 'PASS' if is_zero(lhs - v0sq_sym * D(rr ** 2 * D(Kf, rr), rr) / rr ** 2) else 'FAIL',
           'd_t of the first eq. + second eq.: d_t^2 K = a_N b_N r^-2 d_r(r^2 d_r K) = v0^2 Laplacian K')


def S32_asymptotic_solution():
    """Eq. (asymptotic), L824-836: the printed (trK, a_r) with arbitrary F solves both equations of
    the linearised system (L821-822) exactly; the static Coulomb field a_r = C/r^2 is also a solution (L1045-1047)."""
    header('S32', 826, 'Eq. (asymptotic): outgoing solution', r"""
\trK=\frac{f(u)}{r},\qquad
a_r=-\frac{b_N}{v_0}\,\frac{f(u)}{r}-b_N\,\frac{F(u)}{r^2},\\
u=t-\frac{r}{v_0},\quad F'=f,                                                     (L828-830, gathered)
""")
    Ff = sp.Function('F')
    uu_ = t - rr / v0_ab
    Ksol = D(Ff(uu_), t) / rr                      # f = F'
    asol = -bN_ / v0_ab * D(Ff(uu_), t) / rr - bN_ * Ff(uu_) / rr ** 2
    e1 = sp.simplify(D(Ksol, t) - aN_ * D(rr ** 2 * asol, rr) / rr ** 2)
    e2 = sp.simplify(D(asol, t) - bN_ * D(Ksol, rr))
    report('S32', 828, 'PASS' if e1 == 0 and e2 == 0 else 'FAIL', 'substituting (f = F\', arbitrary F) solves BOTH linearised equations (L821-822) exactly')
    Cst = sp.symbols('C')
    report('S32', 1045, 'PASS' if D(rr ** 2 * Cst / rr ** 2, rr) == 0 else 'FAIL',
           'static Coulomb solution trK = 0, a_r = C/r^2 also solves the linearised system: behind the pulse Theta^2 a_r relaxes to the same value on\n'
           '            every sphere, "the F(u)/r^2 term of Eq. (asymptotic)" (L1045-1047)')
    return Ksol, asol


def S33_wout_asymptotic(Ksol, asol):
    """Eq. (wout-asym), L834-847: the 1/r^2 part of a_r builds up as the pulse goes by and remains
    behind it (F = int f); w_out carries the O(1/r) part, r w_out -> -2 b_N f/v0; w_in has no 1/r part
    and changes while the pulse goes by, d_t w_in = -b_N trK/r."""
    header('S33', 834, 'Eq. (wout-asym): the 1/r^2 part; w_out, w_in; waveform; w_in changes while the pulse goes by', r"""
where the profile $f$ is set by the source and $F$ is its integral.
Besides the $1/r$ falloff of any spherical wave, $a_r$ has a $1/r^2$
part, which builds up as the pulse goes by and remains behind it.               (L834-836)
\wout=-\frac{2b_N}{v_0}\,\frac{f(u)}{r}-b_N\,\frac{F(u)}{r^2},
\qquad
\win=-b_N\,\frac{F(u)}{r^2}.                                                      (L839-841)
The radiation is carried by $\wout$, and $r\,\wout\to-2b_Nf(u)/v_0$ is
the spin-0 waveform we extract. The ingoing variable $\win$ has no $1/r$
part, but it is not zero: it changes while the pulse goes by.                   (L844-846)
""")
    Ff = sp.Function('F')
    uu_ = t - rr / v0_ab
    wout = -bN_ / v0_ab * Ksol + asol
    win = bN_ / v0_ab * Ksol + asol
    okw = is_zero(wout - (-2 * bN_ / v0_ab * D(Ff(uu_), t) / rr - bN_ * Ff(uu_) / rr ** 2)) and is_zero(win + bN_ * Ff(uu_) / rr ** 2)
    report('S33', 839, 'PASS' if okw else 'FAIL', 'w_out carries the O(1/r) part, w_in = O(1/r^2) (no 1/r part): r w_out -> -2 b_N f(u)/v0 is the waveform')
    rate = sp.simplify(D(win, t) - (-bN_ * Ksol / rr))
    bump = lambda x: mpmath.e ** (-1 / (1 - x ** 2)) if abs(x) < 1 else 0
    Fv = [float(mpmath.quad(bump, [-1, min(u_, 1)])) for u_ in (-0.5, 0.0, 0.5, 1.5, 2.0, 5.0)]
    grows = Fv[0] < Fv[1] < Fv[2] < Fv[3] and max(Fv[3:]) - min(Fv[3:]) < 1e-12
    report('S33', 845, 'PASS' if rate == 0 and grows else 'FAIL',
           'd_t w_in = -b_N F\'(u)/r^2 = -b_N trK/r (F\' = f): w_in, like the 1/r^2 part of a_r, changes while the pulse goes by; for a\n'
           f'            positive pulse on |u| < 1, F = int_-inf^u f grows ({Fv[0]:.3f} -> {Fv[2]:.3f}) during it and stays at {Fv[3]:.4f} after it\n'
           '            ("builds up as the pulse goes by and remains behind it"; a static C/r^2 a_r solves the linearised system, S32)')


def S34_retarded_time():
    """Eq. (retarded), L911-921 (in the Extraction paragraph of Sec. V.C; kept here in the order of the
    derivations): u = t - Theta_*/v0 is the flat-space u of Eq. (asymptotic) with r -> Theta_*, the
    Schwarzschild tortoise coordinate, dTheta_*/dTheta = 1/(1 - 2M/Theta); along u = const the wave moves at
    dTheta/dt = v0 (1 - 2M/Theta) < v0 (the "slowing of the wave near the black hole")."""
    header('S34', 911, 'Eq. (retarded): tortoise coordinate', r"""
To compare spheres we use the retarded time                                       (L911)
u=t-\frac{\Theta_*}{v_0},\qquad
\Theta_*=\Theta+2\MADM\ln\!\left(\frac{\Theta}{2\MADM}-1\right),                 (L913-914)
the flat-space $u$ of Eq.~\eqref{eq:asymptotic} with $r$ replaced by
the Schwarzschild tortoise coordinate $\Theta_*$, which accounts
approximately for the slowing of the wave near the black hole.                  (L917-919)
""")
    Mm = sp.symbols('M', positive=True)
    Tstar = Th_s + 2 * Mm * sp.log(Th_s / (2 * Mm) - 1)
    vv = sp.symbols('v0', positive=True)
    speed = sp.simplify(vv / D(Tstar, Th_s))                       # dTheta/dt at u = const
    report('S34', 914, 'PASS' if sp.simplify(D(Tstar, Th_s) - 1 / (1 - 2 * Mm / Th_s)) == 0 and is_zero(speed - vv * (1 - 2 * Mm / Th_s)) else 'FAIL',
           f'dTheta_*/dTheta = 1/(1 - 2M/Theta) (Schwarzschild tortoise); at u = t - Theta_*/v0 = const, dTheta/dt = {speed} < v0: the\n'
           '            slowing near the black hole; Theta_* -> Theta + O(M ln Theta) far out (the flat u of Eq. asymptotic with r -> Theta_*);\n'
           '            the paper flags it as an approximation (L917-919, L920-921, L963-964)')


def S35_boundary_conditions():
    """L849-884, Eq. (bc-win): parities at r = 0 are preserved by the evolution; the entering fields at r_max
    are w_in and P + d_r psi; w_in -> 0 at infinity for outgoing radiation, but is O(1/r^2) at r_max; the
    outgoing-wave condition, Eq. (bc-win), is d_t w_in = -b_N trK/r with d_t -> N^-1(d_t - beta d_r) (the
    derivative along the aether) and r -> Theta, exact for a linear outgoing wave (constant N, beta), and is
    what the code imposes (compute_RHS, outer boundary, collapse.cpp:L892-929;
    the projector algebra of the code update is kept as INFO); P + psi' = -g/r^2 for psi = g(u)/r; for the
    linear spin-0 system d_t w_in = -b_N trK/r <=> (d_t + v0 d_r)(r trK) = 0; the scalar condition (Theta psi
    constant along outgoing rays, imposed through d_t P, equivalent for outgoing waves: INFO code check); the outflow
    gate (INFO code check; the projection and the outflow condition are not stated in the paper)."""
    header('S35', 849, 'Boundary conditions: parities; entering fields; Eq. (bc-win); scalar; code checks (scalar update, outflow)', r"""
At $r=0$ we impose regularity: $\psi$, $P$,
$\trK$ and $Q$ are even in $r$, and $\Theta$, $a_r$ and $\beta$ are odd.
At the outer boundary a condition is needed only for the fields that
enter the grid there, $\win$ for the spin-0 pair and $P+\partial_r\psi$
for the scalar; the outgoing fields simply leave. Physically, no
radiation should come in from infinity. At finite radius, however, a
purely outgoing wave has a nonzero entering field, one power of $r$
below the radiation: $\win=-b_NF(u)/r^2$ by Eq.~\eqref{eq:wout-asym}, and
$P+\partial_r\psi=-g(u)/r^2$ for a scalar wave $\psi=g(u)/r$. Setting
these fields to zero at $r_{\max}$, their value at infinity, therefore
reflects part of the wave, by the ratio of the two falloffs, i.e.\ about
the pulse width over $r_{\max}$. Instead we impose non-reflecting conditions,
which require the fields at $r_{\max}$ to behave as for an outgoing wave: $r\trK=f(u)$ and $r\psi=g(u)$ must be constant along
the outgoing characteristics, of speed $v_0$ and $1$. For the spin-0
pair this is equivalent to $\partial_t\win=-b_N\trK/r$, the form we
implement~\cite{Sarbach:2012pr}; on the collapse background, with the
time derivative at fixed areal radius divided by the lapse and $\Theta$
for $r$, it reads                                                                 (L849-866)
(\partial_t-\beta\partial_r)\win=-Nb_N\,\frac{\trK}{\Theta}.                    (L868)
For the scalar we impose that $\Theta\psi$ be constant along the outgoing
light rays; for a linear wave on flat space this is equivalent to
$\partial_t(P+\partial_r\psi)=-P/r$, the analogue of
Eq.~\eqref{eq:bc-win}, and the entering field again changes at a rate set
by its $1/r^2$ part.                                                              (L871-875)
""")
    # (a) parity preservation
    cs = sp.symbols('e0:12')
    ev_ = lambda o: cs[o] + cs[o + 1] * rr ** 2 + cs[o + 2] * rr ** 4
    od_ = lambda o: cs[o] * rr + cs[o + 1] * rr ** 3
    Th_, a__, b__ = rr + od_(0) * rr ** 2, od_(2), od_(4)
    ps_, P__, K__, Q__, N__ = ev_(6), ev_(9), ev_(6) ** 2 + 1, rr ** 2 * ev_(9), 1 + rr ** 2 * ev_(6)
    vals = {F['Th']: Th_, F['Thr']: D(Th_, rr), F['psr']: D(ps_, rr), F['psrr']: D(ps_, rr, 2), F['P']: P__, F['K']: K__,
            F['Kr']: D(K__, rr), F['a_r']: a__, F['dr_a_r']: D(a__, rr), F['Q']: Q__, Nsym: N__, bsym: b__}
    Prr = sp.symbols('Prr')
    full = {k: (ADV[k] + RHS[k]) for k in RHS}
    full['P'] = full['P'] + bsym * Prr
    vals[Prr] = D(P__, rr)
    par = {'psi': 1, 'P': 1, 'trK': 1, 'a_r': -1, 'Theta': -1}
    okp = all(sp.simplify(full[k].subs(vals) - par[k] * full[k].subs(vals).subs(rr, -rr)) == 0 for k in par)
    report('S35', 849, 'PASS' if okp else 'FAIL',
           'with psi, P, trK, Q, N even and Theta, a_r, beta odd, the RHS of psi, P, trK are even and those of a_r, Theta odd: the\n'
           '            evolution preserves the stated parities (beta odd follows from beta(0) = 0, beta\' = N(Q+K/3) even)')
    # (b), (c) projector algebra for one entering mode (lambda_+)
    vv = sp.symbols('v0', positive=True)
    Mv = Mab.subs(aN_, vv ** 2 / bN_)
    aNv = vv ** 2 / bN_
    ep = sp.Matrix([aNv / vv, 1]); em_ = sp.Matrix([-aNv / vv, 1])
    tp = sp.Matrix([[bN_ / vv, 1]]); tm = sp.Matrix([[-bN_ / vv, 1]])
    eig_ok = (sp.simplify(Mv * ep - (bsym + Nsym * vv) * ep) == sp.zeros(2, 1) and sp.simplify(tp * Mv - (bsym + Nsym * vv) * tp) == sp.zeros(1, 2)
              and sp.simplify(Mv * em_ - (bsym - Nsym * vv) * em_) == sp.zeros(2, 1))
    nrm = sp.simplify((tp * ep)[0])
    eI, thI = ep / nrm, tp                     # normalised pair, theta^I.u = w_in
    thO, eO = tm / sp.simplify((tm * em_)[0]), em_
    dual = sp.simplify(sp.Matrix([[(thI * eI)[0], (thI * eO)[0]], [(thO * eI)[0], (thO * eO)[0]]])) == sp.eye(2)
    Pp_ = sp.eye(2) - eI * thI
    R1, R2, Bs, G_ = sp.symbols('R1 R2 B G')
    Rv = sp.Matrix([R1, R2])
    new = Pp_ * Rv + eI * Bs
    okP = (sp.simplify(Pp_ * Pp_ - Pp_) == sp.zeros(2, 2) and sp.simplify(Pp_ * eI) == sp.zeros(2, 1) and sp.simplify(thI * Pp_) == sp.zeros(1, 2)
           and sp.simplify(thO * Pp_ - thO) == sp.zeros(1, 2) and sp.simplify((thI * new)[0] - Bs) == 0
           and sp.simplify((thO * new)[0] - (thO * Rv)[0]) == 0)
    report('S35', 852, 'INFO' if eig_ok and dual and okP and nrm == 2 else 'FAIL',
           f'(code check) e_+ = (a_N/v0, 1), theta^+ = (b_N/v0, 1) belong to lambda_+ = beta + N v0 (theta^+.e_+ = a_N b_N/v0^2 + 1 = {nrm}); with e_I = e_+/2,\n'
           '            theta^J.e_I = delta^J_I; P = I - e_I theta^I is a projector (P^2 = P), P e_I = 0, theta^I P = 0, theta^out P = theta^out, so\n'
           '            d_t(theta^I.u) = B^I and the outgoing equation is unchanged (projector algebra of the code update; the projector is no\n'
           '            longer printed); w_in = theta^+.u moves along rdot = -lambda_+ (S28) and enters at r_max iff lambda_+ > 0: "a condition is\n'
           '            needed only for the fields that enter the grid there, w_in for the spin-0 pair" (L851-852)')
    code = Rv - ep * ((tp * Rv)[0] - G_) / (tp * ep)[0]
    report('S35', 868, 'PASS' if sp.simplify(code - new.subs(Bs, G_)) == sp.zeros(2, 1) else 'FAIL',
           'release/collapse.cpp (compute_RHS, spin-0 pair, collapse.cpp:L901-915): d_t u = R - e_+ (theta^+.R - G)/(theta^+.e_+) = P R + e_I G: it sets d_t w_in = G and leaves the\n'
           '            outgoing w_out equation unchanged (spin-0 pair), i.e. it imposes the stated condition on the entering field only')
    # (d) spin-0 boundary data: exact for the outgoing linear wave in flat space with constant N, beta
    Nc, bc_ = sp.symbols('N_c beta_c', positive=True)
    Fp = sp.Function('F')
    Tt, Rr = Nc * t, r + bc_ * t                       # T = N t, R = r + beta t: ds^2 = -dT^2 + dR^2
    uret = Tt - Rr / vv
    s_arg = sp.Symbol('s')
    fu = D(Fp(s_arg), s_arg).subs(s_arg, uret)          # f(u) = F'(u)
    Ks = fu / Rr
    as_ = -bN_ / vv * fu / Rr - bN_ * Fp(uret) / Rr ** 2
    # linearised code equations on this background (Theta = R, Theta' = 1, N, beta constants)
    eqK = sp.simplify(D(Ks, t) - (bc_ * D(Ks, r) + Nc * aNv * (D(as_, r) + 2 * as_ / Rr)))
    eqA = sp.simplify(D(as_, t) - (bc_ * D(as_, r) + Nc * bN_ * D(Ks, r)))
    w_in = bN_ / vv * Ks + as_
    rate = sp.simplify(D(w_in, t) - bc_ * D(w_in, r) + Nc * bN_ * Ks / Rr)
    ff_ = sp.Function('f'); FF_ = sp.Function('F')
    win_far = -bN_ * FF_(t - rr / vv) / rr ** 2
    ratio_fall = sp.simplify((bN_ * FF_(t - rr / vv) / rr ** 2) / (2 * bN_ / vv * ff_(t - rr / vv) / rr))
    along = sp.simplify(sum(n_up[a] * D(sp.Function('w')(t, r), XS[a]) for a in range(2)) - (D(sp.Function('w')(t, r), t) - bf * D(sp.Function('w')(t, r), r)) / Nf)
    Bq, Aq = sp.symbols('B_q A_q', positive=True)
    dtTh = bsym * Bq + Nsym * Aq                                          # d_t Theta|_r = beta Theta' + N A (code Theta equation)
    fixTh_shift = sp.simplify(dtTh / Bq)                                   # d_t|_Theta = d_t|_r - (d_t Theta/Theta') d_r
    report('S35', 865, 'INFO',
           f'note: "a time derivative at fixed areal radius becomes d_t - beta d_r" (L817-818) is stated for flat space (Theta = r + c(t), exact,\n'
           f'            S30); on the collapse background ("with the time derivative at fixed areal radius divided by the lapse", L864-865) d_t|_Theta = d_t - ({fixTh_shift}) d_r, A = u(Theta), B = Theta\', which reduces to d_t - beta d_r\n'
           '            only where A = 0 and B = 1; Eq. (bc-win) uses N^-1(d_t - beta d_r) = u^a d_a, the derivative along the aether (below)')
    gsc = sp.Function('g')
    psi_w = gsc(t - rr) / rr                                                # N = 1, beta = 0: d_t psi = N P + beta psi' -> P = d_t psi
    P_w = sp.solve(sp.Eq(D(psi_w, t), (ADV['psi'] + RHS['psi']).subs({Nsym: 1, bsym: 0}).subs(F['P'], sp.Symbol('Pw'))), sp.Symbol('Pw'))[0]
    sc_in = sp.simplify(P_w + D(psi_w, rr))
    sc_out = sp.simplify(P_w - D(psi_w, rr))
    report('S35', 857, 'PASS' if sp.limit(win_far.subs(FF_, sp.atan), rr, sp.oo) == 0 and is_zero(ratio_fall * rr - vv * FF_(t - rr / vv) / (2 * ff_(t - rr / vv)))
           and is_zero(sc_in + gsc(t - rr) / rr ** 2) and is_zero(sc_out - 2 * D(psi_w * rr, t) / rr - gsc(t - rr) / rr ** 2)
           and sp.limit(sc_in.subs(gsc, sp.atan), rr, sp.oo) == 0 else 'FAIL',
           f'w_in = -b_N F(u)/r^2 -> 0 at large radius for bounded F (tested with F = arctan) but is nonzero at finite r, one power of r below\n'
           f'            w_out: |w_in|/|w_out| = v0 F/(2 f r) ~ (pulse width)/r ("the ratio of the two falloffs", L855-859; S36: 3.48% x 100/150 = 2.32%\n'
           f'            vs the derived 2.33% at r_max = 150); scalar psi = g(t - r)/r (P = d_t psi on N = 1, beta = 0): P + d_r psi = {sc_in}\n'
           f'            ("P + d_r psi = -g(u)/r^2", L857), one power below P - d_r psi = 2 g\'(u)/r + g(u)/r^2; both entering fields -> 0 at infinity\n'
           '            ("Setting these fields to zero at r_max, their value at infinity", L857-858)')
    # spin-0: d_t w_in = -b_N trK/r  <=>  (d_t + v0 d_r)(r trK) = 0 for the linear system (a_N b_N = v0^2)
    Kl, al_ = sp.Function('K')(t, rr), sp.Function('a')(t, rr)
    dtK = aNv * D(rr ** 2 * al_, rr) / rr ** 2                             # d_t K = a_N r^-2 d_r(r^2 a_r), a_N = v0^2/b_N
    dta = bN_ * D(Kl, rr)                                                   # d_t a_r = b_N d_r K
    win_rate = bN_ / vv * dtK + dta + bN_ * Kl / rr                         # d_t w_in + b_N K/r
    char_rK = rr * dtK + vv * D(rr * Kl, rr)                                # (d_t + v0 d_r)(r K)
    equiv = sp.simplify(win_rate - bN_ / (vv * rr) * char_rK)
    report('S35', 863, 'PASS' if equiv == 0 and along == 0 else 'FAIL',
           'linear system d_t trK = a_N r^-2 d_r(r^2 a_r), d_t a_r = b_N d_r trK with a_N b_N = v0^2: d_t w_in + b_N trK/r =\n'
           '            (b_N/(v0 r)) (d_t + v0 d_r)(r trK) identically, so (b_N != 0) "r trK = f(u) constant along the outgoing characteristics, of speed\n'
           '            v0" is equivalent to d_t w_in = -b_N trK/r ("For the spin-0 pair this is equivalent to", L862-863); N^-1(d_t - beta d_r) = u^a d_a,\n'
           '            the derivative along the aether (L864-865)')
    report('S35', 868, 'PASS' if eqK == 0 and eqA == 0 and rate == 0 else 'FAIL',
           'Eq. (bc-win), (d_t - beta d_r) w_in = -N b_N trK/Theta, is EXACT for the outgoing linear wave in flat space with constant N, beta (T = N t,\n'
           '            R = r + beta t, d_t - beta d_r = N d_T, Theta = R); B^I = d_t w_in = beta d_r w_in - N b_N trK/Theta is the code\'s\n'
           '            G = bt*dw_in - NN*bN*trK/Th_b (collapse.cpp:L911; bt = NR + R_t, the shift relative to the moving grid)')
    # (e) scalar: Theta psi constant along outgoing light rays
    cc = D(Tf, r) / Tf + Q_c / 2 - K_c / 3
    lplus_Tpsi = sum((n_up[a] + e_up[a]) * D(Tf * pf, XS[a]) for a in range(4)).subs(Lf, 1).doit()
    ok_s = is_zero(lplus_Tpsi / Tf - (P_c + D(pf, r) + cc * pf))
    gfl = sp.Function('g')
    psi_out = gfl(t - r) / r
    ok_flat = sp.simplify(D(psi_out, t) + D(psi_out, r) + psi_out / r) == 0
    report('S35', 871, 'PASS' if ok_s and ok_flat else 'FAIL',
           '(u + e)(Theta psi) = 0  <=>  P + psi\' + c psi = 0, c = Theta\'/Theta + Q/2 - trK/3 (the code\'s c_b, collapse.cpp:L922-923); exact for\n'
           '            psi = g(u)/r in flat space ("r psi = g(u) must be constant along the outgoing characteristics, of speed ... 1", L861-862; L871-872); the entering combination is P + psi\' (S29), so its rate is d_t(P + psi\') = -d_t(c psi)')
    # flat-space form of the scalar condition: X = P + psi' + psi/r = r^-1 (d_t + d_r)(r psi), with P = d_t psi
    psF = sp.Function('psi')(t, rr)
    P_F = D(psF, t)
    Xc = P_F + D(psF, rr) + psF / rr
    okX = is_zero(Xc - D(rr * psF, t) / rr - D(rr * psF, rr) / rr)
    rate_id = sp.simplify(D(Xc, t) - (D(P_F + D(psF, rr), t) + P_F / rr))      # d_t X = d_t(P + psi') + P/r identically
    psg = gfl(t - rr) / rr
    Pg = D(psg, t)
    ent_g = sp.simplify(Pg + D(psg, rr))
    rate_g = sp.simplify(D(Pg + D(psg, rr), t) + Pg / rr)
    wave_g = sp.simplify(D(psg, t, 2) - D(rr ** 2 * D(psg, rr), rr) / rr ** 2)
    report('S35', 873, 'PASS' if okX and rate_id == 0 and is_zero(ent_g + gfl(t - rr) / rr ** 2) and rate_g == 0 and wave_g == 0
           and is_zero(D(ent_g, t) + Pg / rr) else 'FAIL',
           'flat space, P = d_t psi: X = P + psi\' + psi/r = r^-1 (d_t + d_r)(r psi), and d_t X = d_t(P + psi\') + P/r identically, so\n'
           '            d_t(P + d_r psi) = -P/r  <=>  d_t X = 0  <=>  (d_t + d_r)(r psi) = 0 once X = 0 at r_max before the wave arrives (quiescent\n'
           '            initial data; the equivalence is through the time derivative, whereas for the spin-0 pair it is pointwise, S35 check at main.tex L863);\n'
           f'            psi = g(t - r)/r solves the flat wave equation, has P + psi\' = {ent_g} and d_t(P + psi\') = -g\'/r^2 = -P/r: "the entering\n'
           '            field again changes at a rate set by its 1/r^2 part" (L873-875); the code imposes d_t(P + psi\' + c psi) = 0, c -> 1/r (S35 check at main.tex L871)')
    Rp, Rph = sp.symbols('R_P R_phi')
    e_sc, th_sc = sp.Matrix([1, 1]) / 2, sp.Matrix([[1, 1]])            # (P, phi) ordering, entering pair for beta + N
    proj = (sp.eye(2) - e_sc * th_sc) * sp.Matrix([Rp, Rph]) + e_sc * Bs
    codesc = sp.Matrix([Bs - Rph, Rph])                                  # code: d_t P = B - d_t psi', psi evolved unchanged
    d_out = sp.simplify((sp.Matrix([[1, -1]]) * (codesc - proj))[0])
    d_in = sp.simplify((th_sc * (codesc - proj))[0])
    out_wave = sp.simplify(d_out.subs(Bs, Rp + Rph))                      # exactly outgoing: B = R_P + R_phi
    report('S35', 875, 'INFO' if d_in == 0 and d_out != 0 and out_wave == 0 else 'FAIL',
           f'(code check) scalar: the code (collapse.cpp:L927-928) imposes the condition through d_t P only; the entering combination obeys d_t(P + psi\') = B,\n'
           f'            and d_t(P - psi\') differs from the projector form by {d_out}, which vanishes when B = R_P + R_psi\' (outgoing wave): for an\n'
           f'            outgoing wave the outgoing scalar equation is unchanged, as for the spin-0 pair (the projection is not stated in the paper)')
    # (f) outflow
    bre = sp.Symbol('beta_re', real=True)
    lamp_s0, lamp_sc = bre + 2 * sp.Rational(5, 4), bre + 2                # N = 2, v0 = 5/4 (any positive values)
    ok_o = (sp.solve_univariate_inequality(lamp_s0 < 0, bre, relational=False) == sp.Interval.open(-sp.oo, -sp.Rational(5, 2))
            and sp.solve_univariate_inequality(lamp_sc < 0, bre, relational=False) == sp.Interval.open(-sp.oo, -2)
            and sp.simplify((bsym - Nsym * vv) - (bsym + Nsym * vv)) == -2 * Nsym * vv)
    report('S35', 875, 'INFO' if ok_o else 'FAIL',
           '(code check; the outflow condition is not stated in the paper) lambda_- < lambda_+ = beta + N v0 < 0  <=>  beta < -v0 N\n'
           '            (spin-0); beta + N < 0  <=>  beta < -N (scalar): no characteristic enters, and the code gate (ICN_step, collapse.cpp:L1001-1030) switches\n'
           '            the conditions off')


def S36_bc_tests():
    """L857-860, L877-882: the boundary tests.  With w_in set to zero (held at 0) the entering field is frozen while the true one
    follows the Coulomb part -b_N F(u)/r^2 of the outgoing wave; the reflection relative to the
    incident w_out is |b_N| F(inf)/r_max^2 / max|w_out|, independent of v0 and b_N, ~ pulse width
    / r_max.  For the paper's pulse it is 3.48% (r_max = 100) and 2.33% (r_max = 150)."""
    header('S36', 877, 'Boundary-condition tests: reflection with and without the outgoing-wave condition', r"""
Setting
these fields to zero at $r_{\max}$, their value at infinity, therefore
reflects part of the wave, by the ratio of the two falloffs, i.e.\ about
the pulse width over $r_{\max}$.                                                  (L857-860)
These conditions reflect very little. A spin-0 pulse
$\trK\propto\exp[-(r-50)^2/16]$ crossing a flat-space boundary at
$r_{\max}=100$ and $150$, for $v_0=1$ and $1.267$, is reflected into
$\win$ below $10^{-5}$ of the incident $\wout$, against $3.5\%$
($r_{\max}=100$) and $2.3\%$ ($r_{\max}=150$) when $\win$ is set to zero.
A scalar pulse is reflected with $0.03\%$ of its amplitude.                      (L877-882)
""")
    how('(derived) data trK = g(r) = exp[-(r-50)^2/16], a_r = 0 at t = 0: the outgoing half is r trK = G(r - v0 t)/2,\n'
        'G = r g, i.e. f(u) = G(-v0 u)/2 and F(u) = int f du (Eq. asymptotic); ahead of the wave a_r = 0 fixes the static\n'
        'part to zero.  At r_max the true entering field is w_in = -b_N F(u)/r^2 (Eq. wout-asym); B^I = 0 holds it at 0, so the\n'
        'reflected w_in is b_N F(u)/r^2, at most |b_N| F(inf)/r^2; the incident w_out = -(b_N/v0)[2 f/r + v0 F/r^2].  The ratio\n'
        'is independent of v0 and b_N (F scales as 1/v0).')
    c0, w0 = dval('probe_bc')[0], dval('probe_bc')[1]
    s = np.linspace(0, 200, 2000001)
    ds = s[1] - s[0]
    G = s * np.exp(-(s - c0) ** 2 / w0)
    f = 0.5 * G
    Fh = np.cumsum(f[::-1])[::-1] * ds                    # v0 F(u) as a function of sigma = -v0 u
    ratios = []
    for Rm in dval('probe_bc')[2:]:
        ratios.append(100 * Fh[0] / Rm ** 2 / (2 * f / Rm + Fh / Rm ** 2).max())
    p100, p150 = dval('refl_B0_pct')
    report('S36', 880, 'PASS' if round(ratios[0], 1) == p100 else 'FAIL',
           f'r_max = 100: derived reflection {ratios[0]:.2f}% (printed {p100}%); ratio x r_max/100 = {ratios[0]:.2f}, {ratios[1] * 150 / 100:.2f} ~ const:\n'
           '            "by the ratio of the two falloffs, i.e. about the pulse width over r_max" (L858-860)')
    report('S36', 881, 'PASS' if round(ratios[1], 1) == p150 else 'FAIL',
           f'r_max = 150: derived reflection {ratios[1]:.2f}% (printed {p150}%; the B^I = 0 test of release/VALIDATION.md gives {BC_VALIDATION_UNPRINTED["prefix"][0][1]}%)')
    val_spin0, val_prefix, val_scalar = (BC_VALIDATION_UNPRINTED[k] for k in ('spin0', 'prefix', 'scalar'))
    ok_val = (max(val_spin0) < dval('refl_new') and all(abs(a - p100) <= 0.05 + 1e-9 and abs(b - p150) <= 0.05 + 1e-9 for a, b in val_prefix)
              and abs(100 * val_scalar - dval('refl_scalar_pct')) <= 0.005 + 1e-9)
    report('S36', 882, 'DATA' if ok_val else 'FAIL',
           f'release VALIDATION record: spin-0 reflection {min(val_spin0):.1e}-{max(val_spin0):.1e} < {dval("refl_new")} (printed "below 1e-5");\n'
           f'            pre-fix code {val_prefix[0][0]}/{val_prefix[1][0]}% (R = 100) and {val_prefix[0][1]}% (R = 150) = printed {p100}%, {p150}%; '
           f'scalar {val_scalar:.1e} = {100 * val_scalar:.3f}% -> printed {dval("refl_scalar_pct")}%')


def S37_extraction(v0s):
    """L886-894 (opening paragraph of the Results), L896-910 with Eq. (wbar) (L900-906), L911-921 and Fig. 3 caption
    (L926-932): the two runs (v0 = 1 and v0 = 1.267, S15) and the two steps; Eq. (wbar) (normalised kernel, a high-pass
    delta w_out) and its implementation in release/analysis (WBAR_CODE_UNPRINTED); filter widths; u = t - Theta_* in the
    caption is Eq. (retarded) for v0 = 1.

    HOW: the two speeds from S15.  Eq. (wbar) symbolically: on any interval [a, b] of t' (the integrals are written
    without limits; on data they run over the available dumps) the ratio returns a constant exactly, so delta of a
    constant is 0, and it is linear in w_out; on the whole line it also returns a + b t unchanged and t^2 + sigma^2 for
    t^2, and delta w_out has the transfer function 1 - exp(-sigma^2 omega^2/2) (a high-pass filter at omega ~ 1/sigma).
    The code: gaussian_running_mean transcribed (numpy) and compared on a test series with Eq. (wbar) integrated over the
    data range; the truncation at 3 sigma drops a kernel weight 1 - erf(3/sqrt2) = 0.27%; the record on the figure data."""
    header('S37', 886, 'Results: the two runs; Eq. (wbar); extraction of the pulse', r"""
We analyse two collapses of the same initial shell: one with the
baseline couplings, for which $v_0=1$ and the spin-0 and scalar
signals travel together, and one with $c_2=0$, for which $v_0=1.267$ and
the two signals separate (Sec.~\ref{sec:method}).                              (L886-889)
We first extract the spin-0 waveform from the outgoing characteristic
variable $\wout$ and check that it falls off as $1/\Theta$. We then
measure the speed of the signal in two independent ways: with a linear
probe pulse on the collapse background, and with the burst in the
Misner--Sharp luminosity, used as a clock between spheres.                      (L890-894)
To isolate the pulse, we split $\wout$ on each sphere into a slowly
varying background $\bar w_{\rm out}$ and a fast remainder
$\delta\wout$, smoothing $\wout$ in time with a Gaussian kernel:                (L897-899)
  \bar w_{\rm out}(\Theta,t)=
  \frac{\int\wout(\Theta,t')\,e^{-(t-t')^2/2\sigma^2}\,dt'}
       {\int e^{-(t-t')^2/2\sigma^2}\,dt'},
  \label{eq:wbar}                                                               (L901-904)
and $\delta\wout=\wout-\bar w_{\rm out}$.                                       (L906)
With a width $\sigma=10$, the remainder isolates the
pulse $\Theta\,\delta\wout$, whose amplitude depends on $\sigma$ but is
the same on every sphere for any $\sigma\in[5,15]$                              (L907-909)
Fig. 3 caption (abridged): baseline collapse ($v_0=1$) ... retarded time $u=t-\Theta_*$   (L927-931)
""")
    v0b = math.sqrt(float(v0sq_sym.subs(BASE)))
    report('S37', 887, 'PASS' if abs(v0b - 1) < 1e-7 and round(v0s, 3) == dval('v0_super') and v0s > 1 else 'FAIL',
           f'baseline v0 = {v0b:.9f} = 1, the speed of the scalar (light), so the two signals travel together; c2 = 0: v0 = {v0s:.6f}\n'
           f'            -> {dval("v0_super")} > 1, so the spin-0 signal outruns the scalar and the two separate (S15); the two steps are S37 (extraction,\n'
           f'            1/Theta falloff S33) and S38 (probe pulse and luminosity-burst clock)')
    # Eq. (wbar), symbolically
    t, tp, a, b, c0, c1_, sg, om = sp.symbols('t tp a b c0 c1_ sigma omega', real=True)
    sgp = sp.Symbol('sigma', positive=True)
    ker = sp.exp(-(t - tp)**2 / (2 * sgp**2))
    norm_ab = sp.integrate(ker, (tp, a, b))
    wbar = lambda f, lo, hi: sp.integrate(f * ker, (tp, lo, hi)) / sp.integrate(ker, (tp, lo, hi))
    const_ab = sp.simplify(sp.integrate(c0 * ker, (tp, a, b)) / norm_ab - c0)
    inf = (-sp.oo, sp.oo)
    lin = sp.simplify(wbar(c0 + c1_ * tp, *inf) - (c0 + c1_ * t))
    quad = sp.simplify(wbar(tp**2, *inf) - (t**2 + sgp**2))
    # transfer function: w = exp(i omega t') -> w_bar = exp(-sigma^2 omega^2/2) exp(i omega t)
    x = sp.Symbol('x', real=True)   # t' = t + x
    kx = sp.exp(-x**2 / (2 * sgp**2))
    tr = sp.simplify(sp.integrate(sp.expand_trig(sp.cos(om * (t + x))) * kx, (x, -sp.oo, sp.oo)) / sp.integrate(kx, (x, -sp.oo, sp.oo))
                     - sp.exp(-sgp**2 * om**2 / 2) * sp.cos(om * t))
    ok_sym = const_ab == 0 and lin == 0 and quad == 0 and tr == 0
    report('S37', 902, 'PASS' if ok_sym else 'FAIL',
           f'Eq. (wbar) is a normalised average: w_bar[c] - c = {const_ab} on ANY interval [a, b] of t\' (so delta w_out = 0 for a constant,\n'
           f'            also where the kernel is cut by the ends of the data), linear in w_out; on the whole line w_bar[a + b t] - (a + b t) = {lin},\n'
           f'            w_bar[t^2] - (t^2 + sigma^2) = {quad}, and w_bar[cos(omega t)] = exp(-sigma^2 omega^2/2) cos(omega t) (residual {tr}): delta w_out\n'
           f'            = (1 - exp(-sigma^2 omega^2/2)) w_out, a high-pass filter that keeps the narrow pulse and removes the broad hump')
    # the code: transcription of spin0_common.gaussian_running_mean (truncate = 3, weights exp(..) * np.gradient(t), renormalised)
    def grm(tt, y, sigma, truncate=3.0):
        out = np.empty_like(y)
        dtt = np.gradient(tt)
        for i in range(len(tt)):
            s_ = np.abs(tt - tt[i]) <= truncate * sigma
            w_ = np.exp(-0.5 * ((tt[s_] - tt[i]) / sigma) ** 2) * dtt[s_]
            out[i] = np.dot(w_, y[s_]) / w_.sum()
        return out, tt <= tt[-1] - truncate * sigma
    W = WBAR_CODE_UNPRINTED
    sig, dts = W['sigma'], W['dt']
    tt = np.arange(dts, 300 + dts / 2, dts)
    y = 0.5 * np.exp(-((tt - 120) / 40) ** 2) + 0.07 * np.exp(-((tt - 150) / 2) ** 2) * np.cos(tt - 150)   # broad hump + narrow pulse
    yb, full = grm(tt, y, sig, W['truncate'])
    yc, _ = grm(tt, np.full_like(tt, 3.7), sig, W['truncate'])
    trap = lambda f: np.sum(0.5 * (f[1:] + f[:-1]) * np.diff(tt))
    ye = np.array([trap(np.exp(-0.5 * ((tt - ti) / sig) ** 2) * y) / trap(np.exp(-0.5 * ((tt - ti) / sig) ** 2)) for ti in tt])
    dmax = np.max(np.abs(yb - ye))
    cut = 1 - math.erf(W['truncate'] / math.sqrt(2))
    ok_code = np.max(np.abs(yc - 3.7)) < 1e-12 and dmax < cut * np.max(np.abs(y)) and cut < 0.003
    report('S37', 901, 'PASS' if ok_code else 'FAIL',
           f'(code check) gaussian_running_mean (release/analysis/spin0_common.py) is Eq. (wbar) as a quadrature on the dumps, with the\n'
           f'            kernel renormalised where the ends of the data cut it (= the ratio of Eq. (wbar) over the data range) and truncated at\n'
           f'            |t - t\'| <= {W["truncate"]:g} sigma (dropped kernel weight 1 - erf(3/sqrt2) = {100 * cut:.2f}%): constant 3.7 returned to {np.max(np.abs(yc - 3.7)):.0e}\n'
           f'            at every t incl. the ends; test series (hump + pulse of amplitude 0.07, dt = {dts}, sigma = {sig:g}): max |w_bar_code - w_bar_eq| = {dmax:.1e}\n'
           f'            < {100 * cut:.2f}% of max |w| = {cut * np.max(np.abs(y)):.1e} (the dropped weight times max |w|), i.e. {100 * dmax / 0.07:.1f}% of the pulse amplitude')
    d, pu = W['maxdiff'], W['pulse']
    lo, hi = dval('sigma_range')
    ok_rec = (W['sigma'] == dval('sigma_filter') and max(d) <= 0.01 * min(pu) and lo <= dval('sigma_filter') <= hi)
    report('S37', 908, 'DATA' if ok_rec else 'FAIL',
           f'sigma = {dval("sigma_filter")} lies in [{lo}, {hi}]; it is the default --sigma of fig_wout_background.py, the value stored in\n'
           f'            data/wout_background.npz and used by the release recipe (no --sigma), and the scan 5-15 is its --sigma_scan; on the figure data\n'
           f'            max |Theta (w_bar_code - w_bar_eq)| = {min(d):.1e}-{max(d):.1e} on the seven spheres, <= {100 * max(d) / min(pu):.1f}% of the pulse amplitude\n'
           f'            {min(pu):.3f}-{max(pu):.3f}; the background is plotted only where the 3-sigma window is inside the data (t <= 269.5 of 300);\n'
           f'            the caption\'s u = t - Theta_* is Eq. (retarded) with v0 = 1 (baseline run)')


def S38_propagation_speed(v0s):
    """L936-967, Fig. probe caption (L971-977) and Fig. lum caption (L983-991): (i) the probe: by Sec. principal
    (S28) the outgoing spin-0 characteristic has (rdot + beta)/N = v0 exactly, whatever N and beta, so the probe speed
    needs no approximation of the background; the medians 0.01% / 0.02%, the early scatter (t <~ 40) and the bound 3e-4
    afterwards against the record (probe.log, data/probe_speed.npz).  (ii) The luminosity clock: L = -dm/dt|_Theta is the
    L of Eq. (lum) (S40); a signal at speed v in tortoise distance gives lag = Delta Theta_*/v, slope 1/v; the record's
    Theta_* is Eq. (retarded) with M_ADM = 3.7084; the fit gives v = 1.2406 -> 1.241, 2.07% below v0 -> "2%"; the scalar
    speed 1 (slope 1) is excluded; the burst coincides with the narrow pulse of w_out; the burst's u on every sphere
    (u = t_arr - Theta_*/v0) against the printed "u ~ 10"."""
    header('S38', 936, 'Propagation speed: probe pulse and luminosity clock', r"""
To test the spin-0 propagation, we add to the initial data a linear probe pulse of $\trK$ (amplitude
$10^{-6}$, centred at $r=150$), isolate it as the difference with the run
without it, and measure the speed of the peak of its $\wout$ relative to
the \ae ther, $v=(\dot r_{\rm pk}+\beta)/N$.                                     (L936-939)
As shown in Sec.~\ref{sec:principal}, the outgoing spin-0
characteristic moves at exactly $v_0$ relative to the \ae ther, measured in its
proper time and proper radial distance, so
this measurement requires no approximation of the background. Figure~\ref{fig:probe} shows $(v-v_0)/v_0$ during the
evolution. The scatter at early times, $t\lesssim40$, comes from the
splitting of the initial $\trK$ pulse, which has both $\wout$ and
$\win$ components, into its outgoing and ingoing parts; afterwards, the
deviation stays below $3\times10^{-4}$. The median speed agrees with $v_0$ to $0.01\%$ for $v_0=1$ and to $0.02\%$ for $v_0=1.267$.   (L940-947)
To carry out an independent measurement of the
speed of the nonlinear signal emitted by the collapse, without
using the characteristic variables, we use the Misner--Sharp
luminosity through a sphere of fixed areal radius,
$L=-dm/dt|_\Theta$ (its decomposition into \ae ther and scalar fluxes is
given in Sec.~\ref{sec:flux}).
For $v_0=1.267$, $L$ shows on every sphere a sharp burst at
$u\simeq10$--$15$, followed by a broad oscillation (Fig.~\ref{fig:lum}(a));
the burst coincides with the narrow pulse of $\wout$. We take the
arrival time of the burst on each of the spheres $\Theta=120$--$400$
and plot its lag with respect to the innermost sphere against
the tortoise distance $\Delta\Theta_*$ between the spheres
(Fig.~\ref{fig:lum}(b)). A signal moving at speed $v$ in tortoise
distance gives a straight line of slope $1/v$. The arrival times lie on
a line, and a linear fit gives $v=1.241$, $2\%$ below $v_0=1.267$, which
we attribute to the Schwarzschild tortoise approximation in
Eq.~\eqref{eq:retarded}.                                                         (L949-965)
The burst is therefore
carried by the spin-0 mode, not by the scalar field, which propagates at the speed of light.   (L966-967)
Fig. probe caption (abridged): the medians of the measured speeds, which agree
with $v_0$ to $0.01\%$ and $0.02\%$, respectively.                               (L976-977)
Fig. lum caption (abridged): spheres of constant $\Theta$, against retarded time $u=t-\Theta_*/v_0$,
for $v_0=1.267$. It is dominated by the \ae ther Coulomb term ... A linear fit (blue) gives the
speed $1.241$, $2\%$ below $v_0=1.267$ (black, dashed) ... The dotted line
is the speed of the scalar, $v=1$.                                               (L983-991)
""")
    how('(i) probe: the outgoing eigenvalue of the spin-0 principal part (S28, Mab) gives rdot = -beta + N v0, so (rdot + beta)/N\n'
        'is v0 identically in N and beta (gamma_rr = 1: proper radial distance = dr, aether proper time = N dt); the numbers\n'
        'against PROBE_MEASURED_UNPRINTED.  (ii) clock: lag = t_arr - t_arr0 fitted against Theta_* - Theta_*0 (numpy polyfit,\n'
        'transcribing fig_luminosity_clock.py); Theta_* recomputed from Eq. (retarded); u of the burst = t_arr - Theta_*/v0.')
    # (i) the probe
    lam = [e for e in Mab.eigenvals()]
    rd_out = [-e for e in lam if is_zero(-e - (-bsym + Nsym * v0_ab))]
    rel = sp.simplify((rd_out[0] + bsym) / Nsym) if rd_out else None
    report('S38', 940, 'PASS' if rel is not None and is_zero(rel - v0_ab) and not rel.has(Nsym) and not rel.has(bsym) else 'FAIL',
           f'outgoing spin-0 characteristic (S28): rdot = -beta + N v0, (rdot + beta)/N = {rel} = v0 exactly, independent of N and beta, in\n'
           '            proper radial distance (gamma_rr = 1) per aether proper time (N dt), so v = (rdot_pk + beta)/N needs no approximation\n'
           '            of the background (unlike u = t - Theta_*/v0 of Eq. retarded, S34)')
    P = PROBE_MEASURED_UNPRINTED
    ok_med = (tuple(round(x, 2) for x in P['med_pct']) == dval('probe_dev') == dval('probe_dev_text'))
    report('S38', 947, 'DATA' if ok_med else 'FAIL',
           f'medians (record, probe.log): (median - v0)/v0 = +{P["med_pct"][0]:.4f}% (v0 = 1), +{P["med_pct"][1]:.4f}% (v0 = 1.267) -> '
           f'{dval("probe_dev_text")[0]}% and {dval("probe_dev_text")[1]}% (body L947 and\n'
           f'            caption L977; probe amplitude {dval("probe_amp")[0]}, centred at r = {dval("probe_amp")[1]})')
    tb = dval('probe_scatter_t')
    ok_bd = max(P['late_max']) < dval('probe_bound') and all(e > l for e, l in zip(P['early_max'], P['late_max']))
    report('S38', 947, 'DATA' if ok_bd else 'FAIL',
           f'data/probe_speed.npz ({P["n"][0]} / {P["n"][1]} samples from t = {P["t_first"]:g}): for t >= {tb} max |v - v0|/v0 = {P["late_max"][0]:.2e} (v0 = 1),\n'
           f'            {P["late_max"][1]:.2e} (v0 = 1.267) < {dval("probe_bound"):.0e} ("stays below 3e-4"; the v0 = 1 maximum is within 1.5% of the bound);\n'
           f'            for t < {tb}: {P["early_max"][0]:.2e} / {P["early_max"][1]:.2e}, larger (the early scatter; its cause, the w_out/w_in splitting of the\n'
           '            initial trK pulse, is the interpretation of the text)')
    # (ii) the luminosity clock
    C = LUM_CLOCK_UNPRINTED
    Th = np.array(C['Theta'], float); Ts = np.array(C['Theta_star']); ta = np.array(C['t_arr'])
    Ts_eq = Th + 2 * C['MADM'] * np.log(Th / (2 * C['MADM']) - 1)
    dT, lag = Ts - Ts[0], ta - ta[0]
    slope, icpt = np.polyfit(dT, lag, 1)
    vfit = 1 / slope
    rms = float(np.sqrt(np.mean((lag - (icpt + slope * dT)) ** 2)))
    dev = 100 * (vfit / v0s - 1)
    ok_fit = (np.max(np.abs(Ts_eq - Ts)) < 0.006 and abs(vfit - C['v_fit']) < 1e-4 and round(vfit, 3) == dval('lag_speed')
              and round(-dev) == dval('lag_dev_pct') and (C['Theta'][0], C['Theta'][-1]) == dval('lum_spheres') and rms < 0.01 * lag[-1])
    report('S38', dline('lag_speed'), 'DATA' if ok_fit else 'FAIL',
           f'clock record (lum_clock.log, spheres Theta = {C["Theta"][0]}-{C["Theta"][-1]}): Theta_* = Eq. (retarded) with M_ADM = {C["MADM"]} to\n'
           f'            {np.max(np.abs(Ts_eq - Ts)):.3f}; lag = c + Delta Theta_*/v has slope 1/v; fit: slope {slope:.4f}, v = {vfit:.4f} -> {dval("lag_speed")},\n'
           f'            {dev:.2f}% vs v0 = {v0s:.4f} -> "{dval("lag_dev_pct")}%" (also the Fig. lum caption, L989); rms residual {rms:.2f} over a lag of\n'
           f'            {lag[-1]:.1f} ("lie on a line"); the scalar speed 1 (slope 1, the dotted line) misses the last sphere by\n'
           f'            {dT[-1] * (1 - slope):.0f} time units, so the burst is the spin-0 signal')
    rd_sc = sp.simplify((-(bsym - Nsym) + bsym) / Nsym)
    report('S38', 967, 'PASS' if rd_sc == 1 else 'FAIL',
           f'scalar: eigenvalues beta +/- N (S29), outgoing rdot = -beta + N, (rdot + beta)/N = {rd_sc}: the speed of light relative to the\n'
           '            aether, as printed ("which propagates at the speed of light"); v0 = 1.267 != 1 (S15)')
    u_b = ta - Ts / v0s
    W = C['wout_pulse_u']
    common = [th for th in W if th in C['Theta']]
    gap = max(abs(W[th] - u_b[C['Theta'].index(th)]) for th in common)
    report('S38', 957, 'DATA' if gap < 0.5 else 'FAIL',
           f'"the burst coincides with the narrow pulse of w_out": steepest point of the narrow pulse of Theta w_out (data/invariant_wout_dm.npz)\n'
           f'            at u = {", ".join(f"{W[th]:.2f} (Theta = {th})" for th in W)}; burst (min of L) at u = '
           f'{", ".join(f"{u_b[C["Theta"].index(th)]:.2f}" for th in common)} on Theta = {", ".join(map(str, common))}: within {gap:.2f}')
    lo_u, hi_u = dval('lum_burst_u')
    ok_u = round(min(u_b)) == lo_u and round(max(u_b)) == hi_u and all(lo_u - 0.5 <= x < hi_u + 0.5 for x in u_b)
    report('S38', dline('lum_burst_u'), 'DATA' if ok_u else 'CAVEAT',
           f'"on every sphere a sharp burst at u ~ {lo_u}-{hi_u}": u = t_arr - Theta_*/v0 of the burst (min of L) = '
           f'{", ".join(f"{x:.1f}" for x in u_b)}\n'
           f'            on Theta = {", ".join(map(str, C["Theta"]))}: from {u_b[0]:.1f} on the innermost to {u_b[-1]:.1f} on the outermost sphere, i.e. 10-15 to\n'
           f'            the nearest unit; the drift {u_b[-1] - u_b[0]:.1f} is Delta Theta_* (1/v - 1/v0) = {dT[-1] * (slope - 1 / v0s):.1f} of the 2% lag itself, so the range\n'
           '            "u ~ 10-15" (not a single u ~ 10) is the right statement')


def S39_misner_sharp_identity():
    """L1004: "the time derivative of the Misner-Sharp mass, evaluated at fixed areal radius, equals the flux of the
    stress-energy tensor through the sphere".  The identity behind it (cited in the paper, not printed), d_a m = (Theta^2/2)(G_a^b - delta_a^b G^(2)) d_b Theta with G^(2) = G^t_t + G^r_r, is
    proven exactly for a general 2-metric, with a negative control; at fixed Theta (d_b Theta = delta_b^R)
    it gives dm/dt|_Theta = (Theta^2/2) G^R_T = (Theta^2/2) T^R_T (S40).  Slow; runs in a worker."""
    header('S39', 997, 'Opening of "Mass conservation and radiated energy"; Misner-Sharp flux identity (exact)', r"""
We now turn to the energy carried by the wave. We first
check that, at every time, the change of the Misner--Sharp mass on each
sphere is accounted for by the flux of the stress-energy tensor
through it. We then show that most of this flux is
the \ae ther Coulomb field, which is linear in the wave amplitude, and
remove it to obtain the radiated energy.                                         (L997-1002)
The time derivative of the Misner--Sharp mass at fixed areal radius equals the flux of the
stress-energy tensor through the sphere~\cite{Hayward:1994bu}.                    (L1004, part)
""")
    how('(derived) fully general 2-metric g_AB(t,r) (3 free functions) + Theta(t,r)^2 dOmega^2, m = (Theta/2)(1 - (d Theta)^2)\n'
        '(= Eq. (ms), S17); G^a_b computed from scratch; both components a = t, r checked symbolically.  Negative\n'
        'control: with the 4-d trace G instead of the 2-d G^(2) the identity must fail.')
    show_slow('_job_ms', 'S39', 1004, 'symbolic proof of the Misner-Sharp identity')
    report('S39', 997, 'INFO', 'the two steps of the opening paragraph (L997-1002): the flux identity, Eq. (lum), is S39-S40 (exact) and\n'
           '            S42 (residuals on the run); the Coulomb part, linear in the amplitude, and its removal are S43-S46 (Eqs. dm-coulomb, mc,\n'
           '            Mae), the radiated energy S47-S48')


def S40_luminosity():
    """Eq. (lum), L1004-1015: in (T, R) = (t, Theta) coordinates the identity gives
    dm/dt|_Theta = (Theta^2/2) G^R_T; with T_ab built from (rho, j, p) in the code frame and
    G_ab = T_ab this is -(F_ae + F_psi) with the printed F_X; scalar projections checked.  The wording of
    L1004 ("Using the field equations G_ab = T_ab, one obtains the luminosity L") is tested: the identity plus G_ab = T_ab is
    sufficient, and the full tensor equation is needed (L contains p_X; the evolution equations and the
    Hamiltonian constraint alone lack d_t Q)."""
    header('S40', 1004, 'Eq. (lum): luminosity through fixed-Theta spheres', r"""
Using the field equations $G_{ab}=T_{ab}$, one obtains the luminosity $L$, split into \ae ther and scalar-field
contributions as                                                                  (L1004, part)
L\equiv-\left.\frac{dm}{dt}\right|_\Theta=F_{\text{\ae}}+F_\psi ,\\
F_X=\frac{N\Theta^2}{2}\left[A\,(\rho_X+p_X)+\frac{A^2+B^2}{B}\,j_X\right],        (L1007-1008, gathered)
with $F>0$ outward, $A$, $B$, $\rho_X$ and $j_X$ as in
Sec.~\ref{sec:mass} for each component $X$ of $T_{ab}$, and
$p_X=T_{ab}e^ae^b$ the radial pressure measured by the
\ae ther.                                                                         (L1012-1015)
""")
    how('(derived) in coordinates (T, R) = (t, Theta) we have d_b Theta = delta_b^R and\n'
        'delta_T^R = 0, so S39 gives dm/dt|_Theta = (Theta^2/2) G^R_T.  We build T_ab in the code\n'
        'coordinates from (rho, j, p) and the code\'s u_a = (-N,0), e_a = (beta,1), transform to (T,R)\n'
        'using d_t Theta = beta B + N A (the code\'s Theta equation), and compare with F_X.')
    Tt_s, Tr_s = sp.symbols('Theta_t Theta_r')
    ul = sp.Matrix([-Nn, 0]); el = sp.Matrix([bsym, 1])
    gi2d = sp.Matrix([[-1 / Nn ** 2, bsym / Nn ** 2], [bsym / Nn ** 2, 1 - bsym ** 2 / Nn ** 2]])
    uu2 = gi2d * ul; ee2 = gi2d * el
    okframe = sp.simplify((ul.T * uu2)[0] + 1) == 0 and sp.simplify((el.T * ee2)[0] - 1) == 0 and sp.simplify((ul.T * ee2)[0]) == 0
    Tll = rho * ul * ul.T + jj * (ul * el.T + el * ul.T) + pr_ * el * el.T
    chk = [sp.simplify((uu2.T * Tll * uu2)[0] - rho), sp.simplify(-(uu2.T * Tll * ee2)[0] - jj), sp.simplify((ee2.T * Tll * ee2)[0] - pr_)]
    Tmix = gi2d * Tll                                           # T^a_b
    dxdT = sp.Matrix([1, -Tt_s / Tr_s])                          # d x^a / dT at fixed R
    dRdx = sp.Matrix([Tt_s, Tr_s])                              # d R / d x^b
    TRT = (dRdx.T * Tmix * dxdT)[0]
    TRT = TRT.subs(Tt_s, bsym * Tr_s + Nn * AA).subs(Tr_s, BB)
    Lum = -Th_s ** 2 / 2 * TRT
    Fpaper = Nn * Th_s ** 2 / 2 * (AA * (rho + pr_) + (AA ** 2 + BB ** 2) / BB * jj)
    ok_lum = okframe and all(c == 0 for c in chk) and is_zero(Lum - Fpaper)
    report('S40', 1007, 'PASS' if ok_lum else 'FAIL',
           'L = -dm/dt|_Theta = -(Theta^2/2) T^R_T == F_ae + F_psi with the printed F_X, in the normalisation G_ab = T_ab',
           'the prefactor N Theta^2/2 requires G_ab = T_ab, which the paper states where psi is introduced (L315)')
    psi_t_ = bsym * psp + Nn * Pp
    dpsi = sp.Matrix([psi_t_, psp])
    kin2 = (dpsi.T * gi2d * dpsi)[0]
    Tps = dpsi * dpsi.T - sp.Matrix([[-Nn ** 2 + bsym ** 2, bsym], [bsym, 1]]) * kin2 / 2
    rp = sp.simplify((uu2.T * Tps * uu2)[0]); jp = sp.simplify(-(uu2.T * Tps * ee2)[0]); pp_ = sp.simplify((ee2.T * Tps * ee2)[0])
    report('S40', 1014, 'PASS' if is_zero(rp - rho_psi) and is_zero(pp_ - rho_psi) and is_zero(jp + Pp * psp) else 'FAIL',
           f'scalar (G_ab = T_ab, T^psi = d psi d psi - g (d psi)^2/2): rho_psi = {rp}, p_psi = {pp_}, j_psi = {jp}')
    eps_n = sp.symbols('epsilon', positive=True)        # null dust T_ab = eps k_a k_b, k = u +/- e: rho = p = eps, j = +/- eps
    F_out = sp.simplify(Fpaper.subs({rho: eps_n, pr_: eps_n, jj: eps_n}))
    F_in = sp.simplify(Fpaper.subs({rho: eps_n, pr_: eps_n, jj: -eps_n}))
    ok_sign = (is_zero(F_out - Nn * Th_s ** 2 * eps_n * (AA + BB) ** 2 / (2 * BB))
               and is_zero(F_in + Nn * Th_s ** 2 * eps_n * (AA - BB) ** 2 / (2 * BB)))
    report('S40', 1012, 'PASS' if ok_sign else 'FAIL',
           '"F > 0 outward": outgoing null dust (k = u + e) gives F = N Theta^2 eps (A+B)^2/(2B) >= 0, ingoing (k = u - e) gives\n'
           '            -N Theta^2 eps (A-B)^2/(2B) <= 0 (B > 0 outside the horizon); L = -dm/dt > 0 when energy leaves the sphere')
    how('(derived) the wording of L1004.  Sufficiency: the identity of S39 at fixed Theta, dm/dt|_Theta = (Theta^2/2) G^R_T,\n'
        'with G_ab = T_ab is the eq:lum check above.  Necessity of the full tensor equation: m of Eq. (ms) with theta_pm of Eq. (expansions) depends on Theta, Theta\', trK and Q;\n'
        'd_t m at fixed Theta therefore needs d_t Q.  Eq. (evolution) evolves psi, Theta, P, trK, a_r only (Q is fixed on each\n'
        'slice by Eq. (radial)), and Eq. (hamiltonian) has no time derivative.  Independently: F_X contains the radial pressure\n'
        'p_X (dF/dp = N Theta^2 A/2), while the only pressure-carrying printed equation, the trK line, carries rho + S with\n'
        'S = p_r + 2 p_perp (checked on the scalar: rho_psi + p_r + 2 p_perp = 2P^2, its source in the trK line, whereas\n'
        'rho_psi + p_r = P^2 + psi\'^2 is what F_psi needs).')
    Thq, Thpq, Qq, Kq, Pq, phq = sp.symbols('Theta Theta_p Q K P phi')
    Aq = Thq * (3 * Qq - 2 * Kq) / 6
    m_ms = Thq / 2 * (1 - (Thpq + Aq) * (Thpq - Aq))              # Eq. (ms) with Eq. (expansions)
    mQ = sp.simplify(sp.diff(m_ms, Qq))
    ok_mQ = is_zero(mQ - Thq ** 2 * Aq / 2) and not is_zero(mQ)
    evolved = {'psi', 'Theta', 'P', 'trK', 'a_r'}                 # left-hand sides of Eq. (evolution)
    rho_p, pr_p, pt_p = (Pq ** 2 + phq ** 2) / 2, (Pq ** 2 + phq ** 2) / 2, (Pq ** 2 - phq ** 2) / 2
    ok_S = is_zero(rho_p + pr_p + 2 * pt_p - 2 * Pq ** 2) and not is_zero(2 * Pq ** 2 - 2 * (rho_p + pr_p))
    dFdp = sp.simplify(sp.diff(Fpaper, pr_))
    ok_p = is_zero(dFdp - Nn * Th_s ** 2 * AA / 2)
    confirmed = ok_lum and ok_mQ and 'Q' not in evolved and ok_S and ok_p
    report('S40', 1004, 'PASS' if confirmed else 'FAIL',
           '"Using the field equations G_ab = T_ab, one obtains the luminosity L": the identity of the first sentence (S39) at fixed\n'
           '            Theta with G_ab = T_ab gives -dm/dt|_Theta = F_ae + F_psi exactly (above); the full tensor equation is needed:\n'
           f'            dF_X/dp_X = N Theta^2 A/2 != 0 (the rr equation) and dm/dQ = {mQ} needs d_t Q, absent from Eq. (evolution); the\n'
           '            printed trK line alone carries only rho + p_r + 2 p_perp (scalar: 2P^2), so the evolution equations and the Hamiltonian constraint alone would not suffice')


def S41_precursor_shear():
    """The precursor shear behind "Before the front arrives, the freely falling aether makes E slightly
    negative" (L1120-1121, checked in S48).  Q = -2 M t/((1 - c13) Theta^3) is derived here to O(t) from
    the code equations for any time-symmetric data (P = K = Q = a_r = 0 at t = 0), and with it the rise
    m - M_ADM = c13 M^2 t^2/(2(1-c13)^2 Theta^3) ahead of the front (used by S48).  The paper prints
    neither this shear nor the aether stress tensor."""
    header('S41', 1120, 'The precursor shear ahead of the front (derivation behind S48)', r"""
Before the front arrives,
the freely falling \ae ther makes $E$ slightly negative; ...                   (L1120-1121, abridged)
""")
    how('(derived) the precursor shear.  Time-symmetric data: P = K = Q = a_r = 0, N = 1, beta = 0 at t = 0.\n'
        'The code RHS (S10) give d_t K = d_t a_r = d_t Theta = d_t psi = 0 and d_t P = psi\'\' + 2 Theta\' psi\'/Theta\n'
        'at t = 0, so to O(t) only P changes.  The momentum constraint then gives\n'
        '(1-c13)(Theta^3 Q)\' = -Theta^3 P psi\' = -t Theta^3 psi\'(psi\'\' + 2 Theta\' psi\'/Theta), and the integral\n'
        'from 0 to r equals 2 m(r) t with m\'(r) = Theta^2 Theta\' psi\'^2/4 from the Hamiltonian constraint at t = 0.')
    at0 = {F['P']: 0, F['K']: 0, F['Q']: 0, F['a_r']: 0, F['dr_a_r']: 0, F['Kr']: 0, Nsym: 1, bsym: 0}
    ok0 = all(sp.simplify((ADV[k] + RHS[k]).subs(at0)) == 0 for k in ('trK', 'a_r', 'Theta', 'psi'))
    Pt = sp.simplify((ADV['P'] + RHS['P']).subs(at0))
    report('S41', 1120, 'PASS' if ok0 and is_zero(Pt - (F['psrr'] + 2 * F['Thr'] * F['psr'] / F['Th'])) else 'FAIL',
           f'at t = 0: d_t(trK, a_r, Theta, psi) = 0 and d_t P = {Pt}  (so K = O(t^2): "free fall from rest leaves trK ~ 0")')
    ThF, psF = sp.Function('Theta')(rr), sp.Function('psi')(rr)
    integrand = ThF ** 3 * D(psF, rr) * (D(psF, rr, 2) + 2 * D(ThF, rr) * D(psF, rr) / ThF)
    total_der = D(ThF ** 3 * D(psF, rr) ** 2 / 2, rr)
    mprime = ThF ** 2 * D(ThF, rr) * D(psF, rr) ** 2 / 4
    ok_int = is_zero(integrand - total_der - 2 * mprime)
    Thpp0 = -(D(ThF, rr) ** 2 - 1) / (2 * ThF) - ThF * D(psF, rr) ** 2 / 4
    m0 = ThF / 2 * (1 - D(ThF, rr) ** 2)
    ok_m = is_zero(D(m0, rr).subs(D(ThF, rr, 2), Thpp0) - mprime)
    report('S41', 1120, 'PASS' if ok_int and ok_m else 'FAIL',
           'Theta^3 psi\'(psi\'\' + 2 Theta\' psi\'/Theta) = d_r(Theta^3 psi\'^2/2) + 2 m\'(r), with m\' = Theta^2 Theta\' psi\'^2/4 (Ham. constraint at t = 0);\n'
           '            outside the shell psi\' = 0, hence Q = -2 m(r) t/((1-c13) Theta^3) + O(t^2), with m(r) = M_ADM')
    Msym = sp.Symbol('M', positive=True)
    Qs = -2 * Msym * t / ((1 - c13) * Th_s ** 3)
    Tp = sp.symbols('Theta_p', positive=True)
    dm_t2 = (Tp ** 2 / 2) * (-sp.Rational(3, 4) * c13 * Qs.subs(Th_s, Tp) ** 2)
    m_t2 = sp.simplify(-sp.integrate(dm_t2, (Tp, Th_s, sp.oo)))
    report('S41', 1120, 'PASS' if is_zero(m_t2 - c13 * Msym ** 2 * t ** 2 / (2 * (1 - c13) ** 2 * Th_s ** 3)) else 'FAIL',
           f'm(Theta, t) - M = -int_Theta^inf (Theta^2/2)(-3/4 c13 Q^2) dTheta = {m_t2} + O(t^3): m rises as t^2/Theta^3\n'
           '            ahead of the front for c13 > 0, and not at all in GR (c13 = 0, Birkhoff)')
    return m_t2, Msym


def S42_flux_identity_numbers():
    """L1016-1021 and Fig. fluxid caption (L1026-1033): the residual and the flux magnitudes."""
    header('S42', 1016, 'Flux identity: residuals and flux magnitudes', r"""
For $v_0=1.267$ ... the two sides of Eq.~\eqref{eq:lum} agree to within $0.02\%$ of the
peak flux at all times on every sphere from $\Theta=40$ to $200$
(Fig.~\ref{fig:fluxid}). ... The scalar flux is small, $|F_\psi|\lesssim0.02$, against \ae ther-flux
peaks of $0.9$--$4$, mostly the Coulomb contribution discussed below.            (L1016-1021, abridged)
Fig. fluxid caption (abridged): on the spheres $\Theta=40$, $100$ and $200$ ... in units of $10^{-4}$   (L1027-1033)
""")
    lo, hi = dval('flux_spheres')
    okf = all(lo <= x <= hi for x in dval('fluxid_spheres')) and dval('Fpsi_max') / dval('Fae_peak')[0] < 0.05
    report('S42', 1020, 'DATA' if okf else 'FAIL',
           f'|F_psi| <~ {dval("Fpsi_max")} is <= {100 * dval("Fpsi_max") / dval("Fae_peak")[0]:.1f}% of the aether peaks {dval("Fae_peak")}; the figure spheres '
           f'{dval("fluxid_spheres")} lie in {dval("flux_spheres")}; residual {dval("flux_resid_pct")}% = 2e-4 of the peak, the "units of 1e-4" of the caption')


def S43_coulomb_identity():
    """Eq. (dm-coulomb), L1044-1061: d_r m from m = (Theta/2)(1 - Theta'^2 + ...) with Theta'' from
    the Hamiltonian constraint and Q' from the momentum constraint; the O(eps) part of
    d_Theta m - (c14/2) Theta' d_Theta(Theta^2 a_r) vanishes for ANY Theta' (the Theta' factor
    makes the linear statement exact); dropping it leaves an O(M/Theta) remainder; the Q^2
    coefficient of the quadratic energy is -3/4 c13 (L1124)."""
    header('S43', 1052, 'Eq. (dm-coulomb): Coulomb term from the Hamiltonian constraint', r"""
At large radius, where $\trK,Q\to0$, the
Hamiltonian constraint reduces to                                                 (L1051-1052)
\partial_\Theta m={}&\frac{c_{14}}{2}\,\partial_r\Theta\;
\partial_\Theta\!\left(\Theta^2a_r\right)\\
&+\bigl(\text{terms quadratic in }a_r,\,\trK,\,Q,\,\psi\bigr),                  (L1055-1057, split)
with $\partial_r\Theta=1+O(M/\Theta)$. To leading order the linear term
is a total derivative                                                            (L1061-1062)
... this quadratic energy density is not positive definite, due to shear contributions from the \ae ther stress
tensor, yielding $-\tfrac{3}{4}c_{13}Q^2$, which is negative for $c_{13}>0$.          (L1124, part)
""")
    how('(derived) m = (Theta/2)(1 - Theta\'^2 + Theta^2(3Q-2K)^2/36); d_r m with Theta\'\' from the\n'
        'Hamiltonian constraint and Q\' from the momentum constraint (both shown in S10 to be the\n'
        'field equations); d_Theta m = d_r m / Theta\'.  Fields a_r, K, Q, psi, P scaled by eps; Theta\'\n'
        'kept GENERAL (= 1 + delta).')
    ThF = sp.Function('Theta')(rr)
    aF, KF, QF, psF, PF = [sp.Function(s)(rr) for s in ('a', 'K', 'Q', 'psi', 'P')]
    mF = ThF / 2 * (1 - D(ThF, rr) ** 2 + ThF ** 2 * (3 * QF - 2 * KF) ** 2 / 36)
    Thpp = sp.solve(sp.Eq(D(ThF, rr, 2) + (D(ThF, rr) ** 2 - 1) / (2 * ThF) + c14 * aF * D(ThF, rr)
                          + ThF / 4 * (c14 * (2 * D(aF, rr) + aF ** 2) + PF ** 2 + D(psF, rr) ** 2
                                       + sp.Rational(3, 2) * (1 - c13) * QF ** 2 - DEN / 3 * KF ** 2), 0), D(ThF, rr, 2))[0]
    Qp = sp.solve(sp.Eq((1 - c13) * (D(QF, rr) + 3 * QF * D(ThF, rr) / ThF) - DEN / 3 * D(KF, rr) + PF * D(psF, rr), 0), D(QF, rr))[0]
    dm = D(mF, rr).subs(D(ThF, rr, 2), Thpp).subs(D(QF, rr), Qp)
    dmdT = dm / D(ThF, rr)
    coul_paper = c14 / 2 * D(ThF, rr) * (D(ThF ** 2 * aF, rr) / D(ThF, rr))      # (c14/2) Theta' d_Theta(Theta^2 a_r)
    coul_noThp = c14 / 2 * D(ThF ** 2 * aF, rr) / D(ThF, rr)                     # the same term without the Theta' factor
    delta, ee_ = sp.symbols('delta epsilon')
    S_al, S_ap, S_K, S_Kp, S_Q, S_ps, S_psp, S_P, S_T = sp.symbols('al alp kk kkp qq ps psp pP TH')

    def scal(e):
        e = e.subs({D(aF, rr): ee_ * S_ap, D(KF, rr): ee_ * S_Kp, D(psF, rr): ee_ * S_psp})
        e = e.subs({aF: ee_ * S_al, KF: ee_ * S_K, QF: ee_ * S_Q, PF: ee_ * S_P, psF: ee_ * S_ps})
        return e.subs(D(ThF, rr), 1 + delta).subs(ThF, S_T)
    diffc = scal(dmdT - coul_paper)
    lin1 = sp.simplify(sp.diff(diffc, ee_).subs(ee_, 0))
    ord0 = sp.simplify(diffc.subs(ee_, 0))
    report('S43', 1055, 'PASS' if lin1 == 0 and ord0 == 0 else 'FAIL',
           'O(eps^0) and O(eps^1) parts of d_Theta m - (c14/2) Theta\' d_Theta(Theta^2 a_r) vanish for ANY Theta\' = 1 + delta:\n'
           '            with the Theta\' factor the linear statement is exact (without it a remainder\n'
           '            c14 delta [Theta a_r + Theta^2 a_r\'/2] is left, see the next check)')
    lin_noThp = sp.factor(sp.series(sp.simplify(sp.diff(scal(dmdT - coul_noThp), ee_).subs(ee_, 0)), delta, 0, 2).removeO())
    report('S43', 1061, 'PASS' if sp.simplify(lin_noThp) != 0 and sp.simplify(lin_noThp.subs(delta, 0)) == 0 else 'FAIL',
           f'"to leading order the linear term is a total derivative": dropping Theta\' leaves {lin_noThp} + O(delta^2),\n'
           '            which is O(delta) = O(M/Theta) relative to the Coulomb term (also "up to O(M/Theta) corrections", L1122)')
    q2c = sp.simplify(sp.diff(diffc.subs({delta: 0, S_al: 0, S_ap: 0, S_K: 0, S_Kp: 0, S_P: 0, S_ps: 0, S_psp: 0}), ee_, 2).subs(ee_, 0) / 2)
    rho_Q2 = sp.simplify(q2c / (S_T ** 2 / 2) / S_Q ** 2)
    report('S43', 1124, 'PASS' if is_zero(rho_Q2 + sp.Rational(3, 4) * c13) else 'FAIL',
           f'quadratic energy density (d_Theta m = (Theta^2/2) rho_eff): Q^2 coefficient = {rho_Q2} -> "-3/4 c13 Q^2" (L1124), G_ab = T_ab;\n'
           '            negative for c13 > 0, so the density is not positive definite')


def S44_coulomb_field():
    """L1037-1049 and Fig. coulomb caption (L1087-1094): the acceleration a_r ~ Theta^-2 of a
    static aether around a black hole (Theta^2 a_r -> M), the F(u)/r^2 term (S32)."""
    header('S44', 1037, 'The confusing luminosity; the aether Coulomb field a_r ~ Theta^-2', r"""
However, the luminosity gives a confusing answer
for the radiated energy. As the pulse crosses a sphere, $m$ swings by an
amount that grows with radius, and behind
the pulse it settles above $\MADM$
(Fig.~\ref{fig:coulomb}(b), thin lines), so that $\int L\,du$ is
negative on every sphere. We explain this behaviour below.                       (L1037-1042)
Let us first note that the spin-0 front sets up the acceleration $a_r\propto\Theta^{-2}$ of the
\ae ther around a black hole, absent from the initial data: behind the pulse $\Theta^2a_r$ relaxes to the same value on
every sphere ... the $F(u)/r^2$ term of
Eq.~\eqref{eq:asymptotic}.                                                       (L1044-1048, abridged)
We call such a contribution the \ae ther Coulomb field.                          (L1049)
Fig. coulomb caption (abridged): during the pulse its amplitude
grows with $\Theta$ (roughly as $\Theta^{0.7}$).                                 (L1090-1091)
""")
    Mm = sp.symbols('M', positive=True)
    astat = Mm / (Th_s ** 2 * sp.sqrt(1 - 2 * Mm / Th_s))
    lim = sp.limit(Th_s ** 2 * astat, Th_s, sp.oo)
    report('S44', 1044, 'PASS' if lim == Mm else 'FAIL',
           f'static aether along the Schwarzschild Killing vector (gamma_rr = 1: d_r = sqrt(1-2M/Theta) d_Theta): a_r = M/(Theta^2 sqrt(1-2M/Theta)),\n'
           f'            lim Theta^2 a_r = {lim}: "a_r ~ Theta^-2 around a black hole"; the static C/r^2 solution of the linearised system is S32')
    Ew = ENERGY_TIME_UNPRINTED
    ok_m = (all(a < b for a, b in zip(Ew['m_swing'], Ew['m_swing'][1:])) and min(Ew['m_late']) > 0 and min(Ew['m_min_u120']) > 0
            and max(Ew['intL']) < 0)
    report('S44', 1039, 'DATA' if ok_m else 'CAVEAT',
           f'record (data/aether_coulomb.npz, Theta = {", ".join(map(str, Ew["Theta"]))}): swing of m = {", ".join(map(str, Ew["m_swing"]))} ("grows with radius");\n'
           f'            m - M_ADM averaged over u = 200-250 = {", ".join(map(str, Ew["m_late"]))} > 0, and > 0 on every sphere from u = {max(Ew["m_last_below_u"]):.0f} on\n'
           f'            (min for u >= 120: {", ".join(map(str, Ew["m_min_u120"]))}): "settles above M_ADM"; int L du = -Delta m over the record = '
           f'{", ".join(map(str, Ew["intL"]))} < 0\n'
           f'            ("negative on every sphere"; from t = 0, where m = M_ADM outside the shell, it is -(m - M_ADM) at the end, '
           f'{Ew["intL_from_t0"][0]} ... {Ew["intL_from_t0"][1]})')
    report('S44', 1049, 'INFO', '"We call such a contribution the aether Coulomb field": a name (the static C/r^2 = F/r^2 solution, S32; its\n'
           '            linear contribution to m is Eq. (dm-coulomb), S43); the "Coulomb term" of Fig. lum caption and S48 is this term')
    report('S44', 1091, 'DATA', f'amplitude of Theta^2 a_r ~ Theta^{dval("Theta2ar_exp")} during the pulse: simulation number (the two-point values are not printed)')


def S45_mc():
    """Eq. (mc), L1067-1074: definitions of m_c and E."""
    header('S45', 1067, 'Eq. (mc): corrected mass and radiated energy', r"""
m_c\equiv m-\frac{c_{14}}{2}\,\Theta^2a_r ,
\qquad
E(\Theta,u)\equiv\MADM-m_c(\Theta,u),                                            (L1068-1070)
""")
    report('S45', 1068, 'INFO', 'definition; with S43: d_Theta m_c = quadratic terms + O(M/Theta) x linear; the coefficient c14/2 is that of\n'
           '            Eq. (dm-coulomb) (S43) and of the surface term of Eq. (Mae) (S46)')


def S46_Mae():
    """Eq. (Mae), L1075-1082: the surface term reduces to (c14/2) Theta^2 a_r in spherical symmetry
    (same coefficient as Eq. mc)."""
    header('S46', 1075, 'Eq. (Mae): Einstein-aether energy and its spherical reduction', r"""
M_{\text{\ae}}=\MADM-\frac{c_{14}}{8\pi G}\oint_\infty r^a a_a\,dA ,              (L1078)
which gives $\MADM-\tfrac{c_{14}}{2}\lim_{\Theta\to\infty}\Theta^2a_r$ in
spherical symmetry, positive for a class of solutions when $c_{14}<2$~\cite{Garfinkle:2011iw}.   (L1081-1082)
""")
    ar_s = sp.Symbol('a_r')
    surf = c14 / (8 * sp.pi) * 4 * sp.pi * Th_s ** 2 * ar_s
    report('S46', 1081, 'PASS' if is_zero(surf - c14 / 2 * Th_s ** 2 * ar_s) else 'FAIL',
           '(c14/8pi) x 4 pi Theta^2 a_r = (c14/2) Theta^2 a_r (r^a = unit radial normal, G = 1): same coefficient as Eq. (mc)')
    report('S46', 1082, 'INFO', 'positivity for c14 < 2 is the Garfinkle-Jacobson theorem (literature); with the static limit of S44,\n'
           '            M_ae = (1 - c14/2) M_ADM for a static aether')


def S47_Einf():
    """Eq. (Einf), L1098-1111, Fig. fluxinf caption (L1135) and the conclusion (L1450-1453): the
    arithmetic of the radiated-energy split with the printed numbers, with the M_ADM = 3.708 of
    the c2 = 0 run (L1098)."""
    header('S47', 1098, 'Eq. (Einf): radiated energy and its aether/scalar split', r"""
$E(\Theta)$ differs from the energy radiated to infinity by corrections of order $M/\Theta$, because the spheres lie at
finite radius. We therefore extrapolate to infinity with the two-parameter fit                   (L1099)
$E_\infty+b/\Theta+O(1/\Theta^2)$ over $\Theta=75$--$200$, and the scalar flux gives $E_\psi=0.197$.   (L1100)
E_\infty=0.532\pm0.023=(0.144\pm0.006)\,\MADM,\\
E_{\text{\ae}}=E_\infty-E_\psi=0.335=0.090\,\MADM,                                (L1104-1105, gathered)
so the spin-0 wave carries off $9\%$ of the ADM mass,   (L1109)
For the superluminal case $c_2=0$ with $\MADM=3.708$,                            (L1098)
Fig. fluxinf caption (abridged): $E_\infty=0.53\pm0.02$                          (L1135)
Conclusion: $14\%$ of the ADM mass:
$9\%$ in the spin-0 wave and $5\%$ in the escaping scalar field.                 (L1452-1453)
""")
    how('each derived quantity is computed on the interval of values compatible with the printed inputs\n'
        '(half a unit in their last digit) and accepted if the printed result lies in that interval.')

    def iv(x, nd):
        h = 0.5 * 10 ** (-nd)
        return (x - h, x + h)

    def consistent(printed, nd, lo, hi):
        h = 0.5 * 10 ** (-nd)
        return lo <= printed + h and hi >= printed - h
    (e_inf, e_err), (f_inf_p, f_err_p) = dval('Einf'), dval('Einf_frac')
    e_ae, f_ae_p = dval('E_ae')
    Ei, Ep = iv(e_inf, 3), iv(dval('E_psi'), 3)
    report('S47', 1105, 'DATA' if abs(e_inf - dval('E_psi') - e_ae) < 1e-9 else 'FAIL',
           f'E_ae = E_inf - E_psi = {e_inf} - {dval("E_psi")} = {e_inf - dval("E_psi"):.3f} = printed {e_ae} (exact)')
    Mu = iv(dval('MADM_super'), 3)
    Mp = Mu
    f_u = (Ei[0] / Mu[1], Ei[1] / Mu[0])
    ferr = ((e_err - 0.0005) / Mp[1], (e_err + 0.0005) / Mp[0])
    Eae = (Ei[0] - Ep[1], Ei[1] - Ep[0])
    fae_p = (Eae[0] / Mp[1], Eae[1] / Mp[0])
    ok_rest = consistent(f_err_p, 3, *ferr) and consistent(f_ae_p, 3, *fae_p) and consistent(dval('pct_ae') / 100, 2, *fae_p)
    report('S47', 1104, 'DATA' if consistent(f_inf_p, 3, *f_u) and ok_rest else 'FAIL',
           f'with the M_ADM = {dval("MADM_super")} of the c2 = 0 run (L1098): E_inf/M_ADM = {e_inf / dval("MADM_super"):.4f} (interval {f_u[0]:.4f}-{f_u[1]:.4f}) -> {f_inf_p};\n'
           f'            {e_err}/M = {e_err / dval("MADM_super"):.4f} -> {f_err_p}; {e_ae}/M = {e_ae / dval("MADM_super"):.4f} -> {f_ae_p} = {dval("pct_ae")}%')
    c14_, c9, c5 = dval('concl_split')
    ok_c = (round(100 * f_inf_p) == c14_ and c9 == dval('pct_ae') and c9 + c5 == c14_ and round(100 * dval('E_psi') / dval('MADM_super')) == c5 and dval('Einf_caption') == (round(e_inf, 2), round(e_err, 2)))
    report('S47', 1452, 'DATA' if ok_c else 'FAIL',
           f'conclusion: {c14_}% = round(100 x {f_inf_p}), {c9}% = Eq. (Einf), {c5}% = E_psi/M_ADM = {100 * dval("E_psi") / dval("MADM_super"):.1f}%, '
           f'{c9} + {c5} = {c14_};\n'
           f'            Fig. fluxinf caption {dval("Einf_caption")} = Eq. (Einf) rounded to two decimals')
    report('S47', 1110, 'INFO', f'c14(2-c14) = {float((c14 * (2 - c14)).subs(SUPER)):.2f} > 0 for the c2 = 0 run: positive linear spin-0 energy (literature criterion, S05)')


def S48_energy_in_time(m_t2, Msym):
    """L1113-1124: Fig. fluxinf in time.  The corrected rate L = -d_u m_c against the uncorrected L of Fig. lum(a), E on the
    four spheres and E_inf(u) against Eq. (Einf), on the record; E as the integral of the quadratic energy (S43); before
    the front E < 0, because the freely falling aether raises m above M_ADM as c13 M^2 t^2/(2(1-c13)^2 Theta^3) (S41) while
    a_r = O(t^3) there; the quadratic density is not positive definite (S43)."""
    header('S48', 1113, 'Fig. fluxinf: the corrected rate, E on all spheres, E_inf(u); E before the front; overshoot', r"""
Figure~\ref{fig:fluxinf} follows the corrected mass in
time. Once the Coulomb term is removed, the rate
$L=-\partial_u m_c|_\Theta$ (panel (a)) is more than an order of magnitude
smaller than the uncorrected luminosity of Fig.~\ref{fig:lum}(a), and
the cumulative energy $E$ (panel (b)) is nearly the same on all
spheres once the pulse has passed; its extrapolation to infinity, shown once the pulse has
crossed all spheres, settles on the value of Eq.~\eqref{eq:Einf}.               (L1113-1119)
Before the front arrives,
the freely falling \ae ther makes $E$ slightly negative; $E$ then overshoots its asymptotic value before settling to
the total radiated energy. Integrating Eq.~\eqref{eq:dm-coulomb} inward from a sphere not yet
reached by the wave shows that, up to $O(M/\Theta)$ corrections, $E$ is
the integral of the quadratic and higher-order field energy outside
$\Theta$. The transient overshoot is a consequence of the fact that this quadratic energy density is not positive
definite, due to shear contributions from the \ae ther stress tensor, yielding
$-\tfrac{3}{4}c_{13}Q^2$, which is negative for $c_{13}>0$.                       (L1120-1124)
""")
    Ew = ENERGY_TIME_UNPRINTED
    Lc = Ew['L_clock']
    corr_120 = float(np.interp(120, Ew['Theta'], Ew['L_corr']))
    r_same200 = Lc[200] / Ew['L_corr'][Ew['Theta'].index(200)]
    r_120 = Lc[120] / corr_120
    r_panels = max(Lc.values()) / max(Ew['L_corr'])
    r_inner = [u_ / c_ for u_, c_ in zip(Ew['L_uncorr_same'], Ew['L_corr'])]
    k = dval('corr_rate_factor')
    ok_r = r_same200 > k and r_120 > k and r_panels > k
    report('S48', dline('corr_rate_factor'), 'DATA' if ok_r else 'CAVEAT',
           f'"more than an order of magnitude smaller than the uncorrected luminosity of Fig. lum(a)": corrected max |L| (spin0_flux_infinity.npz)\n'
           f'            = {", ".join(map(str, Ew["L_corr"]))} on Theta = {", ".join(map(str, Ew["Theta"]))}; uncorrected (luminosity_clock.npz, Theta = 120-400)\n'
           f'            = {min(Lc.values())}-{max(Lc.values())}; same sphere Theta = 200: {Lc[200]}/{Ew["L_corr"][3]} = {r_same200:.0f}; Theta = 120 (corrected\n'
           f'            interpolated, {corr_120:.3f}): {r_120:.1f}; panel peaks {max(Lc.values())}/{max(Ew["L_corr"])} = {r_panels:.0f} (smallest Fig. lum(a) peak / largest\n'
           f'            fluxinf peak = {min(Lc.values()) / max(Ew["L_corr"]):.1f}).  Not printed: on the inner spheres of Fig. fluxinf, not shown in Fig. lum(a), the\n'
           f'            same-sphere ratio (uncorrected max |d_u m| from aether_coulomb.npz) is only {", ".join(f"{x:.1f}" for x in r_inner)} on Theta = 75, 100, 150, 200')
    Ea = Ew['E_avg']
    spr = max(Ea) - min(Ea)
    ok_e = spr / np.mean(Ea) < 0.1 and Ew['E_spread_late'] / np.mean(Ea) < 0.1
    report('S48', 1117, 'DATA' if ok_e else 'CAVEAT',
           f'"E (panel (b)) is nearly the same on all spheres once the pulse has passed": spread <= {Ew["E_spread_late"]} for all u >= 110; behind the pulse\n'
           f'            (energy.log, <E> over u = 200-250) E = {", ".join(map(str, Ea))} on Theta = 75 ... 200, spread {spr:.3f} = {100 * spr / np.mean(Ea):.1f}% of the mean\n'
           f'            ({spr / dval("Einf")[1]:.1f} x the quoted +-{dval("Einf")[1]}), ordered as b/Theta with b = {Ew["E_b"]} (the O(M/Theta) term that the extrapolation removes).  While the\n'
           f'            pulse crosses (u = 10-95) the curves differ by up to {Ew["E_spread_cross"][0]} (u = {Ew["E_spread_cross"][1]}), as large as E itself there: excluded by "once the pulse has passed"')
    b0, b1 = Ew['Einf_band_late']
    (e_inf, e_err) = dval('Einf')
    ok_s = e_inf - e_err <= b0 and b1 <= e_inf + e_err
    report('S48', 1119, 'DATA' if ok_s else 'CAVEAT',
           f'"its extrapolation to infinity, shown once the pulse has crossed all spheres, settles on the value of Eq. (Einf)": E_inf(u)\n'
           f'            (spin0_flux_infinity.npz key Einf_u) is plotted for u >= {Ew["Einf_u_fit"]} (fig_energy.py --u_fit), overshoots to {Ew["Einf_overshoot"][0]} at\n'
           f'            u = {Ew["Einf_overshoot"][1]:.0f} and stays in {b0}-{b1} for u = 200-255 (end of the record), inside {e_inf} +- {e_err} (the uncertainty is half\n'
           '            the peak-to-peak of E_inf(u) over u = 150-250, so the band is the quoted error by construction)')
    how('ahead of the front the data are time symmetric: a_r = O(t^3) (d_t a_r ~ K\' = O(t^2)), so m_c = m + O(t^3) and, from\n'
        'S41, E = M_ADM - m_c = -(m - M_ADM) = -c13 M^2 t^2/(2(1-c13)^2 Theta^3) + O(t^3).')
    Eahead = -m_t2
    val = Eahead.subs({c1: sp.Rational(7, 10), c3: 0, Msym: 1, t: 1, Th_s: 10})
    report('S48', 1120, 'PASS' if val < 0 else 'FAIL',
           f'E = {sp.simplify(Eahead)} < 0 ahead of the front for 0 < c13 < 1: "the freely falling aether makes E slightly negative"')
    report('S48', 1122, 'INFO', 'E as the integral of the quadratic energy outside Theta, up to O(M/Theta): S43; the overshoot follows from the\n'
           '            indefinite -3/4 c13 Q^2 term (S43)')


def _pn_objects():
    """Objects of Sec. V.E (sec:pn_comparison) and App. B: the PN-counting variables, the matrix D0 of App. B.2 (printed),
    v_L, s_{-1} = det D0 / minor_11 (elimination of B and K_ae), g1, x*, the truncated
    redefinition T1 = D0^-1 D1 to O(x) with D1 the diagonal of Eq. (M1), and s_1."""
    x, y = sp.symbols('x y')
    cc14, cc2, cc1 = sp.symbols('c14 c2 c1', positive=True)       # c13 = 0 in Sec. V.E (sec:pn_comparison) / App. B (c3 = -c1)
    R = sp.Rational
    D0 = sp.Matrix([[-1 + cc14 / 2 - (2 * cc14 - 1) * x, cc14 / 2 * x, 2 * cc14 * x],
                    [-R(9, 2) * cc2, -1 + (1 + R(3, 2) * cc2) * x, 6 * cc2],
                    [-R(1, 4) * (3 * cc2 + cc14) + cc14 * x, R(1, 4) * (cc2 - cc14) * x, cc2 - cc14 * x]])
    vL2 = cc2 * (2 - cc14) / (cc14 * (2 + 3 * cc2))
    detD0 = sp.factor(D0.det())
    minor = sp.factor(D0[1:, 1:].det())
    s_m1 = sp.simplify(detD0 / minor)
    g1 = cc14 / (2 * vL2)
    xs = 2 * vL2 / (2 - cc14)
    vT2 = cc1 / cc14
    D1 = sp.diag(-(1 - cc14 / 2) * (1 - x), -(1 - x), cc1 * (1 - x / vT2))
    T1 = (D0.inv() * D1).applyfunc(lambda e: sp.simplify(sp.series(e, x, 0, 2).removeO()))
    D1t = (D0 * T1).applyfunc(sp.expand)
    s1 = sp.simplify(D1t.det() / D1t[1:, 1:].det())
    kappa1 = cc14 * (2 - cc14) / (4 * vL2 ** 2)
    return dict(x=x, y=y, c14=cc14, c2=cc2, c1=cc1, D0=D0, vL2=vL2, detD0=detD0, minor=minor, s_m1=s_m1, g1=g1, xs=xs,
                vT2=vT2, D1=D1, T1=T1, D1t=D1t, s1=s1, kappa1=kappa1)


def S49_pn_section(PN):
    """Sec. V.E (sec:pn_comparison), L1140-1243: v_L^2 = v0^2 at c13 = 0; Eq. (pn_first_surviving) on a mode is
    k^2 s_1(x); kappa_1 ~ v_L^-4; Eq. (G_exact) exactly (with s_{-1} of App. B.2); Eq.
    (inverse_propagator) on the mode; Gamma singular at x = v_L^2; at c2 = 0 the spin-0 mode
    propagates only because c13 != 0.  The statements about Refs. Taherasghari-Will are INFO."""
    header('S49', 1140, 'Sec. V.E (sec:pn_comparison): the spin-0 mode in the post-Newtonian expansion', r"""
the PN analysis is
carried out for $c_{13}=0$, for which the spin-0 speed $v_0$ of
Eq.~\eqref{eq:speeds} reduces to the speed $v_L$ of Ref.~\cite{Taherasghari:2023rwn},
$v_L^2=c_2(2-c_{14})/[c_{14}(2+3c_2)]$.                                          (L1159-1162)
\left(1-\frac{c_{14}}{2}\right)\Box N
-\kappa_1\eta^2\partial_t^4\nabla^{-2}N+\cdots=0,                                 (L1173-1174)
where $\Box\equiv\nabla^2-\eta\partial_t^2$, $\eta$ counts PN
orders, and the coefficient $\kappa_1\propto v_L^{-4}$ is given in
Eq.~\eqref{eq:kappa1}.                                                           (L1177-1179)
\Gamma(x)&=\frac{-(1-\frac12c_{14})(1-x)}{s_{-1}(x)}\\
&=\frac{1-x/x_*}{1-x/v_L^2}=1+\frac{c_{14}}2\frac{y}{1-y},                        (L1191-1192)
where $x_*\equiv2v_L^2/(2-c_{14})$.                                              (L1196)
\frac{y}{1-y}
=\frac{\eta\partial_t^2}{v_L^2}
\left(\nabla^2-\frac{\eta}{v_L^2}\partial_t^2\right)^{-1}.                        (L1199-1201)
Moreover, $\Gamma(x)$ is singular
at the physical scalar pole $x=v_L^2$                                            (L1209-1210)
... the run used to measure the radiated energy has $c_2=0$, so that the
spin-0 mode propagates only because $c_{13}\neq0$.                               (L1238-1239)
""")
    x, y, cc14, cc2 = PN['x'], PN['y'], PN['c14'], PN['c2']
    vL2 = PN['vL2']
    v0_c13_0 = v0sq_sym.subs(c3, -c1).subs(c4, cc14 - c1).subs(c2, cc2)
    report('S49', 1162, 'PASS' if is_zero(v0_c13_0 - vL2) else 'FAIL', 'v0^2 of Eq. (speeds) at c13 = 0 (c3 = -c1) == v_L^2 = c2(2-c14)/[c14(2+3c2)] (also Eq. vL_def, L1757)')
    kk, om, eta = sp.symbols('k omega eta', positive=True)
    k1 = sp.Symbol('kappa_1')
    on_mode = (1 - cc14 / 2) * (-kk ** 2 + eta * om ** 2) - k1 * eta ** 2 * om ** 4 * (-1 / kk ** 2)      # d_t -> -i w, nabla^2 -> -k^2
    s1_printed = -(1 - cc14 / 2) * (1 - x) + k1 * x ** 2
    ok_m = is_zero(on_mode - kk ** 2 * s1_printed.subs(x, eta * om ** 2 / kk ** 2))
    k1_form = PN['kappa1']
    report('S49', 1174, 'PASS' if ok_m and is_zero(sp.diff(sp.log(k1_form), cc2) - sp.diff(sp.log(vL2 ** -2), cc2)) else 'FAIL',
           'Eq. (pn_first_surviving) on e^{-i w t} f, nabla^2 f = -k^2 f: k^2 [-(1-c14/2)(1-x) + kappa_1 x^2] = k^2 s_1(x) of Eq. (s1),\n'
           '            x = eta w^2/k^2 (2PN: eta^2); kappa_1 = c14(2-c14)/(4 v_L^4) ~ v_L^-4 (Eq. kappa1, S61)')
    G = -(1 - cc14 / 2) * (1 - x) / PN['s_m1']
    yy = x / vL2
    ok_G = (is_zero(G - (1 - x / PN['xs']) / (1 - x / vL2)) and is_zero(G - (1 + cc14 / 2 * yy / (1 - yy))))
    report('S49', 1191, 'PASS' if ok_G else 'FAIL',
           'Eq. (G_exact): with s_{-1} = det D0/minor_11(D0) (elimination of B, K_ae from the printed D0, S60), -(1-c14/2)(1-x)/s_{-1}\n'
           '            = (1 - x/x*)/(1 - x/v_L^2) = 1 + (c14/2) y/(1-y) EXACTLY (all orders in x), x* = 2 v_L^2/(2-c14), y = x/v_L^2')
    vL = sp.symbols('v_L', positive=True)
    ip = (eta * (-om ** 2) / vL ** 2) / (-kk ** 2 - eta / vL ** 2 * (-om ** 2))
    report('S49', 1200, 'PASS' if is_zero(ip - (y / (1 - y)).subs(y, eta * om ** 2 / (kk ** 2 * vL ** 2))) else 'FAIL',
           'Eq. (inverse_propagator): (eta d_t^2/v_L^2)(nabla^2 - eta d_t^2/v_L^2)^-1 on the mode = y/(1-y), y = eta w^2/(v_L^2 k^2)')
    num_at_pole = sp.simplify((1 - x / PN['xs']).subs(x, vL2))
    report('S49', 1210, 'PASS' if is_zero(num_at_pole - cc14 / 2) else 'FAIL',
           f'at x = v_L^2 the denominator of Gamma vanishes while its numerator is {num_at_pole} != 0: Gamma is singular at the scalar pole')
    v0_c2_0 = sp.factor(v0sq_sym.subs({c2: 0}))
    report('S49', 1239, 'PASS' if is_zero(v0_c2_0.subs(c3, -c1)) and not is_zero(v0_c2_0) else 'FAIL',
           f'c2 = 0: v0^2 = {v0_c2_0}, which vanishes for c13 = 0 (c123 = c13 + c2): the spin-0 mode propagates only because c13 != 0')
    report('S49', 1143, 'INFO', 'statements about Refs. Taherasghari-Will (2.5PN field equations, flux independent of v_L, Eqs. (4.7), (4.13a))\n'
           '            are literature inputs; the algebra quoted from them is checked in S59-S61')


def S50_white_hole_couplings():
    """L1255-1256: c2 = -c1^3/(2 - 4c1 + 3c1^2) = -32/45 = -0.7111 at c1 = 0.8 keeps v0 = 1 (S15)."""
    header('S50', 1255, 'White-hole couplings', r"""
We take $c_1=0.8$, $c_3=c_4=0$ and
$c_2=-c_1^3/(2-4c_1+3c_1^2)=-0.7111$, which keeps $v_0=1$.                      (L1255-1256)
""")
    c2wh = float((-c1 ** 3 / (2 - 4 * c1 + 3 * c1 ** 2)).subs(c1, sp.Rational(4, 5)))
    v0wh = math.sqrt(float(v0sq_sym.subs(WHITE)))
    report('S50', 1256, 'PASS' if round(c2wh, 4) == dval('wh_c')[1] and abs(v0wh - 1) < 1e-12 else 'FAIL',
           f'c2 = {c2wh:.7f} (= -32/45), v0 = {v0wh:.12f}')


def S51_ingoing_data():
    """Eq. (wh-id), L1257-1265: P = psi' + psi/r is the ingoing condition r psi = G(t + r) at N = 1,
    beta = 0.  Code check: with P != 0 the momentum constraint sources Q and Q^2 enters the
    Theta solve; solving Theta once with Q = 0 (a single pass) leaves a resolution-independent
    Hamiltonian-constraint violation at t = 0, while iterating the two solves to a fixed point
    (release/collapse.cpp, set_initial_data) removes it.  Both are reproduced here with a
    direct numerical solve of the constraints for the A0 = 0.10 data."""
    header('S51', 1257, 'Eq. (wh-id): ingoing shell; initial-data constraint solve', r"""
\psi(0,r)=\frac{A_0}{100}\,r^2\exp\!\Bigl[-\Bigl(\frac{r^2-r_0^2}{\sigma^2}\Bigr)^{2}\Bigr],\\
P(0,r)=\partial_r\psi+\frac{\psi}{r},                                             (L1260-1261, gathered)
so that $r\psi$ initially moves inward at the speed of light with $r_0=10$ and $\sigma=4$.   (L1265)
""")
    G_ = sp.Function('G')
    psi_in = G_(t + rr) / rr
    report('S51', 1261, 'PASS' if is_zero(D(psi_in, t) - (D(psi_in, rr) + psi_in / rr)) else 'FAIL',
           'r psi = G(t + r) (ingoing, flat) <=> psi_t = psi\' + psi/r; P = u(psi) = psi_t at N = 1, beta = 0 (leading order: N != 1 inside)')
    how('(code) solve Theta\'\' = -(Theta\'^2-1)/(2Theta) - (Theta/4)[P^2 + psi\'^2 + (3/2)(1-c13)Q^2] (RK4) and\n'
        '(Theta^3 Q)\' = -Theta^3 P psi\'/(1-c13) (K = a_r = 0) for the A0 = 0.10, c1 = 0.8 data, at dr = 4e-3 and\n'
        '2e-3: once with Q = 0 (a single pass), and iterated to a fixed point as release/collapse.cpp\n'
        'does; then evaluate the Hamiltonian-constraint residual Theta*C with finite differences.  The\n'
        'single-solve violation must be independent of dr, the iterated one must shrink with dr.')
    c1w, c3w = dval('wh_c')[0], 0.0                  # white-hole runs: c3 = c4 = 0 (L1255)
    c13w = c1w + c3w
    A0v, r0v, sv = dval('wh_A_wh'), 10.0, 4.0

    def constraint_solve(dr_):
        """Return (max|ThetaC| single Q = 0 solve, max|ThetaC| iterated solve, max|Q|, iterations) at spacing dr_."""
        rg = np.arange(0.0, 40.0 + dr_ / 2, dr_)
        gv = (rg * rg - r0v * r0v) / sv ** 2
        psv = A0v / 100 * rg * rg * np.exp(-gv * gv)
        prv = A0v / 100 * np.exp(-gv * gv) * 2 * rg * (1 - 2 * rg * rg * gv / sv ** 2)
        Pv = prv + np.where(rg > 0, psv / np.maximum(rg, 1e-300), 0.0)
        src_s = Pv ** 2 + prv ** 2

        def solve_theta(Qv):
            qq = 1.5 * (1 - c13w) * Qv ** 2 + src_s

            def rhs(x, y):
                s = np.interp(x, rg, qq)
                return np.array([y[1], -(y[1] ** 2 - 1) / (2 * y[0]) - y[0] / 4 * s])
            Th = np.zeros_like(rg)
            Th[1] = rg[1]
            y = np.array([rg[1], 1.0])
            for i in range(1, len(rg) - 1):
                x = rg[i]
                k1 = rhs(x, y); k2 = rhs(x + dr_ / 2, y + dr_ / 2 * k1)
                k3 = rhs(x + dr_ / 2, y + dr_ / 2 * k2); k4 = rhs(x + dr_, y + dr_ * k3)
                y = y + dr_ / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
                Th[i + 1] = y[0]
            return Th

        def solve_Q(Th):
            f = -Th ** 3 * Pv * prv / (1 - c13w)
            cum = np.concatenate([[0.0], np.cumsum(0.5 * (f[1:] + f[:-1]) * dr_)])
            safe_Th = np.where(Th > 0, Th, 1.0)              # Theta = 0 at r = 0 only: Q(0) = 0
            return np.where(Th > 0, cum / safe_Th ** 3, 0.0)

        def residual(Th, Qv):
            Thpp = (Th[2:] - 2 * Th[1:-1] + Th[:-2]) / dr_ ** 2
            Thp = (Th[2:] - Th[:-2]) / (2 * dr_)
            Tm = Th[1:-1]
            C = Thpp + (Thp ** 2 - 1) / (2 * Tm) + Tm / 4 * (src_s[1:-1] + 1.5 * (1 - c13w) * Qv[1:-1] ** 2)
            return np.abs(Tm * C)[rg[1:-1] > 0.5].max()
        Th0 = solve_theta(np.zeros_like(rg))        # single pass: Theta solved with Q = 0 ...
        Q0 = solve_Q(Th0)                           # ... then Q from the momentum constraint
        res_single = residual(Th0, Q0)
        Th, Qv = Th0, Q0
        for it in range(50):
            Th_new = solve_theta(Qv)
            dmax = np.abs(Th_new - Th).max()
            Th = Th_new
            Qv = solve_Q(Th)
            if dmax < 1e-12:
                break
        return res_single, residual(Th, Qv), np.abs(Qv).max(), it + 1
    out = {h: constraint_solve(h) for h in (4e-3, 2e-3)}
    for h, (ro, rn, qm, nit) in out.items():
        print(f'    dr = {h:.0e}: max|Theta C| = {ro:.4f} (single Q = 0 solve) / {rn:.1e} (iterated, {nit} iterations); max|Q| = {qm:.4f}')
    ro4, rn4 = out[4e-3][:2]
    ro2, rn2 = out[2e-3][:2]
    ok = ro2 > 1e-2 and abs(ro2 / ro4 - 1) < 0.05 and rn2 < rn4 / 3 and rn2 < 1e-4
    report('S51', 1262, 'PASS' if ok else 'FAIL',
           f'single Q = 0 solve (one pass): max|Theta C| = {ro4:.4f} -> {ro2:.4f} when dr halves (resolution independent);\n'
           f'            iterated Theta/Q solve (release/collapse.cpp set_initial_data, collapse.cpp:L1265-1280): {rn4:.1e} -> {rn2:.1e} (second order);\n'
           '            the violation of the single pass is removed by the iteration; time-symmetric data (P = 0 => Q = 0) are not affected')


def S52_theta_minus():
    """L1265-1270: (u - e)(Theta) = -theta_-, so theta_- < 0 means that ingoing rays diverge
    (anti-trapped); the GR statement about energy conditions is a literature input."""
    header('S52', 1265, 'theta_- and the anti-trapped condition', r"""
Besides $\theta_+$, we use the ingoing
expansion $\theta_-$ of Eq.~\eqref{eq:expansions}. Where
$\theta_-<0$ the sphere is anti-trapped: ingoing light rays diverge, as in a
white hole. In GR the energy conditions forbid anti-trapped surfaces from
forming out of regular data                                                      (L1265-1269)
""")
    lminus = sp.simplify(sum((n_up[a] - e_up[a]) * D(Tf, XS[a]) for a in range(4)).subs(Lf, 1).doit())
    thm = B_c - (3 * Q_c - 2 * K_c) * Tf / 6
    report('S52', 1267, 'PASS' if is_zero(lminus + thm) else 'FAIL',
           'the ingoing null normal u - e has (u-e)(Theta) = -theta_-: theta_- < 0 <=> ingoing rays diverge (anti-trapped)')
    report('S52', 1268, 'INFO', 'GR forbids anti-trapped surfaces forming from regular data under the energy conditions: literature')


def S53_white_hole_numbers():
    """L1272-1289: the A0 = 0.12 black hole and the A0 = 0.10 white-hole timeline; internal
    consistency of the quoted numbers."""
    header('S53', 1272, 'White hole or black hole: numbers', r"""
At
$A_0=0.12$ the data are initially untrapped everywhere ... but $\theta_+$ crosses zero by $t\simeq2.5$ and a black hole forms. ...
With $A_0=0.10$, after a brief dip of $\theta_+$ below zero at
$t\simeq6$--$17$, an anti-trapped region appears at $t\simeq20$ and ...
grows, its outer edge moving out with the escaping pulse, until it
splits at $t\simeq36$. The inner part shrinks after $t\simeq40$ and
disappears by $t\simeq44$, while a thin anti-trapped layer stays around
the peak of the escaping pulse and travels outward with it
(Fig.~\ref{fig:wh}(a)). The layer is present on all our grids and
narrows with the pulse, from about $0.9$ at $t=40$ to $0.06$ at $t=80$
on the refined grids, which agree on its width to within $20\%$; its
depth, $\min\theta_-$, still grows with resolution. ... whose
anti-trapped surface first appears at a similar time ($t\simeq17$); their
Fig.~11 also shows trapped surfaces at higher amplitude.                        (L1272-1288, abridged)
""")
    d0, d1, ta, tsplit, tg, tgone = dval('wh_timeline')
    ok = d0 < d1 < ta < tsplit < tg < tgone and abs(ta - dval('agg_at_t')) <= 3 and dval('wh_A_bh') > dval('wh_A_wh') and dval('wh_t_cross') < d0
    report('S53', 1276, 'DATA' if ok else 'FAIL',
           f'A0 = 0.10: dip t = {d0}-{d1} < anti-trapped at t = {ta} < split at {tsplit} < inner part shrinks after {tg} < gone by {tgone}; AGG\'s t ~ {dval("agg_at_t")} is "similar";\n'
           f'            the black hole needs the larger amplitude ({dval("wh_A_bh")} > {dval("wh_A_wh")}), as AGG\'s Fig. 11; theta_+ crosses zero at t ~ {dval("wh_t_cross")}')
    z = WH_SPLIT_ZEROS_UNPRINTED
    ok_sp = (len(z[34]) == 2 and len(z[36]) == 4 and z[36][-1] > z[34][-1]
             and z[44][1] - z[44][0] < 0.2 * (z[40][1] - z[40][0]) and z[44][3] > z[40][3])
    report('S53', 1278, 'DATA' if ok_sp else 'FAIL',
           'theta_- = 0 zeros (r > 1) on box40_N51200:\n            ' + '; '.join(f't = {t}: ' + ', '.join(f'{x:.2f}' for x in v) for t, v in z.items()) + '\n'
           f'            one region at t = 34 (outer edge 14.49 -> 16.82 at t = 36, moving out), two at t = 36 ("splits at t ~ {tsplit}"); the inner part\n'
           f'            is {z[36][1] - z[36][0]:.2f} wide at t = 36 and {z[40][1] - z[40][0]:.2f} at t = 40, then shrinks ("after t ~ {tg}") to a remnant\n            {z[44][1] - z[44][0]:.2f} wide at t = 44 (next dump t = 54: "disappears by t ~ {tgone}"), while the outer layer moves from Theta ~ 21 to 25.7 with the escaping pulse (S54)')
    w40, t40_, w80, t80_, pmax = dval('wh_layer_narrow')
    lw, lm = WH_LAYER_LADDER_UNPRINTED['width'], WH_LAYER_LADDER_UNPRINTED['min_theta_minus']
    pair = {t: max(abs(v[i + 1] / v[i] - 1) for i in range(2)) * 100 for t, v in lw.items()}
    spread = {t: (max(v) / min(v) - 1) * 100 for t, v in lw.items()}
    ok_w = (all(abs(x - w40) <= 0.07 for x in lw[t40_]) and all(abs(x - w80) <= 0.0085 for x in lw[t80_])
            and max(spread.values()) <= pmax and all(v[0] > v[1] > v[2] for v in lm.values()))
    report('S53', 1282, 'DATA' if ok_w else 'FAIL',
           'layer width (refined runs whconv_N12800 / 25600 / 51200): t = 40 ' + ' / '.join(f'{x:.3f}' for x in lw[40])
           + ', t = 64 ' + ' / '.join(f'{x:.3f}' for x in lw[64]) + ', t = 80 ' + ' / '.join(f'{x:.3f}' for x in lw[80]) + '\n'
           f'            ("about {w40} at t = {t40_} to {w80} at t = {t80_}": all within 0.07 resp. 0.0085);\n            spread over the three grids (max/min - 1) '
           + ', '.join(f'{spread[t]:.0f}%' for t in lw) + ' (t = 40, 64, 80; successive grids ' + ', '.join(f'{pair[t]:.0f}%' for t in lw)
           + f'): all within the printed {pmax}%;\n            depth min theta_- grows with resolution: t = 40 '
           + ' / '.join(f'{x:.2f}' for x in lm[40]) + ', t = 64 ' + ' / '.join(f'{x:.2f}' for x in lm[64]) + ', t = 80 '
           + ' / '.join(f'{x:.2f}' for x in lm[80]) + '\n            ("still grows with resolution")')


def S54_sharp_features():
    """L1291-1341, Fig. sharp caption (L1400-1408) and Fig. wh caption (L1423-1426): the scaling exponents are mutually
    consistent (rho ~ P^2 ~ dr^-1); width x height -> 0; 3 + 4 = 7 grids; the e-folding time of the compression; the two
    refinements (fixed factor 8 on the escaping pulse; 8 or 16 on the inner shell, held from t = 20 to 30.5)
    against the release recipe; the outermost theta_- = 0 branch is a layer < 1 unit wide on the |P| peak; the third
    feature (inner ingoing front) and the Fig. sharp caption against the record of fig_sharp_features.py and the
    front tracked on box40_N51200."""
    header('S54', 1291, 'Three sharp features: escaping pulse, unresolved inner shell, inner ingoing front', r"""
On coarse grids the white-hole profiles
exhibit three narrow features, each with a different physical origin.
Figure~\ref{fig:sharp}(a) shows them on one profile, at $t=40$, and
Fig.~\ref{fig:wh}(a) follows their areal radius in time: the escaping
pulse (green), the inner shell (black) and the inner front (yellow). ...         (L1291-1295)
white-hole geometry, its width shrinking roughly as $e^{-0.15t}$. We
resolve it by adding a moving refinement to the coordinate map of
Sec.~\ref{sec:method}. The number of grid points does not change: the
computational coordinate is first reparametrised so that the points
cluster in a region, a few units wide, where the spacing is eight times
finer, and are drawn from a band around it, while the grid elsewhere is
left exactly as it was. This is not adaptive mesh refinement: there are
no nested grids, and the refinement factor is fixed. The centre of the
region follows the peak of $|P|$ through a feedback controller            (L1296-1305)
The compression is visible in the outgoing light rays, which converge on
the pulse: the ray that leaves the inner shell at $t=30$, $3.2$ behind
the peak, is $0.3$ behind it at $t=34$ and $0.1$ at $t=40$
(Fig.~\ref{fig:wh}(a)).                                                        (L1309-1312)
piles up, as an ingoing wave, on a sphere of areal radius
$\Theta\simeq6.9$, and at $t\simeq29.5$ it steepens beyond what any of our
grids resolves (Fig.~\ref{fig:sharp}(c)). Figure~\ref{fig:wh}(b,c) measures the width and
height of the peak of Fig.~\ref{fig:sharp}(c) at $t=30$ on seven
grids: three uniform, and four with a refined region of the same kind,
eight or sixteen times finer. Since the shell forms abruptly, with no
earlier peak to follow, this region is placed at the radius where the
uniform runs show the shell forming and is held there from $t=20$ until
$t=30.5$. ... In each group of runs the width of the
peak shrinks roughly as $dr^{0.8}$, its height grows as $dr^{-0.5}$ and
the energy density as $dr^{-1}$.                                                 (L1315-1326)
A third, weaker feature forms inside the shell. From $t\simeq7$ part of
the ingoing scalar pulse travels as a front, a jump of about $0.3$ in
$P$, along ingoing light rays, $\dot r=-\beta-N$. Where it travels the
lapse drops to $N\simeq0.13$, and the front nearly stalls at
$\Theta\simeq3.2$, close to the inner edge of the anti-trapped region
(Fig.~\ref{fig:wh}(a)). Its steepest slope changes little with resolution at
$t=24$ but doubles with each halving of $dr$ by $t=40$
(Fig.~\ref{fig:sharp}(d)), so it too
steepens beyond our resolution. The lapse then recovers, the front falls
to the centre at $t\simeq44$ and has dispersed by $t\simeq54$.                 (L1332-1341)
Fig. sharp caption: The three sharp features of the white-hole run ($c_1=0.8$,
$A_0=0.10$), shown in $P$. (a) $|P|$ against $\Theta$ at $t=40$ (grey
line in Fig.~\ref{fig:wh}(a)) on the uniform grid with $dr=0.0008$. (b) The escaping pulse at $t=40$--$80$,
while its peak moves from $\Theta\simeq21$ to $62$, on the grid with a
region eight times finer that follows it, against the distance from the
peak. (c) The inner shell at $t=30$, against the distance from its peak,
on the three uniform grids (blue; $dr=0.0031$, $0.0016$, $0.0008$) and on
two grids refined around the shell (red; $F=8$ and $16$). (d) The inner
front at $t=40$ on the three uniform grids.                                     (L1400-1408)
Fig. wh caption (abridged): The ray joins the escaping
pulse by $t\simeq34$. The outermost
branch of purple dots is a thin anti-trapped layer, less than one unit
wide (its two edges overlap at this scale), that travels with the peak of
the escaping pulse.                                                              (L1423-1426)
""")
    w1 = dval('wh_width_exp')
    n7, n3, n4 = dval('wh_ngrids')
    ok = abs(2 * dval('wh_height_exp') - dval('wh_rho_exp')) < 1e-12 and n3 + n4 == n7 and w1 + dval('wh_height_exp') > 0 and dval('wh_fig_t') > dval('wh_t_steep')
    report('S54', 1325, 'DATA' if ok else 'FAIL',
           f'rho ~ P^2 ~ (dr^-0.5)^2 = dr^-1: exponents consistent; width x height ~ dr^{w1 + dval("wh_height_exp"):.2f} -> 0 (no\n'
           f'            resolved structure of fixed integral); {n3} + {n4} = {n7} grids; figure time t = {dval("wh_fig_t")} just after t ~ {dval("wh_t_steep")}; '
           f'compression e-folding 1/{dval("wh_compress_rate")} = {1 / dval("wh_compress_rate"):.1f}')
    fits = WH_EXPONENT_FITS_UNPRINTED
    ok_fit = all(abs(w - w1) <= 0.05 and abs(h - dval('wh_height_exp')) <= 0.05 and abs(r_ - dval('wh_rho_exp')) <= 0.06
                 for w, h, r_ in fits.values())
    report('S54', 1326, 'DATA' if ok_fit else 'FAIL',
           f'"in each group of runs ... roughly as dr^{w1}, dr^{dval("wh_height_exp")}, dr^{dval("wh_rho_exp")}" vs the rerun fits at t = 30: '
           f'width dr^{fits["uniform"][0]:.2f} (uniform) /\n'
           f'            dr^{fits["refined"][0]:.2f} (refined), height dr^{fits["uniform"][1]:.2f} / dr^{fits["refined"][1]:.2f}, '
           f'energy density dr^{fits["uniform"][2]:.2f} / dr^{fits["refined"][2]:.2f}: all within 0.06 of the printed exponents')
    rc = WH_REFINE_RECIPE_UNPRINTED
    t_on, t_off = dval('wh_shell_hold')
    ok_ref = (rc['zoneA_F'] == dval('wh_refine_factor') and tuple(rc['zoneB_F']) == tuple(dval('wh_shell_factor'))
              and (rc['refB_ton'], rc['refB_ttrack']) == (t_on, t_off) and t_on < dval('wh_t_steep') < dval('wh_fig_t') < t_off
              and abs(rc['refB_x0'] - rc['r_form']) < rc['refB_w'])
    report('S54', 1320, 'DATA' if ok_ref else 'FAIL',
           f'refinement factors and hold window vs release/README.md: escaping pulse "eight times finer" = --ref_F={rc["zoneA_F"]} (fixed\n'
           f'            factor, PI controller on the |P| peak); inner shell "eight or sixteen times finer" = --refB_F={rc["zoneB_F"][0]} / {rc["zoneB_F"][1]};\n'
           f'            "held there from t = {t_on} until t = {t_off}" = --refB_ton={rc["refB_ton"]} --refB_ttrack={rc["refB_ttrack"]}, which brackets the\n'
           f'            steepening (t ~ {dval("wh_t_steep")}) and the figure time (t = {dval("wh_fig_t")}); "placed at the radius where the uniform runs show\n'
           f'            the shell forming": --refB_x0={rc["refB_x0"]} vs r ~ {rc["r_form"]}, inside the half-width {rc["refB_w"]}')
    lay = WH_LAYER_MEASURED_UNPRINTED
    ok_lay = all(0 < th_o - th_i < dval('wh_layer_width') and th_i < th_p < th_o for th_i, th_o, th_p in lay.values()) and \
        lay[44][0] > lay[40][1]
    report('S54', 1424, 'DATA' if ok_lay else 'FAIL',
           '"a thin anti-trapped layer, less than one unit wide ... travels with the peak of the escaping pulse": on the\n'
           '            dr = 0.0008 run of Fig. wh(a) (box40_N51200) the theta_- = 0 edges are at\n'
           + ''.join(f'            Theta = {th_i:.2f} / {th_o:.2f} at t = {tt} (width {th_o - th_i:.2f}, max |P| at Theta = {th_p:.2f});\n'
                     for tt, (th_i, th_o, th_p) in lay.items())
           + '            both widths < 1 and the |P| peak inside the layer, which moves outward with it')
    g = RAY_GAP_UNPRINTED
    t0r, gap0, t1r, gap1, t2r, gap2 = dval('wh_ray')
    ok_rg = (round(g[t0r], 1) == gap0 and round(g[t1r], 1) == gap1 and round(g[t2r], 1) == gap2 and g[t0r] > g[t1r] > g[t2r]
             and dval('wh_ray_join') == t1r and t0r == dval('wh_fig_t'))
    report('S54', 1310, 'DATA' if ok_rg else 'FAIL',
           f'"the ray that leaves the inner shell at t = {t0r}, {gap0} behind the peak, is {gap1} behind it at t = {t1r} and {gap2} at t = {t2r}":\n'
           f'            escaping-pulse peak minus ray (fig_white_hole.py, box40_N51200) = ' + ', '.join(f'{v:.3f} at t = {t}' for t, v in g.items())
           + f';\n            {g[t0r]:.3f}, {g[t1r]:.3f}, {g[t2r]:.3f} round to the printed {gap0}, {gap1}, {gap2} and the gap shrinks monotonically; '
           f'caption "joins the escaping pulse by t ~ {dval("wh_ray_join")}" (L1424)')
    # the third feature: the inner ingoing front (L1332-1341), against the record (FRONT_MEASURED_UNPRINTED, box40_N51200
    # and wh_wh_N25600; SHARP_SCRIPT_UNPRINTED, release/analysis/fig_sharp_features.py)
    fm, ss = FRONT_MEASURED_UNPRINTED, SHARP_SCRIPT_UNPRINTED
    t7, jump, N_front, Th_st = dval('wh_front')
    t24, t40 = dval('wh_front_slope_t')
    t_fall, t_disp = dval('wh_front_fall')
    ok_nf = dval('wh_nfeatures') == 3
    jumps = list(fm['jump'].values()) + [e[2] for e in fm['early_wh_wh'].values()]
    ok_j = all(abs(abs(j) - jump) <= 0.06 for j in jumps) and all(j < 0 for j in jumps)
    report('S54', 1333, 'DATA' if ok_nf and ok_j else 'FAIL',
           f'"three narrow features" (escaping pulse, inner shell, inner front: {dval("wh_nfeatures")} = the three panels (b)-(d) of Fig. sharp);\n'
           f'            "a jump of about {jump} in P": P falls across the front by {min(abs(j) for j in fm["jump"].values()):.2f}-{max(abs(j) for j in fm["jump"].values()):.2f} '
           f'at t = 20-30 on box40_N51200 (extremum to extremum within\n'
           f'            0.3 in Theta; {abs(fm["jump_t32_36"][0]):.2f}, {abs(fm["jump_t32_36"][1]):.2f}, {abs(fm["jump_t32_36"][2]):.2f} at t = 32, 34, 36 as it steepens) and by '
           f'{min(abs(e[2]) for e in fm["early_wh_wh"].values()):.2f}-{max(abs(e[2]) for e in fm["early_wh_wh"].values()):.2f} at t = 4-8 on wh_wh_N25600')
    ok_ray = all(abs(fm['rdot'][t] - fm['minus_beta_minus_N'][t]) <= 0.0015 for t in fm['rdot']) and \
        all(abs(a - b) <= 0.015 for a, b in fm['early_rdot'].values())
    report('S54', 1334, 'DATA' if ok_ray else 'FAIL',
           '"along ingoing light rays, rdot = -beta - N" (the ingoing null direction of Eq. (metric), dr/dt = -beta - N): front speed\n'
           '            dr/dt (centred, t +- 2) vs -beta - N at the front on box40_N51200, t = 22-30:\n            '
           + ', '.join(f'{fm["rdot"][t]:.4f}/{fm["minus_beta_minus_N"][t]:.4f}' for t in fm['rdot'])
           + ';\n            on wh_wh_N25600, t = 6, 8, 10: ' + ', '.join(f'{a:.3f}/{b:.3f}' for a, b in fm['early_rdot'].values())
           + ' (all within 0.0015 resp. 0.015)')
    sl = [abs(e[1]) for e in fm['early_wh_wh'].values()]
    ok_t7 = min(fm['early_wh_wh']) < t7 < max(fm['early_wh_wh']) and all(sl[i] < sl[i + 1] for i in range(len(sl) - 1))
    report('S54', 1332, 'DATA' if ok_t7 else 'FAIL',
           f'"From t ~ {t7}": on wh_wh_N25600 (dumps every 2 units) the front is present at t = 4-8, steepest dP/dTheta\n            '
           + ', '.join(f'{e[1]:.2f} (Theta {e[0]:.2f})' for e in fm['early_wh_wh'].values())
           + ' at t = 4, 6, 8\n            (growing steadily: no sharp onset criterion, t ~ 7 is a round value inside the window)')
    fe = FRONT_EDGE_UNPRINTED
    ok_N = (all(abs(n - N_front) <= 0.01 for n in fm['N_front'].values()) and all(abs(th - Th_st) <= 0.06 for th in fm['Theta_stall'].values())
            and all(0 < e - f_ < 0.5 for f_, e in fe.values()))
    report('S54', 1335, 'DATA' if ok_N else 'FAIL',
           f'"the lapse drops to N ~ {N_front}": N at the front {min(fm["N_front"].values()):.3f}-{max(fm["N_front"].values()):.3f} at t = 24-28 (its minimum); '
           f'"nearly stalls at Theta ~ {Th_st}":\n            Theta of the front ' + ', '.join(f'{v:.2f}' for v in fm['Theta_stall'].values())
           + f' at t = 30-36 (dr/dt ~ {fm["rdot"][30]:.3f}); "close to the inner edge of the anti-trapped region": front / inner\n'
           '            theta_- = 0 edge ' + ', '.join(f'{f_:.2f} / {e:.2f} at t = {t}' for t, (f_, e) in fe.items()) + ' (just inside the edge, within 0.5)')
    s24, s40 = ss['slope'][24], ss['slope'][40]
    r24 = [s24[i + 1] / s24[i] for i in range(2)]
    r40 = [s40[i + 1] / s40[i] for i in range(2)]
    ok_s = (t24, t40) == (24, 40) and all(r < 1.3 for r in r24) and all(1.9 <= r <= 2.2 for r in r40) and dval('sharp_front_t') == t40
    report('S54', 1338, 'DATA' if ok_s else 'FAIL',
           f'steepest |dP/dTheta| (fig_sharp_features.py, 1.5 < Theta < 5, dr = 0.0031 / 0.0016 / 0.0008):\n            t = {t24}: '
           + ' / '.join(f'{x}' for x in s24) + f' (ratios {r24[0]:.2f}, {r24[1]:.2f}: "changes little");\n            t = {t40}: '
           + ' / '.join(f'{x}' for x in s40) + f' (ratios {r40[0]:.2f}, {r40[1]:.2f}: "doubles with each halving of dr");\n            '
           f'Fig. sharp(d) at t = {dval("sharp_front_t")}')
    ok_f = t40 < t_fall < t_disp and fm['Theta_t44'] < min(fm['Theta_stall'].values()) and fm['slope_t54'] < 0.01
    report('S54', 1340, 'DATA' if ok_f else 'FAIL',
           f'"falls to the centre at t ~ {t_fall} and has dispersed by t ~ {t_disp}": front at Theta = {fm["Theta_t44"]:.2f} at t = 44 (from ~ 3.5 at t = 40-42;\n'
           f'            on wh_wh_N25600 at Theta = 0.21), steepest slope in 1.5 < Theta < 5 at t = 54: {fm["slope_t54"]} (smooth);\n'
           f'            t ~ {t_fall} is also when the anti-trapped region is gone (S53)')
    tf, dfig = dval('sharp_fig')
    te0, te1, th0, th1, Fe = dval('sharp_esc')
    ts, drs, Fs = dval('sharp_shell')
    ok_c = (tf == dval('sharp_front_t') and dfig == dval('wh_fig_dr') and (te0, te1) == (40, 80)
            and round(ss['esc_peak'][te0]) == th0 and round(ss['esc_peak'][te1]) == th1 and Fe == dval('wh_refine_factor')
            and ts == dval('wh_fig_t') and tuple(drs) == tuple(dval('wh_s3_A03')[1:]) and abs(drs[2] - dfig) < 1e-12
            and tuple(Fs) == tuple(dval('wh_shell_factor')) and len(Fs) <= dval('wh_ngrids')[2])
    report('S54', 1402, 'DATA' if ok_c else 'FAIL',
           f'Fig. sharp caption: (a) t = {tf}, dr = {dfig} = the grid of Fig. wh(a); (b) peak of the escaping pulse at Theta = '
           f'{ss["esc_peak"][te0]} (t = {te0})\n            and {ss["esc_peak"][te1]} (t = {te1}) '
           f'(fig_sharp_features.py on whconv_N51200, --ref_F={Fe} --ref_ton=40) -> "~ {th0} to {th1}"; "eight times finer" = the\n            zone-A factor; '
           f'(c) t = {ts} = Fig. wh(b,c), dr = {drs} = r_max 40 / N = 12800, 25600, 51200 (0.003125,\n'
           f'            0.0015625, 0.00078125), F = {Fs[0]}, {Fs[1]} = two of the four refined grids of Fig. wh(b,c); (d) t = {dval("sharp_front_t")}')


def S55_no_unique_continuation():
    """L1343-1356: the frozen-centre artefact; its freezing time ~ dr^(-1/3) grows without bound
    as dr -> 0 (no number left to check)."""
    header('S55', 1343, 'No horizon, no unique continuation; frozen-centre artefact', r"""
the freezing time moves later on every finer grid ... roughly as $dr^{-1/3}$, with no sign of converging.   (L1355-1356, abridged)
""")
    report('S55', 1356, 'DATA', f'freezing time ~ dr^{dval("wh_freeze_exp"):.3f} -> infinity as dr -> 0 ("no sign of converging"): simulation statement '
           '(the five-point ladder is not printed)')


def S56_other_data():
    """L1358-1375: the other initial-data families; the quoted spacing ladders halve; amplitude
    orderings; the A0 = 0.03 artefact is present / absent / present (non-monotonic)."""
    header('S56', 1358, 'Other initial data', r"""
We repeated the amplitude exploration for
$(r_0,\sigma)=(10,3)$ and $(15,4)$, with the same outcome. Weak data
disperse ($A_0=0.01$, $0.015$ and $0.03$ for $\sigma=3$) ...
for $\sigma=3$, $A_0=0.08$ gives a trapped region that later disappears at
$dr\simeq0.009$, but a persistent black hole at $dr\simeq0.0047$. A trapped region
followed by an anti-trapped one does occur, for $r_0=15$, $A_0=0.028$:
the trapped region ... is gone by $t\simeq20$ and the anti-trapped
region appears at $t\simeq24$, identically on three grids. ... at $A_0=0.03$ it occurs at
$dr=0.0031$, not at $0.0016$, and again at $0.0008$, and at $A_0=0.04$ it
happens later on each finer grid                                                 (L1358-1374, abridged)
""")
    def halves(seq):
        return all(1.85 < seq[i] / seq[i + 1] < 2.15 for i in range(len(seq) - 1))
    A08, d09, d0047 = dval('wh_s3_A08')
    A03, *d03 = dval('wh_s3_A03')
    a028, tg, ta, ng = dval('wh_r15')
    ok = (halves((d09, d0047)) and halves(d03) and max(dval('wh_disperse_s3')) < A08 and A03 in dval('wh_disperse_s3') and tg < ta
          and A03 < dval('wh_s3_A04') < A08)
    report('S56', 1365, 'DATA' if ok else 'FAIL',
           f'spacing ladders halve ({d09}, {d0047}; {d03}); sigma = 3: dispersing {dval("wh_disperse_s3")} < near-threshold {A08}; the A0 = {A03} artefact is\n'
           f'            on a dispersing run, present / absent / present on the halving ladder (non-monotonic: "artefact"); r0 = 15, A0 = {a028}: '
           f'trapped gone by t ~ {tg} before anti-trapped at t ~ {ta}, on {ng} grids')


def S57_agg_comparison():
    """L1377-1395: the FWHM of AGG's shell is 2.5 times ours; their two grids differ by a factor 2,
    like our two coarsest uniform grids (factor-2 ladder ending at dr = 0.0008); Theta ~ 7."""
    header('S57', 1377, 'Comparison with Akhoury-Garfinkle-Gupta', r"""
the sharp features it finds near $t\simeq32$ are not fully resolved with
$40{,}000$ and $80{,}000$ points. ... Those
spacings match our two coarsest uniform grids ... their
shell, $(a_1/r)\,e^{-(r-r_0)^2/\sigma^2}$ with $\sigma=2$, is about $2.5$
times wider than ours, ... and their feature sits near the centre rather than at $\Theta\simeq7$.   (L1385-1393, abridged)
""")
    n1, n2 = dval('agg_N')
    dfin = dval('wh_fig_dr')
    report('S57', 1386, 'DATA' if n2 / n1 == 2 and round(dval('wh_Theta_shell')) == dval('agg_Theta_shell') else 'FAIL',
           f'AGG {n1} and {n2} points: a factor 2, like the two coarsest of three halving uniform grids ending at dr = {dfin} '
           f'({4 * dfin:.4f}, {2 * dfin:.4f}; AGG\'s spacings are not printed); Theta ~ {dval("agg_Theta_shell")} = {dval("wh_Theta_shell")} rounded')
    xs = np.linspace(0.01, 30, 300001)
    ours = xs ** 2 * np.exp(-((xs ** 2 - 100) / 16) ** 2)
    theirs = np.exp(-(xs - 10) ** 2 / dval('agg_sigma') ** 2) / xs

    def fwhm(y):
        idx = np.where(y >= y.max() / 2)[0]
        return xs[idx[-1]] - xs[idx[0]]
    rat = fwhm(theirs) / fwhm(ours)
    report('S57', 1391, 'PASS' if abs(rat - dval('agg_width_ratio')) < 0.2 else 'FAIL',
           f'FWHM ratio of the two profiles = {rat:.2f} (paper: ~2.5; Gaussian-width estimate 2/(sigma^2/2r0) = 2.5)')


def S58_conclusion():
    """L1436-1501: the conclusion's statements against the body of the paper."""
    header('S58', 1436, 'Conclusion', r"""
It is the
outgoing characteristic field of the evolution system, it carries the
whole part of the field that falls off as the inverse areal radius, and on
the dynamical black-hole background it propagates at the speed predicted
by the linear theory, also when that speed is superluminal.                    (L1442-1446)
... the collapse
radiates
$14\%$ of the ADM mass:
$9\%$ in the spin-0 wave and $5\%$ in the escaping scalar field.                (L1450-1453)
... Third, the large-$c_1$
case calls for adaptive mesh refinement, whose resolution, unlike the
fixed factor used here, grows as the shell steepens, to decide whether a
naked singularity forms.                                                       (L1498-1501)
""")
    report('S58', 1443, 'INFO', 'outgoing characteristic field and the whole 1/r part (S29, S33); speed (S38: probe 0.01% / 0.02%);\n'
           '            14% / 9% / 5% (S47); Misner-Sharp loss = stress-tensor fluxes (S42); dominant energy condition (S23, S24);\n'
           '            unresolved steepening of the inner shell (S54); "unlike the fixed factor used here": the refinement factors of\n'
           '            Sec. VI are fixed (8, and 8 or 16, S54); statements about Refs. TW: literature; "whereas the spin-0 speed also enters\n'
           '            the far-zone flux of Refs. Foster, Yagi" (L1458-1460) = the introduction (L188-191, "depends also on the speed v_L ...\n'
           '            i.e. the spin-0 speed v0")')
    report('S58', 1474, 'INFO', '"Inside the horizon the lapse collapses and the fields freeze, as in the GEJ simulations (Sec. method)" (L1474-1476)\n'
           '            = Sec. method L410-412 ("the fields inside that radius therefore freeze ... as in the original GEJ simulations")\n'
           '            and the record of L401-402 (DATA lapse_collapse, S13); "Ref. Bhattacharyya read this freezing as the dynamical\n'
           '            formation of a universal horizon" (L1476-1477) = L419-421; "which does form in the collapse of the version of\n'
           '            Horava gravity that survives the GW170817 and solar-system bounds [Franchini:2021bpt]" (L1477-1479) = the\n'
           '            introduction, L169-172 (same wording and citation); literature, not re-derived')


def S59_app_b_toy():
    """App. B.1, L1614-1668: the single-field example: operator on a mode, the truncated
    non-local redefinition T_p, the telescoping of Eq. (Lprime), p = 1, the order of the first
    propagation term, the roots."""
    header('S59', 1614, 'App. B.1: a single-field example of a non-local field redefinition', r"""
\left(\nabla^2-\frac{\eta}{v^2}\partial_t^2\right)\tilde\phi=0,                   (L1619)
\mathcal{L}=-k^2(1-y),\quad
y\equiv\frac{\eta\omega^2}{v^2k^2}
\;\longleftrightarrow\;
\frac{\eta}{v^2}\partial_t^2\nabla^{-2}.                                         (L1626-1629)
T_p=\sum_{j=0}^{p}y^j.                                                            (L1640)
\mathcal{L}'=\mathcal{L}T_p=-k^2(1-y)\sum_{j=0}^{p}y^j
=-k^2(1-y^{p+1}).                                                                 (L1645-1646)
The propagation correction is therefore absent through $p$PN and first
appears at $(p+1)$PN, through a term proportional to
$\eta^{p+1}\partial_t^{2p+2}(\nabla^{-2})^p$.                                    (L1649-1651)
\mathcal{L}'=-k^2+\frac{\eta^2\omega^4}{v^4k^2}.                                  (L1654)
The physical root $y=1$ is still a root of $\mathcal{L}'$. The polynomial
transformation, however, also introduces the roots $y^{p+1}=1$ with
$y\ne1$                                                                           (L1657-1659)
""")
    kk, om, eta, v = sp.symbols('k omega eta v', positive=True)
    y = sp.Symbol('y')
    Lmode = -kk ** 2 - eta / v ** 2 * (-om ** 2)                              # nabla^2 -> -k^2, d_t^2 -> -w^2
    yv = eta * om ** 2 / (v ** 2 * kk ** 2)
    op = eta / v ** 2 * (-om ** 2) * (-1 / kk ** 2)                           # (eta/v^2) d_t^2 nabla^-2 on the mode
    report('S59', 1626, 'PASS' if is_zero(Lmode - (-kk ** 2 * (1 - yv))) and is_zero(op - yv) else 'FAIL',
           'Eq. (simple_L): on e^{-i w t} f, nabla^2 f = -k^2 f, the operator is -k^2 (1 - y), and (eta/v^2) d_t^2 nabla^-2 -> y')
    report('S59', 1645, 'PASS' if all(sp.expand(-(1 - y) * sum(y ** j for j in range(n + 1)) + (1 - y ** (n + 1))) == 0 for n in range(6)) else 'FAIL',
           'Eq. (Lprime): -k^2 (1-y) sum_{j<=p} y^j = -k^2 (1 - y^{p+1}) (telescoping; p = 0..5 expanded)')
    ok_ord = True
    for n in range(5):
        term = kk ** 2 * yv ** (n + 1)                                           # the correction in L' = -k^2 + k^2 y^{p+1}
        opform = -(eta / v ** 2) ** (n + 1) * (-om ** 2) ** (n + 1) * (-1 / kk ** 2) ** n      # -(eta/v^2)^{p+1} d_t^{2p+2} (nabla^-2)^p
        ok_ord &= is_zero(term - opform)
    report('S59', 1650, 'PASS' if ok_ord else 'FAIL',
           'the correction k^2 y^{p+1} is the operator -(eta/v^2)^{p+1} d_t^{2p+2} (nabla^-2)^p on the mode (p = 0..4): order eta^{p+1}, (p+1)PN')
    report('S59', 1654, 'PASS' if is_zero(-kk ** 2 * (1 - yv ** 2) - (-kk ** 2 + eta ** 2 * om ** 4 / (v ** 4 * kk ** 2))) else 'FAIL',
           'Eq. (Lprime_p1): p = 1 gives -k^2 + eta^2 w^4/(v^4 k^2)')
    ok_r = all(sp.simplify((1 - y ** (n + 1)).subs(y, 1)) == 0 and len([z for z in sp.roots(1 - y ** (n + 1), y) if z != 1]) == n for n in range(1, 6))
    report('S59', 1657, 'PASS' if ok_r else 'FAIL', 'y = 1 is a root of 1 - y^{p+1} for every p, and the other p roots are the y^{p+1} = 1, y != 1 (p = 1..5)')


def S60_app_b_reduction(PN):
    """App. B.2, L1670-1785: the spherical reduction (LTF component; Laplacian identity;
    sin(kr)/r); rows 1 and 3 of D0 from the printed reduced Eqs. (4.5a), (4.6) on a mode, with
    the gauge condition; Box -> -k^2(1 - x); det D0; v_L; s_{-1} and g1; the 1PN term."""
    header('S60', 1670, 'App. B.2: spherical reduction and the linear scalar mode', r"""
\tilde B^{jk}=\delta^{jk}\tilde\gamma_1+\partial_j\partial_k\tilde\gamma_2,     (L1679)
so that $\tilde B^{rr}_{\rm LTF}=\tfrac23(\tilde\gamma_2''-\tilde\gamma_2'/r)$. Since     (L1682)
\left(\nabla^2-\frac{2}{r^2}\right)\partial_r\tilde\kappa
=\partial_r\nabla^2\tilde\kappa,                                                  (L1684-1685)
&(1-\tfrac12c_{14})\nabla^2\tilde N-\eta\ddot{\tilde N}
-\tfrac12c_{14}\eta\nabla^2\tilde B ...
+2\eta c_{14}\nabla^2\partial_t(\tilde\kappa_{\text{\ae}}-\tilde\kappa)=0, ...
&-\tfrac14c_{14}\big[\dot{\tilde N}+\eta\dot{\tilde B}
-4\eta(\ddot{\tilde\kappa}_{\text{\ae}}-\ddot{\tilde\kappa})\big] ...
-\tfrac14c_2\big[3\dot{\tilde N}+4\nabla^2\tilde\kappa_{\text{\ae}}
-\eta\dot{\tilde B}\big]=0,                                                       (L1692-1703, split)
\dot{\tilde N}+\nabla^2\tilde\kappa=0,\qquad
\dot{\tilde\kappa}+\tilde\gamma_1+\nabla^2\tilde\gamma_2=0.                       (L1707-1708)
On a mode satisfying $\nabla^2 f=-k^2f$ (for example,
$f=\sin(kr)/r$) ... x\equiv\frac{\eta\omega^2}{k^2} ... so that $\Box\to-k^2(1-x)$                      (L1712-1721, abridged)
\tilde K=(\omega/k)\hat{\tilde K},\qquad
\tilde K_{\text{\ae}}=(\omega/k)\hat{\tilde K}_{\text{\ae}},\qquad
\tilde B=(\omega/k)^2\hat{\tilde B},                                               (L1724-1726)
the gauge relation gives $\hat{\tilde K}=\tilde N$                                (L1728)
[D_0: the 3x3 matrix of L1740-1747, transcribed in _pn_objects]
\det \mathcal{D}_0=\tfrac12(2-c_{14})c_2(1-x)^2                                   (L1751)
v_L^2=\frac{c_2(2-c_{14})}{c_{14}(2+3c_2)}.                                       (L1757)
&=-\left(1-\tfrac12c_{14}\right)
\left[1-(1+g_1)x+O(x^2)\right],                                                  (L1769-1770)
g_1\equiv\frac{c_{14}}{2v_L^2}.                                                   (L1776)
\left(1-\tfrac12c_{14}\right)g_1x
=\frac{c_{14}(2-c_{14})}{4v_L^2}x.                                                (L1782-1783)
""")
    X3 = sp.symbols('X Y Z', real=True)
    rad = sp.sqrt(X3[0] ** 2 + X3[1] ** 2 + X3[2] ** 2)
    g2 = sp.Function('gamma2')
    Rs = sp.Symbol('R', positive=True)
    Hs = sp.Matrix(3, 3, lambda i, j: sp.diff(g2(rad), X3[i], X3[j]))
    nvec = sp.Matrix(X3) / rad
    pt = {X3[0]: Rs, X3[1]: 0, X3[2]: 0}
    Hrr = sp.simplify((nvec.T * Hs * nvec)[0].subs(pt))
    trH = sp.simplify(Hs.trace().subs(pt))
    stf = sp.simplify(Hrr - trH / 3)
    g2p, g2pp = sp.diff(g2(Rs), Rs), sp.diff(g2(Rs), Rs, 2)
    two3 = sp.simplify(stf - sp.Rational(2, 3) * (g2pp - g2p / Rs))
    report('S60', 1682, 'PASS' if two3 == 0 else 'FAIL',
           'B^rr_LTF = (2/3)(gamma2\'\' - gamma2\'/r) (L1682): the delta^jk gamma1 term is pure trace, and the trace-free rr part of\n'
           '            d_j d_k gamma2(r) is (2/3)(gamma2\'\' - gamma2\'/r) (A^<jk> = A^jk - delta^jk A^ll/3)')
    kap = sp.Function('kappa')
    lap = lambda f_: sp.diff(f_, rr, 2) + 2 / rr * sp.diff(f_, rr)
    ok_l = sp.simplify(lap(sp.diff(kap(rr), rr)) - 2 / rr ** 2 * sp.diff(kap(rr), rr) - sp.diff(lap(kap(rr)), rr)) == 0
    kk = sp.Symbol('k', positive=True)
    ok_s = sp.simplify(lap(sp.sin(kk * rr) / rr) + kk ** 2 * sp.sin(kk * rr) / rr) == 0
    report('S60', 1684, 'PASS' if ok_l and ok_s else 'FAIL',
           '(nabla^2 - 2/r^2) d_r kappa = d_r nabla^2 kappa (vector Laplacian of a gradient); nabla^2 (sin kr/r) = -k^2 sin kr/r')
    x, cc14, cc2 = PN['x'], PN['c14'], PN['c2']
    om, eta = sp.symbols('omega eta', positive=True)
    Nt, Kh, Bh = sp.symbols('N_t Kae_h B_h')
    I = sp.I
    # mode amplitudes: d_t -> -i w, nabla^2 -> -k^2; kappa = -i K/k (so that K = (w/k) K_hat with K_hat = N from the gauge)
    kap_ = -I * om / kk ** 2 * Nt                                  # gauge: N_dot + nabla^2 kappa = 0
    gauge_ok = sp.simplify(-I * om * Nt + (-kk ** 2) * kap_) == 0 and sp.simplify(I * kap_ * kk / (om / kk) - Nt) == 0
    kae_ = -I / kk * (om / kk) * Kh
    Bt = (om / kk) ** 2 * Bh
    dt, lapm = -I * om, -kk ** 2
    e45a = (1 - cc14 / 2) * lapm * Nt - eta * dt ** 2 * Nt - cc14 / 2 * eta * lapm * Bt + 2 * eta * cc14 * lapm * dt * (kae_ - kap_)
    e46 = (-cc14 / 4 * (dt * Nt + eta * dt * Bt - 4 * eta * dt ** 2 * (kae_ - kap_))
           - cc2 / 4 * (3 * dt * Nt + 4 * lapm * kae_ - eta * dt * Bt))
    sub = {eta: x * kk ** 2 / om ** 2}
    row1 = sp.expand(sp.simplify((e45a / kk ** 2).subs(sub)))
    row3 = sp.expand(sp.simplify((e46 / (-I * om)).subs(sub)))
    D0 = PN['D0']
    ok1 = all(sp.simplify(row1.coeff(v_) - D0[0, i]) == 0 for i, v_ in enumerate((Nt, Bh, Kh)))
    ok3 = all(sp.simplify(row3.coeff(v_) - D0[2, i]) == 0 for i, v_ in enumerate((Nt, Bh, Kh)))
    report('S60', 1740, 'PASS' if gauge_ok and ok1 and ok3 else 'FAIL',
           'the printed reduced Eqs. (4.5a) and (4.6) on the mode (kappa = -i K/k, the convention implied by K = (w/k) K_hat with\n'
           '            K_hat = N from the gauge condition N_dot + nabla^2 kappa = 0) give EXACTLY rows 1 and 3 of D0 (x = eta w^2/k^2)')
    report('S60', 1730, 'INFO', 'row 2 of D0, from Eq. (4.5c)_tr of Ref. Taherasghari-Will, is not printed in reduced form: taken as given')
    box = (lapm - eta * dt ** 2).subs(sub)
    report('S60', 1721, 'PASS' if sp.simplify(box + kk ** 2 * (1 - x)) == 0 else 'FAIL', 'Eq. (x_def): Box = nabla^2 - eta d_t^2 -> -k^2 (1 - x)')
    vL2 = PN['vL2']
    report('S60', 1751, 'PASS' if is_zero(D0.det() - sp.Rational(1, 2) * (2 - cc14) * cc2 * (1 - x) ** 2 * (1 - x / vL2)) else 'FAIL',
           f'Eq. (det_M0): det D0 = {PN["detD0"]} = (1/2)(2-c14) c2 (1-x)^2 (1 - x/v_L^2) with v_L^2 of Eq. (vL_def)')
    ser = sp.series(PN['s_m1'], x, 0, 2).removeO()
    report('S60', 1768, 'PASS' if is_zero(ser - (-(1 - cc14 / 2) * (1 - (1 + PN['g1']) * x))) else 'FAIL',
           'Eq. (sminus1): s_{-1} = det D0/minor_11 (B, K_ae eliminated) = -(1 - c14/2)[1 - (1 + g1) x + O(x^2)], g1 = c14/(2 v_L^2)')
    report('S60', 1782, 'PASS' if is_zero((1 - cc14 / 2) * PN['g1'] - cc14 * (2 - cc14) / (4 * vL2)) else 'FAIL',
           'Eq. (pn_1_term): (1 - c14/2) g1 = c14(2-c14)/(4 v_L^2): the v_L-dependent 1PN term beyond the Box N part')
    report('S60', 1760, 'INFO', 'interpretation of the (1-x)^2 factors as gauge/constraint modes and of x = v_L^2 as the physical scalar mode')


def S61_app_b_first_order(PN):
    """App. B.3, L1787-1907: the truncated redefinition T1 = D0^-1 D1 (to O(x)) reproduces the
    printed lapse row (Eq. N_transform) and first column; D0 T1 = diag + O(x^2) (Eq. M1); v_T;
    s_1 (Eq. s1) with kappa_1 (Eq. kappa1, both forms); s_1 = s_{-1}(1 + c14 y/2) + O(x^3); the
    series of Gamma; the operator identity; det D1 = det D0 det T1."""
    header('S61', 1787, 'App. B.3: first-order field redefinition', r"""
\tilde N
&=N+\frac{\eta c_{14}}{2-c_{14}}B
+\frac{\eta}{v_L^2}\dot R
\\
&=\left[1+\frac{c_{14}}{2v_L^2}
\eta\partial_t^2\nabla^{-2}\right]N
\\
&\quad
+x\left[\frac{c_{14}}{2-c_{14}}B
+\frac{2c_1}{v_L^2}K_{\text{\ae}}\right],                                         (L1795-1804)
\left(1,\tfrac32c_{14},
\frac{c_{14}+3c_2}{4c_2}\right)^T,                                                (L1818-1819)
\mathcal{D}_1(x)={}&\mathrm{diag}\Bigl(-(1-\tfrac12c_{14})(1-x),\,-(1-x),\\
&\qquad c_1(1-x/v_T^2)\Bigr)+O(x^2),                                              (L1824-1825)
where $v_T^2=c_1/c_{14}$ is the spin-1 speed for $c_{13}=0$.                    (L1829)
s_1(x)=-(1-\tfrac12c_{14})(1-x)+\kappa_1x^2+O(x^3),                              (L1839)
\kappa_1
&=g_1\frac{c_{14}(2+3c_2)}{2c_2}
\\
&=\frac{c_{14}(2-c_{14})}{4v_L^4}.                                               (L1845-1848)
\Gamma=1+\frac{c_{14}}2(y+y^2+\cdots).                                           (L1864)
s_1=s_{-1}\left(1+\frac{c_{14}}2y\right)+O(x^3).                                 (L1870)
&\left(\nabla^2-\frac{\eta}{v^2}\partial_t^2\right)
\sum_{j=0}^{p}
\left(\frac{\eta}{v^2}\partial_t^2\nabla^{-2}\right)^j
\\
&\qquad=\nabla^2\left[1-
\left(\frac{\eta}{v^2}\partial_t^2\nabla^{-2}\right)^{p+1}
\right],                                                                          (L1876-1882)
\det \mathcal{D}_1=\det \mathcal{D}_0\,\det T_1                                   (L1895)
In summary, in the first-order truncation the explicit $v_L$
dependence of the reduced lapse equation is absent through 1PN and
reappears at 2PN.                                                                (L1902-1905)
""")
    how('(derived) the first-order redefinition is linear in x: T1 = T1^(0) + x T1^(1).  Requiring D0 T1 = D1 + O(x^2) with D1\n'
        'the diagonal of Eq. (M1) fixes T1 uniquely (D0(0) is invertible): T1 = D0^-1 D1 expanded to O(x).  Its first row must be\n'
        'the printed lapse transformation, with eta d_t^2 nabla^-2 -> x on the mode and B, K_ae the rescaled amplitudes (B_hat).')
    x, cc14, cc2, cc1 = PN['x'], PN['c14'], PN['c2'], PN['c1']
    vL2, g1, T1 = PN['vL2'], PN['g1'], PN['T1']
    row_p = [1 + g1 * x, cc14 / (2 - cc14) * x, 2 * cc1 / vL2 * x]
    ok_row = all(is_zero(T1[0, j] - row_p[j]) for j in range(3))
    ok_col = all(is_zero(T1[i, 0].subs(x, 0) - v_) for i, v_ in enumerate((1, sp.Rational(3, 2) * cc14, (cc14 + 3 * cc2) / (4 * cc2))))
    report('S61', 1795, 'PASS' if ok_row and ok_col else 'FAIL',
           'T1 = D0^-1 D1 to O(x): first row (1 + g1 x, x c14/(2-c14), x 2c1/v_L^2) = Eq. (N_transform) [eta d_t^2 nabla^-2 -> x,\n'
           '            c14/(2 v_L^2) = g1]; first column at x^0 = (1, 3c14/2, (c14+3c2)/(4c2)) as printed (it cancels the x^0 response of B, K_ae)')
    off = (PN['D1t'] - PN['D1']).applyfunc(lambda e: sp.simplify(sp.series(e, x, 0, 2).removeO()))
    report('S61', 1824, 'PASS' if off == sp.zeros(3, 3) else 'FAIL', 'Eq. (M1): D1 = D0 T1 = diag(-(1-c14/2)(1-x), -(1-x), c1(1 - x/v_T^2)) + O(x^2)')
    v1_13 = sp.simplify(v1sq_sym.subs(c3, -c1).subs(c4, cc14 - c1).subs(c1, cc1))
    report('S61', 1829, 'PASS' if is_zero(v1_13 - PN['vT2']) else 'FAIL', f'v1^2 of Eq. (speeds) at c13 = 0 = {v1_13} = v_T^2 = c1/c14')
    s1ser = sp.series(PN['s1'], x, 0, 3).removeO()
    k1 = PN['kappa1']
    report('S61', 1839, 'PASS' if is_zero(s1ser - (-(1 - cc14 / 2) * (1 - x) + k1 * x ** 2)) else 'FAIL',
           'Eq. (s1): the reduced lapse operator of D0 T1 (B, K_ae eliminated) is -(1-c14/2)(1-x) + kappa_1 x^2 + O(x^3)')
    report('S61', 1845, 'PASS' if is_zero(g1 * cc14 * (2 + 3 * cc2) / (2 * cc2) - k1) else 'FAIL',
           'Eq. (kappa1): g1 c14(2+3c2)/(2c2) = c14(2-c14)/(4 v_L^4) (the two printed forms agree)')
    G = -(1 - cc14 / 2) * (1 - x) / PN['s_m1']
    yy = x / vL2
    report('S61', 1864, 'PASS' if is_zero(sp.series(G, x, 0, 4).removeO() - (1 + cc14 / 2 * (yy + yy ** 2 + yy ** 3))) else 'FAIL',
           'Gamma = 1 + (c14/2)(y + y^2 + y^3 + ...) (series of Eq. G_exact)')
    report('S61', 1870, 'PASS' if is_zero(sp.series(PN['s1'] - PN['s_m1'] * (1 + cc14 / 2 * yy), x, 0, 3).removeO()) else 'FAIL',
           's_1 = s_{-1}(1 + c14 y/2) + O(x^3): the first-order redefinition keeps the first term of the Gamma series')
    kk, om, eta, v = sp.symbols('k omega eta v', positive=True)
    opsum = lambda n: (-kk ** 2 - eta / v ** 2 * (-om ** 2)) * sum((eta / v ** 2 * (-om ** 2) * (-1 / kk ** 2)) ** j for j in range(n + 1))
    rhs_ = lambda n: -kk ** 2 * (1 - (eta / v ** 2 * (-om ** 2) * (-1 / kk ** 2)) ** (n + 1))
    report('S61', 1876, 'PASS' if all(is_zero(opsum(n) - rhs_(n)) for n in range(5)) else 'FAIL',
           'operator form of Eq. (Lprime) on modes (p = 0..4): (nabla^2 - eta d_t^2/v^2) sum_j (eta d_t^2 nabla^-2/v^2)^j = nabla^2 [1 - (...)^{p+1}]')
    report('S61', 1895, 'PASS' if is_zero(PN['D1t'].det() - PN['D0'].det() * T1.det()) else 'FAIL',
           'det D1 = det D0 det T1 (for the truncated T1 of above): the physical root of det D0 survives where det T1 != 0')
    # the summary: at fixed c14, v_L^2 = c2(2-c14)/(c14(2+3c2)) varies with c2, so "no explicit v_L" = no c2 at fixed c14
    sc_ = [sp.simplify(sp.expand(s1ser).coeff(x, n)) for n in range(3)]
    ok_sum = (all(is_zero(sp.diff(sc_[n], cc2)) for n in (0, 1)) and not is_zero(sp.diff(sc_[2], cc2))
              and not is_zero(sp.diff(PN['s_m1'], cc2)))
    report('S61', 1903, 'PASS' if ok_sum else 'FAIL',
           f'summary: s_1 = {sc_[0]} + ({sc_[1]}) x + kappa_1 x^2: the x^0 (Newtonian) and x^1 (1PN) coefficients do not depend on c2,\n'
           f'            i.e. on v_L at fixed c14, while kappa_1 = c14(2-c14)/(4 v_L^4) does (x^2 = 2PN); s_{{-1}} itself depends on v_L')


# =======================================================================================
# Main
# =======================================================================================


def main():
    print(__doc__.split('CONVENTIONS')[0].rstrip())
    print(f'paper source version: main.tex of {PAPER_SOURCE} (quoted line numbers refer to it);  '
          f'mode: {"QUICK" if QUICK else "FULL"}')
    start_slow_jobs()
    S00_data_table()
    S01_introduction()
    S02_action()
    S03_combinations()
    S04_speeds()
    S05_spin0_energy_sign()
    S06_gw170817()
    acc_lo = S07_lapse()
    S08_metric_and_variables()
    S08b_evolution_radial()
    S08c_bonetti_barausse()
    S09_nabla_u(acc_lo)
    S10_gej_from_action()
    print(f'  [{elapsed()}]')
    S11_hamiltonian()
    S12_constraint_propagation()
    S13_singularity_avoidance()
    S14_discretization_mesh()
    v0s = S15_couplings()
    S16_initial_data()
    S17_expansions_misner_sharp()
    S18_convergence_order()
    S19_regularised_constraint()
    S20_scattering()
    S21_collapse_numbers()
    S22_mass_numbers()
    S23_dm_dTheta()
    S24_mass_ratio()
    S25_evolution_block()
    S26_principal_part()
    S27_aNbN()
    S28_eigenvalues()
    S29_characteristic_variables()
    S30_linear_system()
    S31_wave_equation()
    Ksol, asol = S32_asymptotic_solution()
    S33_wout_asymptotic(Ksol, asol)
    S34_retarded_time()
    S35_boundary_conditions()
    S36_bc_tests()
    S37_extraction(v0s)
    S38_propagation_speed(v0s)
    S39_misner_sharp_identity()
    print(f'  [{elapsed()}]')
    S40_luminosity()
    m_t2, Msym = S41_precursor_shear()
    S42_flux_identity_numbers()
    S43_coulomb_identity()
    S44_coulomb_field()
    S45_mc()
    S46_Mae()
    S47_Einf()
    S48_energy_in_time(m_t2, Msym)
    PN = _pn_objects()
    S49_pn_section(PN)
    S50_white_hole_couplings()
    S51_ingoing_data()
    S52_theta_minus()
    S53_white_hole_numbers()
    S54_sharp_features()
    S55_no_unique_continuation()
    S56_other_data()
    S57_agg_comparison()
    S58_conclusion()
    S59_app_b_toy()
    S60_app_b_reduction(PN)
    S61_app_b_first_order(PN)

    print('\n' + '=' * 92)
    print('SUMMARY')
    print('=' * 92)
    cnt = Counter(s for _, _, s, _ in RESULTS)
    for sec, line, st, desc in RESULTS:
        if st in ('FAIL', 'CAVEAT'):
            print(f'  {st:7s} {sec} L{line}: {desc.splitlines()[0]}')
    print('  counts: ' + ', '.join(f'{k}={v}' for k, v in sorted(cnt.items())))
    print(f'  paper main.tex of {PAPER_SOURCE}; run time {time.time() - T_START:.0f} s ({"QUICK" if QUICK else "FULL"} mode)')
    return 1 if cnt.get('FAIL', 0) else 0


if __name__ == '__main__':
    sys.exit(main())
