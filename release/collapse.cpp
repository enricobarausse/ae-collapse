// =====================================================================================
//  collapse.cpp -- spherical scalar-field collapse in Einstein-aether theory
// =====================================================================================
//
//  Companion code of the paper
//      E. Albertini, M. Vaglio, E. Barausse, [arXiv number to be filled by the authors] (2026)
//  Please cite that paper if you use this code or data produced with it.
//  Copyright (c) 2026 Emma Albertini, Massimo Vaglio, Enrico Barausse.
//  Licence: MIT (see the file LICENSE).
//
//  Build:  g++ -O2 -std=c++11 -o collapse collapse.cpp     (or `make`; `make warn` adds
//          -Wall -Wextra, which must stay silent).  Single file, no dependencies beyond
//          the C++11 standard library.
//  Run:    ./collapse --help          (named options; README.md has the recipes)
//
//  FORMULATION
//  -----------
//  Garfinkle, Eling & Jacobson, Phys. Rev. D 76, 024003 (2007), gr-qc/0703093 (GEJ).
//  A massless scalar psi minimally coupled to Einstein-aether theory with couplings
//  c1, c2 (run time) and c3 = c4 = 0 (compile time).  Spherical symmetry; the aether is
//  hypersurface orthogonal and unit, so it is the future unit normal of the t = const
//  slices ("aether-fixed foliation").  Metric
//
//      ds^2 = -N^2 dt^2 + (dr + beta dt)^2 + Theta(t,r)^2 dOmega^2 ,      (paper eq:metric)
//
//  i.e. gamma_rr = 1: r is proper radial distance along the slice, and the AREAL radius
//  Theta is an evolved field.  Signature (-,+,+,+), u.u = -1, geometric units with the
//  GEJ normalisation of the scalar (G_ab = T^ae_ab + T^psi_ab, rho_psi = (P^2+psi'^2)/2).
//
//  The labels "eq:..." and "sec:..." below are the LaTeX labels of the equations and
//  sections of the paper (equation numbers move between versions; the labels do not).
//
//  Variables (names as in the code and in the output files):
//    Psi   = psi, the scalar field                                 (even in r)
//    P     = u^a d_a psi, its momentum along the aether           (even)
//    trK   = K, trace of the extrinsic curvature                  (even)
//    Acr   = a_r, radial aether acceleration, a_r = d_r ln N      (odd)
//    Theta = areal radius                                         (odd)
//    Q     = shear of K_ij (K^r_r = Q + K/3)                      (even)
//    NN    = lapse N;  NR = shift beta of eq:metric (N^r)          (odd)
//  Five fields (Psi, P, trK, Acr, Theta) are FREELY evolved; Q, N, beta are recovered on
//  every slice (and at every ICN stage) from the GEJ radial relations:
//
//    momentum constraint   Q' = -3 Q Theta'/Theta + [ -P psi' + (2+c13+3c2)/3 K' ]/(1-c13)
//    (eq:radial)           integrated OUTWARD from Q(0)=0 (Heun; first panel from the
//                          regular series Q(r_1) = rhs(r_1) r_1/5);
//    lapse                 (ln N)' = a_r  (GEJ eq. 28, no factor 1/2), trapezoid outward
//    (eq:lapse)            from the centre, then shifted so that N(r_max) = 1;
//    shift (eq:radial)     beta' = N (Q + K/3), trapezoid outward from beta(0) = 0.
//
//  Evolution equations (eq:evolution; d_t at fixed computational coordinate x; bt = beta + R_t is the
//  advection speed relative to the moving mesh, see below; c13=c1+c3, c14=c1+c4,
//  c123=c1+c2+c3):
//    psi_t   = bt psi' + N P
//    P_t     = bt P'   + N [ K P + a_r psi' + psi'' + 2 Theta' psi'/Theta ]
//    K_t     = bt K'   + N K^2/3 + N [ (c14-2)(a_r' + 2 Theta' a_r/Theta + a_r^2)
//                                      + 2 P^2 + 3(1-c13) Q^2 ] / (2+c13+3c2)
//    a_r,t   = bt a_r' + N [ (2K/3 - Q) a_r + c13/(c14(1-c13)) P psi'
//                            - c123/(c14(1-c13)) K' ]
//    Theta_t = bt Theta' + N Theta (Q/2 - K/3)
//  A mode with eigenvalue lambda of  d_t u = lambda d_r u  moves at dr/dt = -lambda.  The
//  spin-0 pair (K, a_r) has characteristic speeds bt +- N v0 with (eq:aNbN)
//    aN = (c14-2)/(2+c13+3c2),  bN = -c123/(c14(1-c13)),  v0 = sqrt(aN bN);
//  the scalar pair (P, psi') has bt +- N.  The default couplings c1 = 0.7,
//  c2 = -0.5119403 give v0 = 1.
//
//  The aether-fixed slicing is singularity avoiding: a_r grows inside a forming black hole,
//  N -> 0 there and the interior freezes.  There is NO excision and no inner boundary:
//  the grid extends to r = 0 throughout.
//
//  Hamiltonian constraint (eq:hamiltonian).  Used to build Theta at t = 0 (RK4 outward from Theta(0)=0,
//  Theta'(0)=1) and, multiplied by Theta to remove the 0/0 at the origin, as the
//  diagnostic ThetaC written to the profiles:
//    Theta C = Theta Theta'' + (Theta'^2-1)/2 + c14 a_r Theta Theta'
//              + Theta^2/4 [ c14(2 a_r' + a_r^2) + P^2 + psi'^2 + 3/2 (1-c13) Q^2
//                            - (2+c13+3c2)/3 K^2 ]  (= 0 in the continuum).
//
//  NUMERICAL METHOD
//  ----------------
//  * Uniform computational grid x_i = i dx, i = 0..N-1, dx = r_max/(N-1).  Second-order
//    centred differences; one-sided second-order stencils at i = N-1; parity at r = 0.
//  * Moving mesh (ALE, on by default):  r = R(x,t) = x + g(x) s(t), g a C^2 quintic taper
//    from 0 at x = 0.4 r_max to 1 at x = r_max; s(t) is chosen every step so that
//    Theta(r_max) stays fixed.  d_r = d_x/R_x, d_r^2 = d_x^2/R_x^2 - R_xx d_x/R_x^3, and
//    every advection term gets beta -> beta + R_t.  With the mesh off (--ale=0) the map is
//    the identity and the code is the plain fixed-grid scheme.
//  * Moving refinement (optional, needs the moving mesh): the computational coordinate is
//    first reparametrised, xi = X(x,t), so that up to two zones of width ~w are F times
//    finer; the taper then acts on xi.  See the block "MOVING REFINEMENT" below.
//  * Time integration: iterated Crank-Nicolson (Euler predictor, up to 12 trapezoidal
//    corrector iterations, tolerance 1e-10 on Psi and Theta), dt = CFL dx with CFL = 0.1
//    (CFL dx min(R_x) with the moving refinement).
//  * Fourth-order Kreiss-Oliger dissipation applied after each step, eps = 10.
//  * Outer boundary (paper sec:infinity): characteristic conditions on the INCOMING modes
//    of the scalar and spin-0 pairs, applied only when those modes actually enter (the
//    "gate"): the scalar is outgoing-spherical (Theta psi constant along outgoing rays) and
//    the incoming spin-0 field is evolved exactly as for the outgoing linear spin-0 wave
//    (eq:bc-win); Psi and Theta keep their interior RHS.
//
//  OUTPUT (in the working directory)
//    prof_tTTT.T.csv     full profile at each requested time (columns
//                        r,Psi,P,trK,Acr,Theta,NN,Q,ThetaC), 17 significant digits; the
//                        name keeps ONE decimal of the requested time (see dump_profile)
//    phi.dat             step, t, dt, Psi(r=0) every step
//    psi_times.dat       index of the optional auxiliary dumps (--dump_every)
//    *_dump_NNNNNN.csv   optional auxiliary per-field dumps (off by default)
//    spheres.bin/.txt    optional sphere recorder (--rec_spheres, off by default): every
//                        rec_every steps, a 13-point stencil of r and of all fields around
//                        each requested areal radius (block "SPHERE RECORDER" below)
//    stdout              progress (every 500 steps), apparent-horizon updates, FIRST NaN
//    stderr              set-up, outer-BC gate changes, refinement events, fatal errors
//
//  See README.md for the command line, the reproduction recipes and the file formats.
// =====================================================================================

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

using namespace std;

// =====================================================================================
//  RUN-TIME PARAMETERS (set from the command line in main)
// =====================================================================================
int    N        = 25600;     // number of grid points (incl. r = 0 and r = r_max)
double rmax     = 300.0;     // outer edge of the computational coordinate x
double t_final  = 0.0;       // stop time (only if have_t_final)
bool   have_t_final = false; // no stop time given: the run stops after NT_DEFAULT steps
vector<double> t_profile;    // times at which a full profile is written (sorted)
double c1       = 0.7;       // aether coupling c1
double c2       = -0.5119403;// aether coupling c2  (c1=0.7, c2=-0.5119403 => v0 = 1)
const double c3 = 0.0, c4 = 0.0;
// coupling combinations of the paper, set once in main (c3 = c4 = 0, so c13 = c14 = c1):
double c13 = 0.0, c14 = 0.0, c123 = 0.0;   // c1+c3, c1+c4, c1+c2+c3
double den_K = 1.0;          // 2 + c13 + 3 c2      (K equation, momentum constraint)
double den_a = 1.0;          // c14 (1 - c13)       (a_r equation)
// main stops with FATAL if |den_K| or |den_a| < 1e-12, so they are never divided by zero
double ingoing  = 0.0;       // != 0: purely ingoing scalar shell, P = psi' + psi/r
int    ale_on   = 1;         // 1: moving mesh (default); 0: fixed grid
int    dump_every = 0;       // auxiliary per-field dumps every k steps (0 = off)
vector<double> rec_theta;    // sphere recorder: areal radii (empty = off)
int    rec_every  = 1;       // sphere recorder: record every k steps
// initial data (the scalar shell of the paper, eq:id and eq:wh-id):
//     psi = A0/100 r^2 exp(-((r^2-r0^2)/sigma^2)^2)
struct ShellConfig {
    double A0    = 0.35;
    double r0    = 10.0;
    double sigma = 4.0;
} cfg;
// optional linear spin-0 probe for run-differencing: trK += probe_amp exp(-((r-probe_r)/sigma)^2)
// at t = 0, where sigma is the SHELL width cfg.sigma (--sigma, 4 by default)
double probe_amp = 0.0;
double probe_r   = 150.0;

// the scalar shell psi(0,r) and its radial derivative at an arbitrary radius
inline double shell_psi(double r) {
    double arg = (r*r - cfg.r0*cfg.r0) / (cfg.sigma*cfg.sigma);
    return r*r * cfg.A0 * exp(-arg * arg) / 100.0;
}
inline double shell_dpsi(double r) {
    double s2 = cfg.sigma * cfg.sigma;
    double g  = (r*r - cfg.r0*cfg.r0) / s2;
    return (cfg.A0 / 100.0) * exp(-g*g) * 2.0*r * (1.0 - 2.0*r*r*g / s2);
}

// =====================================================================================
//  FIXED NUMERICAL PARAMETERS
// =====================================================================================
const double CFL           = 0.1;
const int    NT_DEFAULT    = 60000;  // step budget when no stop time is given
const int    ICN_MAX_ITERS = 12;
const int    ICN_MIN_ITERS = 1;
const double ICN_TOL       = 1e-10;
const double EPS_KO        = 10;     // KO strength; per step sigma = eps (dt/dx)/16
// apparent-horizon search limited to r <= this (the finder is a diagnostic only; the
// analysis scripts find the horizons from the profiles, without this limit)
const double R_AH_SEARCH_MAX = 45.0;
// ingoing initial data: fixed-point iteration of the Theta and Q solves
const int    ID_MAX_ITERS  = 100;
const double ID_TOL        = 1e-13;  // on max |Delta Theta| between iterations

// =====================================================================================
//  GLOBAL STATE
// =====================================================================================
// computational spacing dx = rmax/(N-1): the d1_* / d2_* stencils are d/dx; it is NOT a
// physical dr (on the moving mesh dr = R_x dx; the physical derivatives are dr1_* / dr2_*)
double dx = 0.0;
double dt = 0.0;             // current time step
int    g_step = 0;           // step counter
int    g_horizon = 0;        // grid index of the outermost apparent horizon found so far

// evolved fields
vector<double> Psi, P, trK, Acr, Theta;
// recovered fields (Q: momentum constraint, NN: lapse, NR: shift), DTheta = d_r Theta
vector<double> Q, NN, NR, DTheta;

// per-step output streams (phi.dat, psi_times.dat).  Static, not local to main, so that
// the FATAL stops (exit()) still run their destructors and flush the buffered tail.
static ofstream fphi, ftimes;

// Outer-boundary characteristic gate (see compute_RHS).  Evaluated once per step in
// ICN_step from the beginning-of-step state so that every ICN stage applies the same
// boundary operator.  -1 = not yet evaluated (compute_RHS is never called before
// ICN_step sets the gates).
int    gate_spin0 = -1, gate_scalar = -1;
// time derivative of the recovered Q at the outer edge, used by the scalar outer BC
// (compute_RHS), set in ICN_step: backward difference over the previous step for the
// predictor stage, between the ICN stages for the corrector
double g_Qdot_edge = 0.0, g_Q_prev_edge = 0.0, g_dt_prev = 0.0;
double spin0_aN = 0.0, spin0_bN = 0.0, spin0_v0 = 1.0;   // set in main from c1, c2

// =====================================================================================
//  MOVING MESH (ALE)
// =====================================================================================
// r = R(x,t) = x + g(x) s(t).  g = 0 for x <= ale_ra = 0.4 rmax, the C^2 quintic
// u^3(10 - 15u + 6u^2), u = (x-ale_ra)/ale_L, up to g = 1 at x = rmax.  s(t) is advanced
// with the velocity ale_sdot = -(d Theta/dt)/Theta' evaluated at r_max, which holds
// Theta(r_max) fixed against the infalling gamma_rr = 1 grid.
// The taper MUST be C^2 at x = rmax (g''(rmax) = 0): the second-derivative stencil at
// i = N-1 uses R_xx[N-1] = 0, and the KO edge pin of Theta (apply_KO_scaled) fixes
// Theta_xx = 0 there.  A C^1 taper makes both inconsistent with the map and produces an
// h-independent boundary layer.  ale_build_profile checks this.
double ale_s    = 0.0;       // mesh amplitude s(t)
double ale_sdot = 0.0;       // ds/dt (frozen for the step)
double ale_ra   = 0.0;       // inner edge of the taper
double ale_L    = 1.0;       // taper width rmax - ale_ra
vector<double> ale_g, ale_gp, ale_gpp;   // g, g', g'' on the grid x_i
vector<double> Rx, Rxx, Rt;              // R_x, R_xx, R_t on the grid
vector<double> rgrid;                    // physical radius r_i = R(x_i, t)

// =====================================================================================
//  MOVING REFINEMENT
// =====================================================================================
// The computational coordinate x is first reparametrised, xi = X(x,t), to cluster points
// around a moving centre x_c(t); the taper then acts on xi:
//     r = G(xi,t) = xi + g(xi) s(t),
//     X_x = 1 - a [ b(u) - b(u/lambda)/lambda ],   u = (x - x_c)/w,   b(u) = (1-u^2)^4,
//     xi  = x - a w [ B(u) - B(u/lambda) ],         B(u) = int_{-1}^u b,  B(1) = 256/315.
// The narrow bump compresses (X_x = 1/F at u = 0 for a = (1-1/F)/(1-1/lambda)), the wide
// one (half-width lambda w, lambda = 4) re-stretches by at most 1 + a/lambda; their
// integrals cancel, so xi = x EXACTLY for |u| > lambda: the origin, the taper and the
// boundary velocity are untouched.  X_x >= 1/F everywhere, so the map cannot fold.
// Chain rule: R_x = G_xi X_x, R_xx = G_xixi X_x^2 + G_xi X_xx, R_t = G_t + G_xi X_t with
//     X_t = a xdot_c [b(u) - b(u/lambda)/lambda] - adot w [B(u) - B(u/lambda)].
// compute_RHS only reads (rgrid, Rx, Rxx, Rt), so the PDEs are untouched.  The step is
// dt = CFL dx min(R_x): the cost grows by F, not F^2.
//
// Up to two zones; their displacements ADD, which is exact while their fine cores stay
// disjoint (enforced in ref_control_zone).
//   Zone A  (--ref_F, --ref_ton, --ref_w): follows the outgoing front, i.e. the peak of
//           |P| at Theta > REF_THETA_MIN.
//   Zone B  (--refB_F, --refB_ton, --refB_w): follows an inner feature, the |P| peak
//           selected at acquisition in REF_THETAB_LO < Theta < REF_THETA_MIN, then
//           followed within +-w wherever it goes.  With --refB_x0 > 0 the zone is instead
//           switched on at refB_ton at the FIXED centre x0 (no search), held there until
//           --refB_ttrack, and only then tracks the peak.  This is needed when the feature
//           forms abruptly where no |P| peak sits beforehand (the inner shell of the
//           white-hole runs forms at t ~ 29.5, r ~ 9.96).
// The centre is moved by a PI controller (zero lag for a front at constant speed,
// |xdot_c| <= REF_VMAX); the strength a ramps in over REF_T_RAMP after acquisition and
// ramps out on release (centre leaving the allowed range, peak lost or decayed below 5%
// of its value at acquisition, or a clash with the other zone's core).
// The ICN corrector uses the mesh at t+dt (rgrid, Rx, Rxx, Rt of the advanced mesh) with
// the mesh velocities frozen for the step: freezing the geometry at t is a first-order
// time error when the bump moves by O(1) of its width.
double ref_F    = 0.0;       // zone A factor F at the centre; <= 1 => off
double ref_ton  = 50.0;      // zone A: time from which the tracker may acquire the front
double ref_w    = 5.0;       // zone A: half-width of the compressing bump (x units)
double refB_F   = 0.0;       // zone B factor; <= 1 => off
double refB_ton = 20.0;      // zone B: acquisition (or pre-positioning) time
double refB_w   = 2.0;       // zone B: half-width
double refB_x0  = 0.0;       // zone B: pre-positioned centre (0 = search for the peak)
double refB_ttrack = 0.0;    // zone B: hold at x0 until this time, then track
static const double REF_LAMBDA    = 4.0;   // width ratio of the re-stretching bump
static const double REF_THETA_MIN = 15.0;  // zone A tracks only at Theta > this
static const double REF_THETAB_LO = 1.0;   // zone B selects its feature at Theta > this
static const double REF_T_RAMP    = 5.0;   // ramp time of the refinement strength a
static const double REF_TAU       = 0.5;   // PI controller: proportional time
static const double REF_TAU_I     = 1.0;   //                integral time
static const double REF_VMAX      = 2.0;   // |xdot_c| clamp
static const double REF_P_MIN     = 1e-6;  // |P| peak below this is noise: no acquisition

struct RefZone {
    // configuration (set in set_initial_data)
    double F      = 0.0;     // refinement factor at the centre (<= 1: zone off)
    double ton    = 0.0;     // acquisition (or pre-positioning) time
    double w      = 1.0;     // half-width of the compressing bump (x units)
    double th_lo  = 0.0;     // Theta lower bound of the acquisition window
    bool   inner  = false;   // zone B behaviour (window Theta < REF_THETA_MIN)
    const char* name = "";   // "A" or "B" (log messages)
    double x0     = 0.0;     // pre-positioned centre (0: search for the peak)
    double ttrack = 0.0;     // pre-positioned zone: hold at x0 until this time
    // state
    double xc      = 0.0;    // centre x_c
    double xc_dot  = 0.0;    // its velocity (frozen for the step)
    double vi      = 0.0;    // integral term of the PI controller
    double a       = 0.0;    // strength
    double a_dot   = 0.0;    // its rate (frozen for the step)
    double atarget = 0.0;    // strength the ramp is heading for
    double p_acq   = 0.0;    // |P| peak at acquisition (reference for the 5% release)
    bool acquired = false, released = false;
    bool warned   = false;   // range warning already printed
};
RefZone RZ[2];               // RZ[0] = zone A, RZ[1] = zone B
inline bool ref_on() { return ref_F > 1.0 || refB_F > 1.0; }

// bump profile b(u) = (1-u^2)^4 (C^3 at u = +-1), its derivative and its integral B(u)
inline double ref_b(double u)  { if (u <= -1.0 || u >= 1.0) return 0.0; double v = 1.0-u*u; return v*v*v*v; }
inline double ref_bp(double u) { if (u <= -1.0 || u >= 1.0) return 0.0; double v = 1.0-u*u; return -8.0*u*v*v*v; }
inline double ref_B(double u) {
    if (u <= -1.0) return 0.0;
    if (u >=  1.0) return 256.0/315.0;
    double u2 = u*u;
    return u*(1.0 + u2*(-4.0/3.0 + u2*(1.2 + u2*(-4.0/7.0 + u2/9.0)))) + 128.0/315.0;
}

// C^2 quintic taper: g = 0 for x <= x0, u^3(10-15u+6u^2) for 0<u<1, 1 for x >= x0+L.
inline void quintic_taper(double x, double x0, double L, double& g, double& gp, double& gpp) {
    if (x <= x0) { g = 0.0; gp = 0.0; gpp = 0.0; return; }
    double u = (x - x0)/L;
    if (u >= 1.0) { g = 1.0; gp = 0.0; gpp = 0.0; return; }
    g   = u*u*u*(10.0-15.0*u+6.0*u*u);
    gp  = 30.0*u*u*(1.0-u)*(1.0-u)/L;
    gpp = 60.0*u*(1.0-u)*(1.0-2.0*u)/(L*L);
}

// Tabulate the taper g, g', g'' on the computational grid and check that it is C^2 at the
// outer edge (|g''| L^2 over the last four points must be small; it is 300/(N-1) for the
// quintic and ~6 for a C^1 taper -- the test discriminates for N >~ 286).
void ale_build_profile(int n) {
    ale_g.assign(n,0.0); ale_gp.assign(n,0.0); ale_gpp.assign(n,0.0);
    for (int i=0;i<n;++i)
        quintic_taper(i*dx, ale_ra, ale_L, ale_g[i], ale_gp[i], ale_gpp[i]);
    if (ale_on && n > 4) {
        double gpp_edge = 0.0;
        for (int i = n-4; i < n; ++i)
            gpp_edge = max(gpp_edge, fabs(ale_gpp[i]) * ale_L * ale_L);
        if (!(gpp_edge < 1.0)) {
            fprintf(stderr, "FATAL: ALE taper C^2 check failed: max |g''|*L^2 = %.3g over "
                    "the last 4 points; expected %.3g = 300/(N-1) for the C^2 quintic. "
                    "Either the taper is not C^2 at the outer edge or N < ~286.\n",
                    gpp_edge, 300.0/(double)(n-1));
            exit(EXIT_FAILURE);
        }
    }
}

// Rebuild rgrid, Rx, Rxx, Rt from the current mesh state (ale_s, ale_sdot and, with the
// refinement, the zone centres/strengths and their velocities).
void ale_update(int n) {
    if (ale_on==0){
        for(int i=0;i<n;++i){ rgrid[i]=i*dx; Rx[i]=1.0; Rxx[i]=0.0; Rt[i]=0.0; }
        return;
    }
    for(int i=0;i<n;++i){
        if (ref_on()) {
            // composed map r = G(X(x,t),t)
            const double x = i*dx, lam = REF_LAMBDA;
            double xi = x, Xx = 1.0, Xxx = 0.0, Xt = 0.0;
            for (int z = 0; z < 2; ++z) {
                const RefZone& Z = RZ[z];
                if (!(Z.F > 1.0)) continue;
                const double w = Z.w, u = (x - Z.xc)/w, v = u/lam;
                xi  = xi - Z.a*w*(ref_B(u) - ref_B(v));
                Xx  = Xx - Z.a*(ref_b(u) - ref_b(v)/lam);
                Xxx = Xxx - (Z.a/w)*(ref_bp(u) - ref_bp(v)/(lam*lam));
                Xt  = Xt + Z.a*Z.xc_dot*(ref_b(u) - ref_b(v)/lam)
                         - Z.a_dot*w*(ref_B(u) - ref_B(v));
            }
            double go, gop, gopp;
            quintic_taper(xi, ale_ra, ale_L, go, gop, gopp);
            const double G   = xi + go*ale_s;
            const double Gx  = 1.0 + gop*ale_s;
            const double Gxx = gopp*ale_s;
            const double Gt  = go*ale_sdot;
            rgrid[i] = G;
            Rx[i]    = Gx*Xx;
            Rxx[i]   = Gxx*Xx*Xx + Gx*Xxx;
            Rt[i]    = Gt + Gx*Xt;
        } else {
            rgrid[i]= i*dx + ale_g[i]*ale_s;
            Rx[i]  = 1.0 + ale_gp[i]*ale_s;
            Rxx[i] = ale_gpp[i]*ale_s;
            Rt[i]  = ale_g[i]*ale_sdot;
        }
        // A folded (Rx <= 0) or non-finite map would otherwise show up as a silent
        // blow-up, since every d/dr divides by Rx.
        if (!(Rx[i] > 0.0)) {
            if (!isfinite(Rx[i]))
                fprintf(stderr, "FATAL: ALE map is non-finite at i=%d (Rx=%.6g, s=%.6g, "
                        "sdot=%.6g). A NaN/Inf reached the grid velocity; look upstream "
                        "at Theta/NN/NR near the outer boundary.\n",
                        i, Rx[i], ale_s, ale_sdot);
            else
                fprintf(stderr, "FATAL: ALE mesh folded at i=%d (Rx=%.6g, s=%.6g, "
                        "g'=%.6g).\n", i, Rx[i], ale_s, ale_gp[i]);
            exit(EXIT_FAILURE);
        }
    }
}

// Controller of one refinement zone Z (O = the other zone).  Uses the state at the start
// of the step and sets Z.xc_dot and Z.a_dot, which are then frozen for the step.
void ref_control_zone(RefZone& Z, const RefZone& O, double t, double dt_step, int n) {
    Z.xc_dot = 0.0; Z.a_dot = 0.0;
    if (t < Z.ton) return;
    const double w = Z.w, lam = REF_LAMBDA;
    // allowed range of the centre: the support (lambda w) stays clear of the origin and
    // of the last 20 points
    const double x_lo = 1.0 + lam*w;
    const double x_hi = dx*(double)(n-1) - 20.0*dx - lam*w;
    const double th_lo = Z.th_lo;
    const double th_hi = Z.inner ? REF_THETA_MIN : Theta[n-1] - 10.0;
    // a zone still ramping out keeps compressing the mesh, so it counts until a = 0
    const bool other = (O.F > 1.0) && O.acquired && (!O.released || O.a > 0.0);
    if (Z.x0 > 0.0) {                      // pre-positioned zone: switch on at x0, hold until ttrack
        if (!Z.acquired) {
            if (Z.x0 < x_lo || Z.x0 > x_hi) {
                fprintf(stderr, "FATAL: refinement %s x0=%.4g outside the allowed range [%.4g, %.4g] "
                        "(support must stay inside the grid, clear of the origin)\n",
                        Z.name, Z.x0, x_lo, x_hi);
                exit(EXIT_FAILURE);
            }
            if (other && fabs(Z.x0 - O.xc) < w + O.w + 1.0) {
                fprintf(stderr, "FATAL: refinement %s x0=%.4g overlaps the core of zone %s at x_c=%.4g\n",
                        Z.name, Z.x0, O.name, O.xc);
                exit(EXIT_FAILURE);
            }
            Z.xc = Z.x0; Z.vi = 0.0; Z.acquired = true; Z.p_acq = 0.0;
            Z.atarget = (1.0 - 1.0/Z.F)/(1.0 - 1.0/lam);
            int ix = (int)(Z.xc/dx); if (ix < 0) ix = 0; if (ix > n-1) ix = n-1;
            fprintf(stderr, "refinement %s: pre-positioned at x_c=%.4f (r=%.4f) at t=%.4f, holding "
                    "until t=%.4g (F=%.3g, w=%.3g)\n", Z.name, Z.xc, rgrid[ix], t, Z.ttrack, Z.F, w);
        }
        if (t < Z.ttrack) {
            Z.xc_dot = 0.0;
            if (!Z.released && other && fabs(Z.xc - O.xc) < w + O.w) {   // the other zone moved onto us
                Z.released = true; Z.atarget = 0.0;
                fprintf(stderr, "refinement %s: releasing while holding at t=%.4f (core clash with %s "
                        "at x_c=%.4g)\n", Z.name, t, O.name, O.xc);
            }
            const double rate = (1.0 - 1.0/Z.F)/(1.0 - 1.0/lam)/REF_T_RAMP;
            double da = Z.atarget - Z.a;
            if (da >  rate*dt_step) da =  rate*dt_step;
            if (da < -rate*dt_step) da = -rate*dt_step;
            Z.a_dot = da/dt_step;
            return;
        }
        if (Z.p_acq == 0.0) {              // first tracking step: amplitude reference for release
            double pb0 = 0.0;
            for (int i = 1; i < n-1; ++i)
                if (fabs(dx*(double)i - Z.xc) <= w) pb0 = max(pb0, fabs(P[i]));
            Z.p_acq = pb0;
            fprintf(stderr, "refinement %s: tracking from t=%.4f (|P| peak %.4g)\n", Z.name, t, pb0);
        }
    }
    if (!Z.warned && !(x_hi - x_lo > 2.0*w)) {
        fprintf(stderr, "refinement %s: WARNING allowed range [%.4g, %.4g] too narrow for w=%.3g, "
                "lambda=%.3g; it will never acquire (raise rmax or lower w)\n", Z.name, x_lo, x_hi, w, lam);
        Z.warned = true;
    }
    // peak of |P| in this zone's Theta window: global search at acquisition, then within
    // +-w of x_c.  The Theta window only SELECTS the feature at acquisition; afterwards the
    // zone follows its peak wherever it goes (zone A keeps Theta > REF_THETA_MIN).
    int ib = -1; double pb = -1.0;
    for (int i = 1; i < n-1; ++i) {
        const double x = dx*(double)i;
        if (Z.acquired) {
            if (fabs(x - Z.xc) > w || !(Theta[i] > th_lo && Theta[i] < Theta[n-1] - 10.0)) continue;
        } else {
            if (!(Theta[i] > th_lo && Theta[i] < th_hi)) continue;
            if (other && fabs(x - O.xc) < w + O.w + 1.0) continue;   // never pick the other zone's feature
        }
        if (fabs(P[i]) > pb) { pb = fabs(P[i]); ib = i; }
    }
    if (!Z.acquired) {
        // wait until the peak is inside the allowed range with a margin of w and (if the
        // other zone is active) the two fine cores are disjoint
        if (ib < 0 || pb < REF_P_MIN || dx*(double)ib < x_lo + w || dx*(double)ib > x_hi - w) return;
        if (other && fabs(dx*(double)ib - O.xc) < w + O.w + 1.0) return;
        Z.xc = dx*(double)ib; Z.vi = 0.0; Z.acquired = true; Z.p_acq = pb;
        Z.atarget = (1.0 - 1.0/Z.F)/(1.0 - 1.0/lam);   // X_x = 1/F at the centre
        fprintf(stderr, "refinement %s: acquired the front at t=%.4f  x_c=%.4f  r=%.4f  "
                "Theta=%.4f  |P|=%.4g  (F=%.3g, w=%.3g, lambda=%.3g)\n",
                Z.name, t, Z.xc, rgrid[ib], Theta[ib], pb, Z.F, w, lam);
    }
    if (ib > 0 && ib < n-1 && !Z.released) {
        // sub-grid peak position from a parabola through |P| at ib-1, ib, ib+1
        const double ym = fabs(P[ib-1]), y0 = fabs(P[ib]), yp = fabs(P[ib+1]);
        const double den = ym - 2.0*y0 + yp;
        double off = (den < 0.0) ? 0.5*(ym - yp)/den : 0.0;
        if (off > 0.5) off = 0.5; else if (off < -0.5) off = -0.5;
        const double e = dx*((double)ib + off) - Z.xc;
        Z.vi    += dt_step*e/(REF_TAU*REF_TAU_I);
        if (Z.vi >  REF_VMAX) Z.vi =  REF_VMAX;          // anti-windup
        if (Z.vi < -REF_VMAX) Z.vi = -REF_VMAX;
        Z.xc_dot = e/REF_TAU + Z.vi;
        if (Z.xc_dot >  REF_VMAX) Z.xc_dot =  REF_VMAX;
        if (Z.xc_dot < -REF_VMAX) Z.xc_dot = -REF_VMAX;
    }
    const double xn = Z.xc + dt_step*Z.xc_dot;
    // on a clash, release the zone that is closing in (not simply the one checked first)
    const bool clash = other && fabs(xn - O.xc) < w + O.w
                       && (Z.xc - O.xc)*(Z.xc_dot + O.xc_dot) <= 0.0;
    if (!Z.released && (xn < x_lo || xn > x_hi || ib < 0 || pb < 0.05*Z.p_acq || clash)) {
        Z.released = true; Z.atarget = 0.0;
        fprintf(stderr, "refinement %s: releasing at t=%.4f (x_c=%.4f outside [%.4f, %.4f], or "
                "front lost, or |P|=%.3g < 5%% of %.3g at acquisition, or core clash=%d); "
                "ramping out and freezing x_c\n", Z.name, t, Z.xc, x_lo, x_hi, pb, Z.p_acq, (int)clash);
    }
    if (Z.released) { Z.xc_dot = 0.0; Z.vi = 0.0; }
    const double rate = (1.0 - 1.0/Z.F)/(1.0 - 1.0/lam)/REF_T_RAMP;
    double da = Z.atarget - Z.a;
    if (da >  rate*dt_step) da =  rate*dt_step;
    if (da < -rate*dt_step) da = -rate*dt_step;
    Z.a_dot = da/dt_step;
}

void ref_control(double t, double dt_step, int n) {
    for (int z = 0; z < 2; ++z)
        if (RZ[z].F > 1.0) ref_control_zone(RZ[z], RZ[1-z], t, dt_step, n);
}

// =====================================================================================
//  TIME STEP
// =====================================================================================
// dt = CFL dx, or CFL dx min(R_x) with the refinement (without it R_x = 1 + g' s can be
// slightly below 1 while s < 0, as early in ingoing runs; the fixed-grid value is kept
// there, as in all paper runs).
double compute_dt_CFL_from_state() {
    if (ref_on()) {
        double m = 1.0;
        for (size_t i = 0; i < Rx.size(); ++i) m = min(m, Rx[i]);
        return CFL * dx * m;
    }
    return CFL * dx;
}

// =====================================================================================
//  FINITE DIFFERENCES
// =====================================================================================
// d1_*, d2_*: d/dx, d^2/dx^2 on the uniform computational grid.  "even"/"odd" refer to
// the parity of the field at r = 0 (used at i = 0); i = n-1 uses one-sided stencils.
inline double d1_even(const vector<double>& f, int i) {
    int n = (int)f.size();
    if (i == 0)     return 0.0;
    if (i == n - 1) return -(-3.0*f[i] + 4.0*f[i-1] - f[i-2]) / (2.0*dx);
    return (f[i+1] - f[i-1]) / (2.0*dx);
}

inline double d1_odd(const vector<double>& f, int i) {
    int n = (int)f.size();
    if (i == 0)     return f[1] / dx;
    if (i == n - 1) return -(-3.0*f[i] + 4.0*f[i-1] - f[i-2]) / (2.0*dx);
    return (f[i+1] - f[i-1]) / (2.0*dx);
}

inline double d2_even(const vector<double>& f, int i) {
    int n = (int)f.size();
    if (i == 0)     return 2.0*(f[1] - f[0]) / (dx*dx);
    if (i == n - 1) return (2.0*f[i] - 5.0*f[i-1] + 4.0*f[i-2] - f[i-3]) / (dx*dx);
    return (f[i+1] - 2.0*f[i] + f[i-1]) / (dx*dx);
}

inline double d2_odd(const vector<double>& f, int i) {
    int n = (int)f.size();
    if (i == 0)     return 0.0;
    if (i == n - 1) return (2.0*f[i] - 5.0*f[i-1] + 4.0*f[i-2] - f[i-3]) / (dx*dx);
    return (f[i+1] - 2.0*f[i] + f[i-1]) / (dx*dx);
}

// Physical derivatives on the moving mesh:
//   d_r f = d_x f / R_x,   d_r^2 f = d_x^2 f / R_x^2 - R_xx d_x f / R_x^3.
inline double dr1_even(const vector<double>& f, int i) { return d1_even(f,i)/Rx[i]; }
inline double dr1_odd (const vector<double>& f, int i) { return d1_odd (f,i)/Rx[i]; }
inline double dr2_even(const vector<double>& f, int i) {
    double d1 = d1_even(f,i), d2 = d2_even(f,i), rx = Rx[i];
    return d2/(rx*rx) - Rxx[i]*d1/(rx*rx*rx);
}
inline double dr2_odd(const vector<double>& f, int i) {
    double d1 = d1_odd(f,i), d2 = d2_odd(f,i), rx = Rx[i];
    return d2/(rx*rx) - Rxx[i]*d1/(rx*rx*rx);
}

// =====================================================================================
//  KREISS-OLIGER DISSIPATION
// =====================================================================================
// f -= sigma D4 f with sigma = eps (dt/dx)/16, in index space (a fixed rate per unit time).
// Origin: ghost points from the parity of f.  Outer edge:
//  * wave fields (pin_edge = false): ONE-SIDED interior 4th differences -- point n-2 gets
//    the window centred at n-4, n-1 the window centred at n-3.  They annihilate the
//    grid-scale sawtooth and constrain smooth data only at O(h^4).  (Linear ghost
//    extrapolation would pin f_xx = 0 at the edge and reflect ~1/3 of an outgoing
//    scalar pulse.)
//  * Theta (pin_edge = true): ghosts extrapolated linearly in index space.  Its fixed
//    point is linear extrapolation of Theta (Theta_xx = 0 over the last four points), held
//    by a contractive relaxation.  This is the boundary datum of Theta's slowly ENTERING
//    edge mode on the moving mesh (its eigenvalue at r_max is small and positive); without
//    it the edge is under-determined and the ALE feedback runs away.  Its error is
//    O(Theta_rr h) and converges away; it deflates ThetaC at the last three points, which
//    should not be read.
void apply_KO_scaled(vector<double>& f, double eps, double dt_local,
                     bool is_odd_at_origin, bool pin_edge = false) {
    int n  = (int)f.size();
    if (n < 6)      return;   // the outer closure reaches back to n-6
    if (eps == 0.0) return;

    vector<double> old = f;
    double sigma = eps * (dt_local / dx) * (1.0 / 16.0);

    for (int i = 2; i <= n - 3; ++i) {
        double D4 = old[i-2] - 4.0*old[i-1] + 6.0*old[i]
                  - 4.0*old[i+1] + old[i+2];
        f[i] -= sigma * D4;
    }

    {
        double fm1 = is_odd_at_origin ? -old[1] : old[1];
        double fm2 = is_odd_at_origin ? -old[2] : old[2];
        double D4  = fm2 - 4.0*fm1 + 6.0*old[0] - 4.0*old[1] + old[2];
        f[0] -= sigma * D4;
    }
    {
        double fm1 = is_odd_at_origin ? -old[1] : old[1];
        double D4  = fm1 - 4.0*old[0] + 6.0*old[1] - 4.0*old[2] + old[3];
        f[1] -= sigma * D4;
    }

    if (pin_edge) {
        double fp1 = 2.0*old[n-1] - old[n-2], fp2 = 3.0*old[n-1] - 2.0*old[n-2];
        f[n-2] -= sigma * (old[n-4] - 4.0*old[n-3] + 6.0*old[n-2] - 4.0*old[n-1] + fp1);
        f[n-1] -= sigma * (old[n-3] - 4.0*old[n-2] + 6.0*old[n-1] - 4.0*fp1 + fp2);
    } else {
        f[n-2] -= sigma * (old[n-6] - 4.0*old[n-5] + 6.0*old[n-4] - 4.0*old[n-3] + old[n-2]);
        f[n-1] -= sigma * (old[n-5] - 4.0*old[n-4] + 6.0*old[n-3] - 4.0*old[n-2] + old[n-1]);
    }
}

// =====================================================================================
//  APPARENT-HORIZON FINDER (diagnostic only)
// =====================================================================================
// theta_+ = Theta' + (1/6)(3Q - 2K) Theta  (eq:expansions; outgoing null expansion up to a positive
// factor): < 0 inside a trapped region.  Using theta_+ rather than theta_+ theta_- avoids
// latching onto anti-trapped (theta_- = 0) surfaces.
static double theta_plus(int i,
                         const vector<double>& Theta_s,
                         const vector<double>& trK_s,
                         const vector<double>& Q_s) {
    double th    = Theta_s[i];
    double th_r  = dr1_odd(Theta_s, i);
    double combo = 3.0*Q_s[i] - 2.0*trK_s[i];
    return th_r + (1.0/6.0)*combo*th;
}

// Locate the outermost - -> + crossing of theta_+ at r <= R_AH_SEARCH_MAX and report it
// on stdout when it is at or outside the previous one: the index k-1 of the cell that
// contains the zero, the radius of the zero of the linear interpolant, and the time of
// the state (end of the step).  Diagnostic only; the analysis scripts use the dumps.
static void horizon_finder_step(const vector<double>& Theta_s,
                                const vector<double>& trK_s,
                                const vector<double>& Q_s,
                                int Nx, double t) {
    vector<double> hf(Nx, 0.0);
    for (int i = 0; i < Nx; ++i)
        hf[i] = theta_plus(i, Theta_s, trK_s, Q_s);

    int k_max = Nx - 1;
    while (k_max > 1 && rgrid[k_max] > R_AH_SEARCH_MAX) --k_max;

    int i_ah = 0; double r_root = 0.0;   // cell [i_ah, i_ah+1] that contains the zero
    for (int k = k_max; k > 1; --k) {
        if (hf[k] > 0.0 && hf[k-1] < 0.0) {
            // zero of the linear interpolant between r_{k-1} and r_k
            r_root = rgrid[k-1] - hf[k-1]*(rgrid[k] - rgrid[k-1])/(hf[k] - hf[k-1]);
            i_ah = k - 1;
            break;
        }
    }
    if (i_ah > 0 && i_ah >= g_horizon) {
        g_horizon = i_ah;
        printf("Horizon detected/updated: i=%d  r=%.6f  t=%.6f\n", i_ah, r_root, t);
    }
}

// =====================================================================================
//  RECOVERED FIELDS: Q (momentum constraint), N (lapse), beta (shift)
// =====================================================================================
// Q' = -3 Q Theta'/Theta + [ -P psi' + (2+c13+3c2)/3 K' ]/(1-c13) (eq:radial), outward
// from Q(0) = 0.  The ODE is singular at r = 0 (-3Q/r): the first panel uses the regular
// series Q = q2 r^2, 5 q2 = S1  =>  Q(r_1) = rhs[1] r_1/5;  Heun with the physical step
// afterwards.
void compute_Q(const vector<double>& Psi_s, const vector<double>& P_s,
               const vector<double>& Theta_s, const vector<double>& trK_s,
               vector<double>& Q_s) {
    int n = (int)Psi_s.size();
    Q_s.assign(n, 0.0);

    const double denom_factor = 1.0 - c13;   // non-zero: couplings checked in main
    const double factor = (1.0/3.0) * den_K;

    vector<double> rhs_nonQ(n, 0.0);
    for (int i = 0; i < n; ++i) {
        double dPsi = dr1_even(Psi_s, i);
        double dtrK = dr1_even(trK_s, i);
        rhs_nonQ[i] = (-P_s[i] * dPsi + factor * dtrK) / denom_factor;
    }
    auto hstep = [&](int i) -> double { return rgrid[i+1] - rgrid[i]; };

    auto f_rhs = [&](int i, double Qval) -> double {   // called for i >= 1 only
        double th     = Theta_s[i];
        double dTheta = dr1_odd(Theta_s, i);
        double safe_th = (fabs(th) < 1e-14)
                         ? (th >= 0.0 ? 1e-14 : -1e-14) : th;
        return -3.0 * Qval / safe_th * dTheta + rhs_nonQ[i];
    };

    Q_s[0] = 0.0;
    if (n > 1) Q_s[1] = rhs_nonQ[1] * hstep(0) / 5.0;
    for (int i = 1; i < n - 1; ++i) {
        double h      = hstep(i);
        double k1     = f_rhs(i, Q_s[i]);
        double Q_pred = Q_s[i] + h * k1;
        double k2     = f_rhs(i + 1, Q_pred);
        Q_s[i+1]      = Q_s[i] + 0.5 * h * (k1 + k2);
    }
    Q_s[0] = 0.0;
}

// Lapse from (ln N)' = a_r (GEJ eq. 28; a factor 1/2 here would make free evolution of a_r
// inconsistent).  Trapezoid outward from the centre, then shifted so that ln N(r_max) = 0
// (the equation is linear in ln N, so this is the exact solution of the boundary problem).
// N(0) is then reset from the even-parity extrapolation (4 N_1 - N_2)/3.  A non-finite
// lapse stops the run.
void compute_lapse(const vector<double>& Acr_s, vector<double>& NN_s) {
    int n = (int)Acr_s.size();
    NN_s.assign(n, 1.0);
    if (n == 0) return;

    vector<double> lnNN(n, 0.0);
    for (int i = 0; i < n - 1; ++i)
        lnNN[i+1] = lnNN[i] + 0.5 * (rgrid[i+1] - rgrid[i]) * (Acr_s[i] + Acr_s[i+1]);
    const double lnNN_edge = lnNN[n-1];
    for (int i = 0; i < n; ++i) {
        NN_s[i] = exp(lnNN[i] - lnNN_edge);
        if (!isfinite(NN_s[i])) {
            fprintf(stderr, "FATAL: non-finite lapse at i=%d, r=%.6f (step %d); stopping.\n",
                    i, rgrid[i], g_step);
            exit(EXIT_FAILURE);
        }
    }

    if (n > 2)
        NN_s[0] = (4.0 * NN_s[1] - NN_s[2]) / 3.0;
}

// Shift: beta' = N (Q + K/3), trapezoid outward from beta(0) = 0 (beta is odd).
void compute_shift(const vector<double>& Q_s, const vector<double>& trK_s,
                   const vector<double>& NN_s, vector<double>& NR_s) {
    int n = (int)Q_s.size();
    NR_s.assign(n, 0.0);
    if (n == 0) return;

    auto rhs = [&](int i) -> double {
        return NN_s[i] * (Q_s[i] + trK_s[i] / 3.0);
    };

    NR_s[0] = 0.0;
    for (int i = 0; i < n - 1; ++i) {
        double h  = rgrid[i+1] - rgrid[i];
        double k1 = rhs(i);
        double k2 = rhs(i + 1);
        NR_s[i+1] = NR_s[i] + 0.5 * h * (k1 + k2);
    }
    NR_s[0] = 0.0;
}

// d_r Theta (written only to the auxiliary dumps)
void compute_DTheta(const vector<double>& Theta_s, vector<double>& DTheta_s) {
    int n = (int)Theta_s.size();
    DTheta_s.assign(n, 0.0);
    for (int i = 0; i < n; ++i)
        DTheta_s[i] = dr1_odd(Theta_s, i);
}

// =====================================================================================
//  RIGHT-HAND SIDES
// =====================================================================================
// Interior: the five evolution equations of the file header, with mesh-relative advection
// NR + R_t.  At r = 0 the regular limits Theta' psi'/Theta -> psi''(0) and
// Theta' a_r/Theta -> a_r'(0) are used.
//
// Outer boundary (i = n-1): replace the time derivative of the INCOMING characteristic
// field of each radiative pair by boundary data:
//  * spin-0 pair u = (K, a_r): M = [[bt, N aN],[N bN, bt]] (bt = NR + R_t, the shift
//    relative to the grid), speeds bt +- N v0; right eigenvector e+ = (aN/v0, 1), left
//    eigenvectors th+- = (+-bN/v0, 1), th+.e+ = 2, th-.e+ = 0.  With R the interior RHS,
//        d_t u = R - e+ (th+.R - G)/(th+.e+),   G = bt d_r w_in - N bN K/Theta,
//    so d_t w_in = G for w_in = th+.u = (bN/v0) K + a_r, and d_t w_out = th-.R is
//    unchanged.  G is d_t w_in at fixed r of the outgoing linear solution K = f/r,
//    a_r = -(bN/v0) f/r - bN F/r^2 (F' = f), i.e. (d_t - beta d_r) w_in = -N bN K/Theta
//    in the aether frame; d_r w_in uses the one-sided edge derivative of the scalar
//    condition.
//  * scalar pair: the SPHERICAL outgoing condition, imposed at the RHS level,
//        d_t [ P + psi' + (Theta'/Theta + Q/2 - K/3) psi ] = 0,
//    i.e. Theta psi constant along outgoing rays; d_t psi' is taken as d_r(Psi_rhs).
//  * Psi and Theta (speed bt, non-radiative) keep their interior RHS.
// Each condition is applied only if its mode ENTERS (eigenvalue > 0); the gate is set
// once per step in ICN_step.  On the fixed grid the shift at r_max eventually falls below
// -N, both modes leave and the edge is evolved on the interior (one-sided) RHS.
void compute_RHS(const vector<double>& Psi_field,
                 const vector<double>& P_field,
                 const vector<double>& trK_field,
                 const vector<double>& Acr_field,
                 const vector<double>& Theta_field,
                 const vector<double>& Q_field,
                 const vector<double>& NN_field,
                 const vector<double>& NR_field,
                 vector<double>& Psi_rhs,
                 vector<double>& P_rhs,
                 vector<double>& trK_rhs,
                 vector<double>& Acr_rhs,
                 vector<double>& Theta_rhs) {
    int n = (int)Psi_field.size();

    Psi_rhs.assign(n, 0.0);
    P_rhs.assign(n, 0.0);
    trK_rhs.assign(n, 0.0);
    Acr_rhs.assign(n, 0.0);
    Theta_rhs.assign(n, 0.0);

    // coupling coefficients (den_K, den_a non-zero: couplings checked in main)
    const double denom  = den_K;    // 2 + c13 + 3 c2
    const double denom1 = den_a;    // c14 (1 - c13)
    const double numer1 = c13;
    const double numer2 = c123;

    for (int i = 0; i < n; ++i) {
        double th     = Theta_field[i];
        double psi_r  = dr1_even(Psi_field, i);
        double psi_rr = dr2_even(Psi_field, i);
        double p      = P_field[i];
        double k      = trK_field[i];
        double acr    = Acr_field[i];
        double q      = Q_field[i];
        double NN_i   = NN_field[i];
        double bt_i   = NR_field[i] + Rt[i];   // bt: shift + grid velocity (mesh advection)
        double dTheta = dr1_odd(Theta_field, i);

        // Theta' psi'/Theta  (-> psi''(0) at the origin)
        double ThpPsip_over_Th = 0.0;
        if (i != 0 && fabs(th) > 1e-14)
            ThpPsip_over_Th = dTheta * psi_r / th;
        else if (i == 0)
            ThpPsip_over_Th = psi_rr;

        // Theta' a_r/Theta  (-> a_r'(0) at the origin)
        double ThpAr_over_Th = 0.0;
        if (i == 0)
            ThpAr_over_Th = dr1_odd(Acr_field, i);
        else if (fabs(th) > 1e-14)
            ThpAr_over_Th = dTheta * acr / th;

        Psi_rhs[i] = bt_i * psi_r + NN_i * p;

        P_rhs[i] = bt_i * dr1_even(P_field, i)
                 + NN_i * (k*p + acr*psi_r + psi_rr + 2.0*ThpPsip_over_Th);

        double bracket = (c14 - 2.0) *
                         (dr1_odd(Acr_field, i) + 2.0*ThpAr_over_Th + acr*acr);

        trK_rhs[i] = bt_i * dr1_even(trK_field, i)
                   + NN_i * k*k / 3.0
                   + NN_i * (bracket + 2.0*p*p
                             + 3.0*(1.0 - c13)*q*q) / denom;

        Acr_rhs[i] = bt_i * dr1_odd(Acr_field, i)
                   + NN_i * ((2.0*k/3.0 - q) * acr
                   + (numer1 / denom1) * p * psi_r
                   - (numer2 / denom1) * dr1_even(trK_field, i));

        Theta_rhs[i] = bt_i * dr1_odd(Theta_field, i)
                     + NN_i * Theta_field[i] * (q/2.0 - k/3.0);
    }

    // ---- outer boundary: characteristic conditions on the entering modes ----
    {
        const int i = n - 1;
        auto dr_edge = [&](const vector<double>& f) -> double {   // one-sided d_r at the edge
            return (3.0*f[i] - 4.0*f[i-1] + f[i-2]) / (2.0*dx) / Rx[i];
        };
        double Th_b = Theta_field[i];
        if (fabs(Th_b) < 1e-12) Th_b = (Th_b >= 0 ? 1e-12 : -1e-12);

        // spin-0 pair (trK, Acr): d_t w_in = G with w_in = th+.u (eq:bc-win),
        // G = bt d_r w_in - N bN K/Theta (outgoing linear solution at fixed r)
        const double aN = spin0_aN, bN = spin0_bN;
        const double eK = aN / spin0_v0, tK = bN / spin0_v0;  // e+ = (eK,1), th+ = (tK,1)
        // th+.e+ = aN bN/v0^2 + 1 = 2 analytically; computed, as in the formula above
        const double norm = eK*tK + 1.0;
        const double Kr = trK_rhs[i], Ar = Acr_rhs[i];
        if (gate_spin0 != 0) {
            const double bt = NR_field[i] + Rt[i];
            const double dw_in = tK*dr_edge(trK_field) + dr_edge(Acr_field);
            const double G = bt*dw_in - NN_field[i]*bN*trK_field[i]/Th_b;
            const double c = (tK*Kr + Ar - G) / norm;
            trK_rhs[i] = Kr - eK*c;
            Acr_rhs[i] = Ar - c;
        }

        // scalar pair: d_t [P + psi' + c psi] = 0 with c = Theta'/Theta + Q/2 - K/3, i.e.
        // d_t P = -d_r(Psi_rhs) - c Psi_rhs - psi d_t c, where d_t c uses d_t Theta' =
        // d_r(Theta_rhs), d_t Theta = Theta_rhs, d_t K = trK_rhs (after the spin-0
        // condition above) and d_t Q = g_Qdot_edge (backward difference between the ICN
        // stages, set in ICN_step)
        const double Thp_b = dr1_odd(Theta_field, i);
        const double c_b = Thp_b / Th_b + (Q_field[i] / 2.0 - trK_field[i] / 3.0);
        const double dPsi_rhs = dr_edge(Psi_rhs);
        const double dc_dt = dr_edge(Theta_rhs) / Th_b - Thp_b * Theta_rhs[i] / (Th_b*Th_b)
                           + g_Qdot_edge / 2.0 - trK_rhs[i] / 3.0;
        if (gate_scalar != 0)
            P_rhs[i] = -dPsi_rhs - c_b * Psi_rhs[i] - dc_dt * Psi_field[i];
    }
}

// =====================================================================================
//  ICN STEP
// =====================================================================================
// relative L2 change, the ICN convergence criterion (applied to Psi and Theta)
double rel_change_norm(const vector<double>& a, const vector<double>& b) {
    double num = 0.0, den = 0.0;
    int n = (int)a.size();
    for (int i = 0; i < n; ++i) {
        double d = a[i] - b[i];
        num += d*d;
        den += a[i]*a[i];
    }
    if (den < 1e-30) den = 1e-30;
    return sqrt(num / den);
}

// a non-finite field value ends the run (the output written so far is kept and complete)
static void check_fields_finite(double t) {
    const vector<double>* flds[5] = {&Psi, &P, &trK, &Acr, &Theta};
    const char* names[5] = {"Psi", "P", "trK", "Acr", "Theta"};
    for (int f = 0; f < 5; ++f) {
        const vector<double>& v = *flds[f];
        for (int i = 0; i < (int)v.size(); ++i) {
            if (!isfinite(v[i])) {
                // on stdout too, for logs that capture only stdout (README, "Output")
                printf("FIRST NaN: step=%d  field=%s  i=%d  r=%.6f\n", g_step, names[f], i, rgrid[i]);
                fprintf(stderr, "FATAL: non-finite %s after step %d, i=%d, r=%.6f (t=%.6f); "
                        "stopping.\n", names[f], g_step, i, rgrid[i], t);
                fflush(stdout);
                exit(EXIT_FAILURE);
            }
        }
    }
}

// One step t -> t+dt:
//  1. dt; mesh velocity (ALE: hold Theta(r_max)); refinement controller; outer-BC gate.
//  2. ICN: RHS at t; Euler predictor; the mesh is advanced to t+dt (velocities frozen);
//     up to ICN_MAX_ITERS trapezoidal corrections on that mesh, each re-solving Q, N, beta
//     from the iterate, until the relative change of Psi and Theta is below ICN_TOL.
//  3. KO dissipation, recompute Q, N, beta, d_r Theta, finiteness check, apparent-horizon
//     finder.
void ICN_step(double t) {
    g_step++;
    dt = compute_dt_CFL_from_state();

    // ---- ALE: grid velocity that holds the outer boundary at fixed areal radius ----
    // dTheta/dt|_x = NR Theta' + N Theta (Q/2 - K/3) at r_max; R_t(r_max) = -that/Theta'.
    if (ale_on) {
        int n = N;
        double thp = dr1_odd(Theta, n-1);
        // ale_sdot ~ 1/Theta'(r_max): a degrading edge gradient would run the mesh away.
        // Theta' ~ 1 in the weak field at r_max (0.95-0.99 in the paper's runs); outside
        // (0.5, 1.5) the edge has certainly degraded.
        if (!(thp > 0.5 && thp < 1.5)) {
            fprintf(stderr, "FATAL: dTheta/dr at rmax = %.4f at step %d (healthy: 0.95-0.99). "
                    "The ALE mesh velocity is running away; the outer edge of Theta has "
                    "degraded.\n", thp, g_step);
            exit(EXIT_FAILURE);
        }
        double dThdt = NR[n-1]*thp + NN[n-1]*Theta[n-1]*(Q[n-1]/2.0 - trK[n-1]/3.0);
        ale_sdot = -dThdt / thp;                 // g(rmax)=1 => Rt(rmax)=ale_sdot
        for (int i = 0; i < n; ++i) Rt[i] = ale_g[i]*ale_sdot;   // velocity frozen for the step
        if (ref_on()) {         // refinement: zone velocities, then Rt of the composed map
            ref_control(t, dt, n);
            ale_update(n);
        }
    }

    // ---- outer-boundary characteristic gate, once per step ----
    // spin-0 (+) mode eigenvalue bt + N v0, scalar (+) mode bt + N; a mode enters iff > 0.
    // (bt = shift + grid velocity at the edge, printed as "NRm" in the log lines)
    {
        const int n = N;
        const double bt_b = NR[n-1] + Rt[n-1], N_b = NN[n-1];
        if (!isfinite(bt_b) || !isfinite(N_b)) {
            fprintf(stderr, "FATAL: non-finite outer-edge state at step %d (NRm=%g N=%g); "
                    "the blow-up is upstream of the boundary treatment.\n", g_step, bt_b, N_b);
            exit(EXIT_FAILURE);
        }
        const int s0 = (bt_b + N_b*spin0_v0 > 0.0), sc = (bt_b + N_b > 0.0);
        // both modes of a pair entering would need two boundary data; only one is imposed
        if (bt_b - N_b*spin0_v0 > 0.0 || bt_b - N_b > 0.0) {
            fprintf(stderr, "FATAL: at step %d both characteristics of a pair ENTER through "
                    "the outer edge (NRm=%.4f N=%.4f v0=%.4f): two boundary data would be "
                    "needed and only one is imposed.\n", g_step, bt_b, N_b, spin0_v0);
            exit(EXIT_FAILURE);
        }
        if (s0 != gate_spin0) {
            fprintf(stderr, "outer BC gate: spin-0 condition %s at step %d "
                    "(NRm=%.4f N=%.4f v0=%.4f)\n", s0 ? "ON" : "OFF", g_step, bt_b, N_b, spin0_v0);
            gate_spin0 = s0;
        }
        if (sc != gate_scalar) {
            fprintf(stderr, "outer BC gate: scalar condition %s at step %d "
                    "(NRm=%.4f N=%.4f)\n", sc ? "ON" : "OFF", g_step, bt_b, N_b);
            gate_scalar = sc;
        }
    }

    if (g_step % 500 == 0) {
        if (ref_on())
            printf("t=%.6f  dt=%.6e  hor=%d  s_out=%.6g  A: x_c=%.5g xc_dot=%.4g a=%.4g  "
                   "B: x_c=%.5g xc_dot=%.4g a=%.4g\n", t, dt, g_horizon, ale_s,
                   RZ[0].xc, RZ[0].xc_dot, RZ[0].a, RZ[1].xc, RZ[1].xc_dot, RZ[1].a);
        else
            printf("t=%.6f  dt=%.6e  hor=%d\n", t, dt, g_horizon);
    }

    const vector<double> Psi_n   = Psi;
    const vector<double> P_n     = P;
    const vector<double> trK_n   = trK;
    const vector<double> Acr_n   = Acr;
    const vector<double> Theta_n = Theta;

    vector<double> Q_n, NN_n, NR_n;
    compute_Q    (Psi_n, P_n, Theta_n, trK_n, Q_n);
    compute_lapse(Acr_n, NN_n);
    compute_shift(Q_n, trK_n, NN_n, NR_n);

    // d_t Q at the outer edge for the scalar BC: from the previous step for the predictor
    g_Qdot_edge = (g_dt_prev > 0.0) ? (Q_n[N-1] - g_Q_prev_edge) / g_dt_prev : 0.0;

    vector<double> Psi_rhs_n, P_rhs_n, trK_rhs_n, Acr_rhs_n, Theta_rhs_n;
    compute_RHS(Psi_n, P_n, trK_n, Acr_n, Theta_n, Q_n, NN_n, NR_n,
                Psi_rhs_n, P_rhs_n, trK_rhs_n, Acr_rhs_n, Theta_rhs_n);

    // moving mesh: the corrector sees the mesh at t+dt (velocities frozen for the step),
    // with or without refinement
    if (ale_on) {
        ale_s += dt*ale_sdot;
        for (int z = 0; z < 2; ++z)
            if (RZ[z].F > 1.0) { RZ[z].xc += dt*RZ[z].xc_dot; RZ[z].a += dt*RZ[z].a_dot; }
        ale_update(N);
    }

    vector<double> Psi_k(N), P_k(N), trK_k(N), Acr_k(N), Theta_k(N);
    for (int i = 0; i < N; ++i) {
        Psi_k[i]   = Psi_n[i]   + dt * Psi_rhs_n[i];
        P_k[i]     = P_n[i]     + dt * P_rhs_n[i];
        trK_k[i]   = trK_n[i]   + dt * trK_rhs_n[i];
        Acr_k[i]   = Acr_n[i]   + dt * Acr_rhs_n[i];
        Theta_k[i] = Theta_n[i] + dt * Theta_rhs_n[i];
    }

    for (int iter = 0; iter < ICN_MAX_ITERS; ++iter) {
        vector<double> Q_k, NN_k, NR_k;
        compute_Q    (Psi_k, P_k, Theta_k, trK_k, Q_k);
        compute_lapse(Acr_k, NN_k);
        compute_shift(Q_k, trK_k, NN_k, NR_k);
        g_Qdot_edge = (Q_k[N-1] - Q_n[N-1]) / dt;

        vector<double> Psi_rhs_k, P_rhs_k, trK_rhs_k, Acr_rhs_k, Theta_rhs_k;
        compute_RHS(Psi_k, P_k, trK_k, Acr_k, Theta_k, Q_k, NN_k, NR_k,
                    Psi_rhs_k, P_rhs_k, trK_rhs_k, Acr_rhs_k, Theta_rhs_k);

        vector<double> Psi_next(N), P_next(N), trK_next(N),
                       Acr_next(N), Theta_next(N);
        for (int i = 0; i < N; ++i) {
            Psi_next[i]   = Psi_n[i]   + 0.5*dt*(Psi_rhs_n[i]   + Psi_rhs_k[i]);
            P_next[i]     = P_n[i]     + 0.5*dt*(P_rhs_n[i]     + P_rhs_k[i]);
            trK_next[i]   = trK_n[i]   + 0.5*dt*(trK_rhs_n[i]   + trK_rhs_k[i]);
            Acr_next[i]   = Acr_n[i]   + 0.5*dt*(Acr_rhs_n[i]   + Acr_rhs_k[i]);
            Theta_next[i] = Theta_n[i] + 0.5*dt*(Theta_rhs_n[i] + Theta_rhs_k[i]);
        }

        double cn = max(rel_change_norm(Psi_next, Psi_k),
                        rel_change_norm(Theta_next, Theta_k));

        Psi_k.swap(Psi_next); P_k.swap(P_next); trK_k.swap(trK_next);
        Acr_k.swap(Acr_next); Theta_k.swap(Theta_next);

        if (!isfinite(cn)) break;
        if (cn < ICN_TOL && iter >= ICN_MIN_ITERS) break;
    }

    Psi = Psi_k; P = P_k; trK = trK_k; Acr = Acr_k; Theta = Theta_k;

    apply_KO_scaled(Psi,   EPS_KO, dt, false);
    apply_KO_scaled(P,     EPS_KO, dt, false);
    apply_KO_scaled(trK,   EPS_KO, dt, false);
    apply_KO_scaled(Acr,   EPS_KO, dt, true);
    apply_KO_scaled(Theta, EPS_KO, dt, true, /*pin_edge=*/true);

    compute_Q     (Psi, P, Theta, trK, Q);
    compute_lapse (Acr, NN);
    compute_shift (Q, trK, NN, NR);
    compute_DTheta(Theta, DTheta);
    g_Q_prev_edge = Q_n[N-1]; g_dt_prev = dt;   // start-of-step edge value for the next predictor

    check_fields_finite(t + dt);          // the state just produced, last step included
    horizon_finder_step(Theta, trK, Q, N, t + dt);
}

// =====================================================================================
//  INITIAL DATA
// =====================================================================================
// Theta from the Hamiltonian constraint (a second-order ODE for Theta at t = 0, where
// a_r = trK = 0), RK4 outward from Theta(0) = 0, Theta'(0) = 1.  psi_r and, for the
// ingoing shell, P are evaluated analytically at the RK4 stage radii, and the grid
// functions Q, a_r, K are interpolated linearly to them (otherwise the data are first
// order in the shell).  Time-symmetric data have P = Q = 0 and are unaffected.
void integrate_HamiltonianConstraint_initial() {
    int Nloc = (int)Theta.size();
    if (Nloc < 3) return;

    Theta[0] = 0.0;
    vector<double> dth(Nloc, 0.0);
    dth[0] = 1.0;

    double Theta_dd_at0 = 0.0;   // Theta''(0) = 0 (Theta is odd): first step by Taylor
    Theta[1] = Theta[0] + dth[0] * dx + 0.5 * Theta_dd_at0 * dx * dx;
    dth[1]   = dth[0]   + Theta_dd_at0 * dx;

    // P = psi' + psi/r of the ingoing shell at an arbitrary radius
    auto P_ingoing_of = [&](double r) -> double {
        double s2 = cfg.sigma * cfg.sigma;
        double g  = (r*r - cfg.r0*cfg.r0) / s2;
        return shell_dpsi(r) + (cfg.A0 / 100.0) * r * exp(-g*g);
    };
    // w = 0, 1/2, 1: stage at r_i, r_i + dx/2, r_i + dx (grid functions interpolated; the
    // t = 0 map is the identity, r = x).  At t = 0 a_r = K = 0 (the probe is added after
    // this solve), so of the grid functions only Q actually contributes; the a_r and K
    // terms are kept only to mirror eq:hamiltonian (a_r' is taken at the grid point, not
    // interpolated to the stage radius).
    auto Theta_dd = [&](int i, double w, double r_stage, double th_local, double dth_local) -> double {
        double acr_loc = (1.0 - w) * Acr[i] + w * Acr[i+1];
        double p_loc   = (ingoing != 0.0) ? P_ingoing_of(r_stage) : P[i];
        double psi_r   = shell_dpsi(r_stage);
        double k       = (1.0 - w) * trK[i] + w * trK[i+1];
        double q_loc   = (1.0 - w) * Q[i]   + w * Q[i+1];

        double safe_th = th_local;
        if (fabs(safe_th) < 1e-14)
            safe_th = (th_local >= 0 ? 1e-14 : -1e-14);

        double term1 = ((dth_local*dth_local) - 1.0) / (2.0 * safe_th);
        double term2 = c14 * acr_loc * dth_local;
        double term3 = th_local / 4.0 * (
              c14 * (2.0*d1_odd(Acr, i) + acr_loc*acr_loc)
            + p_loc*p_loc
            + psi_r*psi_r
            + 1.5*(1.0 - c13)*q_loc*q_loc
            - (1.0/3.0)*den_K*k*k
        );
        return -(term1 + term2 + term3);
    };

    for (int i = 1; i < Nloc - 1; ++i) {
        double ri    = rgrid[i];          // stage radii: ri, ri+dx/2, ri+dx
        double th    = Theta[i];
        double dth_i = dth[i];

        double k1_th   = dth_i;
        double k1_dth  = Theta_dd(i, 0.0, ri,          th,  dth_i);
        double th2     = th    + 0.5*dx*k1_th;
        double dth2    = dth_i + 0.5*dx*k1_dth;
        double k2_th   = dth2;
        double k2_dth  = Theta_dd(i, 0.5, ri + 0.5*dx, th2, dth2);
        double th3     = th    + 0.5*dx*k2_th;
        double dth3    = dth_i + 0.5*dx*k2_dth;
        double k3_th   = dth3;
        double k3_dth  = Theta_dd(i, 0.5, ri + 0.5*dx, th3, dth3);
        double th4     = th    + dx*k3_th;
        double dth4    = dth_i + dx*k3_dth;
        double k4_th   = dth4;
        double k4_dth  = Theta_dd(i, 1.0, ri + dx,     th4, dth4);

        Theta[i+1] = th    + (dx/6.0)*(k1_th  + 2*k2_th  + 2*k3_th  + k4_th);
        dth[i+1]   = dth_i + (dx/6.0)*(k1_dth + 2*k2_dth + 2*k3_dth + k4_dth);
    }
}

// Allocate the fields, set up the mesh and the refinement zones (with sanity checks),
// set the scalar shell (time symmetric, or ingoing), solve the Hamiltonian constraint for
// Theta, add the optional trK probe, and compute Q, N, beta.
void set_initial_data() {
    rgrid.assign(N, 0.0);
    Psi.assign(N, 0.0);
    P.assign(N, 0.0);
    trK.assign(N, 0.0);
    Acr.assign(N, 0.0);
    Theta.assign(N, 0.0);
    Q.assign(N, 0.0);
    NN.assign(N, 1.0);
    NR.assign(N, 0.0);
    DTheta.assign(N, 0.0);

    // moving mesh: taper over the outer 60% of the box; identity at t = 0
    Rx.assign(N, 1.0); Rxx.assign(N, 0.0); Rt.assign(N, 0.0);
    ale_ra = 0.4 * rmax;
    ale_L  = rmax - ale_ra;
    ale_s = 0.0; ale_sdot = 0.0;
    RZ[0] = RefZone();                     // zone A: the outgoing front
    RZ[0].F = ref_F;  RZ[0].ton = ref_ton;  RZ[0].w = ref_w;
    RZ[0].th_lo = REF_THETA_MIN;  RZ[0].inner = false;  RZ[0].name = "A";
    RZ[1] = RefZone();                     // zone B: an inner feature
    RZ[1].F = refB_F; RZ[1].ton = refB_ton; RZ[1].w = refB_w;
    RZ[1].th_lo = REF_THETAB_LO;  RZ[1].inner = true;   RZ[1].name = "B";
    RZ[1].x0 = refB_x0;  RZ[1].ttrack = refB_ttrack;
    if (ref_on() && ale_on == 0) {
        fprintf(stderr, "FATAL: the moving refinement requires the moving mesh (ale != 0).\n");
        exit(EXIT_FAILURE);
    }
    for (int z = 0; z < 2; ++z)
        if (RZ[z].F > 1.0 && !(RZ[z].w >= 16.0*dx)) {
            fprintf(stderr, "FATAL: refinement %s half-width w=%.4g < 16 dx=%.4g; raise w or N.\n",
                    RZ[z].name, RZ[z].w, 16.0*dx);
            exit(EXIT_FAILURE);
        }
    ale_build_profile(N);
    ale_update(N);

    for (int i = 0; i < N; ++i)
        rgrid[i] = i * dx;                     // t = 0 map is the identity (ale_update
                                               // already gives it; explicit for clarity)

    for (int i = 0; i < N; ++i)
        Psi[i] = shell_psi(rgrid[i]);

    if (ingoing != 0.0) {
        // purely ingoing shell (arXiv:1608.06970): P = psi_r + psi/r at t = 0 (N=1, beta=0).
        // Set before the Theta solve, whose source contains P^2.  (psi/r is formed from the
        // grid value of Psi, which differs from the analytic form used at the RK4 stage
        // radii only by round-off.)
        for (int i = 0; i < N; ++i) {
            double r = rgrid[i];
            P[i] = shell_dpsi(r) + (r > 1e-12 ? Psi[i]/r : 0.0);
        }
    }

    integrate_HamiltonianConstraint_initial();

    if (ingoing != 0.0) {
        // With P != 0 the momentum constraint sources Q, and Q^2 enters the Hamiltonian
        // constraint, so the two solves are coupled: iterate them to a fixed point.  (The
        // time-symmetric data have P = trK = 0, hence Q = 0 identically, and skip this.)
        // Solving Theta once with Q = 0 would leave a constraint violation at t = 0 that
        // does not shrink with resolution.
        for (int it = 0; it < ID_MAX_ITERS; ++it) {
            vector<double> Theta_prev = Theta;
            compute_Q(Psi, P, Theta, trK, Q);
            integrate_HamiltonianConstraint_initial();
            double dmax = 0.0;
            for (int i = 0; i < N; ++i) dmax = max(dmax, fabs(Theta[i] - Theta_prev[i]));
            if (dmax < ID_TOL) { fprintf(stderr, "initial data: Theta/Q solves converged in %d iterations\n", it+1); break; }
            if (it == ID_MAX_ITERS - 1) fprintf(stderr, "WARNING: initial-data Theta/Q iteration did not converge (max dTheta=%.3e)\n", dmax);
        }
    }

    if (probe_amp != 0.0) {
        // linear spin-0 probe for run-differencing, added after the constraint solve (the
        // README recipe "Probe runs" uses --probe_amp=1e-6, small enough to stay linear;
        // paper fig:probe)
        for (int i = 0; i < N; ++i) {
            double g = (rgrid[i] - probe_r) / cfg.sigma;   // width = the shell sigma
            trK[i] += probe_amp * exp(-g*g);
        }
    }

    compute_Q     (Psi, P, Theta, trK, Q);
    compute_lapse (Acr, NN);
    compute_shift (Q, trK, NN, NR);
    compute_DTheta(Theta, DTheta);
    g_Q_prev_edge = Q[N-1]; g_dt_prev = 0.0;
    check_fields_finite(0.0);
}

// =====================================================================================
//  DIAGNOSTICS AND OUTPUT
// =====================================================================================
// Regularised Hamiltonian constraint Theta*C (eq:hamiltonian, see the file header).  Never
// monitor C itself: (Theta'^2-1)/(2 Theta) is 0/0 at r = 0.
double theta_C(int i) {
    const double den = den_K;
    double th   = Theta[i];
    double dTh  = dr1_odd (Theta, i);
    double d2Th = dr2_odd (Theta, i);
    double dps  = dr1_even(Psi,   i);
    double dar  = dr1_odd (Acr,   i);
    return th*d2Th + (dTh*dTh - 1.0)/2.0 + c14*Acr[i]*th*dTh
         + (th*th/4.0)*( c14*(2.0*dar + Acr[i]*Acr[i]) + P[i]*P[i] + dps*dps
                       + 1.5*(1.0 - c13)*Q[i]*Q[i] - (1.0/3.0)*den*trK[i]*trK[i] );
}

// Full profile, written at the first step with t >= the requested time t_req; the file is
// named after t_req (the actual time differs by less than one step).  NOTE: the name keeps
// ONE decimal (%05.1f), so requested times closer than 0.1 overwrite each other, and e.g.
// 29.25 is written as prof_t029.2.csv (README, "Output"; kept for compatibility with the
// paper's data sets).  Scripts that need the exact time read it from phi.dat.
void dump_profile(double t_req) {
    char fname[256];
    snprintf(fname, sizeof(fname), "prof_t%05.1f.csv", t_req);
    ofstream f(fname);
    f << setprecision(17) << scientific;
    f << "# r,Psi,P,trK,Acr,Theta,NN,Q,ThetaC\n";
    for (int i = 0; i < N; ++i)
        f << rgrid[i] << "," << Psi[i]  << "," << P[i]   << "," << trK[i] << ","
          << Acr[i]   << "," << Theta[i]<< "," << NN[i]  << "," << Q[i]   << ","
          << theta_C(i) << "\n";
}

// auxiliary dump of one field: i, r, value
void dump_field_csv(const string& fieldname, const vector<double>& vec, int dump_id) {
    char fname[256];
    snprintf(fname, sizeof(fname), "%s_dump_%06d.csv", fieldname.c_str(), dump_id);
    ofstream f(fname);
    f << setprecision(17) << scientific;
    f << "# i,r," << fieldname << "\n";
    for (int i = 0; i < (int)vec.size(); ++i)
        f << i << "," << rgrid[i] << "," << vec[i] << "\n";
}

// auxiliary combined dump (all fields incl. the shift NR and d_r Theta)
void dump_combined_csv(int dump_id) {
    char fname[256];
    snprintf(fname, sizeof(fname), "psi_dump_%06d.csv", dump_id);
    ofstream f(fname);
    f << setprecision(17) << scientific;
    f << "# i,r,Psi,P,trK,Acr,Theta,NN,NR,Q,DTheta\n";
    int n = (int)Psi.size();
    for (int i = 0; i < n; ++i) {
        f << i << "," << rgrid[i] << ","
          << Psi[i]    << "," << P[i]     << "," << trK[i]    << ","
          << Acr[i]    << "," << Theta[i] << ","
          << NN[i]     << "," << NR[i]    << "," << Q[i]      << ","
          << DTheta[i] << "\n";
    }
}

// the combined file and the per-field files hold the same data; both are written for
// compatibility with the tools of the development version (off by default)
void dump_auxiliary(int dump_id) {
    dump_combined_csv(dump_id);
    dump_field_csv("Psi",    Psi,    dump_id);
    dump_field_csv("P",      P,      dump_id);
    dump_field_csv("trK",    trK,    dump_id);
    dump_field_csv("Acr",    Acr,    dump_id);
    dump_field_csv("Theta",  Theta,  dump_id);
    dump_field_csv("NN",     NN,     dump_id);
    dump_field_csv("NR",     NR,     dump_id);
    dump_field_csv("Q",      Q,      dump_id);
    dump_field_csv("DTheta", DTheta, dump_id);
}

// =====================================================================================
//  SPHERE RECORDER  (--rec_spheres=Th1,Th2,...  --rec_every=k;  off by default)
// =====================================================================================
// For time derivatives on spheres of fixed areal radius (flux identity, -dm/dt|_Theta)
// at the step cadence instead of the profile cadence.  At t = 0 and after every k-th
// completed step (state, mesh and the Q / N / beta solves all at the new time), for each
// requested Theta_s: locate the OUTERMOST i with Theta[i] <= Theta_s < Theta[i+1] (Theta is
// monotonic outside the frozen centre) and store the 13 points i-6 .. i+6 of
//     r (= rgrid, the physical radius of the moving mesh), Psi, P, trK, Acr, Theta, NN, Q,
//     NR (= beta, the shift)
// so that the analysis can take centred r-derivatives within one record and interpolate
// to Theta_s.  The stencil is clamped to [0, N-1]; if no crossing exists (Theta_s beyond
// the grid or inside the frozen region) the stored index is -1 and the stencil is the
// outermost 13 points.  Output: spheres.bin, little-endian float64 records of
//     step, t, dt,  then per sphere: i, r[13], Psi[13], P[13], trK[13], Acr[13],
//                                    Theta[13], NN[13], Q[13], NR[13]
// i.e. 3 + n_spheres*(1 + 9*13) doubles per record (dt = the step just completed; 0 in the
// t = 0 record); spheres.txt describes the layout.  Writes go through a 4 MB stdio
// buffer.  The recorder only reads the state: the evolution is unchanged with it on.
static const int  REC_HALF = 6, REC_NST = 2*REC_HALF + 1, REC_NF = 9;
static FILE*          rec_fp = nullptr;
static vector<double> rec_buf;

static void rec_open() {
    if (rec_theta.empty()) return;
    if (rec_every < 1) { fprintf(stderr, "FATAL: --rec_every=%d < 1.\n", rec_every); exit(EXIT_FAILURE); }
    if (N < REC_NST) { fprintf(stderr, "FATAL: N < %d with the recorder on.\n", REC_NST); exit(EXIT_FAILURE); }
    rec_fp = fopen("spheres.bin", "wb");
    if (!rec_fp) { fprintf(stderr, "FATAL: cannot open spheres.bin.\n"); exit(EXIT_FAILURE); }
    setvbuf(rec_fp, nullptr, _IOFBF, 1 << 22);
    const int ns = (int)rec_theta.size();
    FILE* h = fopen("spheres.txt", "w");
    if (!h) { fprintf(stderr, "FATAL: cannot open spheres.txt.\n"); exit(EXIT_FAILURE); }
    fprintf(h, "# sphere recorder of collapse.cpp: layout of spheres.bin\n");
    fprintf(h, "format float64_le\n");
    fprintf(h, "n_spheres %d\n", ns);
    fprintf(h, "theta_s");
    for (double th : rec_theta) fprintf(h, " %.17g", th);
    fprintf(h, "\n");
    fprintf(h, "n_stencil %d\n", REC_NST);
    fprintf(h, "fields r Psi P trK Acr Theta NN Q NR\n");
    fprintf(h, "record_doubles %d\n", 3 + ns*(1 + REC_NF*REC_NST));
    fprintf(h, "record step t dt [ i_lo r[13] Psi[13] P[13] trK[13] Acr[13] Theta[13] NN[13] Q[13] NR[13] ] x n_spheres\n");
    fprintf(h, "# i_lo: outermost grid index with Theta[i] <= theta_s < Theta[i+1] (-1: none);\n");
    fprintf(h, "# the stencil holds the grid points i_lo-6 .. i_lo+6 (shifted to stay inside 0..N-1)\n");
    fprintf(h, "rec_every %d\n", rec_every);
    fprintf(h, "N %d\nrmax %.17g\nc1 %.17g\nc2 %.17g\nA0 %.17g\nr0 %.17g\nsigma %.17g\n"
               "ingoing %.17g\nale %d\n",
            N, rmax, c1, c2, cfg.A0, cfg.r0, cfg.sigma, ingoing, ale_on);
    fclose(h);
    rec_buf.assign(3 + ns*(1 + REC_NF*REC_NST), 0.0);
}

static void rec_write(int step, double t, double dt_step) {
    if (!rec_fp) return;
    double* b = rec_buf.data();
    *b++ = step; *b++ = t; *b++ = dt_step;
    const vector<double>* F[REC_NF] = { &rgrid, &Psi, &P, &trK, &Acr, &Theta, &NN, &Q, &NR };
    for (double ths : rec_theta) {
        int ic = -1;
        for (int i = N - 2; i >= 0; --i)
            if (Theta[i] <= ths && ths < Theta[i+1]) { ic = i; break; }
        int lo = (ic < 0 ? N - REC_NST : ic - REC_HALF);
        lo = max(0, min(lo, N - REC_NST));
        *b++ = ic;
        for (int f = 0; f < REC_NF; ++f)
            for (int k = 0; k < REC_NST; ++k) *b++ = (*F[f])[lo + k];
    }
    fwrite(rec_buf.data(), sizeof(double), rec_buf.size(), rec_fp);
}

static void rec_close() { if (rec_fp) { fclose(rec_fp); rec_fp = nullptr; } }

// =====================================================================================
//  COMMAND LINE
// =====================================================================================
static const char* USAGE =
"Usage (named options, any order):\n"
"  ./collapse [--N=25600] [--A0=0.35] [--t_final=T] [--times=t1,t2,...] [--rmax=300]\n"
"             [--c1=0.7] [--c2=-0.5119403] [--ingoing=0] [--ale=1]\n"
"             [--ref_F=0] [--ref_ton=50] [--ref_w=5]\n"
"             [--refB_F=0] [--refB_ton=20] [--refB_w=2] [--refB_x0=0] [--refB_ttrack=0]\n"
"             [--r0=10] [--sigma=4] [--probe_amp=0] [--probe_r=150] [--dump_every=0]\n"
"             [--rec_spheres=Th1,Th2,...] [--rec_every=1]\n"
"  (the environment variable COLLAPSE_DUMP_EVERY, if set, overrides --dump_every)\n"
"The positional form of the development version is also accepted: ./collapse --help-legacy\n"
"See README.md.\n";

static const char* USAGE_LEGACY =
"Legacy positional form of the development version (22 slots, trailing ones optional):\n"
"  ./collapse N A0 t_final times rmax c1 c2 ingoing g_ale ale_inner ref_F ref_ton ref_w\n"
"             refB_F refB_ton refB_w refB_Fmax qform refB_x0 refB_ttrack r0 sigma\n"
"  (ale_inner, refB_Fmax and qform were removed and must be 0 or omitted.)\n"
"See README.md, \"Legacy positional form and mapping\".\n";

static void fatal_usage(const char* msg, const char* arg) {
    fprintf(stderr, "FATAL: %s '%s'.\n%s", msg, arg, USAGE);
    exit(EXIT_FAILURE);
}

// numeric option values: the whole string must parse
static double num_value(const char* opt, const char* v) {
    char* end = nullptr;
    double x = strtod(v, &end);
    if (end == v || *end != '\0') fatal_usage("option needs a numeric value, got", opt);
    return x;
}

// comma-separated list of numbers (empty entries are skipped), sorted ascending
static void parse_list(const string& s, vector<double>& out) {
    size_t pos = 0;
    while (pos <= s.size()) {
        size_t c = s.find(',', pos);
        if (c == string::npos) c = s.size();
        if (c > pos) {
            const string item = s.substr(pos, c - pos);
            out.push_back(num_value(item.c_str(), item.c_str()));
        }
        pos = c + 1;
    }
    sort(out.begin(), out.end());
}
static void parse_times(const string& s) { parse_list(s, t_profile); }
static int int_value(const char* opt, const char* v) {
    char* end = nullptr;
    long x = strtol(v, &end, 10);
    if (end == v || *end != '\0') fatal_usage("option needs an integer value, got", opt);
    return (int)x;
}

// Named options: --key=value.
static void parse_named(int argc, char** argv) {
    for (int a = 1; a < argc; ++a) {
        const char* s = argv[a];
        if (strcmp(s, "-h") == 0 || strcmp(s, "--help") == 0) { fputs(USAGE, stdout); exit(0); }
        if (strcmp(s, "--help-legacy") == 0) { fputs(USAGE_LEGACY, stdout); exit(0); }
        if (strncmp(s, "--", 2) != 0) fatal_usage("named and positional arguments cannot be mixed", s);
        const char* eq = strchr(s, '=');
        if (!eq) fatal_usage("option needs a value (--key=value)", s);
        const string key(s + 2, eq - (s + 2));
        const char* v = eq + 1;
        if      (key == "N")           N            = int_value(s, v);
        else if (key == "A0")          cfg.A0       = num_value(s, v);
        else if (key == "t_final")   { t_final      = num_value(s, v); have_t_final = true; }
        else if (key == "times")       parse_times(v);
        else if (key == "rmax")        rmax         = num_value(s, v);
        else if (key == "c1")          c1           = num_value(s, v);
        else if (key == "c2")          c2           = num_value(s, v);
        else if (key == "ingoing")     ingoing      = num_value(s, v);
        else if (key == "ale")         ale_on       = int_value(s, v);
        else if (key == "ref_F")       ref_F        = num_value(s, v);
        else if (key == "ref_ton")     ref_ton      = num_value(s, v);
        else if (key == "ref_w")       ref_w        = num_value(s, v);
        else if (key == "refB_F")      refB_F       = num_value(s, v);
        else if (key == "refB_ton")    refB_ton     = num_value(s, v);
        else if (key == "refB_w")      refB_w       = num_value(s, v);
        else if (key == "refB_x0")     refB_x0      = num_value(s, v);
        else if (key == "refB_ttrack") refB_ttrack  = num_value(s, v);
        else if (key == "r0")          cfg.r0       = num_value(s, v);
        else if (key == "sigma")       cfg.sigma    = num_value(s, v);
        else if (key == "probe_amp")   probe_amp    = num_value(s, v);
        else if (key == "probe_r")     probe_r      = num_value(s, v);
        else if (key == "dump_every")  dump_every   = int_value(s, v);
        else if (key == "rec_spheres") parse_list(v, rec_theta);
        else if (key == "rec_every")   rec_every    = int_value(s, v);
        else fatal_usage("unknown option", s);
    }
}

// Legacy positional form of the development version (argv[1..22]), kept so that its launch
// lines still run.  Slots 10 (ale_inner), 17 (refB_Fmax) and 18 (qform) belonged to removed
// features and must be 0.
static void parse_positional(int argc, char** argv) {
    if (argc > 23) fatal_usage("too many positional arguments, first extra one is", argv[23]);
    if (argc > 1 && string(argv[1]).find('.') != string::npos)
        fatal_usage("first argument is N (an integer), got", argv[1]);
    if (argc > 1) N         = atoi(argv[1]);
    if (argc > 2) cfg.A0    = atof(argv[2]);
    if (argc > 3) { t_final = atof(argv[3]); have_t_final = true; }
    if (argc > 4) parse_times(argv[4]);
    if (argc > 5) rmax      = atof(argv[5]);
    if (argc > 6) c1        = atof(argv[6]);
    if (argc > 7) c2        = atof(argv[7]);
    if (argc > 8) ingoing   = atof(argv[8]);
    if (argc > 9) ale_on    = atoi(argv[9]);
    if (argc > 10 && atof(argv[10]) != 0.0)
        fatal_usage("the inner taper (slot 10, ale_inner) was removed; got", argv[10]);
    if (argc > 11) ref_F    = atof(argv[11]);
    if (argc > 12) ref_ton  = atof(argv[12]);
    if (argc > 13) ref_w    = atof(argv[13]);
    if (argc > 14) refB_F   = atof(argv[14]);
    if (argc > 15) refB_ton = atof(argv[15]);
    if (argc > 16) refB_w   = atof(argv[16]);
    if (argc > 17 && atof(argv[17]) != 0.0)
        fatal_usage("the adaptive zone-B factor (slot 17, refB_Fmax) was removed; got", argv[17]);
    if (argc > 18 && atoi(argv[18]) != 0)
        fatal_usage("the W-form momentum constraint (slot 18, qform) was removed; got", argv[18]);
    if (argc > 19) refB_x0     = atof(argv[19]);
    if (argc > 20) refB_ttrack = atof(argv[20]);
    if (argc > 21) cfg.r0      = atof(argv[21]);
    if (argc > 22) cfg.sigma   = atof(argv[22]);
}

// =====================================================================================
//  MAIN
// =====================================================================================
int main(int argc, char** argv) {
    bool named = false;
    for (int a = 1; a < argc; ++a)
        if (strncmp(argv[a], "--", 2) == 0 || strcmp(argv[a], "-h") == 0) named = true;
    if (named) parse_named(argc, argv);
    else       parse_positional(argc, argv);
    if (ale_on != 0 && ale_on != 1) {
        fprintf(stderr, "FATAL: --ale (slot 9, g_ale) must be 0 or 1, got %d.\n", ale_on);
        return 1;
    }
    // auxiliary dumps: the environment variable of the development version OVERRIDES
    // --dump_every (kept for its launch scripts); say so, since the precedence is unusual
    if (const char* e = getenv("COLLAPSE_DUMP_EVERY")) {
        dump_every = atoi(e);
        fprintf(stderr, "note: COLLAPSE_DUMP_EVERY=%s sets the auxiliary-dump interval to %d "
                "(it overrides --dump_every)\n", e, dump_every);
    }
    if (N < 10) {
        fprintf(stderr, "FATAL: N=%d.\n%s", N, USAGE);
        return 1;
    }

    {   // coupling combinations, and the spin-0 characteristic speed for these couplings
        // (outer projector and gate)
        c13 = c1 + c3;  c14 = c1 + c4;  c123 = c1 + c2 + c3;
        den_K = 2.0 + c13 + 3.0*c2;
        den_a = c14 * (1.0 - c13);
        if (fabs(den_K) < 1e-12 || fabs(den_a) < 1e-12) {
            fprintf(stderr, "FATAL: degenerate couplings (2+c1+c3+3c2 = %.3g, "
                    "(c1+c4)(1-c1-c3) = %.3g).\n", den_K, den_a);
            return 1;
        }
        spin0_aN = (c14 - 2.0) / den_K;
        spin0_bN = -c123 / den_a;
        if (!(spin0_aN*spin0_bN > 0.0)) {
            fprintf(stderr, "FATAL: aN*bN = %.6g <= 0: the spin-0 sector has no real "
                    "characteristic speed for these couplings.\n", spin0_aN*spin0_bN);
            return 1;
        }
        spin0_v0 = sqrt(spin0_aN*spin0_bN);
        fprintf(stderr, "spin-0 speed v0 = %.6f (aN=%.6f bN=%.6f)\n", spin0_v0, spin0_aN, spin0_bN);
    }

    dx = rmax / (N - 1);

    set_initial_data();

    fphi.open("phi.dat");                      // step, t, dt, Psi(r = 0) every step
    fphi << setprecision(17) << scientific;
    fphi << "# step t dt phi\n";

    ftimes.open("psi_times.dat");              // index of the auxiliary dumps
    ftimes << setprecision(17) << scientific;
    ftimes << "# dump_id step t dt\n";

    double t       = 0.0;
    int    dump_id = 0;

    dt = CFL * dx;

    fphi << 0 << " " << t << " " << dt << " " << Psi[0] << "\n";

    rec_open();
    rec_write(0, t, 0.0);

    if (dump_every > 0) {
        dump_auxiliary(dump_id);
        ftimes << dump_id << " " << 0 << " " << t << " " << dt << "\n";
        dump_id++;
    }

    size_t i_prof = 0;

    // Run to t_final if given (so that every grid of a ladder stops at the same time),
    // otherwise for NT_DEFAULT steps.
    for (int step = 1; ; ++step) {
        if (have_t_final && t >= t_final) break;
        if (!have_t_final && step > NT_DEFAULT) break;
        ICN_step(t);
        t += dt;

        while (i_prof < t_profile.size() && t >= t_profile[i_prof]) {
            dump_profile(t_profile[i_prof]);
            ++i_prof;
        }

        fphi << step << " " << t << " " << dt << " " << Psi[0] << "\n";

        if (dump_every > 0 && step % dump_every == 0) {
            dump_auxiliary(dump_id);
            ftimes << dump_id << " " << step << " " << t << " " << dt << "\n";
            dump_id++;
        }

        if (rec_fp && step % rec_every == 0) rec_write(step, t, dt);
    }

    rec_close();
    return 0;
}
