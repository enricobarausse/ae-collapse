"""Shared loaders and numerical helpers of all the analysis scripts: dump I/O, ALE-grid
derivatives, spheres of fixed areal radius, aether and scalar stress tensors, horizon and
peak finding, norms.  Paper labels are cited as eq:NAME / fig:NAME / sec:NAME (the LaTeX
labels of the paper; equation numbers move between versions, the labels do not).

Input: the profile dumps of the collapse code, one file per requested time,

    <rundir>/prof_tTTT.T.csv     columns  r, Psi, P, trK, Acr, Theta, NN, Q, ThetaC

(r = proper radial distance, Psi = scalar, P = its momentum, trK = trace of the extrinsic
curvature K, Acr = a_r, Theta = areal radius, NN = lapse N, Q = shear of K).  With the
moving mesh (--ale=1, the default) the r column moves from dump to dump and is NOT
uniform; radial derivatives here are `ddr` = np.gradient(f, r, edge_order=2), second order
everywhere on the non-uniform grid (some scripts call np.gradient(f, r) directly with the
default edge_order=1, first order at the two end points only; of these only the grid-edge
mass of fig_conv_collapse.m_at uses an end point).

A dump is written at the first step with t >= t_req and named after t_req with ONE
decimal (so 29.25 -> prof_t029.2.csv), and its real time is late by up to one step.  When
<rundir>/phi.dat (step, t, dt, phi at every step) exists we recover the exact time of every
dump from it; this matters for the time derivatives at the sub-percent level.  Which call
to use:
  * list_dumps(rundir)               all dumps, exact times (requested times that are
                                     multiples of 0.1, i.e. every run of the paper except
                                     the quarter-time dumps of the inner-shell ladder);
  * list_dumps(rundir, exact=False)  all dumps, nominal times of the file names (enough
                                     where the time is only a label);
  * requested_dumps(rundir, t_req)   the dumps requested at the times t_req, opened by
                                     name, with exact times (quarter-time dumps);
  * load_dump(rundir, t_req)         one dump by its requested time.

Conventions (the paper's, and the code's):
  * couplings c1, c2 (c3 = c4 = 0 throughout the paper), c13 = c1 + c3, c14 = c1 + c4,
    c123 = c1 + c2 + c3;
  * a_N = (c14 - 2)/(2 + c13 + 3 c2),  b_N = -c123/(c14 (1 - c13)),  v0^2 = a_N b_N
    [paper eq:aNbN];
  * w_out = -(b_N/v0) trK + a_r,  w_in = +(b_N/v0) trK + a_r     [eq:wout];
  * theta_+- = Theta' +- (1/6)(3Q - 2K) Theta                       [eq:expansions];
  * Misner-Sharp mass m = (Theta/2)(1 - h), h = Theta'^2 - (3Q - 2K)^2 Theta^2/36  [eq:ms];
  * u = t - Theta_*/v0,  Theta_* = Theta + 2M ln(Theta/(2M) - 1)                [eq:retarded];
  * shift beta: beta' = N (Q + K/3), beta(0) = 0 (trapezoid on the non-uniform grid);
  * aether frame u^a, unit radial e^a; rho = T_ab u^a u^b, j = -T_ab u^a e^b,
    p = T_ab e^a e^b; flux through a sphere of fixed areal radius        [eq:lum]
       F_X = (N Theta^2/2) [ A (rho_X + p_X) + ((A^2 + B^2)/B) j_X ],
       A = Theta (Q/2 - K/3),  B = Theta'.

Everything is plain numpy (single-threaded: the thread variables are pinned below before
numpy is imported).  Dumps are processed by `map_dumps`, optionally in a few worker
processes (fork).
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import glob
import re
import multiprocessing as mp
import numpy as np

COLS = ("r", "Psi", "P", "trK", "Acr", "Theta", "NN", "Q", "ThetaC")
# the fields whose self-convergence the convergence scripts measure (paper tab:conv-scatter)
CONV_FIELDS = ("Psi", "P", "trK", "Acr", "Theta", "Q")

# baseline couplings of the paper (GEJ values; v0 = 1)
C1_DEFAULT = 0.7
C2_BASELINE = -0.5119403


# ----------------------------------------------------------------------------------------
# couplings
# ----------------------------------------------------------------------------------------
def couplings(c1=C1_DEFAULT, c2=C2_BASELINE, c3=0.0, c4=0.0):
    """Coupling combinations, the principal-part coefficients a_N, b_N and v0."""
    c13, c14, c123 = c1 + c3, c1 + c4, c1 + c2 + c3
    aN = (c14 - 2.0) / (2.0 + c13 + 3.0 * c2)
    bN = -c123 / (c14 * (1.0 - c13))
    v0sq = aN * bN
    if v0sq <= 0:
        raise ValueError(f"no real spin-0 speed for c1={c1}, c2={c2}")
    return dict(c1=c1, c2=c2, c3=c3, c4=c4, c13=c13, c14=c14, c123=c123,
                aN=aN, bN=bN, v0=float(np.sqrt(v0sq)))


# ----------------------------------------------------------------------------------------
# dumps
# ----------------------------------------------------------------------------------------
_TRE = re.compile(r"prof_t([0-9]+(?:\.[0-9]*)?)\.csv$")


def nominal_time(path):
    """Time in the file name of a dump (the requested time rounded to one decimal)."""
    return float(_TRE.search(os.path.basename(path)).group(1))


def dump_path(rundir, t_req):
    """File name of the dump requested at t_req (format of dump_profile in collapse.cpp)."""
    return os.path.join(rundir, f"prof_t{t_req:05.1f}.csv")


def load_dump(rundir, t_req):
    """Profile dict of the dump requested at t_req (see load_profile)."""
    return load_profile(dump_path(rundir, t_req))


def requested_dumps(rundir, t_req):
    """Dumps requested at the times t_req, opened by name -> (t_req, t, paths), sorted by
    t_req, the requested times without a file skipped.  t is the exact time of each dump
    (first step with t >= t_req, from phi.dat; t_req if phi.dat is unavailable).  Use it
    when requested times closer than 0.1 to a multiple of 0.1 matter (e.g. 29.25, whose
    file prof_t029.2.csv would otherwise be read as t = 29.2)."""
    tr = np.array(sorted(float(x) for x in t_req))
    have = np.array([os.path.exists(dump_path(rundir, x)) for x in tr], bool)
    tr = tr[have]
    paths = [dump_path(rundir, x) for x in tr]
    if not paths:
        raise FileNotFoundError(f"none of the requested dumps in {rundir}")
    tex = exact_dump_times(rundir, tr)
    return tr, (tr.copy() if tex is None else tex), paths


def list_dumps(rundir, tmin=-np.inf, tmax=np.inf, exact=True):
    """Sorted (times, paths) of the profile dumps in rundir with tmin <= t <= tmax.

    exact=True: replace the nominal time t_req (file name) by the time of the first step
    with t >= t_req, read from <rundir>/phi.dat when present (the code dumps at that step).
    The file name carries one decimal, so this is exact only for requested times that are
    multiples of 0.1 (quarter-time dumps, e.g. t_req = 29.25 -> prof_t029.2.csv, come back
    up to 0.05 early); for such dumps use requested_dumps (see the module docstring).
    """
    paths = [p for p in glob.glob(os.path.join(rundir, "prof_t*.csv")) if _TRE.search(p)]
    if not paths:
        raise FileNotFoundError(f"no prof_t*.csv in {rundir}")
    tn = np.array([nominal_time(p) for p in paths])
    o = np.argsort(tn)
    tn = tn[o]
    paths = [paths[i] for i in o]
    sel = (tn >= tmin - 1e-9) & (tn <= tmax + 1e-9)
    tn = tn[sel]
    paths = [p for p, s in zip(paths, sel) if s]
    t = tn.copy()
    if exact:
        tex = exact_dump_times(rundir, tn)
        if tex is not None:
            t = tex
    return t, paths


def exact_dump_times(rundir, t_req):
    """Real times of the dumps requested at t_req, from phi.dat (None if unavailable)."""
    f = os.path.join(rundir, "phi.dat")
    if not os.path.exists(f):
        return None
    try:
        ts = np.loadtxt(f, comments="#", usecols=1)
    except Exception:
        return None
    # the code: t += dt; while t >= t_req: dump  -> first step with t >= t_req
    k = np.searchsorted(ts, np.asarray(t_req, float), side="left")
    if np.any(k >= len(ts)):
        return None
    out = ts[k]
    if np.any(np.abs(out - t_req) > 0.05):       # phi.dat inconsistent with the dumps
        return None
    return out


def load_profile(path):
    """Load one prof_t*.csv dump -> dict of 1-d arrays keyed by COLS."""
    a = np.loadtxt(path, delimiter=",", comments="#", ndmin=2)
    if a.shape[1] < 8:
        raise ValueError(f"{path}: expected >= 8 columns, got {a.shape[1]}")
    d = {k: a[:, i] for i, k in enumerate(COLS[:a.shape[1]])}
    return d


# ----------------------------------------------------------------------------------------
# derivatives, interpolation
# ----------------------------------------------------------------------------------------
def ddr(y, r):
    """d/dr on the (non-uniform) grid r: second-order centred, one-sided at the ends."""
    return np.gradient(y, r, edge_order=2)


def fd_weights(x0, x, m):
    """Fornberg weights for the m-th derivative at x0 from the nodes x (any spacing)."""
    x = np.asarray(x, float)
    n = len(x)
    c = np.zeros((n, m + 1))
    c1, c4 = 1.0, x[0] - x0
    c[0, 0] = 1.0
    for i in range(1, n):
        mn = min(i, m)
        c2, c5, c4 = 1.0, c4, x[i] - x0
        for j in range(i):
            c3 = x[i] - x[j]
            c2 *= c3
            if j == i - 1:
                for k in range(mn, 0, -1):
                    c[i, k] = c1 * (k * c[i - 1, k - 1] - c5 * c[i - 1, k]) / c2
                c[i, 0] = -c1 * c5 * c[i - 1, 0] / c2
            for k in range(mn, 0, -1):
                c[j, k] = (c4 * c[j, k] - k * c[j, k - 1]) / c3
            c[j, 0] = c4 * c[j, 0] / c3
        c1 = c2
    return c[:, m]


def ddt_series(t, y, order=4):
    """Time derivative of a series y(t) (axis 0) on the (nearly uniform) dump times t.

    order=4: 5-point centred stencil (Fornberg weights on the actual, jittered times),
    order=2: 3-point; the stencil becomes one-sided at the ends."""
    t = np.asarray(t, float)
    y = np.asarray(y, float)
    h = order // 2
    n = len(t)
    out = np.empty_like(y)
    for i in range(n):
        lo = min(max(i - h, 0), max(n - 2 * h - 1, 0))
        idx = np.arange(lo, min(lo + 2 * h + 1, n))
        w = fd_weights(t[i], t[idx], 1)
        out[i] = np.tensordot(w, y[idx], axes=(0, 0))
    return out


def lagrange_interp(x, y, xq, npts=4):
    """Piecewise Lagrange interpolation (npts nodes, default cubic) of y(x) at xq.

    x must be strictly increasing.  y may be a 1-d array or a list of 1-d arrays (all
    sampled at x); xq outside [x[0], x[-1]] gives NaN."""
    x = np.asarray(x, float)
    xq = np.atleast_1d(np.asarray(xq, float))
    ys = y if isinstance(y, (list, tuple)) else [y]
    j = np.searchsorted(x, xq) - 1
    lo = np.clip(j - (npts // 2 - 1), 0, len(x) - npts)
    idx = lo[:, None] + np.arange(npts)[None, :]
    X = x[idx]
    w = np.ones_like(X)
    for a in range(npts):
        for b in range(npts):
            if a != b:
                w[:, a] *= (xq - X[:, b]) / (X[:, a] - X[:, b])
    bad = (xq < x[0]) | (xq > x[-1])
    outs = []
    for yy in ys:
        v = (np.asarray(yy)[idx] * w).sum(axis=1)
        v[bad] = np.nan
        outs.append(v)
    return outs if isinstance(y, (list, tuple)) else outs[0]


def exterior_start(Theta):
    """First index of the outer region where Theta(r) increases monotonically to r_max.

    Spheres of fixed areal radius are located in this region only (the frozen interior
    can be non-monotonic)."""
    bad = np.nonzero(np.diff(Theta) <= 0.0)[0]
    return 0 if len(bad) == 0 else int(bad[-1]) + 1


def at_theta(prof_Theta, targets, fields, npts=4):
    """Values of fields (dict name -> array on the r grid) on the spheres Theta = targets.

    Cubic Lagrange interpolation with the areal radius as abscissa, restricted to the
    monotonic exterior.  Returns dict name -> array(len(targets))."""
    i0 = exterior_start(prof_Theta)
    Th = prof_Theta[i0:]
    names = list(fields)
    vals = lagrange_interp(Th, [np.asarray(fields[k])[i0:] for k in names], targets, npts)
    return dict(zip(names, vals))


# ----------------------------------------------------------------------------------------
# derived fields
# ----------------------------------------------------------------------------------------
def beta_shift(p):
    """Shift: beta' = N (Q + K/3) integrated outward from beta(0) = 0 (trapezoid)."""
    f = p["NN"] * (p["Q"] + p["trK"] / 3.0)
    r = p["r"]
    return np.concatenate(([0.0], np.cumsum(0.5 * np.diff(r) * (f[:-1] + f[1:]))))


def theta_pm(p, Thp=None):
    """Expansions theta_+- = Theta' +- (1/6)(3Q - 2K) Theta."""
    if Thp is None:
        Thp = ddr(p["Theta"], p["r"])
    s = (3.0 * p["Q"] - 2.0 * p["trK"]) * p["Theta"] / 6.0
    return Thp + s, Thp - s


def misner_sharp(p, Thp=None):
    """Misner-Sharp mass, paper eq:ms: m = (Theta/2)(1 - h),
    h = Theta'^2 - (3Q - 2K)^2 Theta^2 / 36  ( = theta_+ theta_- )."""
    if Thp is None:
        Thp = ddr(p["Theta"], p["r"])
    Th = p["Theta"]
    h = Thp ** 2 - (3.0 * p["Q"] - 2.0 * p["trK"]) ** 2 * Th ** 2 / 36.0
    return 0.5 * Th * (1.0 - h)


def w_out(p, cp):
    """Outgoing spin-0 characteristic  w_out = -(b_N/v0) trK + a_r   [eq:wout]
    (the ingoing one, w_in = +(b_N/v0) trK + a_r, is not used by the scripts)."""
    return -(cp["bN"] / cp["v0"]) * p["trK"] + p["Acr"]


def m_coulomb(m, Theta, a_r, cp):
    """m_c = m - (c14/2) Theta^2 a_r   [eq:mc]: the aether Coulomb term removed."""
    return m - 0.5 * cp["c14"] * Theta ** 2 * a_r


def tortoise(Theta, M):
    """Schwarzschild tortoise radius Theta_* = Theta + 2M ln(Theta/(2M) - 1)."""
    Theta = np.asarray(Theta, float)
    return Theta + 2.0 * M * np.log(Theta / (2.0 * M) - 1.0)


def retarded_time(t, Theta, M, v0):
    """u = t - Theta_*/v0   [eq:retarded]; broadcasts t (times) against Theta."""
    return np.asarray(t, float) - tortoise(Theta, M) / v0


def estimate_madm(p, frac=(0.6, 0.9)):
    """ADM mass read where m(Theta) is flat at large radius: median of m over
    Theta in [frac0, frac1] * Theta_max on an EARLY dump (before any wave arrives)."""
    m = misner_sharp(p)
    Th = p["Theta"]
    s = (Th > frac[0] * Th.max()) & (Th < frac[1] * Th.max())
    return float(np.median(m[s]))


# ----------------------------------------------------------------------------------------
# stress tensors in the aether frame
# ----------------------------------------------------------------------------------------
def aether_set(p, cp, W, B, Qt, beta=None, d=None, deriv=None):
    """Aether stress tensor projected on the aether frame: (rho_ae, j_ae, p_ae).

    Covariant SET of arXiv:1501.07274 eq. (2.34) (mostly-plus, u.u = -1, kappa = 1:
    G_ab = T^ae_ab + T^psi_ab), with the divergence term in the index order of the arXiv
    LaTeX source (J_{(a}^mu u_b) - J^mu_{(a} u_b) - J_(ab) u^mu).  Reduced symbolically to
    the variables of the code with SymPy (validated: rho
    reproduces the aether terms of GEJ's Hamiltonian constraint, the (u,e) projection the
    momentum constraint, and F_ae + F_psi = -dm/dt|_Theta).  The time derivatives of the
    EVOLVED variables (Theta, psi, P, K, a_r) were eliminated with the evolution equations;
    those of the gauge/algebraic quantities N, beta, Q cannot be, and enter as

        W  = d_t ln N |_r,   B = d_t beta |_r,   Qt = d_t Q |_r     (fixed proper radius r)

    which the caller supplies (finite differences across dumps).  c3 = c4 = 0 assumed.

        rho = c14 (a' + a^2/2 + 2 a Theta'/Theta) - (c13 + 3 c2) K^2/6 - (3/4) c13 Q^2
        j   = -T_us,  T_us = c14 [ (W' - beta a')/N + (Q - 2K/3) a ]
        p   = T_ss = (2 c2/(3 D)) [ (2 - c14)(a' + a^2 + 2 a Theta'/Theta) - 2 P^2 ]
                     - 2 c2 (1 - c13) Q^2 / D  -  (c1 + c2) B'/N^2
                     + (6 c1 K^2 + 36 c1 K Q - 27 c1 Q^2 - 18 c1 a^2 + 10 c2 K^2)/36
                     + [ (c1 + c2)(K W + beta K') + 3 c1 (Q W + beta Q')
                         + 3 c2 (Q W + Qt) ] / (3 N),          D = 2 + c13 + 3 c2.

    Arguments: p = profile dict; W, B, Qt arrays on p["r"]; beta (optional) the shift;
    d (optional) dict of precomputed r-derivatives {Theta, Acr, trK, Q}; deriv (optional)
    the d/dr used for the others (default ddr; flux_identity_rec.py passes its row-wise
    version for the (records, stencil) arrays of the sphere recorder).
    """
    if abs(cp["c3"]) > 0 or abs(cp["c4"]) > 0:
        raise NotImplementedError("aether_set is reduced for c3 = c4 = 0")
    c1, c2, c13, c14 = cp["c1"], cp["c2"], cp["c13"], cp["c14"]
    r = p["r"]
    Th, K, Q, a, N, P = p["Theta"], p["trK"], p["Q"], p["Acr"], p["NN"], p["P"]
    if beta is None:
        beta = beta_shift(p)
    if d is None:
        d = {}
    if deriv is None:
        deriv = ddr
    Thp = d.get("Theta", None)
    Thp = deriv(Th, r) if Thp is None else Thp
    ap = d.get("Acr", None)
    ap = deriv(a, r) if ap is None else ap
    Kp = d.get("trK", None)
    Kp = deriv(K, r) if Kp is None else Kp
    Qp = d.get("Q", None)
    Qp = deriv(Q, r) if Qp is None else Qp
    Wp = deriv(W, r)
    Bp = deriv(B, r)
    D = 2.0 + c13 + 3.0 * c2
    with np.errstate(divide="ignore", invalid="ignore"):
        Lr = a * Thp / Th
        rho = c14 * (ap + 0.5 * a ** 2 + 2.0 * Lr) - (c13 + 3.0 * c2) * K ** 2 / 6.0 \
            - 0.75 * c13 * Q ** 2
        Tus = c14 * ((Wp - beta * ap) / N + (Q - 2.0 * K / 3.0) * a)
        pr = (2.0 * c2 / (3.0 * D)) * ((2.0 - c14) * (ap + a ** 2 + 2.0 * Lr) - 2.0 * P ** 2) \
            - 2.0 * c2 * (1.0 - c13) * Q ** 2 / D \
            - (c1 + c2) * Bp / N ** 2 \
            + (6.0 * c1 * K ** 2 + 36.0 * c1 * K * Q - 27.0 * c1 * Q ** 2
               - 18.0 * c1 * a ** 2 + 10.0 * c2 * K ** 2) / 36.0 \
            + ((c1 + c2) * (K * W + beta * Kp) + 3.0 * c1 * (Q * W + beta * Qp)
               + 3.0 * c2 * (Q * W + Qt)) / (3.0 * N)
    return rho, -Tus, pr


def scalar_set(p, psip=None):
    """Scalar stress tensor in the aether frame: rho = p = (P^2 + psi'^2)/2, j = -P psi'."""
    if psip is None:
        psip = ddr(p["Psi"], p["r"])
    rho = 0.5 * (p["P"] ** 2 + psip ** 2)
    return rho, -p["P"] * psip, rho


def sphere_flux(p, rho, j, pr, Thp=None):
    """Energy flux through a sphere of fixed areal radius, paper eq:lum:
    F = (N Theta^2/2)[A (rho + p) + ((A^2 + B^2)/B) j],  A = Theta(Q/2 - K/3), B = Theta'.
    F > 0 outward; L = -dm/dt|_Theta = F_ae + F_psi."""
    if Thp is None:
        Thp = ddr(p["Theta"], p["r"])
    Th = p["Theta"]
    A = Th * (0.5 * p["Q"] - p["trK"] / 3.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        return 0.5 * p["NN"] * Th ** 2 * (A * (rho + pr) + (A ** 2 + Thp ** 2) / Thp * j)


def fixed_r_time_derivs(profs, times, i_c, fields=("lnN", "beta", "Q")):
    """d/dt at FIXED proper radius r of lnN, beta, Q on the grid of dump i_c.

    profs: list of profile dicts on consecutive dumps (stencil), times: their times,
    i_c: position of the central dump in that list.  On the moving mesh the r grids of
    the dumps differ, so each neighbour is first interpolated (cubic) onto the r grid of
    the central dump; then Fornberg weights on the actual times give the derivative
    (5 dumps -> fourth order, 3 -> second order)."""
    rc = profs[i_c]["r"]
    w = fd_weights(times[i_c], times, 1)
    out = {k: np.zeros_like(rc) for k in fields}
    for k, (p, wk) in enumerate(zip(profs, w)):
        vals = {}
        if "lnN" in fields:
            vals["lnN"] = np.log(np.maximum(p["NN"], 1e-300))
        if "beta" in fields:
            vals["beta"] = beta_shift(p)
        if "Q" in fields:
            vals["Q"] = p["Q"]
        same = (len(p["r"]) == len(rc)) and np.array_equal(p["r"], rc)
        for f in fields:
            v = vals[f] if same else lagrange_interp(p["r"], vals[f], rc)
            out[f] = out[f] + wk * v
    return out


# ----------------------------------------------------------------------------------------
# looping over dumps (optionally in a few processes)
# ----------------------------------------------------------------------------------------
_TASK = {}


def _chunk(bounds):
    lo, hi = bounds
    fn, paths, times, win = _TASK["fn"], _TASK["paths"], _TASK["times"], _TASK["win"]
    n = len(paths)
    cache = {}
    res = []
    for i in range(lo, hi):
        need = range(max(0, i - win), min(n, i + win + 1))
        for k in list(cache):
            if k not in need:
                del cache[k]
        for k in need:
            if k not in cache:
                cache[k] = load_profile(paths[k])
        res.append(fn(i, cache, times))
    return res


def map_dumps(fn, times, paths, nproc=1, window=0):
    """[fn(i, cache, times) for every dump i], in order.

    cache maps dump index -> profile dict for the indices i-window .. i+window (clipped);
    fn must not modify the profiles.  nproc > 1 splits the dumps in contiguous chunks over
    forked worker processes (each single-threaded)."""
    n = len(paths)
    _TASK.update(fn=fn, paths=paths, times=np.asarray(times), win=window)
    if nproc <= 1 or n < 4:
        return _chunk((0, n))
    nch = min(n, 4 * nproc)
    edges = np.linspace(0, n, nch + 1).astype(int)
    ctx = mp.get_context("fork")
    with ctx.Pool(nproc) as pool:
        parts = pool.map(_chunk, list(zip(edges[:-1], edges[1:])))
    return [x for part in parts for x in part]


def sphere_series(rundir, thetas, quantities, cp=None, tmin=-np.inf, tmax=np.inf,
                  nproc=1):
    """Time series of local quantities on spheres of fixed areal radius.

    quantities: function prof -> dict name -> array on the r grid (or a list of names among
    'm', 'wout', 'Acr', 'trK', 'Q', 'NN', 'Psi', 'P', 'mc', 'Fpsi' computed here).
    Returns (t, dict name -> array (n_t, n_theta)).

    If the environment variable SPIN0_CACHE names a directory, results for a list of names
    are cached there (.npz keyed by run directory, spheres, names, couplings and the dump
    list), so that several scripts can share one pass over the dumps."""
    times, paths = list_dumps(rundir, tmin, tmax)
    thetas = np.atleast_1d(np.asarray(thetas, float))
    cfile = None
    if os.environ.get("SPIN0_CACHE") and not callable(quantities):
        import hashlib
        key = repr((os.path.abspath(rundir), [round(x, 6) for x in thetas], list(quantities),
                    None if cp is None else (cp["c1"], cp["c2"]), len(paths), paths[0],
                    paths[-1], os.path.getmtime(paths[-1])))
        cdir = os.environ["SPIN0_CACHE"]
        os.makedirs(cdir, exist_ok=True)
        cfile = os.path.join(cdir, "series_" + hashlib.md5(key.encode()).hexdigest() + ".npz")
        if os.path.exists(cfile):
            z = np.load(cfile)
            return z["t"], {k: z[k] for k in z.files if k != "t"}
    if callable(quantities):
        qfun = quantities
    else:
        names = list(quantities)

        def qfun(p):
            Thp = ddr(p["Theta"], p["r"])
            out = {}
            m = misner_sharp(p, Thp) if ("m" in names or "mc" in names) else None
            for k in names:
                if k == "m":
                    out[k] = m
                elif k == "mc":
                    out[k] = m_coulomb(m, p["Theta"], p["Acr"], cp)
                elif k == "wout":
                    out[k] = w_out(p, cp)
                elif k == "Fpsi":
                    out[k] = sphere_flux(p, *scalar_set(p), Thp=Thp)
                else:
                    out[k] = p[k]
            return out

    def fn(i, cache, _t):
        p = cache[i]
        return at_theta(p["Theta"], thetas, qfun(p))

    rows = map_dumps(fn, times, paths, nproc=nproc)
    keys = list(rows[0])
    out = {k: np.array([row[k] for row in rows]) for k in keys}
    if cfile is not None:
        np.savez(cfile, t=times, **out)
    return times, out


def gaussian_running_mean(t, y, sigma, truncate=3.0):
    """Gaussian running mean of y(t) (axis 0) of width sigma, on the dump times t.

    The kernel is renormalised where it is cut by the ends of the series; the returned mask
    `full` is True where the window [t - truncate sigma, t + truncate sigma] lies inside
    the data at the late end (at the early end the signal is ~0 before the wave arrives,
    so the renormalised mean is used)."""
    t = np.asarray(t, float)
    y = np.asarray(y, float)
    out = np.empty_like(y)
    dt = np.gradient(t)
    for i in range(len(t)):
        s = np.abs(t - t[i]) <= truncate * sigma
        w = np.exp(-0.5 * ((t[s] - t[i]) / sigma) ** 2) * dt[s]
        out[i] = np.tensordot(w, y[s], axes=(0, 0)) / w.sum()
    full = t <= t[-1] - truncate * sigma
    return out, full


def parabolic_peak(x, y, j):
    """Sub-sample location and value of an extremum at index j (quadratic through 3 points,
    any spacing)."""
    if j <= 0 or j >= len(x) - 1:
        return x[j], y[j]
    xs, ys = x[j - 1:j + 2], y[j - 1:j + 2]
    c = np.polyfit(xs - xs[1], ys, 2)
    if c[0] == 0:
        return x[j], y[j]
    xm = -c[1] / (2 * c[0])
    if abs(xm) > max(abs(xs[0] - xs[1]), abs(xs[2] - xs[1])):
        return x[j], y[j]
    return xs[1] + xm, np.polyval(c, xm)


# ----------------------------------------------------------------------------------------
# helpers shared by the convergence, threshold, mass-function, moving-mesh and white-hole
# scripts
# ----------------------------------------------------------------------------------------
def norm_l2(v, r, lo, hi, min_points=0):
    """Discrete L2 norm sqrt(mean v^2) over the grid points with lo <= r <= hi; NaN if
    fewer than min_points points lie in the window."""
    m = (r >= lo) & (r <= hi)
    if m.sum() < min_points:
        return np.nan
    return np.sqrt(np.mean(v[m] ** 2))


def run_amplitude(d):
    """Shell amplitude A0 of a run directory: from --A0= in <d>/cmd.txt (an optional copy
    of the command line, saved by the user; collapse does not write it), otherwise from
    the directory name ("...A<value>...", e.g. bis_N12801_A0.130); NaN if neither."""
    cmd = os.path.join(d, "cmd.txt")
    if os.path.exists(cmd):
        m = re.search(r"--A0=([0-9.eE+-]+)", open(cmd).read())
        if m:
            return float(m.group(1))
    m = re.search(r"A(\d+\.\d+)", os.path.basename(d.rstrip("/")))
    return float(m.group(1)) if m else np.nan


def apparent_horizon(p, rmin=0.5, interpolate=True):
    """Outermost zero crossing of theta_+ at r > rmin -> (index i, r_i, Theta_AH).

    theta_+ = Theta' + (3Q - 2K) Theta/6 with Theta' = ddr(Theta, r) (second order, also at
    the end points); Theta_AH is the zero of the linear interpolant between the points i
    and i+1 (or Theta_i if interpolate=False).
    (-1, nan, nan) when theta_+ does not change sign (no horizon).  The cut is on the
    proper radius r, not on Theta: callers that want a detached horizon apply their own
    Theta cut (see analysis/README.md, "Horizon criteria")."""
    r, Th = p["r"], p["Theta"]
    tp, _ = theta_pm(p, ddr(Th, r))
    sg = np.sign(tp)
    z = np.where((sg[:-1] * sg[1:] < 0) & (r[:-1] > rmin))[0]
    if len(z) == 0:
        return -1, np.nan, np.nan
    i = int(z[-1])
    x = Th[i] + (Th[i + 1] - Th[i]) * tp[i] / (tp[i] - tp[i + 1]) if interpolate else Th[i]
    return i, r[i], x


def fwhm(r, f, j):
    """Full width at half maximum of |f| around the peak at index j (walk outward until
    |f| drops below half the peak) -> (width in r, number of points)."""
    h = abs(f[j])
    a = j
    while a > 0 and abs(f[a - 1]) >= h / 2:
        a -= 1
    b = j
    while b < len(f) - 1 and abs(f[b + 1]) >= h / 2:
        b += 1
    return r[b] - r[a], b - a + 1


def peak_index(p, field="P", rmin=3.0, thmin=-np.inf, thmax=35.0):
    """Index of the maximum of |field| among the points with r > rmin, thmin < Theta < thmax.

    The defaults select the inner shell of the white-hole runs (r > 3 skips the centre,
    Theta < 35 the outgoing front of the r_max = 40 box)."""
    k = np.where((p["r"] > rmin) & (p["Theta"] > thmin) & (p["Theta"] < thmax))[0]
    return int(k[np.argmax(np.abs(p[field][k]))])


def steepest_slope(p, lo=1.5, hi=5.0, field="P"):
    """Steepest |d field/d Theta| (one-sided differences between neighbouring points)
    among the points with lo < Theta < hi -> (index j of the left point, slope).  With
    the defaults: the inner ingoing front of the white-hole runs."""
    Th = p["Theta"]
    i = np.where((Th > lo) & (Th < hi))[0][:-1]
    g = np.abs(np.diff(p[field]))[i] / np.diff(Th)[i]
    return int(i[np.argmax(g)]), g.max()


def add_common_args(ap, c2=True, madm=True, out=None):
    """Arguments shared by the spin-0 scripts: RUNDIR, --c1 [, --c2] [, --madm], --nproc,
    --out (default `out`, the figure name of the script)."""
    ap.add_argument("rundir", help="run directory with prof_t*.csv dumps")
    ap.add_argument("--c1", type=float, default=C1_DEFAULT, help="coupling c1 (default 0.7)")
    if c2:
        ap.add_argument("--c2", type=float, default=C2_BASELINE,
                        help="c2 (-0.5119403: v0 = 1;  0: v0 = 1.267)")
    if madm:
        ap.add_argument("--madm", type=float, default=None,
                        help="ADM mass (default: estimated from m(Theta) at large radius "
                             "on the first dump)")
    ap.add_argument("--nproc", type=int, default=4, help="worker processes (default 4)")
    ap.add_argument("--out", default=out,
                    help=f"output figure (png; a pdf is written next to it; default {out})")
    return ap


def resolve_madm(args):
    """M_ADM: --madm if given, otherwise estimated from the first dump (estimate_madm)."""
    if getattr(args, "madm", None) is not None:
        return args.madm
    t, paths = list_dumps(args.rundir)
    M = estimate_madm(load_profile(paths[0]))
    print(f"[M_ADM estimated from the first dump (t={t[0]:.2f}): {M:.4f}]")
    return M
