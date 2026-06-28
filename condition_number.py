#!/usr/bin/env sage
"""
condition_number.py  —  Step 2: condition number growth of d_0^out
=================================================================
Computes, for each internal degree d = 0,...,max_deg, the condition number
of the outgoing differential d_0^out restricted to degree-d components.

This determines whether the formal preimage \psi in Proposition prop:bridge(ii)
(Paper 1) inherits geometric coefficient growth from \phi, which is required for
the Cauchy-Kovalevskaya argument to extend from real-analytic to Gevrey-class
initial data (Paper 2 / Open Problem A).

Method
------
At cochain level k=0, d_0^out is a block matrix whose (i,j)-block consists of
the degree-d slice of morphism \phi_* restricted to source-degree d.

For each d we extract the degree-d block of d_0^out as a QQ matrix and compute:
  - Exact rank over QQ (via p-adic lifting)
  - Rank over GF(p) for sanity check
  - Numerical singular values over RR (via SVD on the RDF approximation)
  - Condition number κ(d) = \sigma_max / \sigma_min (with \sigma_min taken over nonzero values)
  - Ratio κ(d)/\kappa(d-1) to assess growth rate

Growth interpretation
---------------------
  \kappa(d) ~ C           : constant — analytic closure holds trivially
  \kappa(d) ~ d^\alpha         : polynomial (\alpha > 0) — Gevrey-1 closure, explicit radius
  \kappa(d) ~ C^d         : exponential — closure requires Borel summation with s > 1
  \kappa(d) grows faster  : no straightforward Gevrey extension via this approach

Output
------
  condition_number_output.txt   — full table with singular values, ranks, \kappa(d)
  Printed to stdout during run.

Usage (inside SageMath conda env)
----------------------------------
  conda run -n sage sage condition_number.py > condition_number_output.txt 2>&1

Parameters
----------
  T_MIN, T_MAX, A_MAX  — cochain window (can be small; we only need k=0)
  MAX_DEG              — Verma truncation degree (5 to match paper)
  K_LEVEL              — cochain level at which to extract d_0^out (default 0)
  CX_TYPE              — 'A' or 'B'; run both for completeness

Run-time notes
--------------
  At MAX_DEG=5 and A_MAX=4, the degree-5 block can be large (~tens of thousands
  of rows/cols).  SVD over RDF is O(n²m) which is feasible.  The exact rank
  computation uses p-adic lifting and is exact.  Expect minutes per degree level.
"""

from __future__ import print_function
import sys, os, time, datetime

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

from sage.all import QQ, GF, RDF, matrix, vector, floor, log, sqrt

from verma_modules import load_e44
from de_rham_complex import (
    DeRhamComplexA, DeRhamComplexB,
    CochainGroup, window_nodes,
    MORPHISMS_A, MORPHISMS_B,
    assemble_differential,
    DEFAULT_T_MIN, DEFAULT_T_MAX, DEFAULT_A_MAX, MAX_VERMA_DEG,
    _VERMA_CACHE,
)

# ===========================================================================
# Parameters
# ===========================================================================

T_MIN   = -2      # small window suffices; k=0 is interior
T_MAX   =  2
A_MAX   =  4      # match full production parameters
MAX_DEG =  5      # must match paper (Prop 5.5: sv_deg <= 5)
K_LEVEL =  0      # extract d_0^out at cochain level k=0
_P1, _P2 = 65521, 65537   # primes for GF rank certification

# ===========================================================================
# Utilities
# ===========================================================================

def ts():
    return datetime.datetime.now().strftime('%H:%M:%S')


def rank_certified(M, label=""):
    """
    Compute rank over QQ, certified by two GF primes.
    Falls back to exact padic if primes disagree.
    """
    if M.nrows() == 0 or M.ncols() == 0:
        return 0
    r1 = M.change_ring(GF(_P1)).rank()
    r2 = M.change_ring(GF(_P2)).rank()
    if r1 != r2:
        print(f'  {label}: GF rank mismatch ({r1} vs {r2}), using padic', flush=True)
        return M.rank(algorithm='padic')
    return r1


def condition_number_rdf(M):
    """
    Numerical condition number via SVD over RDF.
    Returns (sigma_max, sigma_min_nonzero, kappa, rank_numerical).
    sigma_min_nonzero is the smallest nonzero singular value.
    """
    if M.nrows() == 0 or M.ncols() == 0:
        return (0.0, 0.0, float('inf'), 0)
    M_rdf = M.change_ring(RDF)
    # SVD: M = U * S * V^T, S diagonal with singular values
    try:
        _, S_matrix, _ = M_rdf.SVD()
    except Exception as e:
        return (0.0, 0.0, float('nan'), 0)

    n_sv = min(M.nrows(), M.ncols())
    sigmas = sorted([abs(S_matrix[i, i]) for i in range(n_sv)], reverse=True)

    if not sigmas or sigmas[0] < 1e-14:
        return (0.0, 0.0, float('inf'), 0)

    sigma_max = sigmas[0]
    # Smallest nonzero (use threshold relative to sigma_max)
    threshold = sigma_max * 1e-10
    nonzero = [s for s in sigmas if s > threshold]
    if not nonzero:
        return (sigma_max, 0.0, float('inf'), 0)

    sigma_min = nonzero[-1]
    kappa = sigma_max / sigma_min if sigma_min > 0 else float('inf')
    rank_num = len(nonzero)
    return (sigma_max, sigma_min, kappa, rank_num)


# ===========================================================================
# Extract degree-d slice of the outgoing matrix
# ===========================================================================

def degree_slice(D, g_src, g_tar, d):
    """
    Extract the columns of D corresponding to source degree d
    and the rows corresponding to target degree sv_deg+d (for all blocks).

    D   : sparse QQ matrix, shape (tar_dim, src_dim)
    g_src, g_tar : CochainGroup objects
    d   : int — source internal degree

    Returns a (sub_tar_rows x sub_src_cols) QQ matrix, or None if empty.
    """
    # Column indices in D corresponding to source degree d
    src_cols = []
    for node in g_src.nodes:
        sl = g_src.basis_slice(node, d)
        src_cols.extend(sl)

    if not src_cols:
        return None

    # Row indices: all rows (we keep the full target slice for the columns)
    # This gives the sub-matrix of D restricted to source-degree-d columns.
    sub = D.matrix_from_columns(src_cols)
    return sub


def degree_block_of_d0out(cx_type, k, d, e44_data, t_min, t_max, a_max, max_deg):
    """
    Assemble d_0^out at cochain level k, then extract the degree-d column block.
    Returns (D_full, D_deg_d, g_src).
    """
    morphisms = MORPHISMS_A if cx_type == 'A' else MORPHISMS_B

    nodes_k = window_nodes(k, morphisms, t_min, t_max, a_max)
    g_src = CochainGroup(k, nodes_k, max_deg=max_deg, e44_data=e44_data)

    # Collect all target levels (k + shift for each morphism)
    shifts = (1, 2) if cx_type == 'A' else (1, 3, 4)
    all_blocks = []
    total_tar_rows = 0
    for shift in shifts:
        k_tar = k + shift
        if k_tar > t_max:
            continue
        nodes_ktar = window_nodes(k_tar, morphisms, t_min, t_max, a_max)
        g_tar = CochainGroup(k_tar, nodes_ktar, max_deg=max_deg, e44_data=e44_data)
        D = assemble_differential(g_src, g_tar, morphisms, e44_data,
                                  t_min, t_max, a_max)
        all_blocks.append((g_tar, D))
        total_tar_rows += D.nrows()

    if not all_blocks:
        return None, None, g_src

    src_dim = g_src.total_dim
    D_full = matrix(QQ, total_tar_rows, src_dim, sparse=True)
    r = 0
    for g_tar, D in all_blocks:
        for (i, j), val in D.dict().items():
            D_full[r + i, j] = val
        r += D.nrows()

    # Degree-d column slice
    src_cols = []
    for node in g_src.nodes:
        sl = g_src.basis_slice(node, d)
        src_cols.extend(sl)

    if not src_cols:
        return D_full, None, g_src

    D_deg_d = D_full.matrix_from_columns(src_cols)
    return D_full, D_deg_d, g_src


# ===========================================================================
# Main computation
# ===========================================================================

def run(cx_type, e44_data):
    print(f'\n{"="*70}')
    print(f'Condition number analysis — Complex {cx_type}')
    print(f'  Parameters: t=[{T_MIN},{T_MAX}], a_max={A_MAX}, '
          f'max_deg={MAX_DEG}, k={K_LEVEL}')
    print(f'{"="*70}\n')

    results = []

    for d in range(0, MAX_DEG + 1):
        print(f'[{ts()}] d={d}: assembling degree-{d} block of d_0^out ...', flush=True)
        t0 = time.time()

        _, D_d, g_src = degree_block_of_d0out(
            cx_type, K_LEVEL, d, e44_data,
            T_MIN, T_MAX, A_MAX, MAX_DEG
        )
        t_assemble = time.time() - t0

        if D_d is None or D_d.ncols() == 0:
            print(f'  d={d}: no source basis at this degree, skipping.')
            results.append({
                'd': d, 'nrows': 0, 'ncols': 0,
                'rank_qq': 0, 'sigma_max': 0.0, 'sigma_min': 0.0,
                'kappa': float('nan'), 'rank_num': 0,
            })
            _VERMA_CACHE.clear()
            continue

        nrows, ncols = D_d.nrows(), D_d.ncols()
        print(f'  d={d}: block shape {nrows}×{ncols} '
              f'(assembled in {t_assemble:.1f}s)', flush=True)

        # Exact rank
        print(f'  d={d}: computing exact rank ...', flush=True)
        t0 = time.time()
        r_qq = rank_certified(D_d, f'd={d}')
        print(f'  d={d}: rank_qq={r_qq} ({time.time()-t0:.1f}s)', flush=True)

        # Numerical condition number via SVD
        print(f'  d={d}: computing SVD over RDF ...', flush=True)
        t0 = time.time()
        sigma_max, sigma_min, kappa, rank_num = condition_number_rdf(D_d)
        print(f'  d={d}: \sigma_max={sigma_max:.4e}, \sigma_min={sigma_min:.4e}, '
              f'\kappa={kappa:.4e}, rank_num={rank_num} ({time.time()-t0:.1f}s)',
              flush=True)

        if r_qq != rank_num:
            print(f'  d={d}: WARNING: exact rank {r_qq} != numerical rank {rank_num} '
                  f'(ill-conditioning or SVD tolerance issue)', flush=True)

        # Free Verma module cache before next degree to prevent accumulation.
        _VERMA_CACHE.clear()

        results.append({
            'd': d, 'nrows': nrows, 'ncols': ncols,
            'rank_qq': r_qq, 'sigma_max': sigma_max, 'sigma_min': sigma_min,
            'kappa': kappa, 'rank_num': rank_num,
        })

    # ── Print summary table ──────────────────────────────────────────────
    print(f'\n{"="*70}')
    print(f'Condition number table — Complex {cx_type}')
    print(f'{"="*70}')
    print(f'  {"d":>3}  {"shape":>14}  {"rank":>6}  '
          f'{"\sigma_max":>12}  {"\sigma_min":>12}  {"\kappa(d)":>12}  {"\kappa ratio":>10}')
    print(f'  {"-"*3}  {"-"*14}  {"-"*6}  {"-"*12}  {"-"*12}  {"-"*12}  {"-"*10}')

    prev_kappa = None
    for r in results:
        d = r['d']
        shape = f'{r["nrows"]}×{r["ncols"]}'
        rank  = r['rank_qq']
        smax  = f'{r["sigma_max"]:.4e}'
        smin  = f'{r["sigma_min"]:.4e}'
        kappa_str = f'{r["kappa"]:.4e}' if r['kappa'] != float('nan') else 'N/A'
        if prev_kappa is not None and prev_kappa > 0 and r['kappa'] not in (float('nan'), float('inf')):
            ratio = r['kappa'] / prev_kappa
            ratio_str = f'{ratio:.4f}'
        else:
            ratio_str = '—'
        print(f'  {d:>3}  {shape:>14}  {rank:>6}  '
              f'{smax:>12}  {smin:>12}  {kappa_str:>12}  {ratio_str:>10}')
        if r['kappa'] not in (float('nan'), float('inf')):
            prev_kappa = r['kappa']

    # ── Growth assessment ────────────────────────────────────────────────
    print(f'\n{"="*70}')
    print('Growth assessment')
    print(f'{"="*70}')
    kappas = [(r['d'], r['kappa']) for r in results
              if r['kappa'] not in (float('nan'), float('inf')) and r['kappa'] > 0
              and r['ncols'] > 0]

    if len(kappas) < 2:
        print('  Insufficient data for growth assessment.')
    else:
        ds     = [float(d) for d, _ in kappas]
        logk   = [float(log(k)) for _, k in kappas]
        # Fit log \kappa(d) ~ \alpha·log(d) (polynomial) and log \kappa(d) ~ β·d (exponential)
        # via least squares using the last half of the data points.
        mid = len(ds) // 2
        ds_fit   = ds[mid:]
        logk_fit = logk[mid:]
        logd_fit = [log(RDF(d)) if d > 0 else RDF(0) for d in ds_fit]

        # Polynomial fit: log κ ≈ \alpha·log(d) + c
        if len(ds_fit) >= 2 and all(d > 0 for d in ds_fit):
            A_poly = matrix(RDF, [[float(ld), 1.0] for ld in logd_fit])
            b_poly = vector(RDF, logk_fit)
            try:
                coeffs_poly = (A_poly.T * A_poly).solve_right(A_poly.T * b_poly)
                alpha = float(coeffs_poly[0])
                print(f'  Polynomial fit:    log κ(d) ≈ {alpha:.3f}·log(d) + c')
            except Exception:
                alpha = None
        else:
            alpha = None

        # Exponential fit: log κ ≈ β·d + c
        A_exp = matrix(RDF, [[float(d), 1.0] for d in ds_fit])
        b_exp = vector(RDF, logk_fit)
        try:
            coeffs_exp = (A_exp.T * A_exp).solve_right(A_exp.T * b_exp)
            beta = float(coeffs_exp[0])
            print(f'  Exponential fit:   log κ(d) ≈ {beta:.3f}·d + c')
        except Exception:
            beta = None

        print()
        if beta is not None and alpha is not None:
            if beta < 0.05:
                verdict = ('POLYNOMIAL or CONSTANT growth  →  '
                           'Gevrey-1 closure holds; '
                           'Paper 2 (Gevrey extension) is accessible.')
            elif beta < 0.5:
                verdict = ('MILD EXPONENTIAL growth  →  '
                           'Gevrey-s closure for small s; '
                           'Borel summation argument needed for Paper 2.')
            else:
                verdict = ('FAST EXPONENTIAL growth  →  '
                           'Gevrey approach may be difficult; '
                           'Sobolev-space route (Open Problem B) preferred.')
            print(f'  Verdict: {verdict}')
        print()

    print(f'[{ts()}] Done.')
    return results


# ===========================================================================
# Entry point
# ===========================================================================

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
        description='Condition number growth of d_0^out (Step 2)')
    parser.add_argument('--cx', choices=['A', 'B', 'both'], default='both',
                        help='Which complex to analyse (default: both)')
    parser.add_argument('--max-deg', type=int, default=MAX_DEG)
    parser.add_argument('--a-max',   type=int, default=A_MAX)
    parser.add_argument('--t-min',   type=int, default=T_MIN)
    parser.add_argument('--t-max',   type=int, default=T_MAX)
    args = parser.parse_args()

    MAX_DEG = args.max_deg
    A_MAX   = args.a_max
    T_MIN   = args.t_min
    T_MAX   = args.t_max

    print(f'[{ts()}] Loading e44_brackets.pkl ...', flush=True)
    e44_data = load_e44()
    print(f'[{ts()}] Loaded.\n', flush=True)

    cx_list = ['A', 'B'] if args.cx == 'both' else [args.cx]
    all_results = {}
    for cx in cx_list:
        all_results[cx] = run(cx, e44_data)

    print('\nAll done.')
