#!/usr/bin/env sage
"""
cohomology_h1_fast.py  --  Fast H^1 computation via n-decomposition
===================================================================
Computes H^{K_LEVEL} for Complex A and/or B by exploiting two reductions:

  1. Tight cochain window
     For H^1, we need only the levels that have edges incident to k=K_LEVEL:
       Complex A (sv_deg in {1,2}): levels K_LEVEL-2 .. K_LEVEL+2 = -1..3
       Complex B (sv_deg in {1,3,4}): levels K_LEVEL-4 .. K_LEVEL+4 = -3..5
     (In production_run.py the full window [-6,6] = 13 levels is used.)

  2. n-invariant decomposition
     Every morphism phi of sv_deg s maps M_src[d] at cochain level k to
     M_tar[d+s] at cochain level k+s, so the bigrading invariant

         n  =  k - d    (cochain level minus PBW degree)

     is preserved by every differential component.  The cochain complex
     decomposes as a direct sum indexed by n, and H^{K_LEVEL} decomposes
     accordingly:

         H^{K_LEVEL} = oplus_n  H^{K_LEVEL}_n

     where in each n-slice at cochain level k, only the sub-basis

         C^k_n  =  oplus_{node at level k}  M_node[ k - n ]

     participates.  For K_LEVEL=1 and max_deg=5, the n-slices at level 1
     have dimensions (approximately):

         n=+1 (d=0): ~300 rows  -- trivial
         n= 0 (d=1): ~2.2K
         n=-1 (d=2): ~9K
         n=-2 (d=3): ~25K
         n=-3 (d=4): ~54K
         n=-4 (d=5): ~100K  -- largest, but GF(p)-feasible

     compared with the full un-decomposed dimension of ~191K.

  3. Minimal phi calls
     For each n-slice, assemble_n_slice() calls spec.phi_func with
     max_source_deg = d_src = k_src - n (the one PBW degree needed), so
     _phi_at_degree is called only ONCE per morphism edge (rather than
     once per degree 0..max_src_deg).

Usage
-----
  conda run -n sage sage cohomology_h1_fast.py [--cx A|B|both] [options]
  conda run -n sage sage cohomology_h1_fast.py --cx A > h1_A.txt 2>&1
  conda run -n sage sage cohomology_h1_fast.py --cx B > h1_B.txt 2>&1

Parameters
----------
  --cx       : A, B, or both (default: both)
  --k        : cochain level for cohomology (default: 1)
  --max-deg  : Verma truncation (default: 5)
  --a-max    : a-label bound (default: 4)
  --t-min    : global window lower bound (default: -6)
  --t-max    : global window upper bound (default: +6)
"""

from __future__ import print_function
import sys, os, time, datetime, argparse, traceback

# Force unbuffered output even when sage wraps sys.stdout.
# sage script.py redirects sys.stdout to a Sage-internal stream whose buffer
# only flushes at process exit; this restores raw OS-level fd-1 output so
# that _log(...) calls reach the file/terminal immediately.
try:
    sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', buffering=1)
    sys.stderr = os.fdopen(sys.stderr.fileno(), 'w', buffering=1)
except Exception:
    pass  # already unbuffered or not a real fd (e.g. pytest)

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

from sage.all import QQ, GF, matrix

from verma_modules import load_e44
from de_rham_complex import (
    CochainGroup, window_nodes,
    MORPHISMS_A, MORPHISMS_B,
    _VERMA_CACHE,
)

# ---------------------------------------------------------------------------
# Progress logger -- writes directly to a file, bypassing Sage's stdout buffer.
# Sage's preparser replaces sys.stdout with a fully-buffered internal stream
# that only flushes at process exit.  Use _log() for all progress output so
# it is visible immediately in the log file.
# ---------------------------------------------------------------------------
_log_fh = None   # set in __main__

def _log(*args, **kwargs):
    """Print to stdout AND flush immediately to _log_fh if set."""
    msg = ' '.join(str(a) for a in args)
    print(msg)   # goes to Sage's buffered stdout (visible at process exit)
    if _log_fh is not None:
        try:
            _log_fh.write(msg + '\n')
            _log_fh.flush()
        except Exception:
            pass

# ---------------------------------------------------------------------------
# Default parameters
# ---------------------------------------------------------------------------
T_MIN   = -6
T_MAX   =  6
A_MAX   =  4
MAX_DEG =  5
K_LEVEL =  1       # which H^k to compute
_P1, _P2 = 65521, 65537   # primes for GF rank certification


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------

def ts():
    return datetime.datetime.now().strftime('%H:%M:%S')


def rank_certified(M, label=''):
    """Rank over QQ via two-prime GF certification."""
    if M.nrows() == 0 or M.ncols() == 0:
        return 0
    r1 = M.change_ring(GF(_P1)).rank()
    r2 = M.change_ring(GF(_P2)).rank()
    if r1 != r2:
        info = f' ({label})' if label else ''
        _log(f'  [WARN] GF rank mismatch{info}: {r1} vs {r2}; '
              f'falling back to padic')
        return M.rank(algorithm='padic')
    return r1


def tight_window_bounds(k_level, morphisms):
    """Return (t_min_tight, t_max_tight) -- the minimal window for H^{k_level}."""
    max_sv = max(spec.sv_deg for spec in morphisms)
    return k_level - max_sv, k_level + max_sv


# ---------------------------------------------------------------------------
# Pre-compute all phi matrices once (reused across n-slices)
# ---------------------------------------------------------------------------

def precompute_phi_matrices(morphisms, groups, e44_data,
                            t_min, t_max, a_max, max_deg):
    """
    Call each phi_func exactly once per edge at max_source_deg = max_deg - sv_deg,
    and cache the resulting mats dict.

    Without this, the n-decomposition loop would call phi_func once per n-slice
    (i.e. max_deg+1 times per edge), each time recomputing all degrees 0..d_src.
    Pre-computing reduces total phi_func calls from O(max_deg * n_edges) to n_edges.

    Returns
    -------
    dict mapping (spec_name, src_nd, tar_nd, phi_args) -> mats_dict or None.
    """
    cache = {}
    # Fiber cache: phi matrices depend only on (src fiber, tar fiber, sv_deg, msd),
    # NOT on t.  All edges phi_X(t, a, ...) with the same non-t params share the
    # same matrices.  Key: (spec.name, src_nd_no_t, tar_nd_no_t, msd).
    #
    # src_nd_no_t / tar_nd_no_t = (a, b, c) tuple (omitting t).
    fiber_cache = {}   # (spec_name, src_abc, tar_abc, msd) -> mats
    total = 0
    skipped = 0
    fiber_hits = 0
    for spec in morphisms:
        msd = max_deg - spec.sv_deg   # max source deg for this sv_deg
        if msd < 0:
            continue
        for src_nd, tar_nd, phi_args in spec.edges(t_min, t_max, a_max):
            if src_nd.t not in groups or tar_nd.t not in groups:
                continue
            if src_nd not in groups[src_nd.t].nodes:
                continue
            if tar_nd not in groups[tar_nd.t].nodes:
                continue
            key = (spec.name, src_nd, tar_nd, phi_args)
            if key in cache:
                continue
            total += 1
            # Check fiber cache first (t-independent matrices)
            src_abc = (src_nd.a, src_nd.b, src_nd.c)
            tar_abc = (tar_nd.a, tar_nd.b, tar_nd.c)
            fiber_key = (spec.name, src_abc, tar_abc, msd)
            if fiber_key in fiber_cache:
                cache[key] = fiber_cache[fiber_key]
                fiber_hits += 1
                _log(f'  -> edge {total}: {spec.name}{phi_args} msd={msd} '
                     f'[fiber cache hit]')
                continue
            _log(f'  -> edge {total}: {spec.name}{phi_args} msd={msd} ...')
            t0 = time.time()
            try:
                _, _, _sv, _phi0, mats = spec.phi_func(
                    *phi_args, e44_data, msd, src_e44_data=None
                )
                cache[key] = mats
                fiber_cache[fiber_key] = mats
            except BaseException as _exc:
                _log(f'  !! edge {total} FAILED: {type(_exc).__name__}: {_exc}')
                _log(traceback.format_exc())
                if isinstance(_exc, KeyboardInterrupt):
                    raise
                # SystemExit from Sage's sig_on() (SIGSEGV/OOM in GMP arithmetic)
                # is treated as a skippable edge failure, not a genuine sys.exit().
                cache[key] = None
                fiber_cache[fiber_key] = None
                skipped += 1
            elapsed = time.time() - t0
            if elapsed > 5.0:
                _log(f'    [{ts()}] {spec.name}{phi_args} msd={msd}: '
                      f'{elapsed:.1f}s')
            _VERMA_CACHE.clear()
    _log(f'  Pre-computed {total - skipped}/{total} edges '
         f'({skipped} skipped/failed, {fiber_hits} fiber-cache hits).')
    return cache


# ---------------------------------------------------------------------------
# n-slice assembly (uses pre-computed phi matrices)
# ---------------------------------------------------------------------------

def assemble_n_slice(k_src, k_tar, d_src, d_tar,
                     g_src, g_tar, morphisms, phi_cache,
                     t_min, t_max, a_max):
    """
    Assemble the n-slice of one differential component C^{k_src} to C^{k_tar}
    restricted to source PBW degree d_src and target PBW degree d_tar.

    n = k_src - d_src = k_tar - d_tar  (preserved by every morphism, since
    sv_deg = d_tar - d_src = k_tar - k_src).

    Phi matrices are looked up from phi_cache (pre-computed once via
    precompute_phi_matrices); no phi_func calls are made here.

    Returns
    -------
    Sparse QQ matrix of shape (tar_n_dim, src_n_dim) where
        src_n_dim = sum_{node in g_src} M_node[d_src]
        tar_n_dim = sum_{node in g_tar} M_node[d_tar]
    """
    sv = d_tar - d_src
    if sv != k_tar - k_src:
        raise ValueError(
            f'sv mismatch: d_tar-d_src={sv} vs k_tar-k_src={k_tar-k_src}')

    src_offs = {}
    src_dim = 0
    for nd in g_src.nodes:
        src_offs[nd] = src_dim
        src_dim += g_src.vermas[nd].dim(d_src)

    tar_offs = {}
    tar_dim = 0
    for nd in g_tar.nodes:
        tar_offs[nd] = tar_dim
        tar_dim += g_tar.vermas[nd].dim(d_tar)

    if src_dim == 0 or tar_dim == 0:
        return matrix(QQ, tar_dim, src_dim, sparse=True)

    entries = {}

    for spec in morphisms:
        if spec.sv_deg != sv:
            continue

        for src_nd, tar_nd, phi_args in spec.edges(t_min, t_max, a_max):
            if src_nd.t != k_src or tar_nd.t != k_tar:
                continue
            if src_nd not in g_src.nodes or tar_nd not in g_tar.nodes:
                continue

            cache_key = (spec.name, src_nd, tar_nd, phi_args)
            mats = phi_cache.get(cache_key)
            if mats is None:
                continue

            block = mats.get(d_src)
            if block is None:
                continue

            nrows_exp = g_tar.vermas[tar_nd].dim(d_tar)
            ncols_exp = g_src.vermas[src_nd].dim(d_src)
            if block.nrows() != nrows_exp or block.ncols() != ncols_exp:
                continue  # fiber mismatch -- skip

            r0 = tar_offs[tar_nd]
            c0 = src_offs[src_nd]
            for (ri, ci), val in block.dict().items():
                key = (r0 + ri, c0 + ci)
                prev = entries.get(key)
                entries[key] = val if prev is None else prev + val

    return matrix(QQ, tar_dim, src_dim, entries, sparse=True)


def stack_horizontal(mats):
    """Horizontally stack a list of sparse QQ matrices with the same nrows."""
    mats = [m for m in mats if m.ncols() > 0]
    if not mats:
        return None
    nrows = mats[0].nrows()
    total_cols = sum(m.ncols() for m in mats)
    entries = {}
    col_offset = 0
    for m in mats:
        for (r, c), val in m.dict().items():
            entries[(r, col_offset + c)] = val
        col_offset += m.ncols()
    return matrix(QQ, nrows, total_cols, entries, sparse=True)


def stack_vertical(mats):
    """Vertically stack a list of sparse QQ matrices with the same ncols."""
    mats = [m for m in mats if m.nrows() > 0]
    if not mats:
        return None
    ncols = mats[0].ncols()
    total_rows = sum(m.nrows() for m in mats)
    entries = {}
    row_offset = 0
    for m in mats:
        for (r, c), val in m.dict().items():
            entries[(row_offset + r, c)] = val
        row_offset += m.nrows()
    return matrix(QQ, total_rows, ncols, entries, sparse=True)


# ---------------------------------------------------------------------------
# Main H^k computation
# ---------------------------------------------------------------------------

def compute_hk(cx_type, e44_data, t_min=T_MIN, t_max=T_MAX,
               a_max=A_MAX, max_deg=MAX_DEG, k_level=K_LEVEL):
    """
    Compute H^{k_level} for the given complex type via n-decomposition.

    Algorithm
    ---------
    For each n in {k_level, k_level-1, ..., k_level-max_deg}:
      d = k_level - n  (PBW degree at level k_level in this n-slice)

      1. Collect all incoming contributions (one per sv_deg s with k_level-s in window):
         D^in_s: C^{k_level-s}_n to C^{k_level}_n   (d-s) to d
         Stack horizontally to form d^in_n.

      2. Collect all outgoing contributions (one per sv_deg s with k_level+s in window):
         D^out_s: C^{k_level}_n to C^{k_level+s}_n   d to (d+s)
         Stack vertically to form d^out_n.

      3. H^k_n = (dim C^{k_level}_n - rank d^out_n) - rank d^in_n

    H^{k_level} = Sigma_n H^k_n.

    Returns
    -------
    int -- dimension of H^{k_level}
    """
    morphisms = MORPHISMS_A if cx_type == 'A' else MORPHISMS_B
    sv_degs   = sorted(set(spec.sv_deg for spec in morphisms))

    t_min_tight, t_max_tight = tight_window_bounds(k_level, morphisms)
    t_min_eff = max(t_min_tight, t_min)
    t_max_eff = min(t_max_tight, t_max)

    _log(f'\n{"="*70}')
    _log(f'H^{k_level}  --  Complex {cx_type}  (n-decomposition)')
    _log(f'  Tight window: t=[{t_min_eff},{t_max_eff}]  '
          f'(full: [{t_min},{t_max}])  a_max={a_max}  max_deg={max_deg}')
    _log(f'  sv_degs: {sv_degs}')
    _log(f'{"="*70}\n')

    # Build CochainGroups for all levels in the tight window
    # (use the GLOBAL window for node enumeration so the node sets are correct)
    _log(f'[{ts()}] Building CochainGroups ...')
    groups = {}
    for k in range(t_min_eff, t_max_eff + 1):
        nodes = window_nodes(k, morphisms, t_min, t_max, a_max)
        g = CochainGroup(k, nodes, max_deg=max_deg, e44_data=e44_data)
        groups[k] = g
        _log(f'  k={k:+d}: {len(nodes):2d} nodes  total_dim={g.total_dim:>9,}')

    if k_level not in groups:
        _log(f'  ERROR: k={k_level} not in groups; check window.')
        return None

    g_k = groups[k_level]

    # -- Pre-compute all phi matrices once --------------------------------
    _log(f'\n[{ts()}] Pre-computing phi matrices ...')
    phi_cache = precompute_phi_matrices(
        morphisms, groups, e44_data, t_min, t_max, a_max, max_deg
    )
    _log(f'[{ts()}] Done.')

    # -- n-decomposition loop ----------------------------------------------
    _log(f'\n[{ts()}] Computing H^{k_level} by n-slices ...\n')

    hdr = (f'  {"n":>4}  {"d":>3}  {"dim C_n":>9}  '
           f'{"d^in shape":>22}  {"rk_in":>7}  '
           f'{"d^out shape":>22}  {"rk_out":>7}  '
           f'{"H^k_n":>7}  time')
    _log(hdr)
    _log('  ' + '-' * (len(hdr) - 2))

    total_h = 0

    for d in range(max_deg + 1):
        n = k_level - d

        # Dimension of C^k_n at PBW degree d
        c_k_n = sum(g_k.vermas[nd].dim(d) for nd in g_k.nodes)
        if c_k_n == 0:
            continue

        t0 = time.time()

        # --- Incoming differentials ---
        d_in_blocks = []
        for sv in sv_degs:
            k_src = k_level - sv
            if k_src not in groups:
                continue
            d_src = d - sv
            if d_src < 0 or d_src > max_deg:
                continue
            g_src = groups[k_src]
            D_in = assemble_n_slice(
                k_src, k_level, d_src, d,
                g_src, g_k, morphisms, phi_cache,
                t_min, t_max, a_max
            )
            if D_in.ncols() > 0:
                d_in_blocks.append(D_in)

        # Stack horizontal to get combined d^in: C^{k-*}_n to C^k_n
        d_in_combined = stack_horizontal(d_in_blocks) if d_in_blocks else None
        in_shape = (f'{d_in_combined.nrows()}x{d_in_combined.ncols()}'
                    if d_in_combined is not None else '--')
        rank_in = rank_certified(d_in_combined, f'n={n} d^in') \
            if d_in_combined is not None else 0

        # --- Outgoing differentials ---
        d_out_blocks = []
        for sv in sv_degs:
            k_tar = k_level + sv
            if k_tar not in groups:
                continue
            d_tar = d + sv
            if d_tar > max_deg:
                continue
            g_tar = groups[k_tar]
            D_out = assemble_n_slice(
                k_level, k_tar, d, d_tar,
                g_k, g_tar, morphisms, phi_cache,
                t_min, t_max, a_max
            )
            if D_out.nrows() > 0:
                d_out_blocks.append(D_out)

        # Stack vertical to get combined d^out: C^k_n to C^{k+*}_n
        d_out_combined = stack_vertical(d_out_blocks) if d_out_blocks else None
        out_shape = (f'{d_out_combined.nrows()}x{d_out_combined.ncols()}'
                     if d_out_combined is not None else '--')
        rank_out = rank_certified(d_out_combined, f'n={n} d^out') \
            if d_out_combined is not None else 0

        h_n = c_k_n - rank_in - rank_out
        total_h += h_n

        elapsed = time.time() - t0
        _log(f'  {n:>+4}  {d:>3}  {c_k_n:>9,}  '
              f'{in_shape:>22}  {rank_in:>7,}  '
              f'{out_shape:>22}  {rank_out:>7,}  '
              f'{h_n:>7,}  {elapsed:.1f}s')

    # -- Summary ----------------------------------------------------------
    _log(f'\n{"="*70}')
    _log(f'  H^{k_level}(Complex {cx_type}) = {total_h}')
    _log(f'{"="*70}')

    if total_h == 0:
        _log(f'\n  RESULT: H^{k_level} = 0  [OK] (Theorem 1.1 holds for Complex {cx_type})')
    else:
        _log(f'\n  RESULT: H^{k_level} = {total_h}  '
              f'(nonzero -- review computation or increase max_deg/a_max)')

    return total_h


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Fast H^k via n-decomposition and tight cochain window')
    parser.add_argument('--cx',      choices=['A', 'B', 'both'], default='both')
    parser.add_argument('--k',       type=int, default=K_LEVEL,
                        help='Cohomology level (default: 1)')
    parser.add_argument('--max-deg', type=int, default=MAX_DEG)
    parser.add_argument('--a-max',   type=int, default=A_MAX)
    parser.add_argument('--t-min',   type=int, default=T_MIN)
    parser.add_argument('--t-max',   type=int, default=T_MAX)
    parser.add_argument('--log',     type=str, default=None,
                        help='Path to progress log file (written immediately, '
                             'bypassing Sage stdout buffering)')
    args = parser.parse_args()

    K_LEVEL = args.k
    MAX_DEG = args.max_deg
    A_MAX   = args.a_max
    T_MIN   = args.t_min
    T_MAX   = args.t_max

    # Open progress log file -- all _log() calls write here immediately.
    # This bypasses Sage's fully-buffered stdout wrapper.
    _log_path = args.log
    if _log_path is None:
        cx_tag = args.cx.replace('/', '_')
        _log_path = os.path.join(_HERE, f'h1_{cx_tag}_k{K_LEVEL}_d{MAX_DEG}_a{A_MAX}.log')
    import cohomology_h1_fast as _self
    _self._log_fh = open(_log_path, 'w', buffering=1)
    _log(f'# Log: {_log_path}')
    _log(f'# Started: {datetime.datetime.now().isoformat()}')

    _log(f'[{ts()}] Loading e44_brackets.pkl ...')
    e44_data = load_e44()
    _log(f'[{ts()}] Loaded.\n')

    cx_list = ['A', 'B'] if args.cx == 'both' else [args.cx]
    results = {}
    for cx in cx_list:
        results[cx] = compute_hk(
            cx, e44_data,
            t_min=T_MIN, t_max=T_MAX,
            a_max=A_MAX, max_deg=MAX_DEG,
            k_level=K_LEVEL
        )

    _log(f'\n[{ts()}] Summary:')
    for cx, h in results.items():
        _log(f'  H^{K_LEVEL}(Complex {cx}) = {h}')
    _log('Done.')
