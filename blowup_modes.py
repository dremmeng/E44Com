#!/usr/bin/env sage
"""
blowup_modes.py  --  Explicit representatives of H^1(Complex A) blowup modes
============================================================================
H^1(Euler) != 0 means there exist cocycles in ker(d^out) that are NOT in
im(d^in).  These are the algebraic blowup modes: forcing terms f in C^1 for
which the linearized Euler equation d^in(u) = f has no solution.

For each n-slice (n = k - d) we compute:
  Z_n = ker(d^out_n)            -- the cocycles at this n-value
  B_n = im(d^in_n)  subseteq Z_n      -- the coboundaries
  H^1_n = Z_n / B_n             -- the blowup modes

Then for each representative vector we decompose by node (t,a,b,c) to read
off which sl_4 representations and t-values carry the obstruction.

n=+1 (d=0) is the cleanest: pure fiber modes W(t=1,a,b,c), no PBW structure.
These are the "soft" blowup modes -- spatial constants that can't be resolved.
n=-d+1 modes are "hard" (degree d polynomial spatial structure).

Usage
-----
  sage blowup_modes.py [--cx A] [--max-n 1] [--log blowup.log]
  sage blowup_modes.py --cx A --max-n -1 --full-reps 3
"""

from __future__ import print_function
import sys, os, time, datetime, argparse, traceback

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

try:
    sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', buffering=1)
    sys.stderr = os.fdopen(sys.stderr.fileno(), 'w', buffering=1)
except Exception:
    pass

from sage.all import QQ, matrix, vector, VectorSpace

from verma_modules import load_e44
from de_rham_complex import (
    CochainGroup, window_nodes,
    MORPHISMS_A, MORPHISMS_B,
)
from cohomology_h1_fast import (
    precompute_phi_matrices,
    assemble_n_slice,
    stack_vertical,
    stack_horizontal,
    rank_certified,
    ts,
)

# ===========================================================================
# Parameters
# ===========================================================================
T_MIN   = -1
T_MAX   =  3
A_MAX   =  4
MAX_DEG =  5
K_LEVEL =  1

_log_fh = None

def _log(*args):
    msg = ' '.join(str(a) for a in args)
    print(msg)
    if _log_fh is not None:
        try:
            _log_fh.write(msg + '\n')
            _log_fh.flush()
        except Exception:
            pass


# ===========================================================================
# Node offset map for a single PBW degree (matches assemble_n_slice ordering)
# ===========================================================================

def node_offsets_at_degree(g, d):
    """
    Returns (offsets_dict, total_dim) where offsets_dict[nd] = start index
    in the c_k_n vector for node nd at PBW degree d.
    Ordering matches assemble_n_slice's src_offs construction.
    """
    offs = {}
    cur = 0
    for nd in g.nodes:
        offs[nd] = cur
        cur += g.vermas[nd].dim(d)
    return offs, cur


# ===========================================================================
# H^1 representative extraction for one n-slice
# ===========================================================================

def h1_representatives(n, g_k, groups, morphisms, phi_cache,
                        t_min, t_max, a_max, max_deg, sv_degs):
    """
    Compute basis vectors for H^1_n = ker(d^out_n) / im(d^in_n).

    Returns
    -------
    reps : list of vectors in QQ^{c_k_n}
    offsets : dict{node: int}  -- nodeto start index in the c_k_n vector
    d_src : int  -- PBW degree of the C^k_n slice  (= K_LEVEL - n)
    """
    d_src = K_LEVEL - n
    if d_src < 0 or d_src > max_deg:
        return [], {}, d_src

    offsets, c_k_n = node_offsets_at_degree(g_k, d_src)
    if c_k_n == 0:
        return [], offsets, d_src

    # -- d^out : C^k_nto \oplus C^{k+sv}_n ------------------------------------
    d_out_blocks = []
    for sv in sv_degs:
        k_tar = K_LEVEL + sv
        if k_tar not in groups:
            continue
        d_tar = d_src + sv
        if d_tar > max_deg:
            continue
        g_tar = groups[k_tar]
        blk = assemble_n_slice(
            K_LEVEL, k_tar, d_src, d_tar,
            g_k, g_tar, morphisms, phi_cache,
            t_min, t_max, a_max
        )
        if blk.nrows() > 0 and blk.ncols() > 0:
            d_out_blocks.append(blk)

    D_out = stack_vertical(d_out_blocks) if d_out_blocks else \
            matrix(QQ, 0, c_k_n, sparse=True)

    # -- d^in : \oplus C^{k-sv}_nto C^k_n -------------------------------------
    d_in_blocks = []
    for sv in sv_degs:
        k_src = K_LEVEL - sv
        if k_src not in groups:
            continue
        d_in_deg = d_src - sv
        if d_in_deg < 0:
            continue
        g_src = groups[k_src]
        blk = assemble_n_slice(
            k_src, K_LEVEL, d_in_deg, d_src,
            g_src, g_k, morphisms, phi_cache,
            t_min, t_max, a_max
        )
        if blk.nrows() > 0 and blk.ncols() > 0:
            d_in_blocks.append(blk)

    D_in = stack_horizontal(d_in_blocks) if d_in_blocks else \
           matrix(QQ, c_k_n, 0, sparse=True)

    # -- Compute H^1_n = ker(D_out) / im(D_in) -----------------------------
    V = VectorSpace(QQ, c_k_n)

    # Kernel of D_out
    if D_out.nrows() > 0:
        Z = D_out.right_kernel()
    else:
        Z = V  # no outgoing differentialto everything is a cocycle

    # Image of D_in (as a subspace of V)
    if D_in.ncols() > 0:
        B_vecs = [V(D_in.column(j)) for j in range(D_in.ncols())
                  if D_in.column(j) != 0]
        B = V.subspace(B_vecs) if B_vecs else V.subspace([])
    else:
        B = V.subspace([])

    # Since d^2 = 0, B subseteq Z; compute quotient
    try:
        # Restrict B to Z: B is already subseteq Z so B cap Z = B
        H1 = Z.quotient(B)
        dim_H1 = H1.dimension()
        reps = [H1.lift(H1.basis()[i]) for i in range(dim_H1)]
    except Exception as e:
        _log(f'  WARNING: quotient failed: {e}')
        reps = list(Z.basis())[:max(0, Z.dimension() - B.dimension())]

    return reps, offsets, d_src


# ===========================================================================
# Decompose a representative vector by node
# ===========================================================================

def decompose_by_node(v, g_k, offsets, d):
    """
    Given a vector v in QQ^{c_k_n}, decompose it by node.

    Returns list of (node, component_vector, norm_sq) sorted by norm_sq desc.
    norm_sq = sum of squares of entries (over QQ, exact).
    """
    parts = []
    for nd in g_k.nodes:
        start = offsets[nd]
        dim_d = g_k.vermas[nd].dim(d)
        if dim_d == 0:
            continue
        comp = v[start: start + dim_d]
        norm_sq = sum(x*x for x in comp)
        if norm_sq != 0:
            parts.append((nd, comp, norm_sq))
    parts.sort(key=lambda x: -x[2])
    return parts


# ===========================================================================
# Summary: node-type distribution across all n-slices
# ===========================================================================

def node_type_counts(reps_by_n, g_k, offsets_by_n, d_by_n, nodes_all):
    """
    For each node, count how many blowup representatives have nonzero
    support on that node, across all n-slices.

    Returns dict: nodeto {n: count}
    """
    counts = {nd: {} for nd in nodes_all}
    for n, reps in reps_by_n.items():
        d = d_by_n[n]
        offs = offsets_by_n[n]
        for v in reps:
            for nd in g_k.nodes:
                start = offs.get(nd)
                if start is None:
                    continue
                dim_d = g_k.vermas[nd].dim(d)
                comp = v[start: start + dim_d]
                if any(x != 0 for x in comp):
                    counts[nd][n] = counts[nd].get(n, 0) + 1
    return counts


# ===========================================================================
# Main
# ===========================================================================

def run(cx_type='A', min_n=None, max_n=1, full_reps=5):
    """
    Compute and display H^1 blowup mode representatives.

    Parameters
    ----------
    cx_type  : 'A' or 'B'
    max_n    : largest n-slice to process (default 1 = most interpretable)
    min_n    : smallest n-slice (default = -(MAX_DEG - K_LEVEL))
    full_reps: number of representatives to display per slice (0 = suppress)
    """
    if min_n is None:
        min_n = -(MAX_DEG - K_LEVEL)

    morphisms = MORPHISMS_A if cx_type == 'A' else MORPHISMS_B
    sv_degs   = sorted(set(spec.sv_deg for spec in morphisms))

    _log(f'\n{"="*70}')
    _log(f'Blowup mode analysis -- Complex {cx_type}')
    _log(f'  k={K_LEVEL}, n in [{min_n}, {max_n}],  max_deg={MAX_DEG}, a_max={A_MAX}')
    _log(f'{"="*70}\n')

    _log(f'[{ts()}] Loading e44_brackets.pkl ...')
    e44_data = load_e44()

    _log(f'[{ts()}] Building CochainGroups ...')
    groups = {}
    for k in range(K_LEVEL - max(sv_degs), K_LEVEL + max(sv_degs) + 1):
        nodes = window_nodes(k, morphisms, T_MIN, T_MAX, A_MAX)
        if nodes:
            groups[k] = CochainGroup(k, nodes, max_deg=MAX_DEG,
                                     e44_data=e44_data)
            _log(f'  k={k:+d}: {len(nodes)} nodes, dim={groups[k].total_dim:,}')

    _log(f'\n[{ts()}] Pre-computing phi matrices ...')
    phi_cache = precompute_phi_matrices(
        morphisms, groups, e44_data, T_MIN, T_MAX, A_MAX, MAX_DEG
    )
    _log(f'[{ts()}] Done.\n')

    g_k = groups[K_LEVEL]

    reps_by_n      = {}
    offsets_by_n   = {}
    d_by_n         = {}
    h1_by_n        = {}

    for n in range(max_n, min_n - 1, -1):
        d = K_LEVEL - n
        if d < 0 or d > MAX_DEG:
            continue

        offsets, c_k_n = node_offsets_at_degree(g_k, d)
        if c_k_n == 0:
            continue

        _log(f'{"-"*70}')
        _log(f'n = {n:+d}  (d={d},  dim C^k_n = {c_k_n:,})')

        t0 = time.time()
        reps, offsets, d_src = h1_representatives(
            n, g_k, groups, morphisms, phi_cache,
            T_MIN, T_MAX, A_MAX, MAX_DEG, sv_degs
        )
        elapsed = time.time() - t0

        h1_dim = len(reps)
        _log(f'  H^1_{n:+d} = {h1_dim}   ({elapsed:.1f}s)')

        reps_by_n[n]    = reps
        offsets_by_n[n] = offsets
        d_by_n[n]       = d_src
        h1_by_n[n]      = h1_dim

        if h1_dim == 0 or full_reps == 0 or c_k_n > 50_000:
            if c_k_n > 50_000 and h1_dim > 0:
                _log(f'  (slice too large for full decomposition; node summary only)')
                _log_node_summary(reps, g_k, offsets, d_src, limit=50)
            continue

        # -- Per-representative display ---------------------------------
        n_show = min(full_reps, h1_dim)
        _log(f'\n  Showing {n_show}/{h1_dim} representatives, decomposed by node:\n')
        _log(f'  {"rep":>4}  {"node (t,a,b,c)":>20}  {"dim_W":>6}  '
             f'{"nnz":>5}  {"norm^2":>20}')
        _log(f'  ' + '-'*60)

        for i, v in enumerate(reps[:n_show]):
            parts = decompose_by_node(v, g_k, offsets, d_src)
            first = True
            for nd, comp, norm_sq in parts:
                nd_str = f't={nd.t}, ({nd.a},{nd.b},{nd.c})'
                nnz = sum(1 for x in comp if x != 0)
                dim_W = len(comp)
                prefix = f'  {i+1:>4}' if first else f'  {"":>4}'
                _log(f'  {prefix}  {nd_str:>20}  {dim_W:>6}  '
                     f'{nnz:>5}  {norm_sq}')
                first = False
            if not parts:
                _log(f'  {i+1:>4}  (zero vector?)')
        _log('')

    # -- Cross-slice summary --------------------------------------------
    _log(f'\n{"="*70}')
    _log(f'H^1 by n-slice summary')
    _log(f'{"="*70}')
    _log(f'  {"n":>4}  {"d":>3}  {"dim C^k_n":>12}  {"H^1_n":>8}  {"frac":>8}')
    _log(f'  ' + '-'*48)
    total_h1 = 0
    for n in range(max_n, min_n - 1, -1):
        d = K_LEVEL - n
        if d < 0 or d > MAX_DEG:
            continue
        _, c_k_n = node_offsets_at_degree(g_k, d)
        h = h1_by_n.get(n, 0)
        total_h1 += h
        frac = f'{h/c_k_n:.4f}' if c_k_n > 0 else '--'
        _log(f'  {n:>+4}  {d:>3}  {c_k_n:>12,}  {h:>8,}  {frac:>8}')
    _log(f'  {"":>4}  {"":>3}  {"":>12}  {total_h1:>8,}  (total)')

    # -- Node contribution table ----------------------------------------
    _log(f'\n{"="*70}')
    _log(f'Node contribution to H^1 (over computed n-slices)')
    _log(f'{"="*70}')
    counts = node_type_counts(reps_by_n, g_k, offsets_by_n, d_by_n, g_k.nodes)
    non_empty = [(nd, cnt) for nd, cnt in counts.items() if cnt]
    non_empty.sort(key=lambda x: -sum(x[1].values()))
    _log(f'  {"node (t,a,b,c)":>22}  {"total reps":>10}  by n-slice')
    _log(f'  ' + '-'*60)
    for nd, cnt in non_empty:
        nd_str = f't={nd.t}, ({nd.a},{nd.b},{nd.c})'
        total = sum(cnt.values())
        by_n = '  '.join(f'n={k:+d}:{v}' for k, v in sorted(cnt.items(), reverse=True))
        _log(f'  {nd_str:>22}  {total:>10,}  {by_n}')

    _log(f'\n[{ts()}] Done.')


def _log_node_summary(reps, g_k, offsets, d, limit=50):
    """For large slices: just show which nodes are hit, no full vector display."""
    node_hit = {}
    for v in reps[:limit]:
        parts = decompose_by_node(v, g_k, offsets, d)
        for nd, _, norm_sq in parts:
            node_hit[nd] = node_hit.get(nd, 0) + 1
    _log(f'  Node hits (first {min(limit,len(reps))} reps):')
    for nd, cnt in sorted(node_hit.items(), key=lambda x: -x[1]):
        _log(f'    t={nd.t}, ({nd.a},{nd.b},{nd.c}): {cnt}')


# ===========================================================================
# Entry point
# ===========================================================================

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Blowup mode representatives of H^1(Euler)')
    parser.add_argument('--cx',       choices=['A', 'B'], default='A')
    parser.add_argument('--max-n',    type=int, default=1,
                        help='Largest n-slice to process (default: 1)')
    parser.add_argument('--min-n',    type=int, default=None,
                        help='Smallest n-slice (default: -(MAX_DEG - K_LEVEL))')
    parser.add_argument('--full-reps', type=int, default=5,
                        help='How many reps to display per slice (default: 5)')
    parser.add_argument('--max-deg',  type=int, default=MAX_DEG)
    parser.add_argument('--a-max',    type=int, default=A_MAX)
    parser.add_argument('--t-min',    type=int, default=T_MIN)
    parser.add_argument('--t-max',    type=int, default=T_MAX)
    parser.add_argument('--log',      type=str, default=None)
    args = parser.parse_args()

    MAX_DEG = args.max_deg
    A_MAX   = args.a_max
    T_MIN   = args.t_min
    T_MAX   = args.t_max

    log_path = args.log
    if log_path is None:
        log_path = os.path.join(
            _HERE, f'blowup_{args.cx}_d{MAX_DEG}_a{A_MAX}.log')
    import blowup_modes as _self
    _self._log_fh = open(log_path, 'w', buffering=1)
    _log(f'# Log: {log_path}')
    _log(f'# Started: {datetime.datetime.now().isoformat()}')

    run(
        cx_type=args.cx,
        max_n=args.max_n,
        min_n=args.min_n,
        full_reps=args.full_reps,
    )
