#!/usr/bin/env sage
"""
production_run.py — Resumable production computation of E(4,4) de Rham cohomology.
====================================================================================

Each differential D_{k→k+s} is saved as an individual pickle checkpoint.
Re-run at any time; completed checkpoints are skipped automatically.

Phases
------
  status      Show which checkpoints exist vs. missing.
  build       Compute and save differential matrices one at a time.
  cohomology  Load saved checkpoints and compute H^k level by level.

Usage (inside SageMath, run from the E44Com directory)
------
  sage production_run.py status
  sage production_run.py build A          # Complex A differentials
  sage production_run.py build B          # Complex B differentials
  sage production_run.py cohomology A     # H^k for Complex A
  sage production_run.py cohomology B     # H^k for Complex B
  sage production_run.py cohomology A --level 0   # H^0 only (split across sessions)

Default parameters (override with flags; see --help):
  T_MIN=-6  T_MAX=6  A_MAX=4  MAX_DEG=5

Checkpoint layout (CHECKPOINT_DIR = ./checkpoints/)
--------------------
  diff_A_k+00_to_k+01.pkl   — D_{0→1} for Complex A
  diff_B_k-03_to_k+00.pkl   — D_{-3→0} for Complex B
  cohomology_A_k+00.pkl     — cohomology result at k=0 for Complex A
  ...

Memory notes (32 GB system)
---------------------------
  At full parameters (a_max=4, max_deg=5), each cochain level has dim ~170 K
  (phat4 fibers for nodes (1,0,0), (0,1,0), (0,0,1)).  Sparse differential
  matrices are manageable (~MB each).  The right_kernel() call in cohomology_at()
  does exact QQ row reduction, which is the long pole — it may run for hours or
  days per level at full scale.  Use --level to spread cohomology across sessions.

  If memory is tight between build sessions, the Verma module cache (_VERMA_CACHE)
  grows as more modules are built.  This is intentional for speed; the cache is
  cleared between separate sage invocations automatically.

Decomposition hint (advanced)
------------------------------
  The complex decomposes by the invariant  n = k - d  (cochain level minus
  internal PBW degree).  Within each n-slice, the sub-complex has dimension
  proportional to one block per node per level, reducing the cohomology
  computation to many small systems.  Not yet implemented here; implement as
  a future optimisation if full-scale right_kernel() is too slow.
"""

from __future__ import print_function
import sys, os, argparse, pickle, time, datetime

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

from sage.all import QQ, matrix, VectorSpace

from verma_modules import load_e44
from de_rham_complex import (
    DeRhamComplexA, DeRhamComplexB,
    CochainGroup, window_nodes,
    MORPHISMS_A, MORPHISMS_B,
    assemble_differential, get_verma, _VERMA_CACHE,
    node_fiber_type,
)
from cohomology import outgoing_matrix, incoming_matrix


# ===========================================================================
# Default production parameters  (edit or override via CLI)
# ===========================================================================

T_MIN   = -6
T_MAX   =  6
A_MAX   =  4
MAX_DEG =  5

CHECKPOINT_DIR = os.path.join(_HERE, 'checkpoints')

# Cochain shifts per complex
SHIFTS = {'A': (1, 2), 'B': (1, 3, 4)}
MORPHISMS = {'A': MORPHISMS_A, 'B': MORPHISMS_B}
CX_CLASS  = {'A': DeRhamComplexA, 'B': DeRhamComplexB}


# ===========================================================================
# Utilities
# ===========================================================================

def ts():
    """Current timestamp string."""
    return datetime.datetime.now().strftime('%H:%M:%S')


def diff_path(cx_type, k_src, k_tar):
    fname = f'diff_{cx_type}_k{k_src:+03d}_to_k{k_tar:+03d}.pkl'
    return os.path.join(CHECKPOINT_DIR, fname)


def coh_path(cx_type, k):
    fname = f'cohomology_{cx_type}_k{k:+03d}.pkl'
    return os.path.join(CHECKPOINT_DIR, fname)


def all_diff_pairs(cx_type, t_min, t_max):
    """Generate all (k_src, k_tar) pairs for a complex type."""
    for s in SHIFTS[cx_type]:
        for k in range(t_min, t_max + 1 - s):
            yield k, k + s


def build_groups(cx_type, t_min, t_max, a_max, max_deg, e44_data):
    """Build CochainGroups for all levels.  Caches Verma modules."""
    morphisms = MORPHISMS[cx_type]
    groups = {}
    for k in range(t_min, t_max + 1):
        nodes_k = window_nodes(k, morphisms, t_min, t_max, a_max)
        groups[k] = CochainGroup(k, nodes_k, max_deg=max_deg, e44_data=e44_data)
    return groups


def reconstruct_complex(cx_type, groups, differentials, t_min, t_max, a_max, max_deg):
    """Build a minimal DeRhamComplexA/B shell from pre-built groups and differentials."""
    cls = CX_CLASS[cx_type]
    cx = object.__new__(cls)
    cx.t_min       = t_min
    cx.t_max       = t_max
    cx.a_max       = a_max
    cx.max_deg     = max_deg
    cx.e44_data    = None
    cx.positions   = list(range(t_min, t_max + 1))
    cx.groups      = groups
    cx.differentials = differentials
    return cx


# ===========================================================================
# Phase: status
# ===========================================================================

def phase_status(t_min, t_max, a_max, max_deg):
    print(f'\nCheckpoint status  (t=[{t_min},{t_max}], a_max={a_max}, max_deg={max_deg})')
    print(f'Directory: {CHECKPOINT_DIR}\n')

    for cx_type in ('A', 'B'):
        n_done = n_miss = 0
        print(f'  Complex {cx_type}  [shifts: {SHIFTS[cx_type]}]')
        for k, k_tar in sorted(all_diff_pairs(cx_type, t_min, t_max)):
            p = diff_path(cx_type, k, k_tar)
            if os.path.isfile(p):
                kb = os.path.getsize(p) // 1024
                print(f'    DONE  D_{cx_type}[{k:+d}→{k_tar:+d}]  {kb:6d} KB')
                n_done += 1
            else:
                print(f'    MISS  D_{cx_type}[{k:+d}→{k_tar:+d}]')
                n_miss += 1
        # Cohomology checkpoints
        for k in range(t_min, t_max + 1):
            p = coh_path(cx_type, k)
            if os.path.isfile(p):
                kb = os.path.getsize(p) // 1024
                print(f'    DONE  H_{cx_type}^{k:+d}           {kb:6d} KB')
        print(f'    — {n_done} differentials done, {n_miss} missing\n')


# ===========================================================================
# Phase: build
# ===========================================================================

def phase_build(cx_type, t_min, t_max, a_max, max_deg, e44_data):
    """Compute and checkpoint every differential matrix for one complex."""
    os.makedirs(CHECKPOINT_DIR, exist_ok=True)
    morphisms = MORPHISMS[cx_type]

    # Which differentials still need computing?
    todo = [(k, k_tar) for k, k_tar in sorted(all_diff_pairs(cx_type, t_min, t_max))
            if not os.path.isfile(diff_path(cx_type, k, k_tar))]
    done = [(k, k_tar) for k, k_tar in sorted(all_diff_pairs(cx_type, t_min, t_max))
            if     os.path.isfile(diff_path(cx_type, k, k_tar))]

    print(f'[{ts()}] Complex {cx_type}: {len(done)} differentials already done, '
          f'{len(todo)} to compute.')

    if not todo:
        print(f'[{ts()}] Nothing to do for Complex {cx_type}.')
        return

    # Build all CochainGroups.
    # This constructs Verma modules (phat4 for nodes in FIBER_TYPE) and caches them.
    print(f'[{ts()}] Building CochainGroups ...')
    t0 = time.time()
    groups = build_groups(cx_type, t_min, t_max, a_max, max_deg, e44_data)
    print(f'[{ts()}] Groups built in {time.time()-t0:.1f}s.  Dimensions:')
    for k in sorted(groups):
        g = groups[k]
        print(f'  k={k:+d}: {len(g.nodes):2d} nodes, dim={g.total_dim:7d}')

    # Compute each missing differential.
    n_total = len(todo)
    for idx, (k, k_tar) in enumerate(todo):
        out_path = diff_path(cx_type, k, k_tar)
        g_src = groups[k]
        g_tar = groups[k_tar]

        print(f'[{ts()}] ({idx+1}/{n_total}) D_{cx_type}[{k:+d}→{k_tar:+d}]  '
              f'src_dim={g_src.total_dim}  tar_dim={g_tar.total_dim} ...',
              flush=True)
        t0 = time.time()

        D = assemble_differential(
            g_src, g_tar, morphisms,
            e44_data, t_min, t_max, a_max
        )

        elapsed = time.time() - t0
        nnz = len(D.dict())

        # Save
        payload = {
            'cx_type': cx_type, 'k_src': k, 'k_tar': k_tar,
            't_min': t_min, 't_max': t_max, 'a_max': a_max, 'max_deg': max_deg,
            'D': D,
        }
        with open(out_path, 'wb') as f:
            pickle.dump(payload, f, protocol=4)

        kb = os.path.getsize(out_path) // 1024
        print(f'[{ts()}]   → {D.nrows()}×{D.ncols()}, {nnz} nnz, '
              f'{kb} KB on disk, {elapsed:.1f}s', flush=True)

    print(f'[{ts()}] Complex {cx_type} build phase complete.')


# ===========================================================================
# Phase: cohomology
# ===========================================================================

def cohomology_at_from_matrices(groups, differentials, k, t_min, t_max):
    """
    Compute H^k given pre-built groups and differential dict.

    Mirrors cohomology.py::cohomology_at() but uses our groups/diffs directly.
    Saves the basis matrices only if dim_H > 0 and dim_H <= 200 (to avoid
    multi-GB dense matrices at large scale; set SAVE_BASES = True to override).

    Returns dict with keys: k, dim_Ck, dim_ker, dim_im, dim_H, im_subset_ker.
    """
    SAVE_BASES = False  # set True to save full cocycle/cobound basis matrices

    dim_k = groups[k].total_dim
    V     = VectorSpace(QQ, dim_k)

    # ── Outgoing ──────────────────────────────────────────────────────────
    out_blocks = [D for (ks, _), D in sorted(differentials.items()) if ks == k]
    if not out_blocks:
        O_k = matrix(QQ, 0, dim_k, sparse=True)
    else:
        total_rows = sum(B.nrows() for B in out_blocks)
        O_k = matrix(QQ, total_rows, dim_k, sparse=True)
        r = 0
        for B in out_blocks:
            for (i, j), val in B.dict().items():
                O_k[r + i, j] = val
            r += B.nrows()

    # ── Incoming ─────────────────────────────────────────────────────────
    in_blocks = [D for (_, kt), D in sorted(differentials.items()) if kt == k]
    if not in_blocks:
        I_k = matrix(QQ, dim_k, 0, sparse=True)
    else:
        total_cols = sum(B.ncols() for B in in_blocks)
        I_k = matrix(QQ, dim_k, total_cols, sparse=True)
        c = 0
        for B in in_blocks:
            for (i, j), val in B.dict().items():
                I_k[i, c + j] = val
            c += B.ncols()

    # ── Kernel / image ────────────────────────────────────────────────────
    print(f'  [k={k:+d}] right_kernel of {O_k.nrows()}×{O_k.ncols()} ...', flush=True)
    t0 = time.time()
    if O_k.nrows() == 0:
        K_space = V
    else:
        K_space = V.subspace(O_k.right_kernel().basis())
    print(f'  [k={k:+d}] ker done ({time.time()-t0:.1f}s), dim_ker={K_space.dimension()}', flush=True)

    print(f'  [k={k:+d}] column_space of {I_k.nrows()}×{I_k.ncols()} ...', flush=True)
    t0 = time.time()
    if I_k.ncols() == 0:
        I_space = V.subspace([])
    else:
        I_space = V.subspace(I_k.column_space().basis())
    print(f'  [k={k:+d}] im done ({time.time()-t0:.1f}s), dim_im={I_space.dimension()}', flush=True)

    print(f'  [k={k:+d}] intersection ...', flush=True)
    t0 = time.time()
    inter = K_space.intersection(I_space)
    print(f'  [k={k:+d}] inter done ({time.time()-t0:.1f}s), dim_inter={inter.dimension()}', flush=True)

    dim_ker = K_space.dimension()
    dim_im  = inter.dimension()
    dim_H   = dim_ker - dim_im
    im_sub  = I_space.is_subspace(K_space)

    result = {
        'k': k, 'dim_Ck': dim_k, 'dim_ker': dim_ker,
        'dim_im': dim_im, 'dim_H': dim_H, 'im_subset_ker': im_sub,
    }
    if SAVE_BASES and dim_H <= 200:
        cocycles_rows = list(K_space.basis())
        cobounds_rows = list(inter.basis())
        result['cocycles'] = (matrix(QQ, cocycles_rows) if cocycles_rows
                              else matrix(QQ, 0, dim_k))
        result['cobounds'] = (matrix(QQ, cobounds_rows) if cobounds_rows
                              else matrix(QQ, 0, dim_k))

    return result


def phase_cohomology(cx_type, t_min, t_max, a_max, max_deg, e44_data,
                     levels=None):
    """
    Load saved differential checkpoints and compute H^k.

    levels : list of int or None — if given, compute only those cochain levels.
    """
    # Check which differentials are available
    missing = [(k, k_tar) for k, k_tar in all_diff_pairs(cx_type, t_min, t_max)
               if not os.path.isfile(diff_path(cx_type, k, k_tar))]
    if missing:
        print(f'[{ts()}] WARNING: {len(missing)} differential(s) not yet computed:')
        for k, k_tar in missing:
            print(f'  D_{cx_type}[{k:+d}→{k_tar:+d}]  MISSING')
        print('These levels will be treated as zero differentials.')

    print(f'[{ts()}] Building CochainGroups for Complex {cx_type} ...', flush=True)
    groups = build_groups(cx_type, t_min, t_max, a_max, max_deg, e44_data)
    print(f'[{ts()}] Groups built.', flush=True)

    # Load available differentials
    print(f'[{ts()}] Loading differential checkpoints ...', flush=True)
    differentials = {}
    for k, k_tar in all_diff_pairs(cx_type, t_min, t_max):
        p = diff_path(cx_type, k, k_tar)
        if os.path.isfile(p):
            with open(p, 'rb') as f:
                payload = pickle.load(f)
            differentials[(k, k_tar)] = payload['D']
            print(f'  loaded D_{cx_type}[{k:+d}→{k_tar:+d}]  '
                  f'{len(payload["D"].dict())} nnz', flush=True)
        else:
            # Missing → zero matrix
            g_src = groups[k]
            g_tar = groups[k_tar]
            differentials[(k, k_tar)] = matrix(
                QQ, g_tar.total_dim, g_src.total_dim, sparse=True)

    print(f'[{ts()}] All differentials loaded.\n', flush=True)

    target_levels = levels if levels is not None else list(range(t_min, t_max + 1))

    print(f'Computing H^k for Complex {cx_type}, '
          f'levels: {target_levels}\n')
    print(f"  {'k':>4}  {'dim C^k':>9}  {'dim ker':>9}  {'dim im':>9}  "
          f"{'dim H^k':>9}  {'d²=0':>5}")
    print('  ' + '-' * 55)

    all_results = {}
    for k in target_levels:
        # Skip if already done
        p = coh_path(cx_type, k)
        if os.path.isfile(p):
            with open(p, 'rb') as f:
                r = pickle.load(f)
            print(f"  {k:>4}  {r['dim_Ck']:>9}  {r['dim_ker']:>9}  "
                  f"{r['dim_im']:>9}  {r['dim_H']:>9}  "
                  f"{'✓' if r['im_subset_ker'] else '✗':>5}  [cached]")
            all_results[k] = r
            continue

        print(f'[{ts()}] Computing H^{k} ...')
        t0 = time.time()

        r = cohomology_at_from_matrices(groups, differentials, k, t_min, t_max)

        elapsed = time.time() - t0
        all_results[k] = r

        # Save checkpoint
        with open(p, 'wb') as f:
            pickle.dump(r, f, protocol=4)

        d2_flag = '✓' if r['im_subset_ker'] else '✗'
        print(f"  {k:>4}  {r['dim_Ck']:>9}  {r['dim_ker']:>9}  "
              f"{r['dim_im']:>9}  {r['dim_H']:>9}  {d2_flag:>5}  "
              f"[{elapsed:.0f}s]", flush=True)

    # Summary
    print('\n  Summary:')
    print(f"  {'k':>4}  {'dim H^k':>9}")
    for k in sorted(all_results):
        r = all_results[k]
        print(f"  {k:>4}  {r['dim_H']:>9}")

    return all_results


# ===========================================================================
# CLI
# ===========================================================================

def parse_args():
    p = argparse.ArgumentParser(
        description='Production E(4,4) de Rham cohomology computation.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument('phase', choices=['status', 'build', 'cohomology'],
                   help='Which phase to run.')
    p.add_argument('complex', nargs='?', choices=['A', 'B'],
                   help='Which complex (required for build and cohomology).')
    p.add_argument('--t-min',   type=int, default=T_MIN)
    p.add_argument('--t-max',   type=int, default=T_MAX)
    p.add_argument('--a-max',   type=int, default=A_MAX)
    p.add_argument('--max-deg', type=int, default=MAX_DEG)
    p.add_argument('--level', type=int, nargs='+', metavar='K',
                   help='(cohomology phase only) compute H^k only at these level(s).')
    p.add_argument('--checkpoint-dir', default=CHECKPOINT_DIR,
                   help=f'Directory for pickle files (default: {CHECKPOINT_DIR})')
    return p.parse_args()


def main():
    global CHECKPOINT_DIR
    args = parse_args()
    CHECKPOINT_DIR = args.checkpoint_dir

    t_min   = args.t_min
    t_max   = args.t_max
    a_max   = args.a_max
    max_deg = args.max_deg

    if args.phase == 'status':
        phase_status(t_min, t_max, a_max, max_deg)
        return

    if args.complex is None:
        print('ERROR: --complex A or B required for build/cohomology phases.')
        sys.exit(1)

    cx_type = args.complex

    print(f'[{ts()}] Loading e44_data ...')
    e44_data = load_e44()
    print(f'[{ts()}] e44_data loaded.')
    print(f'Parameters: t=[{t_min},{t_max}], a_max={a_max}, max_deg={max_deg}')
    print(f'Complex: {cx_type}')
    print(f'Checkpoint dir: {CHECKPOINT_DIR}\n')

    if args.phase == 'build':
        phase_build(cx_type, t_min, t_max, a_max, max_deg, e44_data)

    elif args.phase == 'cohomology':
        levels = args.level  # None → all levels
        phase_cohomology(cx_type, t_min, t_max, a_max, max_deg, e44_data,
                         levels=levels)


if __name__ == '__main__':
    main()
