#!/usr/bin/env sage
"""
borel_nu_h1_run.py -- Production H^1 for Borel-reduced Complex B at generic nu.

Model used for generic nu:
  Remove phi_1E from MORPHISMS_B (same proxy used in t_deform_h1.py), then
  apply the Borel 4->3 basis filter:
    alpha_4 = 0 and d_4 absent from wedge set.

For k=1, computes
  H^1 = ker(d_out_red) / im(d_in_red)
with rank-nullity, and supports per-differential checkpoints so reruns resume.
"""

from __future__ import print_function

import argparse
import datetime
import os
import pickle
import time

from sage.all import GF, QQ, matrix

import de_rham_complex as drc


DEFAULT_T_MIN = -6
DEFAULT_T_MAX = 6
DEFAULT_A_MAX = 4
DEFAULT_MAX_DEG = 5
DEFAULT_K = 1


def ts():
    return datetime.datetime.now().strftime('%H:%M:%S')


def keep_basis_state(mon, _w_idx, _node, _d):
    """Borel 4->3 filter on PBW state."""
    alpha, sset = mon
    return alpha[3] == 0 and (3 not in sset)


def reduced_index_data(group):
    """Return kept full indices and map full->reduced for one group."""
    kept_globals = []
    for node in group.nodes:
        V = group.vermas[node]
        for d in range(group.max_deg + 1):
            sl = group.basis_slice(node, d)
            basis_d = V.basis(d)
            for local_idx, (mon, w_idx) in enumerate(basis_d):
                if keep_basis_state(mon, w_idx, node, d):
                    kept_globals.append(sl.start + local_idx)
    full_to_red = {gidx: ridx for ridx, gidx in enumerate(kept_globals)}
    return kept_globals, full_to_red


def projection_and_inclusion(full_dim, kept_globals):
    """Build P: full->reduced and I: reduced->full sparse matrices."""
    m = len(kept_globals)
    p_entries = {(r, g): QQ(1) for r, g in enumerate(kept_globals)}
    i_entries = {(g, r): QQ(1) for r, g in enumerate(kept_globals)}
    P = matrix(QQ, m, full_dim, p_entries, sparse=True)
    I = matrix(QQ, full_dim, m, i_entries, sparse=True)
    return P, I


def assemble_full_safe(src_group, tar_group, morphisms, e44_data, t_min, t_max, a_max):
    """Assemble full differential while skipping failing edges."""
    entries = {}
    for spec in morphisms:
        for src_node, tar_node, phi_args in spec.edges(t_min, t_max, a_max):
            if src_node.t != src_group.k or tar_node.t != tar_group.k:
                continue
            if src_node not in src_group.nodes or tar_node not in tar_group.nodes:
                continue

            msd = drc.src_max_deg(spec.sv_deg, src_group.max_deg)
            src_data = e44_data if drc.node_fiber_type(src_node) == 'phat4' else None
            try:
                block_full = drc.get_morphism_matrix(
                    spec,
                    src_node,
                    tar_node,
                    phi_args,
                    src_group,
                    tar_group,
                    e44_data,
                    msd,
                    src_e44_data=src_data,
                )
            except BaseException:
                continue

            for (i_full, j_full), val in block_full.dict().items():
                key = (i_full, j_full)
                prev = entries.get(key, QQ(0))
                entries[key] = prev + val

    return matrix(QQ, tar_group.total_dim, src_group.total_dim, entries, sparse=True)


def rank_fast(M, label=''):
    """GF(p) rank with two-prime consistency check and QQ fallback."""
    if M.nrows() == 0 or M.ncols() == 0:
        return 0
    if max(M.nrows(), M.ncols()) <= 2000:
        try:
            return M.rank(algorithm='padic')
        except TypeError:
            return M.rank()
    p1 = 65521
    p2 = 65537
    r1 = M.change_ring(GF(p1)).rank()
    r2 = M.change_ring(GF(p2)).rank()
    if r1 != r2:
        print(f'  WARNING {label}: GF rank mismatch ({r1} vs {r2}), using QQ padic', flush=True)
        try:
            return M.rank(algorithm='padic')
        except TypeError:
            return M.rank()
    return r1


def run_tag(args):
    """Stable identifier for a production window/level."""
    return (
        f'tmin{args.t_min:+03d}_tmax{args.t_max:+03d}_'
        f'a{args.a_max}_deg{args.max_deg}_k{args.level:+03d}'
    )


def diff_cache_path(checkpoint_dir, tag, k_src, k_tar):
    return os.path.join(checkpoint_dir, f'diff_borel_nu_{tag}_k{k_src:+03d}_to_k{k_tar:+03d}.pkl')


def result_cache_path(checkpoint_dir, tag, k):
    return os.path.join(checkpoint_dir, f'cohomology_borel_nu_{tag}_k{k:+03d}.pkl')


def linear_cache_path(checkpoint_dir, tag, k):
    return os.path.join(checkpoint_dir, f'linear_data_borel_nu_{tag}_k{k:+03d}.pkl')


def needed_pairs_for_level(k):
    """Complex B shifts are 1,3,4."""
    incoming = [(k - 1, k), (k - 3, k), (k - 4, k)]
    outgoing = [(k, k + 1), (k, k + 3), (k, k + 4)]
    return incoming + outgoing


def run(args):
    os.makedirs(args.checkpoint_dir, exist_ok=True)
    tag = run_tag(args)

    print(f'[{ts()}] Loading e44_data ...', flush=True)
    e44_data = drc.load_e44()
    print(f'[{ts()}] e44_data loaded.', flush=True)

    morph_nu = [m for m in drc.MORPHISMS_B if m.enumerate_fn.__name__ != '_enum_1E']
    print(f'[{ts()}] Using Complex B generic-nu morphisms (phi_1E removed).', flush=True)
    print(f'[{ts()}] Window: t=[{args.t_min},{args.t_max}], a_max={args.a_max}, max_deg={args.max_deg}, k={args.level}', flush=True)
    print(f'[{ts()}] Checkpoint dir: {args.checkpoint_dir}', flush=True)
    print(f'[{ts()}] Run tag: {tag}', flush=True)

    # Build only levels touched by incoming/outgoing maps at k.
    pairs = [
        (ks, kt) for (ks, kt) in needed_pairs_for_level(args.level)
        if args.t_min <= ks <= args.t_max and args.t_min <= kt <= args.t_max
    ]
    needed_levels = sorted(set([x for p in pairs for x in p]))
    if not pairs:
        raise SystemExit('No valid incoming/outgoing pairs in this window for requested level.')

    print(f'[{ts()}] Building CochainGroup objects for levels: {needed_levels}', flush=True)
    groups = {}
    red_data = {}
    for k in needed_levels:
        nodes = drc.window_nodes(k, morph_nu, args.t_min, args.t_max, args.a_max)
        g = drc.CochainGroup(k, nodes, max_deg=args.max_deg, e44_data=e44_data)
        kept, fmap = reduced_index_data(g)
        P, I = projection_and_inclusion(g.total_dim, kept)
        groups[k] = g
        red_data[k] = {
            'kept': kept,
            'full_to_red': fmap,
            'P': P,
            'I': I,
            'dim_red': len(kept),
        }
        print(f'  [k={k:+d}] full_dim={g.total_dim}, red_dim={len(kept)}', flush=True)

    # Assemble/load reduced differentials.
    red_diffs = {}
    for k_src, k_tar in pairs:
        cp = diff_cache_path(args.checkpoint_dir, tag, k_src, k_tar)
        if os.path.isfile(cp):
            with open(cp, 'rb') as f:
                payload = pickle.load(f)
            if (
                payload.get('t_min') == args.t_min and
                payload.get('t_max') == args.t_max and
                payload.get('a_max') == args.a_max and
                payload.get('max_deg') == args.max_deg and
                payload.get('level') == args.level and
                payload.get('model') == 'B_generic_nu_borel_reduced'
            ):
                D_red = payload['D_red']
                red_diffs[(k_src, k_tar)] = D_red
                print(f'[{ts()}] loaded D_red[{k_src:+d}->{k_tar:+d}] nnz={len(D_red.dict())}', flush=True)
                continue
            print(f'[{ts()}] ignoring stale cache for D_red[{k_src:+d}->{k_tar:+d}]', flush=True)

        g_src = groups[k_src]
        g_tar = groups[k_tar]
        P_tar = red_data[k_tar]['P']
        I_src = red_data[k_src]['I']

        print(f'[{ts()}] assembling full D[{k_src:+d}->{k_tar:+d}] ...', flush=True)
        t0 = time.time()
        D_full = assemble_full_safe(g_src, g_tar, morph_nu, e44_data,
                                    args.t_min, args.t_max, args.a_max)
        print(f'[{ts()}] projecting to reduced D_red[{k_src:+d}->{k_tar:+d}] ...', flush=True)
        D_red = P_tar * D_full * I_src
        elapsed = time.time() - t0

        payload = {
            'k_src': k_src,
            'k_tar': k_tar,
            'D_red': D_red,
            'nnz_full': len(D_full.dict()),
            'nnz_red': len(D_red.dict()),
            't_min': args.t_min,
            't_max': args.t_max,
            'a_max': args.a_max,
            'max_deg': args.max_deg,
            'level': args.level,
            'model': 'B_generic_nu_borel_reduced',
        }
        with open(cp, 'wb') as f:
            pickle.dump(payload, f, protocol=4)

        red_diffs[(k_src, k_tar)] = D_red
        print(f'[{ts()}] saved D_red[{k_src:+d}->{k_tar:+d}] nnz={len(D_red.dict())} ({elapsed:.1f}s)', flush=True)

    # Build O_k and I_k in reduced coordinates.
    k = args.level
    dim_k = red_data[k]['dim_red']

    out_blocks = [D for (ks, _), D in sorted(red_diffs.items()) if ks == k]
    if not out_blocks:
        O_k = matrix(QQ, 0, dim_k, sparse=True)
    else:
        O_k = matrix(QQ, sum(B.nrows() for B in out_blocks), dim_k, sparse=True)
        r = 0
        for B in out_blocks:
            for (i, j), val in B.dict().items():
                O_k[r + i, j] = val
            r += B.nrows()

    in_blocks = [D for (_, kt), D in sorted(red_diffs.items()) if kt == k]
    if not in_blocks:
        I_k = matrix(QQ, dim_k, 0, sparse=True)
    else:
        I_k = matrix(QQ, dim_k, sum(B.ncols() for B in in_blocks), sparse=True)
        c = 0
        for B in in_blocks:
            for (i, j), val in B.dict().items():
                I_k[i, c + j] = val
            c += B.ncols()

    print(f'[{ts()}] rank of outgoing O_k {O_k.nrows()}x{O_k.ncols()} ...', flush=True)
    t0 = time.time()
    rank_out = rank_fast(O_k, f'O_k[{k:+d}]')
    dim_ker = dim_k - rank_out
    print(f'[{ts()}] rank_out={rank_out}, dim_ker={dim_ker} ({time.time()-t0:.1f}s)', flush=True)

    print(f'[{ts()}] rank of incoming I_k {I_k.nrows()}x{I_k.ncols()} ...', flush=True)
    t0 = time.time()
    dim_im = rank_fast(I_k, f'I_k[{k:+d}]')
    print(f'[{ts()}] dim_im={dim_im} ({time.time()-t0:.1f}s)', flush=True)

    dim_h = dim_ker - dim_im
    result = {
        'model': 'B_generic_nu_borel_reduced',
        'k': k,
        't_min': args.t_min,
        't_max': args.t_max,
        'a_max': args.a_max,
        'max_deg': args.max_deg,
        'dim_Ck': dim_k,
        'dim_ker': dim_ker,
        'dim_im': dim_im,
        'dim_H': dim_h,
        'im_subset_ker': True,
    }

    out_res = result_cache_path(args.checkpoint_dir, tag, k)
    with open(out_res, 'wb') as f:
        pickle.dump(result, f, protocol=4)
    print(f'[{ts()}] saved cohomology result -> {out_res}', flush=True)

    if args.save_linear_data:
        gk = groups[k]
        layout = [
            {
                'node': {'t': nd.t, 'a': nd.a, 'b': nd.b, 'c': nd.c},
                'd': d,
                'offset_full': gk.offsets[(nd, d)],
                'dim': gk.vermas[nd].dim(d),
            }
            for nd in gk.nodes
            for d in range(gk.max_deg + 1)
        ]
        lin = {
            'model': 'B_generic_nu_borel_reduced',
            'k': k,
            'O_k': O_k,
            'I_k': I_k,
            'basis_layout_full': layout,
            'kept_full_indices_k': red_data[k]['kept'],
            'dim_Ck_reduced': dim_k,
            'dim_O_rows': O_k.nrows(),
            'dim_I_cols': I_k.ncols(),
        }
        out_lin = linear_cache_path(args.checkpoint_dir, tag, k)
        with open(out_lin, 'wb') as f:
            pickle.dump(lin, f, protocol=4)
        print(f'[{ts()}] saved linear data -> {out_lin}', flush=True)

    print('\nSummary')
    print('-------')
    print(f'k={k:+d}')
    print(f'dim C^k (reduced) = {dim_k}')
    print(f'dim ker          = {dim_ker}')
    print(f'dim im           = {dim_im}')
    print(f'dim H^k          = {dim_h}')


def parse_args():
    p = argparse.ArgumentParser(
        description='Production H^1 for Borel-reduced Complex B at generic nu (no phi_1E).'
    )
    p.add_argument('--t-min', type=int, default=DEFAULT_T_MIN)
    p.add_argument('--t-max', type=int, default=DEFAULT_T_MAX)
    p.add_argument('--a-max', type=int, default=DEFAULT_A_MAX)
    p.add_argument('--max-deg', type=int, default=DEFAULT_MAX_DEG)
    p.add_argument('--level', type=int, default=DEFAULT_K,
                   help='Cochain level k to compute (default 1).')
    p.add_argument('--checkpoint-dir',
                   default=os.path.join(os.path.dirname(os.path.abspath(__file__)), 'checkpoints'))
    p.add_argument('--save-linear-data', action='store_true')
    return p.parse_args()


def main():
    args = parse_args()
    run(args)


if __name__ == '__main__':
    main()
