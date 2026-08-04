"""
borel_chain_map_check.py -- Step 1 executable checker for Borel 4->3 chain map.

Defines a concrete reduction map R_B by keeping only PBW basis states with:
  alpha_4 = 0 and d_4 not present in the wedge set S.

Then checks, for each tested (k -> k+s), that two reduced differentials agree:
  (B1) D_red = P_{k+s} * D_full * I_k
  (B2) D_red built directly from phi blocks, filtered to reduced coordinates.

And enforces the actual chain-map identity:
    P_{k+s} * D_full  =  D_red * P_k.

If all defects are zero matrices, Step 1 chain-map commutativity is certified
for the tested window and conventions.
"""

import argparse

from sage.all import QQ, matrix

import de_rham_complex as drc


def _log(msg):
    """Print with flush so long runs stream progress to logs/terminal."""
    print(msg, flush=True)


def keep_basis_state(mon, w_idx, node, d):
    """Borel 4->3 filter on PBW monomial state (fiber-unrestricted)."""
    alpha, sset = mon
    return alpha[3] == 0 and (3 not in sset)


def reduced_index_data(group):
    """
    Build reduced-index metadata for one CochainGroup.

    Returns:
      kept_globals: list of kept full-basis global indices, in reduced order
      full_to_red:  dict mapping full global index -> reduced global index
    """
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
    """
    Build projection P: full -> reduced and inclusion I: reduced -> full.

    P has shape (m, n), I has shape (n, m), where n=full_dim, m=len(kept).
    """
    m = len(kept_globals)
    p_entries = {(r, g): QQ(1) for r, g in enumerate(kept_globals)}
    i_entries = {(g, r): QQ(1) for r, g in enumerate(kept_globals)}

    P = matrix(QQ, m, full_dim, p_entries, sparse=True)
    I = matrix(QQ, full_dim, m, i_entries, sparse=True)
    return P, I


def assemble_reduced_direct(
    src_group,
    tar_group,
    morphisms,
    e44_data,
    t_min,
    t_max,
    a_max,
    src_full_to_red,
    tar_full_to_red,
    src_red_dim,
    tar_red_dim,
):
    """
    Assemble D_red directly from phi blocks, filtering to reduced coordinates.
    """
    entries = {}

    for spec in morphisms:
        for src_node, tar_node, phi_args in spec.edges(t_min, t_max, a_max):
            if src_node.t != src_group.k or tar_node.t != tar_group.k:
                continue
            if src_node not in src_group.nodes or tar_node not in tar_group.nodes:
                continue

            msd = drc.src_max_deg(spec.sv_deg, src_group.max_deg)
            src_data = e44_data if drc.node_fiber_type(src_node) == "phat4" else None

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
                i_red = tar_full_to_red.get(i_full)
                j_red = src_full_to_red.get(j_full)
                if i_red is None or j_red is None:
                    continue

                key = (i_red, j_red)
                prev = entries.get(key, QQ(0))
                entries[key] = prev + val

    return matrix(QQ, tar_red_dim, src_red_dim, entries, sparse=True)


def assemble_full_safe(src_group, tar_group, morphisms, e44_data, t_min, t_max, a_max):
    """
    Assemble the full differential, skipping any edge whose phi builder fails.

    This mirrors the production H^1 pipeline in cohomology_h1_fast.py, where
    failed/skipped edges are treated as absent rather than aborting the run.
    """
    entries = {}

    for spec in morphisms:
        for src_node, tar_node, phi_args in spec.edges(t_min, t_max, a_max):
            if src_node.t != src_group.k or tar_node.t != tar_group.k:
                continue
            if src_node not in src_group.nodes or tar_node not in tar_group.nodes:
                continue

            msd = drc.src_max_deg(spec.sv_deg, src_group.max_deg)
            src_data = e44_data if drc.node_fiber_type(src_node) == "phat4" else None

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


def run_check(cx, t_min, t_max, a_max, max_deg):
    e44_data = drc.load_e44()

    if cx == "A":
        morphisms = drc.MORPHISMS_A
    elif cx == "B":
        morphisms = drc.MORPHISMS_B
    else:
        morphisms = drc.MORPHISMS_A + drc.MORPHISMS_B

    shifts = sorted(set(spec.sv_deg for spec in morphisms))

    _log("=" * 72)
    _log(f"Step 1 chain-map check | cx={cx} | t=[{t_min},{t_max}] a_max={a_max} max_deg={max_deg}")
    _log("=" * 72)

    total_pairs = 0
    passed_pairs = 0

    for k_src in range(t_min, t_max + 1):
        _log(f"[k={k_src:+d}] building source group and reduced indices...")
        nodes_src = drc.window_nodes(k_src, morphisms, t_min, t_max, a_max)
        if not nodes_src:
            continue

        src_group = drc.CochainGroup(k_src, nodes_src, max_deg=max_deg, e44_data=e44_data)
        kept_src, src_map = reduced_index_data(src_group)
        P_src, I_src = projection_and_inclusion(src_group.total_dim, kept_src)

        # Internal consistency of P and I for source level
        id_src = P_src * I_src
        if not id_src.is_one():
            _log(f"[FAIL] k={k_src}: P_k * I_k is not identity on reduced source")
            return False

        for s in shifts:
            k_tar = k_src + s
            if k_tar < t_min or k_tar > t_max:
                continue

            nodes_tar = drc.window_nodes(k_tar, morphisms, t_min, t_max, a_max)
            if not nodes_tar:
                continue

            _log(f"  [pair k={k_src:+d}, s={s}] building target group k'={k_tar:+d}...")

            tar_group = drc.CochainGroup(k_tar, nodes_tar, max_deg=max_deg, e44_data=e44_data)
            kept_tar, tar_map = reduced_index_data(tar_group)
            P_tar, I_tar = projection_and_inclusion(tar_group.total_dim, kept_tar)

            id_tar = P_tar * I_tar
            if not id_tar.is_one():
                _log(f"[FAIL] k={k_tar}: P_k * I_k is not identity on reduced target")
                return False

            _log(f"  [pair k={k_src:+d}, s={s}] assembling full differential...")
            D_full = assemble_full_safe(
                src_group,
                tar_group,
                morphisms,
                e44_data,
                t_min,
                t_max,
                a_max,
            )

            # Construction B1
            _log(f"  [pair k={k_src:+d}, s={s}] projecting full differential...")
            D_red_from_projection = P_tar * D_full * I_src

            # Construction B2
            _log(f"  [pair k={k_src:+d}, s={s}] assembling reduced differential directly...")
            D_red_direct = assemble_reduced_direct(
                src_group,
                tar_group,
                morphisms,
                e44_data,
                t_min,
                t_max,
                a_max,
                src_map,
                tar_map,
                len(kept_src),
                len(kept_tar),
            )

            defect = D_red_from_projection - D_red_direct
            chain_defect = (P_tar * D_full) - (D_red_direct * P_src)
            total_pairs += 1

            nnz_defect = len(defect.dict())
            nnz_chain_defect = len(chain_defect.dict())
            nnz_full = len(D_full.dict())
            nnz_red = len(D_red_direct.dict())

            _log(
                f"k={k_src:>+2d} -> k={k_tar:>+2d} (s={s}): "
                f"full {D_full.nrows()}x{D_full.ncols()} nnz={nnz_full}, "
                f"red {D_red_direct.nrows()}x{D_red_direct.ncols()} nnz={nnz_red}, "
                f"build_defect_nnz={nnz_defect}, chain_defect_nnz={nnz_chain_defect}"
            )

            if defect.is_zero() and chain_defect.is_zero():
                passed_pairs += 1
            else:
                if not defect.is_zero():
                    _log("  [FAIL] Reduced-construction mismatch for this (k,s) pair.")
                if not chain_defect.is_zero():
                    _log("  [FAIL] Chain-map defect nonzero: P*D != D_red*P.")
                return False

    _log("-" * 72)
    _log(f"Passed pairs: {passed_pairs}/{total_pairs}")

    if total_pairs == 0:
        _log("[WARN] No (k,s) pairs were tested in this window.")
        return False

    _log("[OK] Step 1 chain-map check passed on tested window.")
    return True


def main():
    parser = argparse.ArgumentParser(description="Step 1 Borel chain-map checker")
    parser.add_argument("--cx", choices=["A", "B", "both"], default="A")
    parser.add_argument("--t-min", type=int, default=-1)
    parser.add_argument("--t-max", type=int, default=1)
    parser.add_argument("--a-max", type=int, default=1)
    parser.add_argument("--max-deg", type=int, default=1)
    args = parser.parse_args()

    ok = run_check(args.cx, args.t_min, args.t_max, args.a_max, args.max_deg)
    raise SystemExit(0 if ok else 1)


if __name__ == "__main__":
    main()
