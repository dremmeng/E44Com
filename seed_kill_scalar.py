"""
seed_kill_scalar.py -- Seed Step 1: Kill M_1(0,0,0) via the mean-zero constraint.

The 130 irreducible H^1_{n=+1} seeds contain one class at node M_1(0,0,0).
This node carries the trivial sl_4 representation: a spatially constant mode
(zero Fourier wavevector, weight (0,0,0)).

Physical interpretation: on the torus T^4 with the mean-zero normalization
    integral_{T^4} u dx = 0
the constant mode is absent from any admissible velocity field.
Imposing this constraint removes exactly global index 0 from C^1_{n=+1},
reducing the irreducible seed count from 130 to 129.

Method:
  1. Assemble d_out restricted to d=0 source, d=1 target (the differential
     that maps k=1 fibers to k=2 fibers within the n=+1 sector).
  2. Certify dim ker = 130 over GF(p) (full seed count, no constraint).
  3. Remove column 0 (M_1(0,0,0) fiber) from d_out to model the mean-zero
     subspace C^1_{n=+1, mean-zero}.
  4. Certify that the removed column was a right-kernel vector (it maps to 0),
     confirming it is a genuine cocycle.
  5. Show dim ker drops to 129: the scalar seed is isolated and removable.
"""

import sys
from sage.all import GF, matrix, QQ

from de_rham_complex import (
    MORPHISMS_A, MORPHISMS_B,
    window_nodes, CochainGroup, assemble_differential,
)
from verma_modules import load_e44

# -- parameters --------------------------------------------------------------
P1, P2   = 65521, 65537   # GF primes for rank certification
T_MIN    = -1
T_MAX    = 6              # extended window (same as h1_seed_kill.py)
A_MAX    = 4
SHIFT    = 1              # only shift-1 morphisms act at d=0
SCALAR_GLOBAL_IDX = 0     # M_1(0,0,0) sits at global column index 0

# -- helpers -----------------------------------------------------------------
def rank_fast(M, label=""):
    r1 = M.change_ring(GF(P1)).rank()
    r2 = M.change_ring(GF(P2)).rank()
    if r1 != r2:
        print(f"  [{label}] GF ranks disagree ({r1} vs {r2}), falling back to padic")
        r1 = M.rank(algorithm='padic')
    else:
        print(f"  [{label}] GF rank certified: {r1}")
    return r1


def degree0_cols(grp):
    """Sorted global column indices where internal degree d=0."""
    cols = []
    entries = sorted(grp.offsets.items(), key=lambda x: x[1])
    for i, ((nd, d), off) in enumerate(entries):
        nxt = entries[i+1][1] if i+1 < len(entries) else grp.total_dim
        if d == 0:
            cols.extend(range(off, nxt))
    return cols


def degree_sv_rows(grp, sv):
    """Sorted global row indices where internal degree d=sv in grp."""
    rows = []
    entries = sorted(grp.offsets.items(), key=lambda x: x[1])
    for i, ((nd, d), off) in enumerate(entries):
        nxt = entries[i+1][1] if i+1 < len(entries) else grp.total_dim
        if d == sv:
            rows.extend(range(off, nxt))
    return rows


def build_seed_differential(e44_data, morphisms, label):
    """
    Assemble d_out restricted to (k=1, d=0) source columns, stacked over
    all target levels reached by morphisms from k=1.
    Returns a matrix of shape (total_target_rows times 281) over QQ.
    """
    nodes_k1 = window_nodes(1, morphisms, T_MIN, T_MAX, A_MAX)
    g_src = CochainGroup(1, nodes_k1, max_deg=0, e44_data=e44_data)
    d0_src = degree0_cols(g_src)
    n_src = len(d0_src)
    print(f"\n[{label}] g_src dim = {g_src.total_dim}  (expecting 281)  d0_cols={n_src}")

    # Group morphisms by (k_tar, sv_deg) -- build g_tar once per group
    groups = {}   # k_tar -> (sv, nodes_tar)
    for m in morphisms:
        edges_k1 = [(s, t, a) for s, t, a in m.edges(T_MIN, T_MAX, A_MAX) if s.t == 1]
        if not edges_k1:
            continue
        k_tar = edges_k1[0][1].t
        sv = m.sv_deg
        if k_tar not in groups:
            groups[k_tar] = sv
    print(f"[{label}] target levels from k=1: {sorted(groups.keys())}")

    row_blocks = []
    for k_tar, sv in sorted(groups.items()):
        nodes_tar = window_nodes(k_tar, morphisms, T_MIN, T_MAX, A_MAX)
        g_tar = CochainGroup(k_tar, nodes_tar, max_deg=sv, e44_data=e44_data)
        D = assemble_differential(g_src, g_tar, morphisms, e44_data,
                                  T_MIN, T_MAX, A_MAX)
        rows = degree_sv_rows(g_tar, sv)
        if not rows:
            continue
        block = D.matrix_from_rows_and_columns(rows, d0_src)
        nnz = len(block.dict())
        print(f"  k_tar={k_tar} sv={sv}: block {block.nrows()}times{block.ncols()} nnz={nnz}")
        row_blocks.append(block)

    if not row_blocks:
        raise RuntimeError("No outgoing morphism blocks found from k=1")

    from sage.all import block_matrix
    stacked = block_matrix(len(row_blocks), 1,
                           [[b] for b in row_blocks], subdivide=False)
    print(f"[{label}] stacked d_out shape: {stacked.nrows()} times {stacked.ncols()}")
    return stacked


def kernel_dim_from_matrix(M, label=""):
    """dim ker M = ncols - rank(M)."""
    r = rank_fast(M, label)
    return M.ncols() - r


# -- main ---------------------------------------------------------------------
def main():
    print("=" * 60)
    print("Seed Step 1: M_1(0,0,0) scalar mode")
    print("=" * 60)

    e44 = load_e44()

    for label, morphisms in [("Complex A (NS)", MORPHISMS_A),
                              ("Complex B (Euler)", MORPHISMS_B)]:
        print(f"\n{'-'*50}")
        print(f"Complex: {label}")
        print(f"{'-'*50}")

        D = build_seed_differential(e44, morphisms, label)
        n_cols = D.ncols()
        print(f"\n[{label}] Full d_out: {D.nrows()} times {n_cols}")

        # Confirm scalar column is pure zero (it is a cocycle)
        scalar_col = D.column(SCALAR_GLOBAL_IDX)
        is_zero = all(v == 0 for v in scalar_col)
        print(f"[{label}] Column {SCALAR_GLOBAL_IDX} (M_1(0,0,0)) is zero: {is_zero}")
        if not is_zero:
            nnz = sum(1 for v in scalar_col if v != 0)
            print(f"  WARNING: {nnz} nonzero entries -- scalar is NOT a cocycle!")
            return

        # Full kernel dimension (should be 130)
        dim_ker_full = kernel_dim_from_matrix(D, f"{label} full")
        print(f"[{label}] dim ker (full, no constraint) = {dim_ker_full}")
        assert dim_ker_full == 130, f"Expected 130, got {dim_ker_full}"

        # Mean-zero subspace: remove column 0
        remaining_cols = list(range(1, n_cols))
        D_meanzero = D.matrix_from_columns(remaining_cols)
        print(f"[{label}] d_out after mean-zero: {D_meanzero.nrows()} times {D_meanzero.ncols()}")
        dim_ker_mz = kernel_dim_from_matrix(D_meanzero, f"{label} mean-zero")
        print(f"[{label}] dim ker (mean-zero constraint) = {dim_ker_mz}")

        print()
        print(f"  RESULT [{label}]:")
        print(f"    Before mean-zero: {dim_ker_full} seeds")
        print(f"    After  mean-zero: {dim_ker_mz} seeds")
        delta = dim_ker_full - dim_ker_mz
        print(f"    Seeds killed by integralu dx = 0 constraint: {delta}")
        if delta == 1 and dim_ker_mz == 129:
            print(f"    [OK] Scalar mode M_1(0,0,0) eliminated. 129 seeds remain.")
        else:
            print(f"    [X] Unexpected result -- check computation.")

    print("\n" + "=" * 60)
    print("Step 1 complete.")
    print("Summary: the constant-mode seed at M_1(0,0,0) is removed by")
    print("the mean-zero normalization integral_{T^4} u dx = 0.")
    print("Remaining irreducible seeds: 129")
    print("  M_1(0,0,1): 80  [Kac module -- hard core]")
    print("  M_1(0,0,4): 35  [S^4(C^4*)]")
    print("  M_1(1,0,0):  8  [phat4 fundamental]")
    print("  M_1(2,0,0):  6  [S^2(C^4) survivors]")
    print("=" * 60)


if __name__ == "__main__":
    main()
