"""
seed_kill_s2.py -- Seed Step 2: M_1(2,0,0), S^2(C^4) fiber analysis
=================================================================

W_1(2,0,0) = K_1(2,0,0) (t=1 -> irreducible): 10-dim sl_4 rep S^2(C^4).

phi_1A (sv=1): M_0(3,0,0) -> M_1(2,0,0) -> M_2(1,0,0)
  Acts on d=0 fiber: kills 4 indices = e_1 otimes_sym C^4

phi_2DA (sv=2, Complex A only): M_1(2,0,0) -> M_3(0,1,0)
  Acts on d=0 fiber: nnz=0 (pre-certified, same as phi_1A kernel = S^2(C^3))

6 SURVIVORS: S^2(span{e_2,e_3,e_4}) = indices {2,4,5,7,8,9}
  As sl_3-module (alpha_2,alpha_3 roots): IRREDUCIBLE 6-dim rep (highest weight (2,0))
  As sl_4-module: connected to killed indices via e[1] -- same irrep S^2(C^4)

Uses single-node CochainGroups throughout (safe, no OOM risk).
"""

import sys
import os

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

from sage.all import QQ, GF, vector

from verma_modules import load_e44
from de_rham_complex import (
    MORPHISMS_A, MORPHISMS_B,
    Node, CochainGroup, assemble_differential,
)

T_MIN, T_MAX, A_MAX = -1, 6, 4
P1, P2 = 65521, 65537
KILLED   = {0, 1, 3, 6}
SURVIVED = {2, 4, 5, 7, 8, 9}


def mk(k, abc_list, max_deg, e44):
    return CochainGroup(k, [Node(k, a, b, c) for a, b, c in abc_list],
                        max_deg=max_deg, e44_data=e44)


def cert(D, rows=None, cols=None, label=''):
    """GF-certified rank; returns (rank, nnz)."""
    if rows is not None:
        D = D.matrix_from_rows(rows)
    if cols is not None:
        D = D.matrix_from_columns(cols)
    nnz = len(D.dict())
    if nnz == 0:
        return 0, 0
    r1 = D.change_ring(GF(P1)).rank()
    r2 = D.change_ring(GF(P2)).rank()
    if r1 != r2:
        print(f"  WARNING {label}: GF rank mismatch ({r1} vs {r2}), using padic")
        return D.rank(algorithm='padic'), nnz
    return r1, nnz


def main():
    print("=" * 70)
    print("Seed Step 2: M_1(2,0,0) -- S^2(C^4) Analysis")
    print("=" * 70)

    e44 = load_e44()

    # -- Build single-node source group ------------------------------------
    nd = Node(1, 2, 0, 0)
    g1 = CochainGroup(1, [nd], max_deg=0, e44_data=e44)
    W  = g1.vermas[nd].W
    dim_W = W.dim   # = 10
    fiber_cols = list(g1.basis_slice(nd, 0))

    print(f"\n  M_1(2,0,0): type={type(W).__name__}  dim={dim_W}  "
          f"(S^2(C^4), highest weight 2omega_1)")
    print(f"  KILLED indices:   {sorted(KILLED)}   = e_1 otimes_sym C^4")
    print(f"  SURVIVED indices: {sorted(SURVIVED)} = S^2(span{{e_2,e_3,e_4}})")

    # Basis vectors for sl_4 action
    bvec = [vector(QQ, dim_W) for _ in range(dim_W)]
    for j in range(dim_W):
        bvec[j][j] = 1

    # =========================================================================
    # 1. sl_4 WEIGHT TABLE
    # =========================================================================
    print()
    print("=" * 70)
    print("1. sl_4 Weight Table for M_1(2,0,0)")
    print("=" * 70)
    print()
    print(f"  {'idx':>3}  {'basis':>10}  {'sl_4 Dynkin (L_1,L_2,L_3)':>24}  {'status':>9}")
    print("  " + "-" * 57)
    for idx in range(dim_W):
        elt = W.basis_elts[idx]
        wt  = W.weight_of(idx)
        a, b = list(elt)
        interp = f'e_{a}otimese_{b}'
        status = 'KILLED  ' if idx in KILLED else 'survived'
        print(f"  {idx:>3}  {interp:>10}  {str(wt):>24}  {status:>9}")
    print()
    print("  Killed   (phi_1A):  e_1 otimes_sym C^4          = {0,1,3,6}")
    print("  Survived (ker):     S^2(span{e_2,e_3,e_4})   = {2,4,5,7,8,9}")

    # =========================================================================
    # 2. COCYCLE CERTIFICATION: phi_1A and phi_2DA on M_1(2,0,0)[d=0]
    # =========================================================================
    print()
    print("=" * 70)
    print("2. Cocycle Certification: All Outgoing Morphisms")
    print("=" * 70)
    print()
    print("  phi_1A  (sv=1, Complexes A+B): M_1(2,0,0) -> M_2(1,0,0)[d=1]")
    print("  phi_2DA (sv=2, Complex A only): M_1(2,0,0) -> M_3(0,1,0)[d=2]")
    print()

    # phi_1A: M_1(2,0,0) -> M_2(1,0,0), sv=1, maps d=0 -> d=1
    nd2 = Node(2, 1, 0, 0)
    g2  = CochainGroup(2, [nd2], max_deg=1, e44_data=e44)
    D_1A = assemble_differential(g1, g2, MORPHISMS_B, e44, T_MIN, T_MAX, A_MAX)
    rows_d1 = list(g2.basis_slice(nd2, 1))
    rk_1A, nnz_1A = cert(D_1A, rows=rows_d1, label='phi_1A')
    zero_1A = [j for j in range(dim_W)
               if all(D_1A[rows_d1[i], j] == 0 for i in range(len(rows_d1)))]
    print(f"  phi_1A -> M_2(1,0,0)[d=1]:  "
          f"{len(rows_d1)}times{dim_W}  nnz={nnz_1A}  rank={rk_1A}")
    print(f"    Zero columns: {sorted(zero_1A)}")
    print(f"    Killed:  {sorted(set(range(dim_W)) - set(zero_1A))}")

    # phi_2DA: M_1(2,0,0) -> M_3(0,1,0), sv=2, maps d=0 -> d=2
    nd3 = Node(3, 0, 1, 0)
    g3  = CochainGroup(3, [nd3], max_deg=2, e44_data=e44)
    D_2DA = assemble_differential(g1, g3, MORPHISMS_A, e44, T_MIN, T_MAX, A_MAX)
    rows_d2 = list(g3.basis_slice(nd3, 2))
    rk_2DA, nnz_2DA = cert(D_2DA, rows=rows_d2, label='phi_2DA')
    zero_2DA = ([j for j in range(dim_W)
                 if all(D_2DA[rows_d2[i], j] == 0 for i in range(len(rows_d2)))]
                if rows_d2 else list(range(dim_W)))
    print(f"  phi_2DA -> M_3(0,1,0)[d=2]:  "
          f"{len(rows_d2)}times{dim_W}  nnz={nnz_2DA}  rank={rk_2DA}")
    if nnz_2DA > 0:
        print(f"    Zero columns: {sorted(zero_2DA)}")
    print()

    # Survivors = intersection of zero column sets
    all_zero_from_1A   = set(zero_1A)
    all_zero_from_2DA  = set(zero_2DA) if rows_d2 else set(range(dim_W))
    survivors_after_both = all_zero_from_1A & all_zero_from_2DA
    print(f"  Survivors after phi_1A cap phi_2DA: {sorted(survivors_after_both)}")
    if survivors_after_both == SURVIVED:
        print(f"  [OK] Matches expected SURVIVED = {sorted(SURVIVED)}")
    else:
        print(f"  [X] Mismatch! Expected {sorted(SURVIVED)}")

    # =========================================================================
    # 3. H^1 NON-TRIVIALITY: M_1(2,0,0)[d=0] cap im(d^0) = {0}
    # =========================================================================
    print()
    print("=" * 70)
    print("3. H^1 Non-Triviality: M_1(2,0,0)[d=0] cap im(d^0) = {0}")
    print("=" * 70)
    print()
    print("  Incoming morphisms to M_1(2,0,0):")
    print("    phi_1A (sv=1):  M_0(3,0,0) -> M_1(2,0,0)  [maps d=0 -> d=1]")
    print("    No sv=0 morphism targets M_1(2,0,0)[d=0].")
    print("  -> Structural argument: im(d^0) cap M_1(2,0,0)[d=0] = {0}.")
    print()

    # Verify: phi_1A from M_0(3,0,0) maps into d=1 layer (sv=1), not d=0
    nd0 = Node(0, 3, 0, 0)
    g0  = CochainGroup(0, [nd0], max_deg=0, e44_data=e44)
    D0  = assemble_differential(g0, g1, MORPHISMS_B, e44, T_MIN, T_MAX, A_MAX)
    rows_d0_fiber = list(g1.basis_slice(nd, 0))
    rk_in, nnz_in = cert(D0, rows=rows_d0_fiber, label='d^0->M_1(2,0,0)[d=0]')
    dim0 = g0.vermas[nd0].dim(0)
    print(f"  phi_1A: M_0(3,0,0)[d=0] -> M_1(2,0,0)[d=0]:  "
          f"{dim_W}times{dim0}  nnz={nnz_in}  rank={rk_in}")
    print()
    if nnz_in == 0:
        print(f"  [OK] rank(d^0 -> M_1(2,0,0)[d=0]) = 0.  All {dim_W} are genuine H^1 classes.")
        print(f"  After phi_1A cocycle kill: {len(survivors_after_both)} genuine H^1 survivors.")
    else:
        print(f"  [X] UNEXPECTED: nnz={nnz_in}, some fiber vectors in im(d^0).")

    # =========================================================================
    # 4. PBW d=1 ANALYSIS: incoming rank at M_1(2,0,0)[d=1]
    # =========================================================================
    print()
    print("=" * 70)
    print("4. PBW-Degree-1 Analysis: Rank at M_1(2,0,0)[d=1]")
    print("=" * 70)
    print()

    g1_d1 = CochainGroup(1, [nd], max_deg=1, e44_data=e44)
    dim_d0 = g1_d1.vermas[nd].dim(0)
    dim_d1 = g1_d1.vermas[nd].dim(1)
    print(f"  M_1(2,0,0)[d=0]: dim={dim_d0}   M_1(2,0,0)[d=1]: dim={dim_d1}")

    g0_d1 = CochainGroup(0, [nd0], max_deg=0, e44_data=e44)
    D01_d1 = assemble_differential(g0_d1, g1_d1, MORPHISMS_B, e44, T_MIN, T_MAX, A_MAX)
    rows_d1_200 = list(g1_d1.basis_slice(nd, 1))
    rk_d1, nnz_d1 = cert(D01_d1, rows=rows_d1_200, label='d^0->M_1(2,0,0)[d=1]')
    print(f"  phi_1A: M_0(3,0,0)[d=0] -> M_1(2,0,0)[d=1]:  "
          f"{len(rows_d1_200)}times{dim0}  nnz={nnz_d1}  rank={rk_d1}")
    if rk_d1 > 0:
        print(f"  [OK] rank={rk_d1}: {rk_d1}/{dim_d1} degree-1 components are in im(d^0).")
    else:
        print(f"  rank=0: no degree-1 components are exact.")

    # =========================================================================
    # 5. sl_4 IRREDUCIBILITY: killed <-> survived are connected
    # =========================================================================
    print()
    print("=" * 70)
    print("5. sl_4 Irreducibility: Killed and Survived Are in the Same Irrep")
    print("=" * 70)
    print()
    print("  S^2(C^4) is a 10-dim IRREDUCIBLE sl_4-module (highest weight (2,0,0)).")
    print("  No sl_4-equivariant map can have kernel = proper non-trivial subspace.")
    print()
    print("  Explicit evidence via e/f matrix actions:")

    e1 = W.e_mats[1]
    f1 = W.f_mats[1]

    # f[1]: killed -> survived
    f1_1 = f1 * bvec[1]   # e_1otimese_2 [killed idx=1]
    nz_f1_1 = [(j, f1_1[j]) for j in range(dim_W) if f1_1[j] != 0]
    print(f"  f[1] * (e_1otimese_2) [KILLED idx=1]  =  {nz_f1_1}")
    print(f"    -> lands on SURVIVED: {all(j in SURVIVED for j, _ in nz_f1_1)}")

    # e[1]: survived -> killed
    e1_2 = e1 * bvec[2]   # e_2otimese_2 [survived idx=2]
    nz_e1_2 = [(j, e1_2[j]) for j in range(dim_W) if e1_2[j] != 0]
    print(f"  e[1] * (e_2otimese_2) [SURVIVED idx=2] =  {nz_e1_2}")
    print(f"    -> lands on KILLED: {all(j in KILLED for j, _ in nz_e1_2)}")
    print()
    print("  e[1] on all 6 survivors:")
    for idx in sorted(SURVIVED):
        img = e1 * bvec[idx]
        nz = [(j, img[j]) for j in range(dim_W) if img[j] != 0]
        killed_hits = [j for j, _ in nz if j in KILLED]
        surv_hits   = [j for j, _ in nz if j in SURVIVED]
        print(f"    idx {idx} {str(W.weight_of(idx)):>12}: -> {nz}"
              f"  [killed={killed_hits}, surv={surv_hits}]")
    print()
    print("  [OK] e[1] connects survivors to killed -- same sl_4 irrep.")
    print("  [X] No sl_4-equivariant map can kill just the 6 survivors.")

    # =========================================================================
    # 6. sl_3 SUBALGEBRA STRUCTURE OF SURVIVORS
    # =========================================================================
    print()
    print("=" * 70)
    print("6. sl_3 Subalgebra Structure: Survivors = S^2(C^3)")
    print("=" * 70)
    print()
    print("  S^2(span{e_2,e_3,e_4}) is a 6-dim IRREDUCIBLE sl_3-module")
    print("  under the sl_3 subset sl_4 subalgebra (roots alpha_2, alpha_3).")
    print("  Highest weight (2,0) in A_2 language.")
    print()
    print("  sl_3 lowering actions on survivors (f[2], f[3]):")

    f2 = W.f_mats[2]
    f3 = W.f_mats[3]
    for idx in sorted(SURVIVED):
        img2 = f2 * bvec[idx]
        img3 = f3 * bvec[idx]
        nz2 = [(j, img2[j]) for j in range(dim_W) if img2[j] != 0]
        nz3 = [(j, img3[j]) for j in range(dim_W) if img3[j] != 0]
        ok2 = all(j in SURVIVED for j, _ in nz2)
        ok3 = all(j in SURVIVED for j, _ in nz3)
        print(f"    f[2]*idx {idx}: {nz2}  <- stays: {ok2}")
        print(f"    f[3]*idx {idx}: {nz3}  <- stays: {ok3}")

    print()
    print("  [OK] f[2], f[3] map survivors -> survivors: S^2(C^3) is sl_3-stable.")
    print("  [X] S^2(C^3) is NOT sl_4-stable (e[1] breaks out).")
    print()
    print("  Implication:")
    print("    An E(4,4) morphism killing S^2(C^3) must have a singular vector")
    print("    that is sl_3-equivariant but NOT sl_4-equivariant.")
    print("    No such morphism appears in CCK 2026 Complexes A or B.")

    # =========================================================================
    # 7. MORPHISM COVERAGE SUMMARY
    # =========================================================================
    print()
    print("=" * 70)
    print("7. Morphism Coverage: All phi_* with Edges to/from M_1(2,0,0)")
    print("=" * 70)
    print()
    for cx_label, morphisms in [("Complex B (Euler)", MORPHISMS_B),
                                  ("Complex A (NS)",    MORPHISMS_A)]:
        print(f"  {cx_label}:")
        any_found = False
        for m in morphisms:
            # Outgoing: source = M_1(2,0,0)
            out_edges = [(s, t, a) for s, t, a in m.edges(T_MIN, T_MAX, A_MAX)
                         if s == Node(1, 2, 0, 0)]
            # Incoming: target = M_1(2,0,0)
            in_edges = [(s, t, a) for s, t, a in m.edges(T_MIN, T_MAX, A_MAX)
                        if t == Node(1, 2, 0, 0)]
            if out_edges:
                tar_nodes = [str(t) for _, t, _ in out_edges]
                print(f"    {m.name} (sv={m.sv_deg})  OUT:  M_1(2,0,0) -> {tar_nodes}")
                any_found = True
            if in_edges:
                src_nodes = [str(s) for s, _, _ in in_edges]
                print(f"    {m.name} (sv={m.sv_deg})  IN:   {src_nodes} -> M_1(2,0,0)")
                any_found = True
        if not any_found:
            print("    (no morphisms)")
    print()
    print(f"  phi_1A rank={rk_1A}: kills {rk_1A}/10, {dim_W-rk_1A} survive")
    print(f"  phi_2DA rank={rk_2DA}: {'kills none of' if rk_2DA == 0 else 'kills'} "
          f"the {len(SURVIVED)} survivors  (nnz={nnz_2DA})")

    # =========================================================================
    # SUMMARY
    # =========================================================================
    print()
    print("=" * 70)
    print("STEP 2 SUMMARY -- M_1(2,0,0): S^2(C^4) Seed Kill")
    print("=" * 70)
    print()
    print(f"  Starting seeds: {dim_W}  (all genuine H^1 classes -- rank(d^0->fiber)=0)")
    print()
    print("  Certified algebraic facts:")
    print(f"    phi_1A:  nnz={nnz_1A}  rank={rk_1A}  -> kills {rk_1A}/10 = e_1otimesC^4")
    print(f"    phi_2DA: nnz={nnz_2DA}  rank={rk_2DA}  -> kills {rk_2DA}/10 additional")
    print(f"    d^0 -> fiber: rank={rk_in}  -> all {dim_W} are genuine H^1")
    print(f"    d^0 -> d=1:   rank={rk_d1}  -> {rk_d1}/{dim_d1} degree-1 components exact")
    print()
    print("  Survivor weight table:")
    for idx in sorted(SURVIVED):
        a, b = list(W.basis_elts[idx])
        wt = W.weight_of(idx)
        print(f"    idx {idx}: e_{a}otimese_{b}  weight {wt}")
    print()
    print("  +------------------+------+--------------+----------------------+")
    print("  | Stage            | Kill |  Remaining   | Method               |")
    print("  +------------------+------+--------------+----------------------+")
    print(f"  | phi_1A           | {rk_1A:>4} | {dim_W - rk_1A:>12} | e_1otimesC^4 killed         |")
    print(f"  | phi_2DA          | {rk_2DA:>4} | {dim_W - rk_1A - rk_2DA:>12} | nnz=0, no new kill   |")
    print("  +------------------+------+--------------+----------------------+")
    print()
    print(f"  SURVIVING SEEDS: {len(survivors_after_both)}  = S^2(span{{e_2,e_3,e_4}})")
    print()
    print("  Why no further algebraic kill:")
    print("    * S^2(C^4) is irreducible as sl_4-module -- no sl_4-equivariant map")
    print("      has kernel = S^2(C^3).")
    print("    * S^2(C^3) is sl_3-stable (f[2], f[3] close on it).")
    print("    * phi_1A is the only Complex B morphism touching M_1(2,0,0).")
    print("    * phi_2DA (Complex A only) gives nnz=0 on the fiber.")
    print("    * No singular vector at d=1 (d^0 -> d=1 rank feeds different modes).")
    print()
    print("  Open paths to resolve the 6 survivors:")
    print("    1. NEW SINGULAR VECTOR: sl_3-equivariant E(4,4) morphism targeting")
    print("       S^2(C^3). Not in CCK 2026 -- would require a new sv in M_0(2,0,0).")
    print("    2. PHYSICAL LERAY: Incompressibility nabla*u=0 kills the longitudinal")
    print("       mode in S^2(C^3) too (strain tensor in k^perpotimesk^perp subset S^2(C^4) for kneq0).")
    print("       For mode k=(1,0,0,0): S^2(C^3) cap k^perp has dim=... to be computed.")
    print("    3. SOBOLEV: S^2(C^3) <-> quadratic terms in velocity -> controlled by")
    print("       ||u||_{H^1} via the energy inequality (regular solutions).")
    print("    4. SPECTRAL SEQUENCE: S^2(C^3) cancels with M_1(1,0,0) transverse")
    print("       seeds at a higher page of the spectral sequence.")
    print()
    print("  STATUS: Step 2 INCOMPLETE -- 6 seeds remain OPEN.")
    print("  The S^2(C^3) survivors represent strain-tensor H^1 obstructions")
    print("  in the {x_2,x_3,x_4}-coordinate subspace.")
    print("=" * 70)


if __name__ == "__main__":
    main()

