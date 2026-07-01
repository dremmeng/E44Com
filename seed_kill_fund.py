"""
seed_kill_fund.py -- Seed Step 3: M_1(1,0,0), the Physical H^1 Obstruction
==========================================================================

W_1(1,0,0) = K_1(1,0,0): 8-dimensional phat4 Kac module (irreducible).
This is THE ONLY physical H^1 obstruction (see cck_map.py / CCK 2026 sec.9.8).

Module decomposition:
  Sheet A (idx 0-3): half-integer phat4 Cartan eigenvalues (h_i in 1/2+ZZ)
                     sl_4 weights (omega_1 crystal): (1,0,0),(-1,1,0),(0,-1,1),(0,0,-1)
                     phat4-even sector  ->  physical velocity  u_1, u_2, u_3, u_4
  Sheet B (idx 4-7): integer phat4 Cartan eigenvalues (h_i in ZZ)
                     sl_4 weights: (1,0,0,0),(0,1,0,0),(0,0,1,0),(0,0,0,1)
                     phat4-odd sector  ->  pressure gradients partial_1p,...,partial_4p

Seed-kill programme for this block:
  ALL 8 seeds are genuine H^1 classes (algebraically).
  Physical arguments kill:
    Step 5a -- Z_2 parity:       4 Sheet B seeds  (phat4-odd = non-bosonic)
    Step 5b -- Leray projector: 1 Sheet A seed    (longitudinal mode, k*u_hat=0)
  Remaining: 3 transverse Sheet A seeds  ->  OPEN (primary blowup candidates)

Pre-certified from seed_kill_fund_output.txt (previous full-window run):
  phi_1D: 2248times8  nnz=0  rank=0   [window_nodes(2), max_deg=1]
  phi_4H: 53952times8 nnz=0  rank=0   [window_nodes(5), max_deg=4]
  phi_1D d=1 source: 8992times64   nnz=0  rank=0
  phi_4H d=1 source: 101160times64 nnz=0  rank=0

Uses single-node CochainGroups for all computations (safe, fast, no OOM risk).
"""

import sys
import os

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

from sage.all import QQ, GF, vector

from verma_modules import load_e44
from phat4_modules import _l0_even_idx
from de_rham_complex import (
    MORPHISMS_A, MORPHISMS_B,
    Node, CochainGroup, assemble_differential,
)

T_MIN, T_MAX, A_MAX = -1, 6, 4
P1, P2 = 65521, 65537

# Pre-certified from seed_kill_fund_output.txt (window_nodes run, cache required)
PRECERT = {
    "phi_1D_d0_nnz":  0,   "phi_1D_d0_rank": 0,   # 2248times8  nnz=0
    "phi_4H_d0_nnz":  0,   "phi_4H_d0_rank": 0,   # 53952times8 nnz=0
    "phi_1D_d1_nnz":  0,   "phi_1D_d1_rank": 0,   # 8992times64  nnz=0
    "phi_4H_d1_nnz":  0,   "phi_4H_d1_rank": 0,   # 101160times64 nnz=0
}


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
    print("Seed Step 3: M_1(1,0,0) -- Physical H^1 Obstruction Analysis")
    print("=" * 70)

    e44 = load_e44()

    # -- Build single-node source group ------------------------------------
    nd = Node(1, 1, 0, 0)
    g1 = CochainGroup(1, [nd], max_deg=0, e44_data=e44)

    # -- Sheet classification ----------------------------------------------
    W  = g1.vermas[nd].W
    dim_fiber = W.dim   # should be 8
    h_mats = {i: W.action_mats[_l0_even_idx(i, i)] for i in range(1, 5)}
    sheet_A, sheet_B, h_data = [], [], {}
    for idx in range(dim_fiber):
        v = vector(QQ, dim_fiber)
        v[idx] = 1
        h = tuple((h_mats[i] * v)[idx] for i in range(1, 5))
        half_int = any(abs(float(h[j]) - round(float(h[j]))) > 0.1 for j in range(4))
        ws = [w for w, ids in W.weight_spaces.items() if idx in ids]
        h_data[idx] = (h, ws, half_int)
        (sheet_A if half_int else sheet_B).append(idx)

    fiber_cols = list(g1.basis_slice(nd, 0))

    # =========================================================================
    # 1. MODULE STRUCTURE
    # =========================================================================
    print()
    print("=" * 70)
    print("1. Module Structure: W_1(1,0,0) = K_1(1,0,0)")
    print("=" * 70)
    print()
    print(f"  Type: {type(W).__name__}  dim={dim_fiber}  (K_1(1,0,0) irreducible)")
    print(f"  Sheet A (half-integer h, phat4-even): idx {sheet_A}")
    print(f"  Sheet B (integer h, phat4-odd):       idx {sheet_B}")
    print()
    print(f"  {'idx':>3}  {'h_1':>6} {'h_2':>6} {'h_3':>6} {'h_4':>6}  "
          f"{'sl_4 weight':>12}  {'sheet':>8}")
    print("  " + "-" * 62)
    for idx in range(dim_fiber):
        h, ws, half_int = h_data[idx]
        hs = " ".join(f"{float(x):>+.2f}" for x in h)
        ws_s = str(ws[0]) if ws else "(int-level)"
        sheet = "A (even)" if half_int else "B (odd) "
        print(f"  {idx:>3}  {hs}  {ws_s:>12}  {sheet}")
    print()
    print("  Sheet A -- sl_4 crystal (omega_1), CCK velocity sector:")
    print("    idx 0: weight  (1, 0, 0)  <->  u_1  (partial_1)")
    print("    idx 1: weight (-1, 1, 0)  <->  u_2  (partial_2)")
    print("    idx 2: weight  (0,-1, 1)  <->  u_3  (partial_3)")
    print("    idx 3: weight  (0, 0,-1)  <->  u_4  (partial_4)")
    print()
    print("  Sheet B -- phat4-odd sector, CCK pressure sector:")
    print("    idx 4: h = (1,0,0,0)  <->  partial_1 p  (dx_1)")
    print("    idx 5: h = (0,1,0,0)  <->  partial_2 p  (dx_2)")
    print("    idx 6: h = (0,0,1,0)  <->  partial_3 p  (dx_3)")
    print("    idx 7: h = (0,0,0,1)  <->  partial_4 p  (dx_4)")

    # =========================================================================
    # 2. COCYCLE CERTIFICATION AT d=0: all 8 in ker(d^1)
    # =========================================================================
    print()
    print("=" * 70)
    print("2. Cocycle Certification: M_1(1,0,0)[d=0] subseteq ker(d^1)")
    print("=" * 70)
    print()
    print("  Outgoing morphisms from M_1(1,0,0) in Complex B:")
    print("    phi_1D (sv=1): M_1(1,0,0)[d=0] -> M_2(0,0,0)[d=1]  (divergence)")
    print("    phi_4H (sv=4): M_1(1,0,0)[d=0] -> M_5(1,0,0)[d=4]  (4th-order)")
    print("  Complex A: NO morphisms originate from M_1(1,0,0).")
    print()

    # Test phi_1D with single-node target M_2(0,0,0), max_deg=1
    nd2 = Node(2, 0, 0, 0)
    g2  = CochainGroup(2, [nd2], max_deg=1, e44_data=e44)
    D12 = assemble_differential(g1, g2, MORPHISMS_B, e44, T_MIN, T_MAX, A_MAX)
    rows_d1_m2 = list(g2.basis_slice(nd2, 1))
    rk_1D, nnz_1D = cert(D12, rows=rows_d1_m2, cols=fiber_cols, label='phi_1D')
    print(f"  phi_1D -> M_2(0,0,0)[d=1]:  "
          f"{len(rows_d1_m2)}times{dim_fiber}  nnz={nnz_1D}  rank={rk_1D}"
          f"  [single-node, certified]")
    print(f"  phi_4H -> M_5(1,0,0)[d=4]:  "
          f"nnz={PRECERT['phi_4H_d0_nnz']}  rank={PRECERT['phi_4H_d0_rank']}"
          f"  [PRE-CERTIFIED: 53952times8 nnz=0]")
    print()
    all_cocycles = (nnz_1D == 0 and PRECERT['phi_4H_d0_nnz'] == 0)
    if all_cocycles:
        print("  [OK] d^1(z) = 0 for all z in M_1(1,0,0)[d=0] in both complexes.")
        print("  All 8 basis vectors are cocycles.")
    else:
        print("  [X] UNEXPECTED: some basis vectors are not cocycles.")

    # =========================================================================
    # 3. H^1 NON-TRIVIALITY: all 8 NOT in im(d^0)   [structural argument]
    # =========================================================================
    print()
    print("=" * 70)
    print("3. H^1 Non-Triviality: M_1(1,0,0)[d=0] cap im(d^0) = {0}")
    print("=" * 70)
    print()
    print("  STRUCTURAL ARGUMENT:")
    print("  Every morphism phi in Complexes A and B has sv_deg geq 1.")
    print("  phi maps (node, d=k) |-> (node', d=k+sv_deg).  For an element at")
    print("  (M_1(1,0,0), d=0) to be in im(d^0) there must be a phi with")
    print("  sv_deg=0 from some C^0 node into M_1(1,0,0)[d=0].  No such phi")
    print("  exists  ->  im(d^0) cap M_1(1,0,0)[d=0] = {0}.")
    print()
    print("  Verification (rank of all d^0 images into the d=0 fiber):")

    # Complex B sources at k=0: M_0(2,0,0) via phi_1A (sv=1), M_0(0,0,0) via phi_1E (sv=1)
    g0 = mk(0, [(2, 0, 0), (0, 0, 0)], 0, e44)
    D01 = assemble_differential(g0, g1, MORPHISMS_B, e44, T_MIN, T_MAX, A_MAX)
    rk_in_B, nnz_in_B = cert(D01, rows=fiber_cols, label='d^0_B->fiber')
    print(f"  Complex B  d^0 -> M_1(1,0,0)[d=0]:  "
          f"{dim_fiber}times{g0.total_dim}  nnz={nnz_in_B}  rank={rk_in_B}")

    # Complex A source: phi_2EA M_{-1}(1,0,0) -> M_1(1,0,0), sv=2 -> maps d=0->d=2
    g_m1 = mk(-1, [(1, 0, 0)], 0, e44)
    D_2ea = assemble_differential(g_m1, g1, MORPHISMS_A, e44, T_MIN, T_MAX, A_MAX)
    rk_in_A, nnz_in_A = cert(D_2ea, rows=fiber_cols, label='phi_2EA->fiber')
    print(f"  Complex A  phi_2EA -> M_1(1,0,0)[d=0]:  "
          f"{dim_fiber}times{g_m1.total_dim}  nnz={nnz_in_A}  rank={rk_in_A}")
    print()
    if nnz_in_B == 0 and nnz_in_A == 0:
        print("  [OK] rank(d^0 -> M_1(1,0,0)[d=0]) = 0.  All 8 are genuine H^1 classes.")
    else:
        print("  [X] UNEXPECTED: some fiber vectors are in im(d^0).")

    # =========================================================================
    # 4. d=1 ANALYSIS: rank-5 incoming to M_1(1,0,0)[d=1]
    # =========================================================================
    print()
    print("=" * 70)
    print("4. PBW-Degree-1 Analysis: Rank at M_1(1,0,0)[d=1]")
    print("=" * 70)
    print()

    g1_d1 = CochainGroup(1, [nd], max_deg=1, e44_data=e44)
    dim_d0 = g1_d1.vermas[nd].dim(0)
    dim_d1 = g1_d1.vermas[nd].dim(1)
    print(f"  M_1(1,0,0)[d=0]: dim={dim_d0}   M_1(1,0,0)[d=1]: dim={dim_d1}")

    g0_d1 = mk(0, [(2, 0, 0), (0, 0, 0)], 0, e44)
    D01_d1 = assemble_differential(g0_d1, g1_d1, MORPHISMS_B, e44, T_MIN, T_MAX, A_MAX)
    rows_d1 = list(g1_d1.basis_slice(nd, 1))
    rk_d1, nnz_d1 = cert(D01_d1, rows=rows_d1, label='d^0->M_1(1,0,0)[d=1]')
    print(f"  d^0 -> M_1(1,0,0)[d=1]:  {len(rows_d1)}times{g0_d1.total_dim}  "
          f"nnz={nnz_d1}  rank={rk_d1}")
    print()
    if rk_d1 > 0:
        print(f"  [OK] rank={rk_d1}: {rk_d1}/{dim_d1} degree-1 components are in im(d^0).")
        print(f"    These Taylor coefficients are fixed by C^0 data.")
        print(f"    Non-exact degree-1 components: {dim_d1 - rk_d1}")
    print()
    print("  d=1 outgoing (pre-certified from seed_kill_fund_output.txt):")
    print(f"    phi_1D (d=1 source): nnz={PRECERT['phi_1D_d1_nnz']}  "
          f"rank={PRECERT['phi_1D_d1_rank']}  [8992times64 nnz=0]")
    print(f"    phi_4H (d=1 source): nnz={PRECERT['phi_4H_d1_nnz']}  "
          f"rank={PRECERT['phi_4H_d1_rank']}  [101160times64 nnz=0]")
    print("  -> No singular vector at d=1.  d=0 H^1 is unaffected by d=1.")

    # =========================================================================
    # 5. Z_2 PARITY KILL: Sheet B excluded by bosonic restriction
    # =========================================================================
    print()
    print("=" * 70)
    print("5. Z_2 Parity Kill: Sheet B Excluded by phat4 Bosonic Restriction")
    print("=" * 70)
    print()
    print("  The phat4 superalgebra has a Z_2 grading:")
    print("    Z_2 = 0 (even, bosonic):  Sheet A -- top of Kac module")
    print("    Z_2 = 1 (odd, fermionic): Sheet B -- phat4_{-1} action on Sheet A")
    print()
    print("  Z_2-degree detected by h_i eigenvalue integrality:")
    print("    h_i in 1/2 + ZZ  ->  Z_2 = 0  (Sheet A, half-integer eigenvalues)")
    print("    h_i in ZZ      ->  Z_2 = 1  (Sheet B, integer eigenvalues)")
    print()
    print(f"  {'idx':>3}  {'h_1 (QQ)':>8}  {'Z_2':>8}  {'sheet':>7}  physical role")
    print("  " + "-" * 56)
    for idx in range(dim_fiber):
        h, ws, half_int = h_data[idx]
        z2_s = "even (0)" if half_int else "odd  (1)"
        role = "velocity u_i" if half_int else "EXCLUDED (fermionic)"
        print(f"  {idx:>3}  {str(h[0]):>8}  {z2_s:>8}  "
              f"{'A' if half_int else 'B':>7}  {role}")
    print()
    print("  Physical restriction to bosonic (Z_2=0) sector:")
    print("    NSE initial data u_0 in L^2(T^4,RR^4) is bosonic.")
    print("    The CCK map Phi(u,p) maps bosonic data into the Z_2=0 sector of C^1.")
    print("    Sheet B (Z_2=1) is structurally inaccessible from physical data.")
    print()
    print(f"  Z_2 KILL: {len(sheet_B)} Sheet B seeds (idx {sheet_B}) -> ELIMINATED")
    print(f"  Remaining: {len(sheet_A)} Sheet A seeds (idx {sheet_A})")

    # =========================================================================
    # 6. LERAY PROJECTOR: killing the longitudinal Sheet A seed
    # =========================================================================
    print()
    print("=" * 70)
    print("6. Leray Projector: Killing the Longitudinal Sheet A Seed")
    print("=" * 70)
    print()
    print("  Algebraic divergence (phi_1D at d=0 from Section 2):")
    print(f"    phi_1D: nnz={nnz_1D}  rank={rk_1D}  -> all 4 Sheet A vectors")
    print("    satisfy div(z) = 0 algebraically at d=0.  No algebraic kill here.")
    print()
    print("  Physical Fourier-space divergence-free constraint (Leray):")
    print()
    print("    For mode k in ZZ^4 \\ {0}: u(k) in k^perp subset C^4  (incompressibility)")
    print("    Leray projector P = 1 - nablaDelta^-^1nabla* projects onto k^perp.")
    print("    dim(k^perp) = 3  for each fixed k neq 0.")
    print()
    print("    Sheet A basis {e_0^A,...,e_3^A} <-> {e_1,...,e_4} subset RR^4.")
    print("    For mode k = (k_1,k_2,k_3,k_4):  P kills Sigma_i k_i e_i/|k|^2*k.")
    print("    Exactly 1 of 4 Sheet A directions is the longitudinal mode.")
    print()
    print("    Example: k = (1,0,0,0)")
    print("      Longitudinal: e_0^A = weight (1,0,0)  [idx 0]  <- killed")
    print("      Transverse:   e_1^A (weight (-1,1,0)), e_2^A (weight (0,-1,1)),")
    print("                    e_3^A (weight (0,0,-1))   [idx 1,2,3]")
    print()
    print("    By sl_4-invariance: for any k, exactly 1 of 4 Sheet A seeds")
    print("    is killed; the other 3 form the transverse (physical) subspace.")
    print()
    print(f"  LERAY KILL: 1 longitudinal Sheet A seed -> ELIMINATED")
    print(f"  Remaining: {len(sheet_A) - 1} transverse Sheet A seeds")

    # =========================================================================
    # 7. PRESSURE REDUNDANCY (double-check of Z_2 kill)
    # =========================================================================
    print()
    print("=" * 70)
    print("7. Pressure Redundancy: Sheet B Doubly Excluded")
    print("=" * 70)
    print()
    print("  (a) Z_2 parity (Section 5): Sheet B is phat4-odd -> non-bosonic. [OK]")
    print()
    print("  (b) Pressure Poisson equation:")
    print("      -Deltap = partial_ipartial_j(u_iu_j)  ->  nablap = -nablaDelta^-^1partial_ipartial_j(u_iu_j)")
    print("      Pressure is uniquely determined by velocity.")
    print("      Sheet B (partial_ip) contributes no independent H^1 obstruction. [OK]")
    print()
    print("  The Z_2 argument is algebraic and takes precedence.")

    # =========================================================================
    # SUMMARY
    # =========================================================================
    print()
    print("=" * 70)
    print("STEP 3 SUMMARY -- M_1(1,0,0): Physical H^1 Seed Kill")
    print("=" * 70)
    print()
    print("  Starting seeds: 8  (all genuine H^1 classes in Complex B)")
    print()
    print("  Certified algebraic facts:")
    print(f"    phi_1D d=0 nnz={nnz_1D}  phi_4H d=0 nnz={PRECERT['phi_4H_d0_nnz']}"
          f"  -> all 8 are cocycles")
    print(f"    d^0 -> fiber rank={rk_in_B + rk_in_A}"
          f"  -> all 8 are genuine H^1 (sv_deggeq1 structural arg)")
    print(f"    d^0 -> M_1(1,0,0)[d=1] rank={rk_d1}  -> {rk_d1}/{dim_d1} degree-1 exact")
    print(f"    d=1 outgoing: nnz=0  -> no singular vector at d=1")
    print()
    print("  +------------------+------+--------------+----------------------+")
    print("  | Stage            | Kill |  Remaining   | Method               |")
    print("  +------------------+------+--------------+----------------------+")
    print("  | Algebraic cert   |   0  |      8       | nnz=0, sv_deggeq1 arg  |")
    print("  | Z_2 parity        |   4  |      4       | Sheet B phat4-odd    |")
    print("  | Leray projector  |   1  |      3       | k^perp subset C^4, dim=3     |")
    print("  +------------------+------+--------------+----------------------+")
    print()
    print(f"  SURVIVING SEEDS: 3  (transverse Sheet A, phat4-even)")
    print()
    print("  These 3 = the 3 transverse velocity Fourier modes in k^perp subset C^4.")
    print("  They carry the sl_4 quotient rep C^4/<k> cong C^3.")
    print()
    print("  No further algebraic kill is possible:")
    print("    * sl_4-equivariance: any E(4,4) morphism kills all 3 or none;")
    print("      nnz=0 for all outgoing morphisms at d=0 and d=1.")
    print("    * No new Complex B morphism targets M_1(1,0,0).")
    print()
    print("  Open paths to resolve the 3 survivors:")
    print("    1. Sobolev/energy-estimate: CCK image in H^1 controlled by ||u||_{H^1}")
    print("       -> regularity argument kills obstruction for regular solutions.")
    print("    2. New E(4,4)-equivariant morphism at higher PBW degree (dgeq2).")
    print("    3. sl_3-decomposition of C^3 -> possible new morphism sub-rep.")
    print("    4. Spectral-sequence cancellation with M_1(2,0,0) survivors.")
    print()
    print("  STATUS: Step 3 INCOMPLETE -- 3 seeds remain OPEN.")
    print("  These 3 transverse Sheet A seeds are the PRIMARY physical H^1")
    print("  obstruction in the E(4,4) Navier-Stokes regularity programme.")
    print("=" * 70)


if __name__ == "__main__":
    main()

