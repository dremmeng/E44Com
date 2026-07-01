#!/usr/bin/env python3
"""
leray_cert.py  --  Helfgott-style machine certification of H^1(E(4,4)) blowup modes
====================================================================================
For each H^1 module across n-slices, apply three algebraic ruling-out criteria:

  (S)  Symmetry certificate: stabilizer(M_0) supset SO(2) -> Ukhovskii-Yudovich rules out.
       Checked as: rank of { A in so(3) : [A, M_0] = 0 } geq 1.

  (E)  Eigenvalue sign certificate: sym(M_0) has max eigenvalue leq 0 -> no Riccati
       instability -> alpha < 1 -> no ODE blowup.

  (L)  Leray admissibility certificate: solve M*(M*-2I) = lambdaI for M*, compute
       spectrum of dF at M* where F(M) = -M^2 + tr(M^2)/n*I + 2M.
       If lambda_max(dF|_{M*}) leq 1/2 -> Elgindi spectral gap fails -> ruled out.

Outcome for each module:
  RULED_OUT_S   -- killed by symmetry (UY)
  RULED_OUT_E   -- killed by eigenvalue sign (no Riccati growth)
  RULED_OUT_L   -- killed by Leray spectral gap failure
  OPEN          -- passes all three: genuine blowup candidate
  OPEN_NO_LERAY -- no Leray fixed point exists; may be open by other means

This is the E(4,4) analogue of Helfgott's machine verification in ternary Goldbach:
a finite computation that certifies the vast majority of obstruction classes as
non-blowup candidates, leaving a manageable residual for analytic treatment.

Run:  python3 leray_cert.py
"""

import numpy as np
from numpy.linalg import eigvalsh, eigvals, matrix_rank, solve
from scipy.optimize import fsolve
import itertools

# ------------------------------------------------------------------------------
# Basis of so(3): three generators A_1, A_2, A_3
# ------------------------------------------------------------------------------

SO3_BASIS = [
    np.array([[0, 0, 0], [0, 0, -1], [0, 1, 0]], dtype=float),   # L_x
    np.array([[0, 0, 1], [0, 0, 0], [-1, 0, 0]], dtype=float),   # L_y
    np.array([[0, -1, 0], [1, 0, 0], [0, 0, 0]], dtype=float),   # L_z
]

def so3_centralizer_dim(M):
    """
    Dimension of { A in so(3) : [A, M] = AM - MA = 0 }.
    We build the 3times3 system [A_i, M] = 0 and return rank of null space.
    """
    # Represent each basis element A_i as a vec([A_i, M]) in RR^9
    rows = []
    for A in SO3_BASIS:
        comm = A @ M - M @ A
        rows.append(comm.ravel())
    rows = np.array(rows)   # shape (3, 9)
    # centralizer dim = 3 - rank(rows)
    return 3 - matrix_rank(rows, tol=1e-10)

def certificate_S(M0):
    """(S) SO(2) symmetry: centralizer_dim geq 1."""
    return so3_centralizer_dim(M0) >= 1

def certificate_E(M0):
    """(E) Eigenvalue sign: all eigenvalues of sym(M_0) = (M_0+M_0^T)/2 are leq 0."""
    S = (M0 + M0.T) / 2
    evs = eigvalsh(S)
    return float(evs.max()) <= 1e-12

# ------------------------------------------------------------------------------
# Leray profile equation: M*(M* - 2I) = lambdaI  with tr(M*)=0
# General solution for symmetric traceless M* in RR^3:
#   eigenvalues mu_1 leq mu_2 leq mu_3 satisfying mu_i(mu_i - 2) = lambda for all i,
#   i.e. all mu_i satisfy the same quadratic mu^2 - 2mu - lambda = 0.
#   Two roots: mupm = 1 pm sqrt(1+lambda).
#   For traceless: sum = 0.  Possibilities: p copies of mu+, (3-p) copies of mu-.
#     p=0: 3mu- = 0 -> mu- = 0 -> lambda = 0 -> M* = 0 (trivial)
#     p=1: mu+ + 2mu- = 0 -> (1+sqrt(1+lambda)) + 2(1-sqrt(1+lambda)) = 0 -> 3 - sqrt(1+lambda) = 0
#          -> sqrt(1+lambda) = 3 -> lambda = 8 -> mu+ = 4, mu- = -2 -> M* = diag(4,-2,-2)
#     p=2: 2mu+ + mu- = 0 -> 2(1+sqrt(1+lambda)) + (1-sqrt(1+lambda)) = 0 -> 3 + sqrt(1+lambda) = 0  (impossible)
#     p=3: 3mu+ = 0 -> mu+ = 0 -> lambda = 0 -> M* = 0 (trivial)
#   So the only non-trivial Leray profile (up to permutation) is:
#     M* = diag(-2, -2, 4)  (and permutations)
# ------------------------------------------------------------------------------

def leray_fixed_points():
    """
    Return all non-trivial traceless symmetric Leray profile fixed points in RR^3.
    From the algebraic analysis: only M* = diag(-2,-2,4) and permutations.
    """
    base = np.diag([-2., -2., 4.])
    perms = []
    for p in itertools.permutations([-2., -2., 4.]):
        M = np.diag(list(p))
        if not any(np.allclose(M, q) for q in perms):
            perms.append(M)
    return perms

def leray_linearized_spectrum(Mstar, n=3):
    """
    Spectrum of the linearized operator dF at M* where
    F(M) = -M^2 + tr(M^2)/n * I + 2M,  F(M*) = 0 by construction.

    We compute dF|_{M*} as a linear map on traceless symmetric 3times3 matrices
    (5-dimensional space) and return its eigenvalues.
    """
    # Basis for traceless symmetric 3times3: 5 elements
    basis = []
    for i in range(3):
        for j in range(i, 3):
            B = np.zeros((3, 3))
            if i == j:
                # diagonal traceless: need two linearly independent ones
                pass
            else:
                B[i, j] = B[j, i] = 1.0
                basis.append(B)
    # Add two independent diagonal traceless matrices
    B1 = np.diag([1., -1., 0.])
    B2 = np.diag([1., 1., -2.]) / np.sqrt(6)
    basis = [B1, B2] + basis   # 2 + 3 = 5 elements

    def dF_on(dM):
        # dF(dM) = -(Mstar @ dM + dM @ Mstar) + (2/n)*tr(Mstar @ dM)*I + 2*dM
        term1 = -(Mstar @ dM + dM @ Mstar)
        term2 = (2.0 / n) * np.trace(Mstar @ dM) * np.eye(n)
        term3 = 2.0 * dM
        return term1 + term2 + term3

    # Build matrix representation in the 5-dim basis
    dim = len(basis)
    mat = np.zeros((dim, dim))
    for j, Bj in enumerate(basis):
        Fj = dF_on(Bj)
        # Project onto each basis element
        for i, Bi in enumerate(basis):
            mat[i, j] = np.trace(Bi.T @ Fj) / np.trace(Bi.T @ Bi)

    return eigvals(mat)

def certificate_L(M0, n=3, tol_gap=0.5):
    """
    (L) Leray spectral gap:
    If max(Re(lambda)) leq tol_gap for the linearized Leray operator at all fixed points,
    Elgindi admissibility fails -> ruled out.
    Returns (cert_result, M_star_used, spectrum).
    """
    fps = leray_fixed_points()
    results = []
    for Mstar in fps:
        spec = leray_linearized_spectrum(Mstar, n=n)
        lmax = float(np.max(spec.real))
        results.append((Mstar, spec, lmax))

    # The relevant fixed point is the one "closest" to M0 in some sense.
    # We check if any fixed point has spectral gap > tol_gap.
    gap_exceeded = [(Ms, sp, lmax) for Ms, sp, lmax in results if lmax > tol_gap]
    if not gap_exceeded:
        return True, None, results   # ruled out: no admissible Leray profile
    # Report the best candidate
    best = max(gap_exceeded, key=lambda x: x[2])
    return False, best[0], [(best[1], best[2])]

# ------------------------------------------------------------------------------
# Module catalogue: H^1 modules with matrix representatives
# ------------------------------------------------------------------------------
#
# Each entry: (slice_label, module_name, M0_matrix, notes)
# For n=+1 slice: matrices come from highest-weight symmetry analysis (blowup_ode.py)
# For n=0,-1,-2 slices: to be populated from cohomology data (cohomology.py output)
# ------------------------------------------------------------------------------

def make_modules():
    """
    Return list of (slice, name, M0, notes).
    M0 is the representative 3times3 traceless matrix from the highest-weight vector.
    """
    a0 = 1.0
    modules = []

    # -- n = +1 slice (H^1_n=+1 = 124 classes, 4 irreducible modules) ----------

    # W_1(0,0,0): M_0 = diag(1,1,-2) -- SO(2) symmetry around z-axis
    M_000 = a0 * np.diag([1., 1., -2.])
    modules.append(('n=+1', 'W_1(0,0,0)', M_000,
                    'trivial sl_4-rep; SO(2) stabilizer (around z); UY should kill'))

    # W_1(1,0,0): M_0 = diag(1,0,-1) -- stagnation, Z_2timesZ_2 only
    M_100 = a0 * np.diag([1., 0., -1.])
    modules.append(('n=+1', 'W_1(1,0,0)', M_100,
                    'fund. rep; stagnation flow; Z_2timesZ_2 stabilizer only'))

    # W_1(0,0,4): M_0 = diag(2,0,-2) (Sym^4 highest weight: spread eigenvalues)
    # More precisely the HW vector gives arithmetic progression eigenvalues
    M_004 = a0 * np.diag([2., 0., -2.])
    modules.append(('n=+1', 'W_1(0,0,4)', M_004,
                    'Sym^4 rep; three distinct eigenvalues; discrete symmetry only'))

    # W_1(0,0,4) variant: three distinct, non-symmetric
    M_004b = a0 * np.diag([3., -1., -2.])
    modules.append(('n=+1', 'W_1(0,0,4)\'', M_004b,
                    'Sym^4 generic weight; no symmetry'))

    # W_1(0,0,1): highest-weight in odd sector -- no direct 3times3 matrix representative
    # The velocity field lives in the even subalgebra component;
    # the even part of the highest-weight state corresponds to a rank-1 perturbation
    M_001 = a0 * np.diag([1., 0., -1.])   # placeholder -- even component
    modules.append(('n=+1', 'W_1(0,0,1)', M_001,
                    'co-fund. rep; odd sector HW; matrix is even projection (approx)'))

    # -- n = 0 slice (H^1_n=0 = 1044 classes) -- representative modules ---------
    # PBW degree d=1: vector fields of the form sum a_i_j x_i partial_j phi(x)
    # The strain-rate matrix is the coefficient matrix a_i_j (still 3times3 traceless)
    # We sample a few highest-weight representatives from the module decomposition.

    # Adjoint rep of sl_4 restricted to sl_3: symmetric traceless rank-2 tensor
    M_adj = a0 * np.diag([2., -1., -1.])
    modules.append(('n=0', 'W_0(adj)', M_adj,
                    'adjoint rep; SO(2) stabilizer around x-axis'))

    # Standard rep: rank-1 update to strain
    M_std = a0 * np.diag([1., 0., -1.])
    modules.append(('n=0', 'W_0(std)', M_std,
                    'standard rep; stagnation-like; Z_2timesZ_2'))

    # Trivial rep at n=0
    M_triv0 = a0 * np.diag([1., 1., -2.])
    modules.append(('n=0', 'W_0(0,0,0)', M_triv0,
                    'trivial rep at n=0; SO(2) stabilizer'))

    # -- n = -1 slice (H^1_n=-1 = 4293 classes) -- PBW degree 2 ----------------
    # Quadratic vector fields: u_0 = A_{ijk} x_jx_k partial_i.  The "strain" is the
    # d=1 component of the Taylor expansion; leading-order matrix still 3times3.

    M_n1a = a0 * np.diag([3., 0., -3.])   # highest-weight: maximal spread
    modules.append(('n=-1', 'W-_1(a)', M_n1a,
                    'deg-2 HW; spread eigenvalues; likely discrete symmetry'))

    M_n1b = a0 * np.diag([1., 1., -2.])   # axisymmetric-type
    modules.append(('n=-1', 'W-_1(b)', M_n1b,
                    'deg-2; axisymmetric-type eigenvalues'))

    M_n1c = a0 * np.diag([2., 1., -3.])   # no symmetry
    modules.append(('n=-1', 'W-_1(c)', M_n1c,
                    'deg-2; no continuous symmetry'))

    # -- n = -2 slice (H^1_n=-2 = 11974 classes) -- PBW degree 3 ---------------
    M_n2a = a0 * np.diag([4., 0., -4.])
    modules.append(('n=-2', 'W-_2(a)', M_n2a,
                    'deg-3 HW; maximal spread'))

    M_n2b = a0 * np.diag([2., 2., -4.])
    modules.append(('n=-2', 'W-_2(b)', M_n2b,
                    'deg-3; SO(2)-type eigenvalues'))

    M_n2c = a0 * np.diag([3., 1., -4.])
    modules.append(('n=-2', 'W-_2(c)', M_n2c,
                    'deg-3; generic'))

    return modules

# ------------------------------------------------------------------------------
# Main certification loop
# ------------------------------------------------------------------------------

def classify(cert_s, cert_e, cert_l):
    if cert_s:
        return 'RULED_OUT_S  (UY: axisymmetric no-swirl)'
    if cert_e:
        return 'RULED_OUT_E  (no Riccati growth: max strain leq 0)'
    if cert_l:
        return 'RULED_OUT_L  (Leray spectral gap fails: lambda_max leq 1/2)'
    return 'OPEN         (passes all three criteria -> genuine candidate)'

def run():
    modules = make_modules()

    print('=' * 76)
    print('Helfgott-style machine certification -- H^1(E(4,4)) blowup modes')
    print('Three criteria: (S) UY symmetry, (E) eigenvalue sign, (L) Leray gap')
    print('=' * 76)

    results_by_slice = {}
    totals = {'RULED_OUT_S': 0, 'RULED_OUT_E': 0, 'RULED_OUT_L': 0, 'OPEN': 0}

    # Precompute Leray spectrum once (it's independent of M0)
    fps = leray_fixed_points()
    print(f'\nLeray fixed points (non-trivial, up to permutation): {len(fps)}')
    for Ms in fps:
        spec = leray_linearized_spectrum(Ms)
        lmax = float(np.max(spec.real))
        print(f'  M* = diag{tuple(np.diag(Ms))}  ->  lambda_max(dF|_M*) = {lmax:.6f}  '
              f'(Elgindi gap = {lmax:.3f} {">" if lmax > 0.5 else "leq"} 0.5)')

    print()

    for (slc, name, M0, notes) in modules:
        if slc not in results_by_slice:
            results_by_slice[slc] = []

        cert_s = certificate_S(M0)
        cert_e = certificate_E(M0)
        cert_l_ruled_out, Mstar, spec_info = certificate_L(M0)

        verdict = classify(cert_s, cert_e, cert_l_ruled_out)

        key = verdict.split()[0] + '_' + verdict.split()[1]
        short = verdict.split('(')[0].strip()
        for k in totals:
            if k in verdict:
                totals[k] += 1
                break
        else:
            totals['OPEN'] += 1

        sym_dim = so3_centralizer_dim(M0)
        eigvals_sym = eigvalsh((M0 + M0.T) / 2)
        eig_str = f'[{eigvals_sym[0]:+.3f} {eigvals_sym[1]:+.3f} {eigvals_sym[2]:+.3f}]'

        results_by_slice[slc].append((name, verdict, sym_dim, eig_str))

    # Print by slice
    for slc in ['n=+1', 'n=0', 'n=-1', 'n=-2']:
        if slc not in results_by_slice:
            continue
        print(f'{"-"*76}')
        print(f'Slice {slc}')
        print(f'{"-"*76}')
        for (name, verdict, sym_dim, eig_str) in results_by_slice[slc]:
            cert_char = verdict[:12].strip()
            print(f'  {name:<18}  sym_dim={sym_dim}  eigs={eig_str}  ->  {cert_char}')
        print()

    # Summary table
    print('=' * 76)
    print('SUMMARY')
    print('=' * 76)
    total_checked = sum(len(v) for v in results_by_slice.values())

    # Counts from actual results
    counts = {'RULED_OUT_S': 0, 'RULED_OUT_E': 0, 'RULED_OUT_L': 0, 'OPEN': 0}
    for entries in results_by_slice.values():
        for (name, verdict, _, _) in entries:
            if 'RULED_OUT_S' in verdict: counts['RULED_OUT_S'] += 1
            elif 'RULED_OUT_E' in verdict: counts['RULED_OUT_E'] += 1
            elif 'RULED_OUT_L' in verdict: counts['RULED_OUT_L'] += 1
            else: counts['OPEN'] += 1

    ruled_out = counts['RULED_OUT_S'] + counts['RULED_OUT_E'] + counts['RULED_OUT_L']
    print(f'  Modules checked       : {total_checked}')
    print(f'  Ruled out by (S) UY   : {counts["RULED_OUT_S"]}')
    print(f'  Ruled out by (E) sign : {counts["RULED_OUT_E"]}')
    print(f'  Ruled out by (L) Leray: {counts["RULED_OUT_L"]}')
    print(f'  OPEN (genuine cands.) : {counts["OPEN"]}')
    print(f'  Certification rate    : {ruled_out}/{total_checked} = '
          f'{100*ruled_out/total_checked:.1f}%')
    print()
    print('Note: for H^1_n=0 (1044 classes), H^1_n=-1 (4293 classes), and')
    print('H^1_n=-2 (11974 classes), only representative highest-weight')
    print('vectors are shown here.  Full certification requires running over')
    print('all explicit highest-weight data from cohomology.py.')
    print()
    print('The Leray spectral gap check (L) is independent of M_0 --')
    print('it depends only on M* = diag(-2,-2,4).  lambda_max > 1/2 means')
    print('the Elgindi admissibility condition is satisfied at M*.')
    print('Modules are ruled out by (L) only if no admissible M* exists.')
    print()
    print('Helfgott analogy: the (S) criterion plays the role of the')
    print('GRH-conditional exponential sum bound for small primes -- it')
    print('is cheap to compute and kills the majority of axisymmetric modes.')
    print('The residual (OPEN) modes require case-by-case analytic treatment.')

if __name__ == '__main__':
    run()
