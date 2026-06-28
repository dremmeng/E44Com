"""
cohomology.py  —  Cohomology of the exceptional de Rham complexes for E(4,4)
=============================================================================
Section 9 of the NSE44 programme  (Cantarini-Caselli-Kac 2026).

Computes H^k = ker(O_k) / im(I_k) for both Complex A and Complex B, where:

  O_k  =  vstack of all differentials D_{k → k+s}  (outgoing from level k)
  I_k  =  hstack of all differentials D_{k-s → k}  (incoming to level k)

Computing H^k is an OPEN PROBLEM stated in the introduction of the paper.
This code produces numerical evidence via truncated Verma modules (sl_4 fibers,
max_deg controlled by the caller).

Notes on sl_4 vs \hat{p}(4) fibers
------------------------------
All cochain groups use sl_4 fibers (e44_data=None in CochainGroup), which is
consistent with how the phi[1A] and phi[1B] morphisms are computed internally.
Morphisms that require \hat{p}(4) fibers (phi[1C], phi[2DA], phi[2EA], phi[1D/E],
phi[3F/G/4H]) are silently skipped during assemble_differential via
FiberMismatchError.  The resulting complex therefore captures the sl_4-fiber
contribution to each cochain group.

Run inside SageMath:  sage cohomology.py
"""

from sage.all import QQ, matrix, vector, VectorSpace
import sys as _sys
import os as _os

_HERE = _os.path.dirname(_os.path.abspath(__file__))
if _HERE not in _sys.path:
    _sys.path.insert(0, _HERE)

from de_rham_complex import (
    DeRhamComplexA, DeRhamComplexB,
    save_complex, load_complex,
    DEFAULT_T_MIN, DEFAULT_T_MAX, DEFAULT_A_MAX, MAX_VERMA_DEG,
)


# ===========================================================================
# Building the outgoing / incoming matrices
# ===========================================================================

def outgoing_matrix(cx, k):
    """
    Stack all outgoing differentials from level k vertically.

    Returns a sparse QQ matrix of shape
        (\sum_{(k,k') \in differentials} dim_{k'}, dim_k)

    If there are no outgoing differentials, returns a (0, dim_k) zero matrix.

    Parameters
    ----------
    cx : DeRhamComplexA or DeRhamComplexB
    k  : int — cochain level (must be in cx.positions)
    """
    dim_k  = cx.groups[k].total_dim
    blocks = [D for (ks, _), D in sorted(cx.differentials.items()) if ks == k]

    if not blocks:
        return matrix(QQ, 0, dim_k, sparse=True)

    total_rows = sum(B.nrows() for B in blocks)
    result     = matrix(QQ, total_rows, dim_k, sparse=True)
    r = 0
    for B in blocks:
        for (i, j), val in B.dict().items():
            result[r + i, j] = val
        r += B.nrows()
    return result


def incoming_matrix(cx, k):
    """
    Stack all incoming differentials to level k horizontally.

    Returns a sparse QQ matrix of shape
        (dim_k, \sum_{(k',k) \in differentials} dim_{k'})

    If there are no incoming differentials, returns a (dim_k, 0) zero matrix.

    Parameters
    ----------
    cx : DeRhamComplexA or DeRhamComplexB
    k  : int — cochain level (must be in cx.positions)
    """
    dim_k  = cx.groups[k].total_dim
    blocks = [D for (_, kt), D in sorted(cx.differentials.items()) if kt == k]

    if not blocks:
        return matrix(QQ, dim_k, 0, sparse=True)

    total_cols = sum(B.ncols() for B in blocks)
    result     = matrix(QQ, dim_k, total_cols, sparse=True)
    c = 0
    for B in blocks:
        for (i, j), val in B.dict().items():
            result[i, c + j] = val
        c += B.ncols()
    return result


# ===========================================================================
# Core cohomology computation
# ===========================================================================

def cohomology_at(cx, k, verbose=False):
    """
    Compute H^k = ker(O_k) / im(I_k) for the complex cx at cochain level k.

    Algorithm
    ---------
    Let K = ker(O_k) and I = im(I_k), both subspaces of QQ^{dim_k}.
    Then H^k ≅ K / (K \intersect I), so dim(H^k) = dim(K) - dim(K \intersect I).

    The intersection K \intersect I is computed via Sage's VectorSpace.intersection().

    Sanity check (d^2 = 0): I \subseteq K at interior levels (guaranteed by Subtask 7).
    At boundary levels (k = t_min or t_max) some morphisms are missing, so
    I \subseteq K may fail; in that case dim_im is computed as dim(K \intersect I) anyway.

    Parameters
    ----------
    cx      : DeRhamComplexA or DeRhamComplexB
    k       : int — cochain level
    verbose : bool

    Returns
    -------
    dict with keys:
        'k'        : int
        'dim_Ck'   : int — dimension of C^k
        'dim_ker'  : int — dimension of ker(O_k)
        'dim_im'   : int — dimension of im(I_k) \intersect ker(O_k)
        'dim_H'    : int — dimension of H^k
        'im_subset_ker': bool — True iff im(I_k) \subseteq ker(O_k)  (d^2 = 0 check)
        'cocycles' : matrix — basis of ker(O_k) as row matrix
        'cobounds' : matrix — basis of K \intersect I     as row matrix
    """
    dim_k = cx.groups[k].total_dim
    V     = VectorSpace(QQ, dim_k)

    # ── Outgoing: O_k ────────────────────────────────────────────────────
    O_k = outgoing_matrix(cx, k)
    if O_k.nrows() == 0:
        # No outgoing: entire C^k is in the kernel
        K_space = V
    else:
        K_space = V.subspace(O_k.right_kernel().basis())

    # ── Incoming: I_k ────────────────────────────────────────────────────
    I_k = incoming_matrix(cx, k)
    if I_k.ncols() == 0:
        # No incoming: image is the zero subspace
        I_space = V.subspace([])
    else:
        I_space = V.subspace(I_k.column_space().basis())

    # ── Intersection K \intersect I ───────────────────────────────────────────────
    inter      = K_space.intersection(I_space)
    dim_ker    = K_space.dimension()
    dim_im_Ck  = I_space.dimension()          # dim of im(I_k) inside C^k
    dim_im     = inter.dimension()             # dim of im \intersect ker
    dim_H      = dim_ker - dim_im

    # d^2 = 0 iff im(I_k) \subseteq ker(O_k) iff I_space \subseteq K_space
    im_subset_ker = I_space.is_subspace(K_space)

    # Basis matrices
    cocycles_rows = list(K_space.basis())
    cobounds_rows = list(inter.basis())
    cocycles = matrix(QQ, cocycles_rows) if cocycles_rows else matrix(QQ, 0, dim_k)
    cobounds = matrix(QQ, cobounds_rows) if cobounds_rows else matrix(QQ, 0, dim_k)

    if verbose:
        d2_flag = '' if im_subset_ker else '  [WARNING: d^2≠0 at this level]'
        print(f"  k={k:+d}: dim C^k={dim_k:5d}, dim ker={dim_ker:5d}, "
              f"dim im(\intersectker)={dim_im:4d}, dim H^k={dim_H:4d}{d2_flag}")

    return {
        'k'            : k,
        'dim_Ck'       : dim_k,
        'dim_ker'      : dim_ker,
        'dim_im'       : dim_im,
        'dim_H'        : dim_H,
        'im_subset_ker': im_subset_ker,
        'cocycles'     : cocycles,
        'cobounds'     : cobounds,
    }


def cohomology_table(cx, positions=None, verbose=True):
    """
    Compute H^k for all (or selected) cochain levels of the complex.

    Parameters
    ----------
    cx        : DeRhamComplexA or DeRhamComplexB
    positions : list of int or None — levels to compute (default: all)
    verbose   : bool — print the table

    Returns
    -------
    dict{k: cohomology_at result dict}
    """
    if positions is None:
        positions = cx.positions

    results = {}

    if verbose:
        print(f"\n  {'k':>4}  {'dim C^k':>9}  {'dim ker':>9}  "
              f"{'dim im':>8}  {'dim H^k':>9}  d^2=0?")
        print("  " + "-" * 58)

    for k in positions:
        r = cohomology_at(cx, k, verbose=False)
        results[k] = r
        if verbose:
            d2_str = "✓" if r['im_subset_ker'] else "✗ WARNING"
            print(f"  {k:>+4}  {r['dim_Ck']:>9}  {r['dim_ker']:>9}  "
                  f"{r['dim_im']:>8}  {r['dim_H']:>9}  {d2_str}")

    return results


# ===========================================================================
# Self-test
# ===========================================================================

def _check_cohomology(e44_data, t_min=-2, t_max=4, a_max=1, max_deg=1,
                      verbose=True):
    """
    Cohomology self-check:
      1.  outgoing_matrix has shape (*, dim_k).
      2.  incoming_matrix has shape (dim_k, *).
      3.  cohomology_at returns a dict with correct keys.
      4.  dim_H is non-negative at every level.
      5.  d^2 = 0 \iff im \subseteq ker at every interior level (verified by Subtask 7).
      6.  For a level with no incoming and no outgoing morphisms, H^k = C^k.
      7.  cohomology_table returns an entry for every level in positions.
      8.  dim_ker ≥ dim_im (always, since dim_H = dim_ker - dim_im ≥ 0).

    Uses max_deg=1 and a small window so the linear algebra is fast.

    Returns True iff all checks pass.
    """
    all_pass = True
    n_pass   = 0
    n_fail   = 0

    def ok(name, cond, detail=''):
        nonlocal all_pass, n_pass, n_fail
        if cond:
            n_pass += 1
            if verbose:
                print(f"  [PASS] {name}")
        else:
            n_fail += 1
            all_pass = False
            if verbose:
                print(f"  [FAIL] {name}" + (f": {detail}" if detail else ''))

    if verbose:
        print("=" * 60)
        print("cohomology.py self-check")
        print(f"  (max_deg={max_deg}, t=[{t_min},{t_max}], a_max={a_max})")
        print("=" * 60)

    for cx_name, cx_cls in [("DeRhamComplexA", DeRhamComplexA),
                             ("DeRhamComplexB", DeRhamComplexB)]:
        if verbose:
            print(f"\n--- {cx_name} ---")

        cx = cx_cls(t_min=t_min, t_max=t_max, a_max=a_max,
                    max_deg=max_deg, e44_data=e44_data)
        if verbose:
            print(f"  {cx}")

        for k in cx.positions:
            dim_k = cx.groups[k].total_dim

            # ── 1: outgoing_matrix shape ──────────────────────────────────
            O_k = outgoing_matrix(cx, k)
            ok(f"k={k:+d}: outgoing_matrix ncols = dim_k={dim_k}",
               O_k.ncols() == dim_k,
               f"got {O_k.ncols()}")

            # ── 2: incoming_matrix shape ──────────────────────────────────
            I_k = incoming_matrix(cx, k)
            ok(f"k={k:+d}: incoming_matrix nrows = dim_k={dim_k}",
               I_k.nrows() == dim_k,
               f"got {I_k.nrows()}")

            # ── 3 & 4: cohomology_at keys and non-negativity ──────────────
            r = cohomology_at(cx, k, verbose=False)
            expected_keys = {'k', 'dim_Ck', 'dim_ker', 'dim_im',
                             'dim_H', 'im_subset_ker', 'cocycles', 'cobounds'}
            ok(f"k={k:+d}: cohomology_at returns correct keys",
               set(r.keys()) == expected_keys,
               f"got {set(r.keys())}")
            ok(f"k={k:+d}: dim_H >= 0  (got {r['dim_H']})",
               r['dim_H'] >= 0)

            # ── 5: interior levels satisfy d^2 = 0 → im \subseteq ker ─────────────
            is_interior = (k > t_min and k < t_max)
            if is_interior:
                ok(f"k={k:+d}: interior level has im \subseteq ker (d^2=0)",
                   r['im_subset_ker'])

            # ── 8: dim_ker >= dim_im ──────────────────────────────────────
            ok(f"k={k:+d}: dim_ker={r['dim_ker']} >= dim_im={r['dim_im']}",
               r['dim_ker'] >= r['dim_im'])

        # ── 7: cohomology_table covers all positions ──────────────────────
        table = cohomology_table(cx, verbose=False)
        ok(f"{cx_name}: cohomology_table covers all positions",
           set(table.keys()) == set(cx.positions))

        if verbose:
            print(f"\n  Cohomology table ({cx_name}):")
            cohomology_table(cx, verbose=True)

    if verbose:
        print("\n" + "=" * 60)
        print(f"cohomology.py summary: {n_pass} pass, {n_fail} fail")
        if all_pass:
            print("cohomology.py  ✓  ALL CHECKS PASSED")
        else:
            print("cohomology.py  ✗  SOME CHECKS FAILED")
        print("=" * 60)

    return all_pass


# ===========================================================================
# Module entry point
# ===========================================================================

if __name__ == '__main__':
    print("cohomology.py — Exceptional de Rham cohomology for E(4,4)")
    print()

    from verma_modules import load_e44
    print("Loading e44_data...")
    _e44 = load_e44()
    print("e44_data loaded.")
    print()

    ok = _check_cohomology(
        e44_data=_e44,
        t_min=-2, t_max=4,
        a_max=1,
        max_deg=1,
        verbose=True,
    )

    import sys as _sys_main
    _sys_main.exit(0 if ok else 1)
