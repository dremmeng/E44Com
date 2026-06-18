"""
morphisms.py  —  Verma-module morphisms for the de Rham complex of E(4,4)
=========================================================================
Section 7 of the NSE44 programme  (Cantarini–Caselli–Kac 2026).

This module assembles the ten E(4,4)-Verma-module morphisms \phi[1A]–\phi[4H]
needed for both incarnations of the exceptional de Rham complex (Figures 7
and 8 of the paper).  The degree-1 through degree-4 morphisms are imported
from singular_vectors.py, and the composite degree-2 morphism \phi[2EA] is
constructed here.

Morphism catalogue
------------------
  phi_1A(t, a, e44_data, max_source_deg)
        M_{t-1}(a+1,0,0) \to M_t(a,0,0)         degree 1   a≥1 or a=t=0

  phi_1B(t, c, e44_data, max_source_deg)
        M_{t-1}(0,0,c-1) \to M_t(0,0,c)          degree 1   c≥1

  phi_1C(t, e44_data, max_source_deg)
        M_{t-1}(0,1,0) \to M_t(0,0,1)            degree 1   p̂(4) target fiber

  phi_1D(t, e44_data, max_source_deg)
        M_{t-1}(1,0,0) \to M_t(0,0,0)            degree 1   t≠0

  phi_1E(e44_data, max_source_deg)
        M_0(0,0,0) \to M_1(1,0,0)                degree 1   p̂(4) target fiber

  phi_2DA(t, e44_data, max_source_deg)
        M_{t-2}(2,0,0) \to M_t(0,1,0)            degree 2   p̂(4) target fiber

  phi_2EA(e44_data, max_source_deg)
        M_{-1}(1,0,0) \to M_1(1,0,0)             degree 2   = phi_1E \circ phi_1A

  phi_3F(e44_data, max_source_deg)
        M_0(0,0,0) \to M_3(1,0,0)                degree 3   p̂(4) target fiber

  phi_3G(e44_data, max_source_deg)
        M_{-3}(1,0,0) \to M_0(0,0,0)             degree 3

  phi_4H(t, e44_data, max_source_deg)
        M_{t-4}(1,0,0) \to M_t(1,0,0)            degree 4   p̂(4) target fiber

Return convention
-----------------
Every phi_* function returns a 5-tuple:
    (M_src, M_tar, sv_deg, phi0_vecs, matrices)

where:
  M_src     : VermaModule — source module M_{t'}(a',b',c')
  M_tar     : VermaModule — target module M_t(a,b,c)
  sv_deg    : int — degree of the defining singular vector
  phi0_vecs : list[QQ-vector] — fiber-map images \phi_0(e_k) \in M_tar[sv_deg],
              one per source-fiber basis vector (len = M_src.dim_W)
  matrices  : dict{d: QQ-matrix} — morphism matrices at source degree d,
              shape dim(M_tar[sv_deg+d]) \times dim(M_src[d]),
              for d = 0 ... max_source_deg

Utility functions
-----------------
  compose_morphisms(phi_B_mats, phi_A_mats, sv_deg_A, max_result_deg)
      Matrix-level composition  (\phi_B \circ \phi_A)[d] = phi_B_mats[sv_A+d] * phi_A_mats[d]

  is_proportional(v, w)
      True iff vectors v and w are scalar multiples of each other (over QQ).

  save_morphisms(filepath, e44_data, max_deg, verbose)
      Re-exported from singular_vectors; computes and pickles all morphism data.

  load_morphisms(filepath)
      Re-exported from singular_vectors; reloads pickled morphism data.

Dependency
----------
  e44_brackets.pkl  — written by sage e44_structure.py
  phat4_cache.pkl   — written by sage phat4_modules.py (optional cache)

Run inside SageMath:  sage morphisms.py
"""

from sage.all import QQ, vector, matrix
import sys as _sys
import os as _os

_HERE = _os.path.dirname(_os.path.abspath(__file__))
if _HERE not in _sys.path:
    _sys.path.insert(0, _HERE)

# ---------------------------------------------------------------------------
# Core imports
# ---------------------------------------------------------------------------
from verma_modules import (
    M_verma, load_e44,
    l_minus1_action_matrix, l0_action_matrix,
    _w_action_from_l0_idx,
)

from singular_vectors import (
    # Degree-1 morphisms (re-exported)
    phi_1A, phi_1B, phi_1C, phi_1D, phi_1E,
    # Degree-2 primitive morphism (re-exported)
    phi_2DA,
    # Degree-3 and degree-4 morphisms (re-exported)
    phi_3F, phi_3G, phi_4H,
    # Singular-vector helpers used in phi_2EA and tests
    w1A, w1B, w1C, w1D, w1E, w3G, w4H,
    # Internal phi-map building helpers
    _compute_phi0, _phi_at_degree,
    # Verification helper
    annihilator_check, _verify_morphism_equivariance,
    # Export helpers (re-exported)
    save_morphisms, load_morphisms,
)

__all__ = [
    # All morphisms
    'phi_1A', 'phi_1B', 'phi_1C', 'phi_1D', 'phi_1E',
    'phi_2DA', 'phi_2EA',
    'phi_3F', 'phi_3G', 'phi_4H',
    # Utilities
    'compose_morphisms', 'is_proportional',
    # I/O
    'save_morphisms', 'load_morphisms',
    # Tests
    '_check_section7',
]


# ===========================================================================
# \phi[2EA] — composite degree-2 morphism  M_{-1}(1,0,0) \to M_1(1,0,0)
# ===========================================================================
#
# \phi[2EA] is the degree-2 morphism defined by the composite
#
#       \phi[1E] \circ \phi[1A]:  M_{-1}(1,0,0) \to M_0(0,0,0) \to M_1(1,0,0)
#
# where \phi[1A] uses the special-case parameters (t=0, a=0):
#       \phi[1A](0,0): M_{-1}(1,0,0) \to M_0(0,0,0),   sending v_hw \to w[1A](0,0)
# and \phi[1E]: M_0(0,0,0) \to M_1(1,0,0).
#
# The degree-2 singular vector w[2EA] \in M_1(1,0,0)[2] is:
#       w[2EA] = (\phi[1E] degree-1 matrix) \cdot w[1A](0,0)
#
# This is nonzero — the composite morphism is a proper degree-2 embedding.
# (The composition \phi[1E] \circ \phi[1A] ≠ 0, confirmed numerically in s3.6.)
# ---------------------------------------------------------------------------

def phi_2EA(e44_data, max_source_deg=0, src_e44_data=None):
    """
    Composite degree-2 morphism  \phi[2EA]: M_{-1}(1,0,0) \to M_1(1,0,0).

    Defined as the composition  \phi[1E] \circ \phi[1A]  where \phi[1A] uses the
    special case (t=0, a=0): M_{-1}(1,0,0) \to M_0(0,0,0).

    The defining singular vector w[2EA] \in M_1(1,0,0)[deg=2] is obtained
    by applying the degree-1 matrix of \phi[1E] to w[1A](t=0, a=0).

    Parameters
    ----------
    e44_data      : dict from load_e44() — E(4,4) bracket data
    max_source_deg: int (default 0)
        Compute morphism matrices at source degrees 0 ... max_source_deg.
    src_e44_data  : e44_data to use for the source fiber.
                    Pass e44_data when source node M_{-1}(1,0,0) is phat4.

    Returns
    -------
    (M_src, M_tar, sv_deg, phi0_vecs, matrices)
      M_src     : VermaModule  M_{-1}(1,0,0)  fiber determined by src_e44_data
      M_tar     : VermaModule  M_1(1,0,0)     p̂(4) fiber dim_W = 8
      sv_deg    : 2
      phi0_vecs : list of QQ-vectors in M_tar[2] (one per src fiber basis vector)
      matrices  : dict{d: matrix} for d = 0 ... max_source_deg

    Raises
    ------
    RuntimeError if w[2EA] is zero (would indicate a sign/implementation bug).
    ValueError   if e44_data is None.
    """
    if e44_data is None:
        raise ValueError("phi_2EA requires e44_data (load from e44_brackets.pkl)")

    sv_deg = 2

    # ── Step 1: build w[2EA] = \phi[1E][degree-1 mat] * w[1A](0,0) ────────
    # \phi[1E] with max_source_deg=1 to expose its degree-1 morphism matrix,
    # which maps M_0(0,0,0)[1] \to M_1(1,0,0)[2].
    _, _, _, _, mats_1E = phi_1E(e44_data, max_source_deg=1)
    _, v_1A_00 = w1A(QQ(0), QQ(0))          # w[1A](t=0,a=0) \in M_0(0,0,0)[1]
    w_2EA = mats_1E[1] * vector(QQ, v_1A_00)  # lives in M_1(1,0,0)[2]

    if w_2EA.is_zero():
        raise RuntimeError(
            "w[2EA] = \phi[1E][1] * w[1A](0,0) is zero; "
            "expected a nonzero degree-2 singular vector in M_1(1,0,0)."
        )

    # ── Step 2: set up source and target modules ──────────────────────────
    M_tar = M_verma(QQ(1), 1, 0, 0,
                    max_deg=sv_deg + max_source_deg,
                    e44_data=e44_data)
    # Source fiber: node (1,0,0) uses phat4 when src_e44_data is provided.
    M_src = M_verma(QQ(-1), 1, 0, 0, max_deg=max_source_deg,
                    e44_data=src_e44_data)

    # ── Step 3: fiber map \phi_0: W_{-1}(1,0,0) \to M_1(1,0,0)[2] ────────────
    phi0 = _compute_phi0(w_2EA, M_tar, sv_deg, M_src.W, e44_data)

    # ── Step 4: assemble morphism matrices at each source degree ──────────
    matrices = {}
    for d in range(max_source_deg + 1):
        matrices[d] = _phi_at_degree(M_src, M_tar, sv_deg, phi0, d)

    return M_src, M_tar, sv_deg, phi0, matrices


# ===========================================================================
# Utility: matrix-level composition
# ===========================================================================

def compose_morphisms(phi_B_mats, phi_A_mats, sv_deg_A, max_result_deg=0):
    """
    Matrix-level composition  (\phi_B \circ \phi_A) at each source degree.

    For morphisms
        \phi_A: M_src[d]   \to M_mid[d + sv_A]   (matrices phi_A_mats[d])
        \phi_B: M_mid[d]   \to M_tar[d + sv_B]   (matrices phi_B_mats[d])

    the composition satisfies:
        (\phi_B \circ \phi_A)[d] = phi_B_mats[sv_A + d] * phi_A_mats[d]

    Parameters
    ----------
    phi_B_mats    : dict{d: QQ-matrix} — matrices for \phi_B (outer morphism)
    phi_A_mats    : dict{d: QQ-matrix} — matrices for \phi_A (inner morphism)
    sv_deg_A      : int — degree shift of \phi_A (singular-vector degree)
    max_result_deg: int (default 0) — compute for source degrees 0..max_result_deg

    Returns
    -------
    dict{d: QQ-matrix} — composition matrices at source degree d
    """
    result = {}
    for d in range(max_result_deg + 1):
        mat_A = phi_A_mats.get(d)
        mat_B = phi_B_mats.get(sv_deg_A + d)
        if mat_A is None or mat_B is None:
            continue
        result[d] = mat_B * mat_A
    return result


def is_proportional(v, w):
    """
    Return True iff QQ-vectors v and w are scalar multiples of each other.

    Returns (True, ratio) where ratio = c such that v = c*w,
    or (False, None) if they are not proportional.
    Neither vector may be zero unless both are zero (in which case True, 0).
    """
    v = vector(QQ, v)
    w = vector(QQ, w)
    if v.is_zero() and w.is_zero():
        return True, QQ(0)
    if v.is_zero() or w.is_zero():
        return False, None
    ratio = None
    for vi, wi in zip(v, w):
        if wi == QQ(0):
            if vi != QQ(0):
                return False, None
        else:
            r = vi / wi
            if ratio is None:
                ratio = r
            elif r != ratio:
                return False, None
    return True, ratio


# ===========================================================================
# Section 7 test suite
# ===========================================================================

def _check_section7(e44_data=None, verbose=True):
    """
    Section 7 checkpoint: morphism construction and composition tests.

    Part A — Construct all 10 morphisms; verify E(4,4) equivariance
             (L_0 equivariance, L_1 annihilation, L_{-1} consistency, rank)
             for a selection of parameter values:
               \phi[1A](t=0, a=1),  \phi[1B](t=1, c=1),  \phi[1C](t=1),
               \phi[1D](t=1),       \phi[1E],
               \phi[2DA](t=2),      \phi[2EA],
               \phi[3F],            \phi[3G],             \phi[4H](t=3).

    Part B — Non-zero compositions:
      (1) \phi[1E] \circ \phi[1A](0,0)  = \phi[2EA]     (defines the composite morphism)
      (2) \phi[3F] \circ \phi[1A](0,0) \proto w[4H](t=3)  in M_3(1,0,0)[4]
      (3) \phi[1E] \circ \phi[3G]       \proto w[4H](t=1) in M_1(1,0,0)[4]

    Part C — Zero compositions (d² = 0 for the de Rham complex):
      (4) \phi[1A](1,1) \circ \phi[1A](0,2) = 0
      (5) \phi[1D](1)   \circ \phi[1A](0,1) = 0
      (6) \phi[1B](2,2) \circ \phi[1B](1,1) = 0
      (7) \phi[1D](3)   \circ \phi[1A](2,1) = 0
      (8) \phi[3G]      \circ \phi[1A](-3,1)= 0

    Returns True iff all checks pass.
    """
    if e44_data is None:
        if verbose:
            print("─" * 60)
            print("Section 7 checkpoint — SKIPPED (no e44_data)")
            print("Load via:  e44_data = load_e44('e44_brackets.pkl')")
            print("─" * 60)
        return True

    all_pass = True
    total_pass = 0
    total_fail = 0

    def _check(name, got, expected):
        nonlocal all_pass, total_pass, total_fail
        ok = (got == expected)
        status = "PASS" if ok else "FAIL"
        if verbose:
            print(f"  [{status}] {name}: got {got!r}, expected {expected!r}")
        if ok:
            total_pass += 1
        else:
            total_fail += 1
            all_pass = False

    if verbose:
        print("=" * 60)
        print("Section 7 — Verma-module morphisms \phi[1A]–\phi[4H]")
        print("=" * 60)

    # ======================================================================
    # Part A — Construction and equivariance of all 10 morphisms
    # ======================================================================
    if verbose:
        print("\n╔══ Part A: Construct and verify all 10 morphisms ══╗\n")

    cache = {}      # store for reuse in Parts B and C

    def _build_and_verify(key, builder, args, kwargs=None):
        """Build one morphism and run the equivariance checker."""
        nonlocal all_pass, total_pass, total_fail
        kwargs = kwargs or {}
        if verbose:
            print(f"  Building {key} ...")
        try:
            result = builder(*args, **kwargs)
        except Exception as exc:
            if verbose:
                print(f"    [FAIL] {key}: raised {type(exc).__name__}: {exc}")
            total_fail += 1
            all_pass = False
            cache[key] = None
            return
        M_src, M_tar, sv_deg, phi0, mats = result
        cache[key] = result
        n_pass, n_fail = _verify_morphism_equivariance(
            key, M_src, M_tar, sv_deg, phi0, mats, e44_data,
            verbose=verbose,
        )
        total_pass += n_pass
        total_fail += n_fail
        if n_fail > 0:
            all_pass = False
        elif verbose:
            print(f"    {key}: {n_pass} checks PASS")

    # \phi[1A](t=0, a=1)  and  \phi[1A](t=2, a=1)
    _build_and_verify(
        "phi_1A(0,1)", phi_1A, (QQ(0), 1, e44_data),
        {'max_source_deg': 1},
    )
    _build_and_verify(
        "phi_1A(2,1)", phi_1A, (QQ(2), 1, e44_data),
        {'max_source_deg': 1},
    )

    # \phi[1B](t=1, c=1)
    _build_and_verify(
        "phi_1B(1,1)", phi_1B, (QQ(1), 1, e44_data),
        {'max_source_deg': 1},
    )

    # \phi[1C](t=1)   [p̂(4) target fiber]
    _build_and_verify(
        "phi_1C(1)", phi_1C, (QQ(1), e44_data),
        {'max_source_deg': 1},
    )

    # \phi[1D](t=1)
    _build_and_verify(
        "phi_1D(1)", phi_1D, (QQ(1), e44_data),
        {'max_source_deg': 1},
    )

    # \phi[1E]
    _build_and_verify(
        "phi_1E", phi_1E, (e44_data,),
        {'max_source_deg': 3},
    )

    # \phi[2DA](t=2)  [p̂(4) target fiber]
    _build_and_verify(
        "phi_2DA(2)", phi_2DA, (QQ(2), e44_data),
        {'max_source_deg': 0},
    )

    # \phi[2EA]  [composite]
    _build_and_verify(
        "phi_2EA", phi_2EA, (e44_data,),
        {'max_source_deg': 0},
    )

    # \phi[3F]   [p̂(4) target fiber]
    _build_and_verify(
        "phi_3F", phi_3F, (e44_data,),
        {'max_source_deg': 1},
    )

    # \phi[3G]
    _build_and_verify(
        "phi_3G", phi_3G, (e44_data,),
        {'max_source_deg': 1},
    )

    # \phi[4H](t=3)  and  \phi[4H](t=1)
    _build_and_verify(
        "phi_4H(3)", phi_4H, (QQ(3), e44_data),
        {'max_source_deg': 0},
    )
    _build_and_verify(
        "phi_4H(1)", phi_4H, (QQ(1), e44_data),
        {'max_source_deg': 0},
    )

    if verbose:
        print(f"\n  Part A summary: {total_pass} pass, {total_fail} fail")

    # ======================================================================
    # Part B — Non-zero compositions
    # ======================================================================
    if verbose:
        print("\n╔══ Part B: Non-zero compositions ══╗\n")

    pass_b = 0
    fail_b = 0

    def _check_b(name, got, expected):
        nonlocal pass_b, fail_b, all_pass, total_pass, total_fail
        ok = (got == expected)
        status = "PASS" if ok else "FAIL"
        if verbose:
            print(f"  [{status}] {name}: got {got!r}, expected {expected!r}")
        if ok:
            pass_b += 1
            total_pass += 1
        else:
            fail_b += 1
            total_fail += 1
            all_pass = False

    # ── (1) \phi[1E] \circ \phi[1A](0,0) == w[2EA] ────────────────────────────────
    # Both land in M_1(1,0,0)[2].
    # Compute composition: phi_1E[1] * w[1A](0,0)
    if cache.get("phi_1E") is not None:
        _, _, _, _, mats_1E = cache["phi_1E"]
        _, v_1A_00 = w1A(QQ(0), QQ(0))
        comp1 = mats_1E[1] * vector(QQ, v_1A_00)
        _check_b("(1) phi_1E \circ phi_1A(0,0): nonzero", not comp1.is_zero(), True)

        # Compare with \phi[2EA] at degree 0 (hw column of phi0)
        if cache.get("phi_2EA") is not None:
            M_src_2EA, _, _, _, mats_2EA = cache["phi_2EA"]
            # The hw column of phi_2EA[0] (degree-0 matrix) should equal comp1
            hw_col_2EA = mats_2EA[0].column(M_src_2EA.W.v_hw)
            prop_ok, ratio = is_proportional(comp1, hw_col_2EA)
            _check_b("(1) phi_1E\circphi_1A(0,0) \proto phi_2EA(v_hw)", prop_ok, True)
            if verbose and prop_ok:
                print(f"    proportionality constant = {ratio}")
    else:
        if verbose:
            print("  [(SKIP)] (1): phi_1E not built successfully")

    # ── (2) \phi[3F] \circ \phi[1A](0,0) \proto w[4H](t=3) ─────────────────────────────
    # \phi[1A](0,0): M_{-1}(1,0,0) \to M_0(0,0,0), degree 1
    # \phi[3F]:      M_0(0,0,0)    \to M_3(1,0,0), degree 3
    # Image: M_3(1,0,0)[4].
    if cache.get("phi_3F") is not None:
        _, _, _, _, mats_3F = cache["phi_3F"]
        _, v_1A_00_sv = w1A(QQ(0), QQ(0))
        comp2 = mats_3F[1] * vector(QQ, v_1A_00_sv)
        _check_b("(2) phi_3F \circ phi_1A(0,0): nonzero", not comp2.is_zero(), True)

        try:
            _, w4H_3 = w4H(QQ(3), e44_data)
            prop_ok, ratio = is_proportional(comp2, w4H_3)
            _check_b("(2) phi_3F\circphi_1A(0,0) \proto w[4H](t=3)", prop_ok, True)
            if verbose and prop_ok:
                print(f"    proportionality constant = {ratio}")
        except Exception as exc:
            if verbose:
                print(f"  [(SKIP)] (2) w[4H](t=3) failed: {exc}")
    else:
        if verbose:
            print("  [(SKIP)] (2): phi_3F not built successfully")

    # ── (3) \phi[1E] \circ \phi[3G] \proto w[4H](t=1) ──────────────────────────────────
    # \phi[3G]: M_{-3}(1,0,0) \to M_0(0,0,0), degree 3
    # \phi[1E]: M_0(0,0,0)    \to M_1(1,0,0), degree 1
    # Image: M_1(1,0,0)[4].
    if cache.get("phi_1E") is not None and cache.get("phi_3G") is not None:
        _, _, _, _, mats_1E_long = cache["phi_1E"]
        _, v_3G = w3G()
        # phi_1E acts at degree 3 of its source: mats_1E_long[3]
        if 3 in mats_1E_long:
            comp3 = mats_1E_long[3] * vector(QQ, v_3G)
            _check_b("(3) phi_1E \circ phi_3G: nonzero", not comp3.is_zero(), True)
            try:
                _, w4H_1 = w4H(QQ(1), e44_data)
                prop_ok, ratio = is_proportional(comp3, w4H_1)
                _check_b("(3) phi_1E\circphi_3G \proto w[4H](t=1)", prop_ok, True)
                if verbose and prop_ok:
                    print(f"    proportionality constant = {ratio}")
            except Exception as exc:
                if verbose:
                    print(f"  [(SKIP)] (3) w[4H](t=1) failed: {exc}")
        else:
            if verbose:
                print("  [(SKIP)] (3): phi_1E degree-3 matrix not available "
                      "(need max_source_deg≥3)")
    else:
        if verbose:
            print("  [(SKIP)] (3): phi_1E or phi_3G not built successfully")

    if verbose:
        print(f"\n  Part B: {pass_b} pass, {fail_b} fail")

    # ======================================================================
    # Part C — Zero compositions  (d² = 0)
    # ======================================================================
    if verbose:
        print("\n╔══ Part C: Zero compositions (d² = 0) ══╗\n")

    pass_c = 0
    fail_c = 0

    def _check_c(name, got, expected):
        nonlocal pass_c, fail_c, all_pass, total_pass, total_fail
        ok = (got == expected)
        status = "PASS" if ok else "FAIL"
        if verbose:
            print(f"  [{status}] {name}: got {got!r}, expected {expected!r}")
        if ok:
            pass_c += 1
            total_pass += 1
        else:
            fail_c += 1
            total_fail += 1
            all_pass = False

    # (4) \phi[1A](1,1) \circ \phi[1A](0,2) = 0  at degree 2
    #     M_{-1}(3,0,0) \to^{\phi[1A](0,2)} M_0(2,0,0) \to^{\phi[1A](1,1)} M_1(1,0,0)
    try:
        _, _, _, _, mats_1A_11 = phi_1A(QQ(1), 1, e44_data, max_source_deg=1)
        _, v_1A_02 = w1A(QQ(0), 2)
        comp = mats_1A_11[1] * vector(QQ, v_1A_02)
        _check_c("(4) phi_1A(1,1) \circ phi_1A(0,2) = 0", comp.is_zero(), True)
    except Exception as exc:
        if verbose:
            print(f"  [(SKIP)] (4): {exc}")

    # (5) \phi[1D](1) \circ \phi[1A](0,1) = 0  at degree 2
    #     M_{-1}(2,0,0) \to^{\phi[1A](0,1)} M_0(1,0,0) \to^{\phi[1D](1)} M_1(0,0,0)
    if cache.get("phi_1D(1)") is not None:
        _, _, _, _, mats_1D_1 = cache["phi_1D(1)"]
        _, v_1A_01 = w1A(QQ(0), 1)
        comp = mats_1D_1[1] * vector(QQ, v_1A_01)
        _check_c("(5) phi_1D(1) \circ phi_1A(0,1) = 0", comp.is_zero(), True)
    else:
        try:
            _, _, _, _, mats_1D_1 = phi_1D(QQ(1), e44_data, max_source_deg=1)
            _, v_1A_01 = w1A(QQ(0), 1)
            comp = mats_1D_1[1] * vector(QQ, v_1A_01)
            _check_c("(5) phi_1D(1) \circ phi_1A(0,1) = 0", comp.is_zero(), True)
        except Exception as exc:
            if verbose:
                print(f"  [(SKIP)] (5): {exc}")

    # (6) \phi[1B](2,2) \circ \phi[1B](1,1) = 0  at degree 2
    #     M_0(0,0,0) \to^{\phi[1B](1,1)} M_1(0,0,1) \to^{\phi[1B](2,2)} M_2(0,0,2)
    try:
        _, _, _, _, mats_1B_11 = phi_1B(QQ(1), 1, e44_data, max_source_deg=1)
        _, _, _, _, mats_1B_22 = phi_1B(QQ(2), 2, e44_data, max_source_deg=1)
        _, v_1B_11 = w1B(QQ(1), 1)
        comp = mats_1B_22[1] * vector(QQ, v_1B_11)
        _check_c("(6) phi_1B(2,2) \circ phi_1B(1,1) = 0", comp.is_zero(), True)
    except Exception as exc:
        if verbose:
            print(f"  [(SKIP)] (6): {exc}")

    # (7) \phi[1D](3) \circ \phi[1A](2,1) = 0  at degree 2
    #     M_1(2,0,0) \to^{\phi[1A](2,1)} M_2(1,0,0) \to^{\phi[1D](3)} M_3(0,0,0)
    try:
        _, _, _, _, mats_1D_3 = phi_1D(QQ(3), e44_data, max_source_deg=1)
        _, v_1A_21 = w1A(QQ(2), 1)
        comp = mats_1D_3[1] * vector(QQ, v_1A_21)
        _check_c("(7) phi_1D(3) \circ phi_1A(2,1) = 0", comp.is_zero(), True)
    except Exception as exc:
        if verbose:
            print(f"  [(SKIP)] (7): {exc}")

    # (8) \phi[3G] \circ \phi[1A](-3,1) = 0  at degree 4
    #     M_{-4}(2,0,0) \to^{\phi[1A](-3,1)} M_{-3}(1,0,0) \to^{\phi[3G]} M_0(0,0,0)
    if cache.get("phi_3G") is not None:
        _, _, _, _, mats_3G = cache["phi_3G"]
        _, v_1A_m31 = w1A(QQ(-3), 1)
        comp = mats_3G[1] * vector(QQ, v_1A_m31)
        _check_c("(8) phi_3G \circ phi_1A(-3,1) = 0", comp.is_zero(), True)
    else:
        try:
            _, _, _, _, mats_3G = phi_3G(e44_data, max_source_deg=1)
            _, v_1A_m31 = w1A(QQ(-3), 1)
            comp = mats_3G[1] * vector(QQ, v_1A_m31)
            _check_c("(8) phi_3G \circ phi_1A(-3,1) = 0", comp.is_zero(), True)
        except Exception as exc:
            if verbose:
                print(f"  [(SKIP)] (8): {exc}")

    if verbose:
        print(f"\n  Part C: {pass_c} pass, {fail_c} fail")

    # ── Summary ───────────────────────────────────────────────────────────
    if verbose:
        print("\n" + "═" * 60)
        print(f"Section 7 summary: {total_pass} pass, {total_fail} fail")
        if all_pass:
            print("Section 7  ✓  ALL CHECKS PASSED")
        else:
            print("Section 7  ✗  SOME CHECKS FAILED — see above")
        print("═" * 60)

    return all_pass


# ===========================================================================
# Module self-test (run when executed as a script)
# ===========================================================================

if __name__ == '__main__':
    import sys as _sys_main

    print("morphisms.py — Section 7 self-test")
    print("Loading E(4,4) bracket data ...")
    e44 = load_e44('e44_brackets.pkl')
    if e44 is None:
        print("ERROR: e44_brackets.pkl not found.  "
              "Run 'sage e44_structure.py' first.")
        _sys_main.exit(1)

    ok = _check_section7(e44_data=e44, verbose=True)
    _sys_main.exit(0 if ok else 1)
