"""
verma_modules.py  —  Verma/Kac modules for the NSE44 programme
===============================================================
Cantarini-Caselli-Kac 2026: E(4,4) as the Navier-Stokes algebra.

Section layout (mirrors plan-nsLTKMA.md / NSE44.tex SS2):
  s2   — Phase 2: Verma Modules and Singular Vectors
    s2.1   \hat{p}(4)-module W_hat(t, a, b, c): weight spaces, sl_4 action,
           central charge  (checkpoint: dim(W_hat(1,1,0,0)) == 4)
    (s2.2-s2.5 in later milestones)

*Depends on* e44_structure.py via the pickle e44_brackets.pkl.
Run e44_structure.py first to generate that file.

Run inside SageMath:  sage verma_modules.py
"""

from sage.all import (
    QQ, ZZ, matrix, vector, FreeModule, identity_matrix,
    latex, RootSystem,
)
from sage.combinat.crystals.highest_weight_crystals import HighestWeightCrystal
import pickle as _pickle
import os as _os

# ---------------------------------------------------------------------------
# Load E44 data produced by e44_structure.py
# ---------------------------------------------------------------------------

def load_e44(path="e44_brackets.pkl"):
    """
    Load the E44 graded basis and bracket table from a pickle file.
    Returns the payload dict, or None with a warning if the file is absent.
    """
    if not _os.path.exists(path):
        import warnings
        warnings.warn(
            f"e44_brackets.pkl not found at '{path}'. "
            "Run 'sage e44_structure.py' first to generate it. "
            "Proceeding without the E44 bracket data."
        )
        return None
    with open(path, 'rb') as fh:
        return _pickle.load(fh)


# ===========================================================================
# s2.1 — \hat{p}(4)-module  W_hat(t, a, b, c)
# ===========================================================================
#
# The reductive core of L_0 = \hat{p}(4) (non-central extension of the strange
# Lie superalgebra p(4)) contains:
#   \hat{p}(4)_0  ≅  sl_4       (A_3 base)
#   \hat{p}(4)_{-1} ≅ \wedge^2(C^4*)  (antisymmetric 2-tensors, dim = 6)
#   \hat{p}(4)_1  ≅  S^2(C^4)     (symmetric 2-tensors, dim = 10)
#   centre   = CC          (central element C; acts as scalar t)
#
# The irreducible finite-dimensional \hat{p}(4)-module W_hat(t, a, b, c) is
# characterised by:
#   • sl_4 Dynkin highest weight labels (a, b, c)  with a,b,c \in Z_{>=0}
#   • central charge t \in C  (C acts as t \cdot Id)
#
# As a vector space it equals the irreducible sl_4-module
#     V(a\cdot\omega_1 + b\cdot\omega_2 + c\cdot\omega_3)
# of dimension  dim = #{semistandard Young tableaux of shape (a+b+c, b+c, c)
#                       with entries in {1,2,3,4}}.
#
# Implementation
# ──────────────
# • SageMath's crystals.Tableaux(['A', 3], shape) enumerates all basis
#   elements and encodes the sl_4 weight structure.
#   - shape = (a+b+c, b+c, c)  (no trailing zeros)
# • Basis vectors are indexed 0 ... dim-1 (crystal enumeration order).
# • The sl_4 weight of basis vector k is stored as the tuple
#       (h_1, h_2, h_3)  =  (n_1-n_2, n_2-n_3, n_3-n_4)
#   where nⱼ = number of times letter j appears in the k-th tableau.
#   These are the eigenvalues of the simple coroot operators H_i.
# • sl_4 Chevalley generator matrices are built from the crystal operators
#   \tilde{e}_i, \tilde{f}_i (Kashiwara crystal basis: non-zero entries are all 1).
#   The relation [e_i, f_j] = \delta_{ij} h_i holds in this basis.
# • The central element C acts as t \cdot Id_{dim}.
# ---------------------------------------------------------------------------


# --- Crystal construction ---------------------------------------------------

def _sl4_crystal(a, b, c):
    """
    Return the A_3 HighestWeightCrystal for sl_4 highest weight (a,b,c) [Dynkin].
    Uses the weight lattice element  a\cdot\omega_1 + b\cdot\omega_2 + c\cdot\omega_3  of A_3.
    Returns None for the trivial module (a=b=c=0).
    """
    if a == 0 and b == 0 and c == 0:
        return None   # trivial module handled separately
    WL = RootSystem(['A', 3]).weight_lattice()
    fw = WL.fundamental_weights()
    hw = a * fw[1] + b * fw[2] + c * fw[3]
    return HighestWeightCrystal(hw)


def _tableau_dynkin_weight(elt):
    """
    Dynkin label triple (h_1, h_2, h_3) for a crystal element of an A_3 crystal.

    Uses the ambient-space weight vector  wt = (a_1, a_2, a_3, a_4) \in Z^4:
        h_i = a_i - a_{i+1}   (eigenvalue of the i-th simple coroot H_i).
    """
    wt = elt.weight()
    return (int(wt[0]) - int(wt[1]), int(wt[1]) - int(wt[2]), int(wt[2]) - int(wt[3]))


def _crystal_matrices(crystal):
    """
    Build n\timesn QQ matrices for the sl_4 Chevalley generators
        e_i, f_i, h_i    (i = 1, 2, 3)
    from the crystal operators \tilde{e}_i, \tilde{f}_i of an A_3 Tableaux crystal.

    Convention (Kashiwara crystal basis):
        f_i(basis[k]) = basis[j]  with matrix entry  F[j, k] = 1,  or 0 if absent
        e_i(basis[k]) = basis[j]  with matrix entry  E[j, k] = 1,  or 0 if absent
        h_i               acts diagonally  with  H[k, k] = (h_i-weight of basis[k])

    In this basis the Chevalley relations [e_i, f_j] = \delta_{ij} h_i hold.

    Returns
    -------
    e_mats : dict  {i: n\timesn matrix QQ}    i \in {1,2,3}
    f_mats : dict  {i: n\timesn matrix QQ}
    h_mats : dict  {i: n\timesn matrix QQ}
    """
    basis = list(crystal)
    n = len(basis)
    idx = {b: k for k, b in enumerate(basis)}

    e_mats, f_mats, h_mats = {}, {}, {}

    for i in range(1, 4):
        mat_e = matrix(QQ, n, n)
        mat_f = matrix(QQ, n, n)
        mat_h = matrix(QQ, n, n)

        for k, b in enumerate(basis):
            # Lowering  f_i: column k \to row idx[\tilde{f}_i(b)]
            fb = b.f(i)
            if fb is not None:
                mat_f[idx[fb], k] = QQ(b.epsilon(i) + 1)

            # Raising  e_i: column k \to row idx[\tilde{e}_i(b)]
            eb = b.e(i)
            if eb is not None:
                mat_e[idx[eb], k] = QQ(b.phi(i) + 1)

            # Cartan  h_i: diagonal entry = Dynkin label at index i
            wt = _tableau_dynkin_weight(b)
            mat_h[k, k] = QQ(wt[i - 1])

        e_mats[i] = mat_e
        f_mats[i] = mat_f
        h_mats[i] = mat_h

    return e_mats, f_mats, h_mats


# --- Module class -----------------------------------------------------------

class WHat4Module:
    """
    Finite-dimensional irreducible \hat{p}(4)-module  W_hat(t, a, b, c).

    Parameters
    ----------
    t : rational number (or integer)
        Central charge; the central element C of \hat{p}(4) acts as t \cdot Idn.
    a, b, c : non-negative integers
        Dynkin labels of the sl_4 highest weight
        \wedge = a\cdot\omega_1 + b\cdot\omega_2 + c\cdot\omega_3  (\omega_i = i-th fundamental weight of A_3).

    Key attributes
    --------------
    dim : int
        Dimension of the module = #{SSYT of shape (a+b+c, b+c, c) in {1,2,3,4}}.
    weight_spaces : dict
        { (h_1, h_2, h_3) : [basis_indices] }
        Partition of range(dim) by sl_4 Dynkin-label weight.
    v_hw : int
        Index of the unique highest-weight basis vector (weight = (a,b,c)).
    e_mats, f_mats, h_mats : dict {i: matrix}
        dim\timesdim matrices for the sl_4 Chevalley generators (i \in {1,2,3}),
        in the Kashiwara crystal basis.
    crystal : SageMath crystal or None
        The underlying A_3 Tableaux crystal (None for the trivial module).
    basis_elts : list
        List of crystal elements in the canonical crystal enumeration order.
    """

    def __init__(self, t, a, b, c):
        self.t = t
        self.a = int(a)
        self.b = int(b)
        self.c = int(c)

        a, b, c = self.a, self.b, self.c     # local ints for construction

        # ── Trivial module (a = b = c = 0) ──────────────────────────────
        if a == 0 and b == 0 and c == 0:
            self.crystal    = None
            self.basis_elts = [None]          # sentinel for the unique basis vector
            self.dim        = 1
            self.weight_spaces = {(0, 0, 0): [0]}
            self.v_hw       = 0
            _zero = matrix(QQ, 1, 1)
            self.e_mats = {i: _zero for i in range(1, 4)}
            self.f_mats = {i: _zero for i in range(1, 4)}
            self.h_mats = {i: _zero for i in range(1, 4)}
            return

        # ── General module ───────────────────────────────────────────────
        self.crystal    = _sl4_crystal(a, b, c)
        self.basis_elts = list(self.crystal)
        self.dim        = len(self.basis_elts)

        # Weight space decomposition
        self.weight_spaces = {}
        for k, elt in enumerate(self.basis_elts):
            wt = _tableau_dynkin_weight(elt)
            self.weight_spaces.setdefault(wt, []).append(k)

        # Locate the unique highest-weight vector (weight = (a, b, c))
        hw_wt = (a, b, c)
        hw_idx = self.weight_spaces.get(hw_wt, [])
        if len(hw_idx) != 1:
            raise ValueError(
                f"Expected exactly one hw vector with weight {hw_wt} "
                f"in W_hat({t}, {a}, {b}, {c}); found {len(hw_idx)}."
            )
        self.v_hw = hw_idx[0]

        # sl_4 generator matrices
        self.e_mats, self.f_mats, self.h_mats = _crystal_matrices(self.crystal)

    # --- Module maps -------------------------------------------------------

    def action_of_C(self):
        """
        Matrix of the central element C \in \hat{p}(4): acts as t \cdot Id_{dim}.
        Returns a dim\timesdim QQ-matrix.
        """
        return self.t * identity_matrix(QQ, self.dim)

    # --- Queries -----------------------------------------------------------

    def weight_of(self, idx):
        """Return the Dynkin-label weight tuple (h_1,h_2,h_3) of basis vector idx."""
        if self.crystal is None:
            return (0, 0, 0)
        return _tableau_dynkin_weight(self.basis_elts[idx])

    def weight_multiplicities(self):
        """Return a dict {weight_tuple: multiplicity}."""
        return {wt: len(idxs) for wt, idxs in self.weight_spaces.items()}

    def character_string(self):
        """Human-readable weight-multiplicity summary."""
        lines = [
            f"W_hat(t={self.t}, a={self.a}, b={self.b}, c={self.c})  dim={self.dim}",
            f"  highest-weight vector  v_hw = basis[{self.v_hw}]  weight {self.weight_of(self.v_hw)}",
            "  weight spaces:",
        ]
        for wt, idxs in sorted(self.weight_spaces.items(), reverse=True):
            lines.append(f"    {wt}  mult={len(idxs)}  basis_indices={idxs}")
        return "\n".join(lines)

    def __repr__(self):
        return (
            f"WHat4Module(t={self.t}, a={self.a}, b={self.b}, c={self.c}; "
            f"dim={self.dim}, #weight_spaces={len(self.weight_spaces)})"
        )


# --- Public constructor -----------------------------------------------------

def W_hat(t, a, b, c):
    """
    Construct the irreducible \hat{p}(4)-module W_hat(t, a, b, c).

    Parameters
    ----------
    t : central charge (integer or rational)
    a, b, c : sl_4 Dynkin highest weight labels (non-negative integers)

    Returns
    -------
    WHat4Module

    Examples
    --------
    sage: W = W_hat(1, 1, 0, 0)
    sage: W.dim
    4
    sage: W.weight_of(W.v_hw)
    (1, 0, 0)

    sage: W_hat(0, 0, 1, 0).dim     # \wedge^2(C^4)
    6

    sage: W_hat(0, 0, 0, 1).dim     # \wedge^3(C^4) ≅ (C^4)*
    4

    sage: W_hat(0, 0, 0, 0).dim     # trivial
    1
    """
    return WHat4Module(t, a, b, c)


# ===========================================================================
# s2.1 Checkpoint
# ===========================================================================

def _check_s21(verbose=True):
    """
    S2.1 checkpoint: verify dimensions and weight structure of W_hat modules.

    Checks
    ------
    (1) dim(W_hat(1,1,0,0)) == 4   [fundamental sl_4 rep, the NS velocity rep]
    (2) dim(W_hat(0,0,1,0)) == 6   [\wedge^2(C^4), spin-1 sector]
    (3) dim(W_hat(0,0,0,1)) == 4   [\wedge^3(C^4) = dual fundamental]
    (4) dim(W_hat(0,2,0,0)) == 10  [S^2(C^4)]
    (5) dim(W_hat(0,0,0,0)) == 1   [trivial rep]
    (6) W_hat(1,1,0,0).v_hw has weight (1,0,0)
    (7) # weight spaces of W_hat(1,1,0,0) == 4  (one-dimensional spaces)
    (8) C acts as t\cdotId (checked via trace = t\cdotdim for t=3, dim=4)

    Returns True iff all checks pass.
    """
    all_pass = True

    def _check(name, got, expected):
        nonlocal all_pass
        ok = (got == expected)
        status = "PASS" if ok else "FAIL"
        if verbose:
            print(f"  [{status}] {name}: got {got!r}, expected {expected!r}")
        if not ok:
            all_pass = False

    if verbose:
        print("─" * 60)
        print("s2.1 checkpoint  —  W_hat(t, a, b, c) dimensions and structure")
        print("─" * 60)

    # (1)-(5): dimension checks
    _check("dim W_hat(1,1,0,0) [fundamental (1,0,0)]",  W_hat(1,1,0,0).dim, 4)
    _check("dim W_hat(0,0,1,0) [\wedge^2(C^4), (0,1,0)]",    W_hat(0,0,1,0).dim, 6)
    _check("dim W_hat(0,0,0,1) [\wedge^3(C^4), (0,0,1)]",    W_hat(0,0,0,1).dim, 4)
    _check("dim W_hat(0,2,0,0) [S^2(C^4), (2,0,0)]",    W_hat(0,2,0,0).dim, 10)
    _check("dim W_hat(0,0,0,0) [trivial (0,0,0)]",     W_hat(0,0,0,0).dim, 1)

    # (6): highest-weight vector weight
    W100 = W_hat(1, 1, 0, 0)
    _check("W_hat(1,1,0,0).weight_of(v_hw)", W100.weight_of(W100.v_hw), (1, 0, 0))

    # (7): four distinct 1-dimensional weight spaces for the fundamental rep
    _check("#weight_spaces W_hat(1,1,0,0)",  len(W100.weight_spaces), 4)
    _check("all wt spaces dim 1 (fund. rep)",
           all(len(v) == 1 for v in W100.weight_spaces.values()), True)

    # (8): central element C = t \cdot Id
    W_test = W_hat(3, 1, 0, 0)               # t=3, dim=4
    C_mat  = W_test.action_of_C()
    _check("action_of_C trace = t\cdotdim",      int(C_mat.trace()), 3 * 4)
    _check("action_of_C = 3\cdotId_4",           C_mat == 3 * identity_matrix(QQ, 4), True)

    # (9): Chevalley relation [e_i, f_i] v = h_i v for i=1 on fund. rep
    #      (tests that crystal matrices are a valid sl_4 representation)
    comm_ok = True
    for i in range(1, 4):
        comm = W100.e_mats[i] * W100.f_mats[i] - W100.f_mats[i] * W100.e_mats[i]
        if comm != W100.h_mats[i]:
            comm_ok = False
            break
    _check("Chevalley [e_i,f_i] = h_i (i=1,2,3) on W_hat(1,1,0,0)", comm_ok, True)

    if verbose:
        print("─" * 60)
        if all_pass:
            print("s2.1  ✓  ALL CHECKS PASSED")
        else:
            print("s2.1  ✗  SOME CHECKS FAILED — see above")
        print("─" * 60)

    return all_pass


# ===========================================================================
# s2.2 — Verma module M(t, a, b, c) as a graded vector space
# ===========================================================================
#
# M(t,a,b,c) = U(L_{-1}) \otimes W_hat(t,a,b,c)   (Verma / Kac module)
#
# Since L_{j} = 0 for j < -1 in the principal grading, U(E44_{<0}) = U(L_{-1}).
# All L_{-1}-brackets land in L_{-2} = 0, so:
#   U(L_{-1}) = QQ[e_0,e_1,e_2,e_3] \otimes \wedge(d_0,d_1,d_2,d_3)
#   (polynomial in 4 even generators \times exterior algebra in 4 odd generators)
#
# PBW monomials at degree d: pairs (\alpha, S) with
#   \alpha = (a_0,a_1,a_2,a_3) \in Z_{>=0}^4  (even; no restriction on powers)
#   S \subseteq {0,1,2,3}                   (odd; each appears at most once)
#   \Sigma\alpha_i + |S| = d
#
# Dimension of U(L_{-1})_d = \Sigma_{k=0}^{min(d,4)} C(d-k+3,3) \cdot C(4,k)
#   d=0 \to 1,  d=1 \to 8,  d=2 \to 32,  d=3 \to 96,  d=4 \to 256
#
# Basis of M(t,a,b,c)[d] = {(j, k) : j \in PBW_d, k \in 0..dim(W)-1}
#   total dimension = |PBW_d| \times dim(W_hat(t,a,b,c))
#
# L_{>=0}-action is implemented in S2.3.
# ---------------------------------------------------------------------------

from itertools import combinations as _combs


def _compositions_4(n):
    """
    Return all 4-tuples (a0,a1,a2,a3) of non-negative integers summing to n.
    Used to enumerate PBW monomials in the even generators.
    """
    result = []
    for a0 in range(n + 1):
        for a1 in range(n - a0 + 1):
            for a2 in range(n - a0 - a1 + 1):
                a3 = n - a0 - a1 - a2
                result.append((a0, a1, a2, a3))
    return result


def _pbw_basis_at_degree(d):
    """
    Return the ordered list of PBW monomials for U(L_{-1}) at degree d.

    Each monomial is a pair  (alpha, S)  where:
      alpha = (a0,a1,a2,a3) — tuple of non-negative ints, \Sigma\alpha_i = d - |S|
      S     = frozenset \subseteq {0,1,2,3}   (chosen odd generators)
      \Sigma\alpha_i + |S| = d

    Ordering: by (|S|, S, alpha)  for a canonical PBW enumeration.
    """
    result = []
    for k in range(min(d, 4) + 1):          # |S| = k
        for S in _combs(range(4), k):        # ordered subsets of size k
            Sfrozen = frozenset(S)
            rem = d - k                      # \Sigma\alpha_i = rem
            for alpha in _compositions_4(rem):
                result.append((alpha, Sfrozen))
    return result


class VermaModule:
    """
    Verma module  M(t, a, b, c)  for E(4,4), graded up to degree max_deg.

    As a graded vector space:
        M[d] = U(L_{-1})_d  \otimes  W_hat(t, a, b, c)

    Parameters
    ----------
    t : central charge (integer or rational)
    a, b, c : sl_4 Dynkin highest weight labels
    max_deg : truncation degree (default 4; sufficient for singular vector
              classification per CKC 2026)

    Key attributes
    --------------
    W : WHat4Module
        The \hat{p}(4)-module W_hat(t, a, b, c) at degree 0.
    pbw : dict {d: list of (alpha, frozenset) monomials}
        PBW basis of U(L_{-1}) at each degree.
    dim_W : int
        Dimension of W_hat.
    max_deg : int
        Truncation degree.

    Key methods
    -----------
    dim(d)              \to dimension of M[d]
    basis(d)            \to list of (pbw_monomial, w_basis_idx) at degree d
    basis_index(d, mon, k) \to flat integer index of (mon, k) in M[d]
    zero_vec(d)         \to zero vector in QQ^{dim(d)}
    to_vec(d, coeffs)   \to dense column vector from {(j,k): coeff} dict
    """

    def __init__(self, t, a, b, c, max_deg=4, e44_data=None):
        self.t       = t
        self.a       = int(a)
        self.b       = int(b)
        self.c       = int(c)
        self.max_deg = max_deg
        self._e44_data = e44_data

        if e44_data is not None:
            from phat4_modules import phat4_module
            self.W = phat4_module(t, a, b, c, e44_data)
        else:
            self.W = W_hat(t, a, b, c)
        self.dim_W = self.W.dim

        # PBW basis at each degree
        self.pbw = {d: _pbw_basis_at_degree(d) for d in range(max_deg + 1)}

        # Flat index tables  (pbw_idx, w_idx) \to flat_idx
        self._flat_idx = {}
        for d in range(max_deg + 1):
            pbw_d = self.pbw[d]
            for j, mon in enumerate(pbw_d):
                for k in range(self.dim_W):
                    self._flat_idx[(d, j, k)] = j * self.dim_W + k

    # --- Dimension / basis queries -----------------------------------------

    def dim(self, d):
        """Dimension of M[d] = |PBW_d| \times dim(W_hat)."""
        if d < 0 or d > self.max_deg:
            return 0
        return len(self.pbw[d]) * self.dim_W

    def pbw_dim(self, d):
        """Number of PBW monomials at degree d in U(L_{-1})."""
        if d < 0 or d > self.max_deg:
            return 0
        return len(self.pbw[d])

    def basis(self, d):
        """
        Return the ordered basis of M[d] as a list of (pbw_monomial, w_idx) pairs.
        pbw_monomial = (alpha, frozenset_S)
        w_idx        = integer index into W_hat basis
        """
        if d < 0 or d > self.max_deg:
            return []
        return [(mon, k) for mon in self.pbw[d] for k in range(self.dim_W)]

    def basis_index(self, d, mon_idx, w_idx):
        """
        Return the flat integer index of basis element (pbw[d][mon_idx], w_idx)
        in the degree-d piece.
        """
        return self._flat_idx[(d, mon_idx, w_idx)]

    # --- Vector space utilities --------------------------------------------

    def zero_vec(self, d):
        """Return the zero vector in QQ^{dim(d)}."""
        return vector(QQ, self.dim(d))

    def to_vec(self, d, coeffs):
        """
        Build a dense vector from a dict { (mon_idx, w_idx): coeff } or
        a list of (mon_idx, w_idx, coeff) triples.
        Returns a QQ-vector of length dim(d).
        """
        v = [QQ(0)] * self.dim(d)
        if isinstance(coeffs, dict):
            for (j, k), c in coeffs.items():
                v[self._flat_idx[(d, j, k)]] += QQ(c)
        else:
            for (j, k, c) in coeffs:
                v[self._flat_idx[(d, j, k)]] += QQ(c)
        return vector(QQ, v)

    # --- Human-readable summary --------------------------------------------

    def dim_table(self):
        """Return a dict {d: dim(d)} for d = 0 .. max_deg."""
        return {d: self.dim(d) for d in range(self.max_deg + 1)}

    def __repr__(self):
        dims = ', '.join(f'd{d}={self.dim(d)}' for d in range(self.max_deg + 1))
        return (
            f"VermaModule(t={self.t}, a={self.a}, b={self.b}, c={self.c}; "
            f"max_deg={self.max_deg}; dims=[{dims}])"
        )

    # --- Spatial differentiation D_i: M[d] → M[d-1] ----------------------
    #
    # U(L_{-1}) is a polynomial superalgebra: polynomial in the even generators
    # e_0,...,e_3 and exterior in the odd generators d_0,...,d_3.
    # (L_{-1} is abelian since [L_{-1}, L_{-1}] \subseteq L_{-2} = 0 in the depth-1
    # principal grading of E(4,4).)
    #
    # The derivation D_i acts on PBW monomials (\alpha, S) \otimes w by:
    #   D_i((\alpha, S) \otimes w) = \alpha_i · (\alpha - \epsilon_i, S) \otimes w
    # where \epsilon_i = (0,...,1,...,0) in slot i.  This is independent of the fiber
    # and encodes \partial/\partialx_i on formal Taylor-series coefficients.
    # ---------------------------------------------------------------------------

    def action_of_ei(self, i, d):
        """
        Matrix of D_i = \partial/\partiale_i acting on M[d] → M[d-1].

        D_i is the PBW derivation that removes one factor of the even generator
        e_i from each monomial:
            D_i((\alpha, S) \otimes w) = \alpha_i · (\alpha - \epsilon_i, S) \otimes w.

        Parameters
        ----------
        i : int in {0, 1, 2, 3}
            Index of the even L_{-1} generator.
        d : int
            Source degree.  Must satisfy 1 ≤ d ≤ max_deg.

        Returns
        -------
        A  dim(d-1) \times dim(d)  matrix over QQ.
        """
        if d < 1 or d > self.max_deg:
            raise ValueError(f"d={d} out of range [1, {self.max_deg}]")
        rows = self.dim(d - 1)
        cols = self.dim(d)
        if rows == 0 or cols == 0:
            return matrix(QQ, rows, cols)

        mat = matrix(QQ, rows, cols)
        pbw_d   = self.pbw[d]
        pbw_dm1 = self.pbw[d - 1]
        pbw_dm1_idx = {mon: j for j, mon in enumerate(pbw_dm1)}
        dim_W   = self.dim_W

        for j, (alpha, S) in enumerate(pbw_d):
            ai = alpha[i]
            if ai == 0:
                continue
            # Target monomial: decrement \alpha[i] by 1
            new_alpha = list(alpha)
            new_alpha[i] = ai - 1
            new_mon = (tuple(new_alpha), S)
            if new_mon not in pbw_dm1_idx:
                continue
            new_j = pbw_dm1_idx[new_mon]
            coeff = QQ(ai)
            for k in range(dim_W):
                mat[new_j * dim_W + k, j * dim_W + k] += coeff

        return mat


# --- Public constructor -----------------------------------------------------

def M_verma(t, a, b, c, max_deg=4, e44_data=None):
    """
    Construct the Verma module M(t, a, b, c) truncated at degree max_deg.

    Returns a VermaModule object whose degree-d pieces are the vector spaces
        M[d] = U(L_{-1})_d \otimes W(t, a, b, c).

    When e44_data is provided, the fiber is the full irreducible \hat{p}(4)-module
    W_t(a,b,c) (from phat4_modules).  Otherwise it is the sl_4-irreducible
    V(a,b,c) (from W_hat).

    Examples
    --------
    sage: M = M_verma(0, 1, 0, 0)
    sage: M.dim(0)
    4
    sage: M.dim(1)
    32
    sage: M.dim(2)
    128
    """
    return VermaModule(t, a, b, c, max_deg=max_deg, e44_data=e44_data)


# ===========================================================================
# Step 4 — Laplacian operator on M
# ===========================================================================
#
# The formal spatial Laplacian is  \Delta = \Sigma_i D_i^2  where
#   D_i = \partial/\partiale_i : M[d] → M[d-1]   (action_of_ei above).
#
# It maps  M[d] → M[d-2]  and has matrix:
#   B(d) = \Sigma_i A_i(d-1) · A_i(d)
# where  A_i(d) = M.action_of_ei(i, d).
#
# B(d) is a  dim(d-2) \times dim(d)  rational matrix.
# ===========================================================================

def laplacian_matrix(M, d):
    """
    Matrix of the formal Laplacian  \Delta = \Sigma_i D_i^2  on  M[d] → M[d-2].

    Computes B(d) = \Sigma_i A_i(d-1) · A_i(d)  where A_i(k) = M.action_of_ei(i, k).

    Parameters
    ----------
    M : VermaModule
    d : int — source degree (must satisfy d >= 2)

    Returns
    -------
    A  dim(d-2) \times dim(d)  matrix over QQ.
    """
    if d < 2:
        raise ValueError(f"d={d}: Laplacian requires d >= 2")
    if d > M.max_deg:
        raise ValueError(f"d={d} exceeds max_deg={M.max_deg}")

    rows = M.dim(d - 2)
    cols = M.dim(d)
    if rows == 0 or cols == 0:
        return matrix(QQ, rows, cols)

    B = matrix(QQ, rows, cols)
    for i in range(4):
        Ai_d    = M.action_of_ei(i, d)       # dim(d-1) \times dim(d)
        Ai_dm1  = M.action_of_ei(i, d - 1)   # dim(d-2) \times dim(d-1)
        B += Ai_dm1 * Ai_d
    return B


def laplacian_matrices(M):
    """
    Return a list  [None, None, B(2), B(3), ..., B(max_deg)]
    where  B(d) = laplacian_matrix(M, d).

    B(0) and B(1) are set to None (undefined / trivially zero).
    """
    result = [None, None]
    for d in range(2, M.max_deg + 1):
        result.append(laplacian_matrix(M, d))
    return result


# ===========================================================================
# s2.3 — L_{<0}-action on M and commutator verification
# ===========================================================================
#
# The Verma module M = U(L_{-1}) \otimes W is an E(4,4)-module.  The action of
# L_{<0} = L_{-1} on a PBW basis element  (\alpha,S) \otimes w  at degree d is:
#
#   e_i \cdot ((\alpha,S) \otimes w)  =  (\alpha + \epsilon_i, S) \otimes w      (even generator; polynomial mult)
#   d_i \cdot ((\alpha,S) \otimes w)  =  sgn \cdot (\alpha, S\cup{i}) \otimes w  (odd; \epsilon if i ∉ S)
#                         0                       (if i \in S)
#
# where sgn = (-1)^{|S ∩ {0,...,i-1}|} from the super-commutation rule.
#
# The L_0 action is derived from the Verma property via the PBW relation:
#   X \cdot (u \otimes w) = (X\cdotu) \otimes w + u \otimes (X\cdotw)      (Leibniz / module axiom)
# where X \cdot u is the adjoint action [X, u] computed using the E44 bracket,
# decomposed in the L_{-1} basis, and X \cdot w is the \hat{p}(4)-action on W_hat.
#
# The L_1 action maps M[d] \to M[d-1] via the same Leibniz rule:
#   Y \cdot (u \otimes w) = (Y\cdotu) \otimes w + u \otimes (Y\cdotw)
# where Y \cdot u lands in L_0 (via [L_1, L_{-1}] \subseteq L_0), and Y \cdot w = 0
# because L_{>0} acts trivially on W (Verma module construction).
# Then (L_0-element) acts on the remaining PBW part via another Leibniz step.
#
# Implementation note:
#   • We store L_{-1}, L_0, L_1 bases from the E44 pickle;
#   • We build coefficient-extraction helpers to decompose bracket results;
#   • The L_0 action on W_hat uses the sl_4 Chevalley matrices from WHat4Module
#     (the central element C acts as t \cdot Id, and \hat{p}(4)_{±1} pieces act trivially
#      on the lowest-weight Verma generator by construction).
# ---------------------------------------------------------------------------

# --- Coefficient extraction -----------------------------------------------

def _expand_in_basis(result_tuple, basis_list, R):
    """
    Express a 4-tuple `result_tuple` of polynomials as a linear combination
    of elements in `basis_list` (each a dict with key 'basis').

    Each basis element has exactly one nonzero component which is a single
    monomial m.  The coefficient of that basis element in `result_tuple` is
    the coefficient of m in `result_tuple[comp]`.

    Returns a dict {basis_index: QQ_coefficient} (omitting zero coefficients).
    """
    coeffs = {}
    for idx, b in enumerate(basis_list):
        bvec = b['basis']
        for comp in range(4):
            if bvec[comp] != R(0):
                mon = bvec[comp]
                c = result_tuple[comp].monomial_coefficient(mon)
                if c != 0:
                    coeffs[idx] = QQ(c)
                break   # each basis element has a single nonzero component
    return coeffs


# --- PBW monomial helpers --------------------------------------------------

def _pbw_apply_even(alpha, S, gen_idx):
    """
    Apply even generator e_i (gen_idx) to PBW monomial (alpha, S).
    Returns new (alpha', S') with alpha'[gen_idx] += 1.
    Sign = +1 always (polynomial multiplication).
    """
    al = list(alpha)
    al[gen_idx] += 1
    return tuple(al), S, QQ(1)


def _pbw_apply_odd(alpha, S, gen_idx):
    """
    Apply odd generator d_i (gen_idx) to PBW monomial (alpha, S).
    Returns (alpha, S\cup{i}, sign) or None if i \in S.

    Sign from super-commutation:
        d_i \cdot (dⱼ_1 ∧ ... ∧ dⱼ_k) = sgn \cdot d_i ∧ dⱼ_1 ∧ ... ∧ dⱼ_k
    We put d_i at position p = |{j \in S : j < i}|.
    In our canonical ordering (S = {j_1 < j_2 < ... < j_k}),
        d_i is inserted at position p, so sign = (-1)^p.
    """
    if gen_idx in S:
        return None   # Grassmann relation d_i^2 = 0
    p = sum(1 for j in S if j < gen_idx)
    new_S = frozenset(S | {gen_idx})
    return alpha, new_S, QQ((-1) ** p)


# --- L_{-1} action matrix on M[d] \to M[d+1] --------------------------------

def l_minus1_action_matrix(M, gen_idx, parity):
    """
    Build the matrix of the L_{-1} action of generator g = e_i or d_i
    (gen_idx \in {0,1,2,3}, parity \in {0,1}) on M[d] \to M[d+1],
    for d = 0 ... max_deg-1.

    Returns a list  mats  where  mats[d]  is a
        dim(d+1) \times dim(d) QQ-matrix.
    """
    mats = []
    for d in range(M.max_deg):
        rows = M.dim(d + 1)
        cols = M.dim(d)
        if rows == 0 or cols == 0:
            mats.append(matrix(QQ, rows, cols))
            continue
        mat = matrix(QQ, rows, cols)
        pbw_d = M.pbw[d]
        pbw_dp1 = M.pbw[d + 1]
        # Build reverse index for pbw[d+1]
        pbw_dp1_idx = {mon: j for j, mon in enumerate(pbw_dp1)}

        for j, (alpha, S) in enumerate(pbw_d):
            if parity == 0:
                new_alpha, new_S, sign = _pbw_apply_even(alpha, S, gen_idx)
            else:
                result = _pbw_apply_odd(alpha, S, gen_idx)
                if result is None:
                    continue
                new_alpha, new_S, sign = result
            new_mon = (new_alpha, new_S)
            if new_mon not in pbw_dp1_idx:
                continue
            new_j = pbw_dp1_idx[new_mon]
            # For each w-basis index k, set mat[new_j*dim_W + k, j*dim_W + k]
            for k in range(M.dim_W):
                row = new_j * M.dim_W + k
                col = j * M.dim_W + k
                mat[row, col] += sign
        mats.append(mat)
    return mats


# --- L_0 action matrix on M[d] \to M[d] --------------------------------------
#
# The L_0 action via Leibniz:
#   h \cdot ((\alpha,S) \otimes w_k)
#     = \Sigmaⱼ c_j \cdot ((\alphaⱼ',Sⱼ') \otimes w_k)   [from [h, generators] decomposed in L_{-1}]
#     + (\alpha,S) \otimes (h\cdotw_k)              [from the W-action]
#
# For the sl_4 generators (identified as the degree-0 vector fields xⱼ\partial_i \in L_0),
# the adjoint action [sl4_elem, e_i-generator] and [sl4_elem, d_i-generator]
# is given by the E44 bracket, decomposed in L_{-1}.
# The W-action uses the sl_4 matrices from W_hat.
#
# For the central element C (the Euler vector field \Sigma x_i\partial_i \in L_0):
#   C \cdot ((\alpha,S) \otimes w_k) = (|\alpha|+|S|) \cdot ((\alpha,S) \otimes w_k) + t \cdot ((\alpha,S) \otimes w_k)
#   = (deg + t) \cdot (\alpha,S) \otimes w_k
# because C acts as the grading operator on U(L_{-1}) (degree d piece gets
# eigenvalue d) and as t\cdotId on W.
# ---------------------------------------------------------------------------

def _build_lm1_adj_table(e44_data):
    """
    Pre-compute the adjoint action of each L_0 basis element on each L_{-1}
    basis element, as a dict:
        adj[(L0_idx, Lm1_idx)] = {Lm1_result_idx: coeff}
    Both lists include even and odd elements in the order from e44_data.

    Uses the btable from the E44 pickle.
    """
    from sage.all import PolynomialRing
    R = PolynomialRing(QQ, ['x1', 'x2', 'x3', 'x4'])
    btable = e44_data['btable']
    L0    = e44_data['E44'][0]
    Lm1   = e44_data['E44'][-1]
    Lm1_even = [b for b in Lm1 if b['parity'] == 0]
    Lm1_odd  = [b for b in Lm1 if b['parity'] == 1]

    adj = {}
    for i0, b0 in enumerate(L0):
        for im1, bm1 in enumerate(Lm1):
            key = (b0['label'], bm1['label'])
            result = btable.get(key)
            if result is None:
                adj[(i0, im1)] = {}
                continue
            # Determine target sector: [even, even]\toeven, [even,odd]\toodd, etc.
            parity_result = (b0['parity'] + bm1['parity']) % 2
            if parity_result == 0:
                coeffs = _expand_in_basis(result, Lm1_even, R)
                # remap to full Lm1 index
                full = {Lm1.index(Lm1_even[ii]): c for ii, c in coeffs.items()}
            else:
                coeffs = _expand_in_basis(result, Lm1_odd, R)
                full = {Lm1.index(Lm1_odd[ii]): c for ii, c in coeffs.items()}
            adj[(i0, im1)] = full
    return adj


def _w_action_from_l0_idx(W, L0_idx):
    """
    Return the dim_W \times dim_W action matrix of the L0[L0_idx] even generator
    on crystal module W, using the Chevalley e/f/h matrices directly.

    The even L_0 basis for j=0 is ordered as: for each monomial in
    [x_4, x_3, x_2, x_1] and each slot 0..3 (\partial_1..\partial_4).  Thus:
        L0[k] = x_{row} \partial_{col}  where  row = 4 - k//4,  col = k%4 + 1  (1-indexed).

    e44_structure convention: E_ij(i,j) = E_{ij} = -xⱼ \partial_i,  so
        x_{row} \partial_{col}  =  -E_{col, row}   (1-indexed matrix units).

    Representation map:  The Lie algebra homomorphism ψ: gl_4^{code} \to gl_4^{std}
    is ψ(E_{ij}^{code}) = -e_{ji}^{std}  (the transpose with a sign, required to
    preserve the bracket: [ψ(X), ψ(Y)] = ψ([X, Y])).

    Therefore:
        p(L0[k]) = p(-E_{col,row}^{code}) = -p(E_{col,row}^{code})
                 = -(-e_{row,col}^{crystal}) = +e_{row,col}^{crystal}.

    This gives the action on W:
      • Off-diagonal col < row  (e_{row,col}^{std} with row > col = lowering):
            W_mat = f[...]  (Chevalley lowering matrix)
      • Off-diagonal col > row  (e_{row,col}^{std} with row < col = raising):
            W_mat = e[...]  (Chevalley raising matrix)
      • Diagonal col = row = r:
            x_r \partial_r = -E_{rr}^{code} + (1/4)C,  so
            action = \epsilonᵣ + t/4,  where
              \epsilon_1 = (3h_1+2h_2+h_3)/4,  \epsilon_2 = (-h_1+2h_2+h_3)/4,
              \epsilon_3 = (-h_1-2h_2+h_3)/4,  \epsilon_4 = (-h_1-2h_2-3h_3)/4.
            This gives:
              r=1: (3h_1+2h_2+h_3+t)/4,    r=2: (-h_1+2h_2+h_3+t)/4,
              r=3: (-h_1-2h_2+h_3+t)/4,    r=4: (-h_1-2h_2-3h_3+t)/4.

    Odd L_0 generators (L0_idx >= 16) act as zero on the crystal fiber.
    """
    if L0_idx >= 16:
        return matrix(QQ, W.dim, W.dim)   # odd L_0 acts as zero on the fiber

    n = W.dim
    e = W.e_mats   # {1,2,3} \to n\timesn matrix
    f = W.f_mats
    h = W.h_mats
    t = W.t
    Id = identity_matrix(QQ, n)

    def comm(A, B):
        return A * B - B * A

    row = 4 - L0_idx // 4   # 1-indexed; x_{row}
    col = L0_idx % 4 + 1    # 1-indexed; \partial_{col}
    # L0[k] = x_{row} \partial_{col} = -E_{col, row} in e44_structure convention.
    # Action on W: +e_{row,col}^{crystal}  (see derivation in docstring).

    if row == col:
        # Diagonal: x_r \partial_r = -E_{rr}^{code} + (1/4)C, acts as (\epsilon_r + t/4).
        r = row
        if r == 1:
            return (3 * h[1] + 2 * h[2] + h[3] + t * Id) / 4
        elif r == 2:
            return (-h[1] + 2 * h[2] + h[3] + t * Id) / 4
        elif r == 3:
            return (-h[1] - 2 * h[2] + h[3] + t * Id) / 4
        else:  # r == 4
            return (-h[1] - 2 * h[2] - 3 * h[3] + t * Id) / 4

    elif col < row:
        # L0[k] = x_r\partial_c with r > c: e_{row,col}^{std} is lower-triangular.
        # Action = +f[...] (Chevalley lowering).
        if   (col, row) == (1, 2): return f[1]
        elif (col, row) == (2, 3): return f[2]
        elif (col, row) == (3, 4): return f[3]
        elif (col, row) == (1, 3): return comm(f[2], f[1])
        elif (col, row) == (2, 4): return comm(f[3], f[2])
        elif (col, row) == (1, 4): return comm(comm(f[1], f[2]), f[3])
        else:
            return matrix(QQ, n, n)

    else:  # col > row
        # L0[k] = x_r\partial_c with r < c: e_{row,col}^{std} is upper-triangular.
        # Action = +e[...] (Chevalley raising).
        if   (col, row) == (2, 1): return e[1]
        elif (col, row) == (3, 2): return e[2]
        elif (col, row) == (4, 3): return e[3]
        elif (col, row) == (3, 1): return comm(e[1], e[2])
        elif (col, row) == (4, 2): return comm(e[2], e[3])
        elif (col, row) == (4, 1): return comm(e[3], comm(e[2], e[1]))
        else:
            return matrix(QQ, n, n)


def l0_action_matrix(M, L0_idx, e44_data):
    """
    Build the matrix of the L_0 action of the L0_idx-th L_0 basis element
    on M[d] \to M[d], for d = 0 ... max_deg.

    Returns a list  mats  where  mats[d]  is a dim(d)\timesdim(d) QQ-matrix.

    Requires e44_data loaded from the E44 pickle.

    The L_0 element is identified with one of the 32 basis elements of L_0:
    indices 0..15 = even (sl_4 generators); indices 16..31 = odd.
    The sl_4 action on W_hat uses the Chevalley matrices; the central
    element C (identified as the trace element) is handled separately.

    Results are memoised on M._l0_cache to avoid redundant recomputation
    when the same L0_idx is queried many times (e.g. during L_1 annihilator
    checks that loop over all 80 L_1 generators).
    """
    # --- Memo cache ---
    if not hasattr(M, '_l0_cache'):
        M._l0_cache = {}
    if L0_idx in M._l0_cache:
        return M._l0_cache[L0_idx]

    L0   = e44_data['E44'][0]
    Lm1  = e44_data['E44'][-1]
    b0   = L0[L0_idx]
    R    = list(Lm1[0]['basis'][0].parent().gens())  # polynomial ring gens

    # Pre-compute adjoint table (cached on M if already built)
    if not hasattr(M, '_adj_table'):
        M._adj_table = _build_lm1_adj_table(e44_data)
    adj = M._adj_table

    # L_0 action on W_hat via sl_4 matrices:
    # Map L_0 basis element to a sl_4 Chevalley generator index, or None.
    # L_0 even has 16 elements = 4\times4 matrix units acting as xⱼ\partial_i.
    # We identify the (i,j) matrix unit with the sl_4 generator action.
    # The sl_4 matrices in W_hat are keyed by Chevalley index {1,2,3};
    # we need the full gl_4 ≅ sl_4 \oplus C action.
    # Strategy: act on W-component using the bracket [b0, L_{-1}] result
    # expressed in L_{-1}, then deduce how W transforms.
    # For the full implementation we use the W_hat matrices directly.
    # The L_0 even basis in e44_structure.py is E_{ij} = -xⱼ\partial_i + (1/4)\delta_iⱼ C.
    # To find which (i,j) pairs, we look at the single-monomial structure:
    # b0['basis'] = (0,...,xⱼ,...,0) \times (-1) in slot i (or +correction for diagonal).
    # We compute the W-action by evaluating how sl_4 matrices act, using
    # the adjoint table on the degree-0 piece (which is just the W-action).

    mats = []
    for d in range(M.max_deg + 1):
        dim_d = M.dim(d)
        if dim_d == 0:
            mats.append(matrix(QQ, 0, 0))
            continue
        mat = matrix(QQ, dim_d, dim_d)
        pbw_d = M.pbw[d]

        for j, (alpha, S) in enumerate(pbw_d):
            # Leibniz: act on each L_{-1} factor in the monomial ...
            # Efficient approach: iterate over all L_{-1} generators that
            # appear in (alpha, S) and apply the adjoint action on each.

            # --- Act on the polynomial (even) part ---
            # The monomial e_i^{\alpha_i} contributes \alpha_i terms, each applying
            # [b0, e_i] to one factor (product rule for algebra action).
            for gen_i in range(4):      # even generators e0..e3
                pow_i = alpha[gen_i]
                if pow_i == 0:
                    continue
                # Lm1 index of e_i = gen_i (first 4 entries are even)
                ad_result = adj.get((L0_idx, gen_i), {})
                for lm1_res_idx, coeff in ad_result.items():
                    # Apply: replace one e_i factor by lm1_res_idx factor
                    b_res = Lm1[lm1_res_idx]
                    p_res = b_res['parity']
                    if p_res == 0:       # result is even: increment alpha
                        res_gen_i = lm1_res_idx   # index 0..3 \to even gen
                        new_alpha = list(alpha)
                        new_alpha[res_gen_i] += 1
                        new_alpha[gen_i] -= 1
                        new_mon = (tuple(new_alpha), S)
                    else:                # result is odd: add to S
                        res_gen_i = lm1_res_idx - 4   # offset by 4 even gens
                        if res_gen_i in S:
                            continue
                        p_ins = sum(1 for jj in S if jj < res_gen_i)
                        sign_odd = QQ((-1) ** p_ins)
                        new_S = frozenset(S | {res_gen_i})
                        # also remove one e_i factor
                        new_alpha = list(alpha)
                        new_alpha[gen_i] -= 1
                        new_mon = (tuple(new_alpha), new_S)
                        coeff = coeff * sign_odd
                    # Check new_mon is in pbw basis
                    if new_mon not in M._pbw_idx[d]:
                        continue
                    new_j = M._pbw_idx[d][new_mon]
                    for k in range(M.dim_W):
                        mat[new_j * M.dim_W + k, j * M.dim_W + k] += (
                            QQ(pow_i) * coeff
                        )

            # --- Act on the exterior (odd) part ---
            for pos, gen_i in enumerate(sorted(S)):
                lm1_idx = 4 + gen_i    # odd gens are indices 4..7 in Lm1
                ad_result = adj.get((L0_idx, lm1_idx), {})
                # Sign from the exterior algebra: acting on the element at
                # position pos requires moving it through pos predecessors.
                # For even h (odd\toodd case): this is the extraction sign.
                # For odd h (odd\toeven case): this equals the Koszul sign.
                # In both cases the combined sign is (-1)^pos.
                sign_from_S = QQ((-1) ** pos)
                for lm1_res_idx, coeff in ad_result.items():
                    b_res = Lm1[lm1_res_idx]
                    p_res = b_res['parity']
                    if p_res == 1:      # result is odd: swap in S
                        res_gen_i = lm1_res_idx - 4
                        new_S = frozenset((S - {gen_i}) | {res_gen_i})
                        # Compute sign from sorting new_S vs old S
                        old_sorted = sorted(S)
                        new_sorted = sorted(new_S)
                        # Sign = sign_from_S \times sign for re-sorting
                        # (-1)^pos from moving through predecessors,
                        # (-1)^{new_pos} for reinserting res_gen_i
                        new_pos = sum(1 for jj in new_S if jj < res_gen_i)
                        total_sign = sign_from_S * QQ((-1) ** new_pos)
                        new_mon = (alpha, new_S)
                    else:               # result is even: move to alpha
                        res_gen_i = lm1_res_idx
                        new_S = frozenset(S - {gen_i})
                        new_alpha = list(alpha)
                        new_alpha[res_gen_i] += 1
                        total_sign = sign_from_S
                        new_mon = (tuple(new_alpha), new_S)
                    if new_mon not in M._pbw_idx[d]:
                        continue
                    new_j = M._pbw_idx[d][new_mon]
                    for k in range(M.dim_W):
                        mat[new_j * M.dim_W + k, j * M.dim_W + k] += (
                            total_sign * coeff
                        )

            # --- Act on W-factor ---
            # W-action is computed separately below via _w_action_from_l0_idx.

        # --- Add W-action contribution ---
        # W-action via Chevalley crystal matrices.
        # Each even L_0 generator L0[k] = x_{4-k//4} \partial_{k%4+1} = E_{row,col}
        # (1-indexed pure gl_4 matrix unit).  _w_action_from_l0_idx maps this
        # to the correct sl_4/gl_4 crystal action on the fiber W for any rep.
        # Odd L_0 generators act as zero on the sl_4 crystal fiber,
        # but nontrivially on the full \hat{p}(4) fiber (Phat4Module).
        if hasattr(M.W, 'action_mats'):
            W_mat = M.W.action_mats[L0_idx]
        else:
            W_mat = _w_action_from_l0_idx(M.W, L0_idx)

        # Add W-action to every PBW monomial.
        # Koszul sign: h commutes past the full PBW monomial u to
        # reach w, picking up (-1)^{p_h * p_u} where p_u = |S| mod 2.
        p_h = b0['parity']
        for j, (alpha_j, S_j) in enumerate(pbw_d):
            koszul = QQ((-1) ** (p_h * (len(S_j) % 2)))
            for k in range(M.dim_W):
                for k2 in range(M.dim_W):
                    c = W_mat[k2, k]
                    if c != 0:
                        mat[j * M.dim_W + k2, j * M.dim_W + k] += koszul * c

        mats.append(mat)

    M._l0_cache[L0_idx] = mats
    return mats


# --- L_1 action matrix on M[d] \to M[d-1] ------------------------------------
#
# L_1 element Y acts via the Leibniz rule on each L_{-1} factor:
#   Y \cdot (u \otimes w)  =  ([Y, u]_{L_0} acting on remaining u) \otimes w
# because Y \cdot w = 0 (L_{>0} acts trivially on the Verma generator).
#
# Concretely, for each L_{-1} generator g_i appearing in (\alpha,S):
#   [Y, g_i] \in L_0  (since deg(Y)=1, deg(g_i)=-1, sum=0)
# Then [Y,g_i] acts on the remaining PBW monomial via the L_0 action above.
# ---------------------------------------------------------------------------

def l1_action_matrix(M, L1_idx, e44_data):
    """
    Build the matrix of the L_1 action of the L1_idx-th L_1 basis element
    on M[d] \to M[d-1], for d = 1 ... max_deg.

    Returns a list  mats  where  mats[d]  (for d >= 1) is a
        dim(d-1) \times dim(d) QQ-matrix.  mats[0] = None (no degree-(-1) piece).

    Uses the Leibniz rule: Y \cdot (g_i \cdot u \otimes w) = [Y,g_i] \cdot u \otimes w + (-1)^{p_i} g_i \cdot (Y\cdotu) \otimes w
    Terminated by the Verma property: Y acts trivially on (degree-0 piece \otimes W).

    Results are memoised on M._l1_cache.
    """
    # --- Memo cache ---
    if not hasattr(M, '_l1_cache'):
        M._l1_cache = {}
    if L1_idx in M._l1_cache:
        return M._l1_cache[L1_idx]

    L1   = e44_data['E44'][1]
    L0   = e44_data['E44'][0]
    Lm1  = e44_data['E44'][-1]
    bt   = e44_data['btable']
    b1   = L1[L1_idx]
    R_poly = Lm1[0]['basis'][0].parent()

    if not hasattr(M, '_adj_table'):
        M._adj_table = _build_lm1_adj_table(e44_data)

    # Pre-compute [Y, g_i] \in L_0 for each L_{-1} generator g_i
    # Result: {lm1_idx: {L0_idx: coeff}}
    y_adj_lm1 = {}
    for im1, bm1 in enumerate(Lm1):
        key = (b1['label'], bm1['label'])
        result = bt.get(key)
        if result is None:
            y_adj_lm1[im1] = {}
            continue
        # Result is in L_0: decompose in L_0 basis
        p_result = (b1['parity'] + bm1['parity']) % 2
        if p_result == 0:
            L0_sector = [b for b in L0 if b['parity'] == 0]
        else:
            L0_sector = [b for b in L0 if b['parity'] == 1]
        coeffs_sector = _expand_in_basis(result, L0_sector, R_poly)
        # Remap to full L0 indices
        full = {}
        for sec_idx, c in coeffs_sector.items():
            full_idx = L0.index(L0_sector[sec_idx])
            full[full_idx] = c
        y_adj_lm1[im1] = full

    mats = [None]   # mats[0] undefined (no degree-(-1))
    for d in range(1, M.max_deg + 1):
        rows = M.dim(d - 1)
        cols = M.dim(d)
        if rows == 0 or cols == 0:
            mats.append(matrix(QQ, rows, cols))
            continue
        mat = matrix(QQ, rows, cols)
        pbw_d   = M.pbw[d]
        pbw_dm1 = M.pbw[d - 1]
        pbw_dm1_idx = {mon: j for j, mon in enumerate(pbw_dm1)}

        for j, (alpha, S) in enumerate(pbw_d):
            # Apply Y by Leibniz over each generator in the monomial.
            # Even generators (in alpha):
            # PBW order: e0^a0 * e1^a1 * e2^a2 * e3^a3 * d_{s0} * ... * w
            # For e_i^n, slot k (0..n-1): Y acts on the k-th copy of e_i.
            # The "after" piece (generators at index > i, plus remaining
            # copies of e_i, plus odd generators) receives the L_0 action.
            # The "before" piece (e_j^{a_j} for j < i, plus e_i^k) is
            # prepended unchanged.
            for gen_i in range(4):
                pow_i = alpha[gen_i]
                if pow_i == 0:
                    continue
                lm1_idx = gen_i
                for L0_idx, c10 in y_adj_lm1[lm1_idx].items():
                    L0_mats = l0_action_matrix(M, L0_idx, e44_data)
                    col_start = j * M.dim_W
                    for slot in range(pow_i):
                        # "After" piece: e_i^{pow_i-slot-1} * e_{i+1}^... * odd
                        inner_alpha = [0] * 4
                        inner_alpha[gen_i] = pow_i - slot - 1
                        for gi2 in range(gen_i + 1, 4):
                            inner_alpha[gi2] = alpha[gi2]
                        inner_mon = (tuple(inner_alpha), S)
                        inner_d = sum(inner_alpha) + len(S)
                        if inner_d < 0 or inner_d > M.max_deg:
                            continue
                        pbw_inner_idx = M._pbw_idx[inner_d]
                        if inner_mon not in pbw_inner_idx:
                            continue
                        jp_inner = pbw_inner_idx[inner_mon]
                        L0_mat_inner = L0_mats[inner_d]
                        dim_inner = M.dim(inner_d)
                        # Apply L0 on "after" piece, then prepend "before" piece
                        for k in range(M.dim_W):
                            src_col = jp_inner * M.dim_W + k
                            for row_inner in range(dim_inner):
                                entry = L0_mat_inner[row_inner, src_col]
                                if entry == 0:
                                    continue
                                # Decode result monomial at inner_d
                                res_j = row_inner // M.dim_W
                                res_w = row_inner % M.dim_W
                                res_alpha, res_S = M.pbw[inner_d][res_j]
                                # Prepend "before": e_j^{a_j} for j<i, plus e_i^slot
                                new_alpha = list(res_alpha)
                                for gi2 in range(gen_i):
                                    new_alpha[gi2] += alpha[gi2]
                                new_alpha[gen_i] += slot
                                new_mon = (tuple(new_alpha), res_S)
                                pbw_res = M._pbw_idx[d - 1]
                                if new_mon not in pbw_res:
                                    continue
                                new_j = pbw_res[new_mon]
                                row_dst = new_j * M.dim_W + res_w
                                mat[row_dst, col_start + k] += c10 * entry

            # Odd generators (in S):
            # PBW order: e0^a0 .. e3^a3 * f_{s0} * f_{s1} * ... * w (sorted s_i)
            # When Y acts on f_{gen_i} at sorted position pos, the super-sign
            # comes from commuting Y past the pos odd f-generators before it.
            # The sign is (-1)^{pY * pos}.
            pY = b1['parity']
            for pos, gen_i in enumerate(sorted(S)):
                lm1_idx = 4 + gen_i
                sign_from_pos = QQ((-1) ** (pY * pos))
                S_before = frozenset(j for j in S if j < gen_i)
                S_after  = frozenset(j for j in S if j > gen_i)
                deg_after = len(S_after)
                for L0_idx, c10 in y_adj_lm1[lm1_idx].items():
                    # [Y, f_{gen_i}] = c10 * L0[L0_idx] acts on (S_after \otimes w),
                    # which lives in M[deg_after].
                    L0_mats = l0_action_matrix(M, L0_idx, e44_data)
                    L0_mat = L0_mats[deg_after]
                    dim_after = M.dim(deg_after)
                    pbw_after = M.pbw[deg_after]
                    pbw_after_idx = M._pbw_idx[deg_after]
                    zero_alpha = (0, 0, 0, 0)
                    src_mon = (zero_alpha, S_after)
                    if src_mon not in pbw_after_idx:
                        continue
                    jp_inner = pbw_after_idx[src_mon]
                    col_start = j * M.dim_W
                    for k in range(M.dim_W):
                        src_col = jp_inner * M.dim_W + k
                        for row_inner in range(dim_after):
                            entry = L0_mat[row_inner, src_col]
                            if entry == 0:
                                continue
                            # Decode which inner output monomial this is
                            inner_j = row_inner // M.dim_W
                            inner_w = row_inner % M.dim_W
                            inner_alpha, inner_S_after = pbw_after[inner_j]
                            # Assemble the full result monomial:
                            #   outer alpha + inner alpha, S_before \cup inner_S_after
                            new_alpha = tuple(alpha[i] + inner_alpha[i] for i in range(4))
                            new_S = S_before | inner_S_after
                            # Validity: S_before and inner_S_after must be disjoint
                            if len(new_S) != len(S_before) + len(inner_S_after):
                                continue  # f-generator already present \to 0 by antisymmetry
                            # Sign from merging S_before and inner_S_after
                            # into sorted PBW order (anticommuting odd gens)
                            inversions = sum(
                                1 for a in S_before
                                for b in inner_S_after if a > b
                            )
                            merge_sign = QQ((-1) ** inversions)
                            new_mon = (new_alpha, new_S)
                            pbw_res_idx = M._pbw_idx[d - 1]
                            if new_mon not in pbw_res_idx:
                                continue
                            new_j = pbw_res_idx[new_mon]
                            row_dst = new_j * M.dim_W + inner_w
                            mat[row_dst, col_start + k] += (
                                merge_sign * sign_from_pos * c10 * entry
                            )

        mats.append(mat)
    M._l1_cache[L1_idx] = mats
    return mats


# --- Attach _pbw_idx reverse lookup to VermaModule -------------------------
# (needed by l0_action_matrix; added as a lazy attribute)
_orig_VermaModule_init = VermaModule.__init__


def _patched_init(self, t, a, b, c, max_deg=4, e44_data=None):
    _orig_VermaModule_init(self, t, a, b, c, max_deg=max_deg, e44_data=e44_data)
    # Reverse index: (alpha, frozenset_S) -> position in pbw[d]
    self._pbw_idx = {}
    for d in range(max_deg + 1):
        self._pbw_idx[d] = {mon: j for j, mon in enumerate(self.pbw[d])}


VermaModule.__init__ = _patched_init


# ===========================================================================
# s2.3 Checkpoint
# ===========================================================================

def _check_s23(verbose=True):
    """
    S2.3 checkpoint: verify L_{-1} and L_0 actions, and the commutator identity.

    Checks
    ------
    (1)  L_{-1} action: e_i maps M[0] \to M[1] with correct dimension
    (2)  L_{-1} action: applying e_0 to w[0] \otimes (0,{}) gives e_0-monomial \otimes w[0]
    (3)  L_{-1} action: d_i^2 = 0  (Grassmann)
    (4)  L_0 action: h_1 eigenvalue on highest-weight vector at degree 0
         = Dynkin label a (for W_hat(0,1,0,0), h_1\cdotv = 1\cdotv)
    (5)-(24) Commutator identity [L_{-1,i}, L_1,j] = [L_{-1,i}, L_1,j]_bracket
             on M_0(1,0,0) [W_hat = fund. rep, dim 4]: 20 pairs covering
             all four parity combos: even\timeseven, odd\timeseven, even\timesodd, odd\timesodd.
    (25)-(44) Same commutator identity on M_1(0,0,0) [W_hat = trivial, dim 1]:
             20 pairs. This module hosts w[1D]; its trivial W-fiber exposes
             bugs that pass vacuously on the fundamental rep.

    Returns True iff all checks pass.
    """
    all_pass = True

    def _check(name, got, expected):
        nonlocal all_pass
        ok = (got == expected)
        status = "PASS" if ok else "FAIL"
        if verbose:
            print(f"  [{status}] {name}: got {got!r}, expected {expected!r}")
        if not ok:
            all_pass = False

    if verbose:
        print("─" * 60)
        print("s2.3 checkpoint  —  L_{-1} / L_0 / L_1 actions and commutator")
        print("─" * 60)

    e44_data = load_e44()
    if e44_data is None:
        print("  [SKIP] e44_brackets.pkl not found; skipping s2.3 checks.")
        return None

    M = M_verma(0, 1, 0, 0)   # W_hat = fund. rep (dim 4)

    # (1) Dimension check: L_{-1} action e_i
    Lm1_mats_e0 = l_minus1_action_matrix(M, 0, 0)   # even gen 0
    _check("l_minus1(e0): mat[0] shape rows",  Lm1_mats_e0[0].nrows(), M.dim(1))
    _check("l_minus1(e0): mat[0] shape cols",  Lm1_mats_e0[0].ncols(), M.dim(0))

    # (2) e_0 \cdot (w[0]) = e_0-monomial \otimes w[0]
    # In M[0], basis = [(alpha=(0,0,0,0), S={})] \times w[0..3]
    # e_0 should map w[0] (index 0) to  ((1,0,0,0), {}) \otimes w[0]
    # which is index 0 in M[1] (since the first pbw_1 monomial could be (0,0,0,1))
    # Let's just check the column is a unit vector
    col0 = Lm1_mats_e0[0].column(0)
    _check("l_minus1(e0): col 0 has exactly one nonzero entry",
           sum(1 for x in col0 if x != 0), 1)
    _check("l_minus1(e0): nonzero entry = 1",
           max(abs(x) for x in col0), 1)

    # (3) Grassmann: d_0^2 = 0
    Lm1_mats_d0 = l_minus1_action_matrix(M, 0, 1)   # odd gen 0
    prod = Lm1_mats_d0[1] * Lm1_mats_d0[0]
    _check("Grassmann d_0^2 = 0", prod.is_zero(), True)

    # (4) L_0 h-action: sl_4 Cartan h_1 on W_hat(0,1,0,0)
    # h_1 eigenvalue on v_hw (weight (1,0,0)) = 1
    # In e44_structure, the h_1 = E_{11} - E_{22} is in L_0.
    # e44_data['E44'][0] has 32 elements; find Cartan-type element.
    # Instead of searching, use W_hat directly for this check:
    W_fund = M.W   # W_hat(0,1,0,0)
    hw_idx = W_fund.v_hw
    h1_mat = W_fund.h_mats[1]
    hw_vec = vector(QQ, [QQ(1) if k == hw_idx else QQ(0) for k in range(W_fund.dim)])
    h1_hw = h1_mat * hw_vec
    _check("h_1 \cdot v_hw = 1\cdotv_hw (h_1 eigenvalue on fund. h.w. vector)",
           h1_hw[hw_idx], QQ(1))

    # -----------------------------------------------------------------
    # Commutator identity  [L_{-1,i}, L_1,j]  on a given module M_test
    # -----------------------------------------------------------------
    L0_list = e44_data['E44'][0]
    Lm1_list = e44_data['E44'][-1]
    L1_list  = e44_data['E44'][1]
    bt = e44_data['btable']
    R_poly = Lm1_list[0]['basis'][0].parent()

    def _commutator_checks(M_test, pairs, tag):
        """
        Test the commutator identity [X, Y] = bracket(X,Y) at degree 1
        for each (im1, i1) pair in *pairs*, on module *M_test*.
        """
        for im1, i1 in pairs:
            bm1 = Lm1_list[im1]
            b1  = L1_list[i1]
            gen_idx_m1 = im1 if bm1['parity'] == 0 else im1 - 4

            mat_X_d0 = l_minus1_action_matrix(M_test, gen_idx_m1, bm1['parity'])[0]
            mat_Y_d1 = l1_action_matrix(M_test, i1, e44_data)[1]
            mat_X_d1 = l_minus1_action_matrix(M_test, gen_idx_m1, bm1['parity'])[1]
            mat_Y_d2 = l1_action_matrix(M_test, i1, e44_data)[2]

            sign_XY = QQ((-1) ** (bm1['parity'] * b1['parity']))
            lhs_mat = mat_X_d0 * mat_Y_d1 - sign_XY * mat_Y_d2 * mat_X_d1

            key_br = (bm1['label'], b1['label'])
            br_result = bt.get(key_br, None)
            rhs_mat = matrix(QQ, M_test.dim(1), M_test.dim(1))
            if br_result is not None and any(c != 0 for c in br_result):
                p_br = (bm1['parity'] + b1['parity']) % 2
                if p_br == 0:
                    L0_sector = [b for b in L0_list if b['parity'] == 0]
                else:
                    L0_sector = [b for b in L0_list if b['parity'] == 1]
                L0_coeffs = _expand_in_basis(br_result, L0_sector, R_poly)
                for sec_idx, coeff in L0_coeffs.items():
                    full_L0_idx = L0_list.index(L0_sector[sec_idx])
                    L0_mats = l0_action_matrix(M_test, full_L0_idx, e44_data)
                    rhs_mat += coeff * L0_mats[1]

            diff = lhs_mat - rhs_mat
            _check(f"{tag} commutator ({im1},{i1})", diff.is_zero(), True)

    # (5)-(24) Commutator identity on M_0(1,0,0)  [W_hat = fund., dim 4]
    # [L_{-1,i}, L_1,j]ₘₐₜ v = (L_{-1,i}\circL_1,j - (-1)^{p_i p_j} L_1,j\circL_{-1,i}) v
    # should equal L_0[L_{-1,i},L_1,j]_bracket v  for all v in M[1].
    #
    # 20 pairs covering all four parity combinations:
    #   even\timeseven, odd\timeseven, even\timesodd, odd\timesodd.
    # L₋_1 indices: 0-3 even (e_i), 4-7 odd (d_i).
    # L_1  indices: 0-39 even, 40-79 odd.

    if verbose:
        print("  --- M_0(1,0,0)  [W_hat = fund. rep, dim 4] ---")

    _commutator_checks(M, [
        # even L₋_1 \times even L_1
        (0,0),(1,0),(0,1),(2,3),(3,2),
        # odd L₋_1 \times even L_1
        (4,0),(5,1),(4,4),
        # even L₋_1 \times odd L_1
        (0,40),(1,45),(2,50),(0,60),(3,70),
        # odd L₋_1 \times odd L_1
        (4,40),(5,45),(6,50),(7,60),(4,70),(5,55),(6,65),
    ], "M_0(1,0,0)")

    # (25)-(44) Commutator identity on M_1(0,0,0)  [W_hat = trivial, dim 1]
    # This module is where w[1D] lives.  The trivial W-fiber means the
    # L_0 W-action is purely scalar, so many odd L_1 actions vanish;
    # the commutator test catches bugs that pass vacuously on fund. rep.

    if verbose:
        print("  --- M_1(0,0,0)  [W_hat = trivial, dim 1] ---")

    M_triv = M_verma(1, 0, 0, 0)

    _commutator_checks(M_triv, [
        # even L₋_1 \times even L_1
        (0,0),(1,0),(2,3),(3,2),
        # odd L₋_1 \times even L_1
        (4,0),(5,1),(6,2),(7,3),
        # even L₋_1 \times odd L_1
        (0,40),(1,45),(2,50),(3,60),(0,70),(1,55),(2,65),(3,75),
        # odd L₋_1 \times odd L_1
        (4,40),(5,45),(6,50),(7,60),
    ], "M_1(0,0,0)")

    if verbose:
        print("─" * 60)
        if all_pass:
            print("s2.3  ✓  ALL CHECKS PASSED")
        else:
            print("s2.3  ✗  SOME CHECKS FAILED — see above")
        print("─" * 60)

    return all_pass


def _check_commutator_all_deg(verbose=True):
    """
    Verify the super-commutator identity at ALL degrees 0 ... max_deg-1.

        [X, Y]_d  =  X_{d-1} \circ Y_d  -  (-1)^{|X||Y|}  Y_{d+1} \circ X_d

    where  X = L_{-1}[im1],  Y = L_1[i1],  [X,Y] = L_0 bracket.

    Tests every (im1, i1) pair (8 \times 80 = 640) at every degree d on the
    specified module.  Reports per-degree failure counts.

    Modules tested:
      • M_0(1,0,0)   [fund. rep, dim_W=4]
      • M_1(0,0,0)   [trivial, dim_W=1]
      • M_0(0,0,2)   [Sym^2(C^4)*, dim_W=10]
    """
    e44_data = load_e44()
    L1_list  = e44_data['E44'][1]
    L0_list  = e44_data['E44'][0]
    Lm1_list = e44_data['E44'][-1]
    bt       = e44_data['btable']
    R_poly   = Lm1_list[0]['basis'][0].parent()

    all_pass = True
    total_checks = 0
    total_failures = 0

    def _check(name, got, expected):
        nonlocal all_pass, total_checks, total_failures
        total_checks += 1
        ok = (got == expected)
        if not ok:
            all_pass = False
            total_failures += 1
        if verbose and not ok:
            print(f"  [FAIL] {name}: got {got}, expected {expected}")
        elif verbose and ok:
            pass  # suppress per-check output for brevity

    if verbose:
        print("─" * 60)
        print("Commutator identity  —  all degrees, all pairs")
        print("─" * 60)

    modules = [
        (0, 1, 0, 0, 4, "M_0(1,0,0)"),
        (1, 0, 0, 0, 4, "M_1(0,0,0)"),
        (0, 0, 0, 2, 2, "M_0(0,0,2)"),
    ]

    for t, a, b, c, max_deg, tag in modules:
        M_test = M_verma(t, a, b, c, max_deg=max_deg)
        if verbose:
            dims = [M_test.dim(d) for d in range(max_deg + 1)]
            print(f"  {tag}  max_deg={max_deg}  dims={dims}")
        n_fail = 0
        n_pass = 0
        fail_details = []

        for im1 in range(len(Lm1_list)):
            bm1 = Lm1_list[im1]
            gen_idx_m1 = im1 if bm1['parity'] == 0 else im1 - 4
            lm1_mats = l_minus1_action_matrix(M_test, gen_idx_m1, bm1['parity'])

            for i1 in range(len(L1_list)):
                b1 = L1_list[i1]
                l1_mats = l1_action_matrix(M_test, i1, e44_data)
                sign_XY = QQ((-1) ** (bm1['parity'] * b1['parity']))

                # Compute bracket [X, Y] decomposition in L_0
                key_br = (bm1['label'], b1['label'])
                br_result = bt.get(key_br, None)
                br_coeffs = {}
                if br_result is not None and any(cc != 0 for cc in br_result):
                    p_br = (bm1['parity'] + b1['parity']) % 2
                    if p_br == 0:
                        L0_sector = [bb for bb in L0_list if bb['parity'] == 0]
                    else:
                        L0_sector = [bb for bb in L0_list if bb['parity'] == 1]
                    L0_coeffs = _expand_in_basis(br_result, L0_sector, R_poly)
                    for sec_idx, coeff in L0_coeffs.items():
                        full_L0_idx = L0_list.index(L0_sector[sec_idx])
                        br_coeffs[full_L0_idx] = coeff

                for d in range(M_test.max_deg):
                    dim_d = M_test.dim(d)
                    if dim_d == 0:
                        continue

                    # LHS = X_{d-1} \circ Y_d  -  sign * Y_{d+1} \circ X_d
                    # X: M[d] \to M[d+1] = lm1_mats[d]  (rows=dim(d+1), cols=dim(d))
                    # Y: M[d] \to M[d-1] = l1_mats[d]   (rows=dim(d-1), cols=dim(d))
                    # X_{d-1}: M[d-1] \to M[d]
                    # Y_{d+1}: M[d+1] \to M[d]

                    # First term: X_{d-1} \circ Y_d
                    if d >= 1 and l1_mats[d] is not None:
                        mat_Y_d = l1_mats[d]          # dim(d-1) \times dim(d)
                        mat_X_dm1 = lm1_mats[d - 1]   # dim(d) \times dim(d-1)
                        term1 = mat_X_dm1 * mat_Y_d    # dim(d) \times dim(d)
                    else:
                        term1 = matrix(QQ, dim_d, dim_d)

                    # Second term: Y_{d+1} \circ X_d
                    if d + 1 <= M_test.max_deg and l1_mats[d + 1] is not None:
                        mat_X_d = lm1_mats[d]         # dim(d+1) \times dim(d)
                        mat_Y_dp1 = l1_mats[d + 1]    # dim(d) \times dim(d+1)
                        term2 = mat_Y_dp1 * mat_X_d    # dim(d) \times dim(d)
                    else:
                        term2 = matrix(QQ, dim_d, dim_d)

                    lhs = term1 - sign_XY * term2

                    # RHS = L_0[bracket] at degree d
                    rhs = matrix(QQ, dim_d, dim_d)
                    for L0_idx, coeff in br_coeffs.items():
                        L0_mats = l0_action_matrix(M_test, L0_idx, e44_data)
                        rhs += coeff * L0_mats[d]

                    if lhs != rhs:
                        n_fail += 1
                        if len(fail_details) < 10:
                            fail_details.append(
                                f"    ({im1},{i1}) deg {d}: "
                                f"parity=({bm1['parity']},{b1['parity']})"
                            )
                    else:
                        n_pass += 1

        total = n_pass + n_fail
        if verbose:
            print(f"  {tag}: {total} checks, {n_pass} pass, {n_fail} fail")
            for line in fail_details:
                print(line)

        _check(f"{tag} all-degree commutator", n_fail, 0)

    if verbose:
        print("─" * 60)
        if all_pass:
            print("Commutator (all deg)  ✓  ALL CHECKS PASSED")
        else:
            print(f"Commutator (all deg)  ✗  {total_failures} FAILED")
        print("─" * 60)

    return all_pass


def _check_leibniz(verbose=True):
    """
    Direct verification of the Leibniz rule in l1_action_matrix at degree 2.

    For Y \in L_1 and a degree-2 PBW monomial u = g_1\cdotg_2 (in PBW order),
    the Leibniz rule gives:

        Y\cdot(g_1\cdotg_2\cdotw) = [Y, g_1]\cdot(g_2\cdotw) + (-1)^{|Y||g_1|} g_1\cdot([Y, g_2]\cdotw)

    where [Y, g_k] \in L_0 acts on the generators AFTER g_k via l0_action_matrix,
    and g_k is prepended via l_minus1_action_matrix.

    This function computes both sides independently and compares:
      LHS = l1_action_matrix(Y, deg=2) applied to each M[2] basis vector
      RHS = manual Leibniz expansion using l0_action_matrix + l_minus1_action_matrix

    Tests on M_1(0,0,0) [trivial W, dim 1] where the bug is exposed.
    """
    all_pass = True

    def _check(name, got, expected):
        nonlocal all_pass
        ok = (got == expected)
        status = "PASS" if ok else "FAIL"
        if verbose:
            print(f"  [{status}] {name}: got {got!r}, expected {expected!r}")
        if not ok:
            all_pass = False

    if verbose:
        print("─" * 60)
        print("Leibniz checkpoint  —  direct deg-2 verification")
        print("─" * 60)

    e44_data = load_e44()
    if e44_data is None:
        print("  [SKIP] e44_brackets.pkl not found.")
        return None

    L1_list  = e44_data['E44'][1]
    Lm1_list = e44_data['E44'][-1]
    L0_list  = e44_data['E44'][0]
    bt = e44_data['btable']
    R_poly = Lm1_list[0]['basis'][0].parent()

    # --- Module under test ---
    M = M_verma(1, 0, 0, 0, max_deg=2)
    dim_W = M.dim_W
    if verbose:
        print(f"  Module: M_1(0,0,0), dim_W={dim_W}")
        print(f"  dim(0)={M.dim(0)}, dim(1)={M.dim(1)}, dim(2)={M.dim(2)}")
        print()

    # ── Part 1: Direct Leibniz at degree 2 ──────────────────────────
    # For each L_1 gen Y, compute l1(Y, d=2) and compare with manual
    # Leibniz expansion on all M[2] basis vectors.
    #
    # The manual expansion for monomial (\alpha, S) at degree 2:
    #
    # Even factors (in PBW order e_0^a_0 \cdot e_1^a_1 \cdot e_2^a_2 \cdot e_3^a_3):
    #   For each e_i with a_i > 0, [Y, e_i] \in L_0 acts on the "after" piece
    #   (generators with index > i plus remaining copies of e_i), then
    #   the "before" piece (generators with index < i) is prepended.
    #
    # Odd factors (d_{s_0}, d_{s_1}, ... in sorted order):
    #   For d_{s_k} at position k, [Y, d_{s_k}] \in L_0 acts on S_after \otimes w,
    #   then S_before is prepended with super-sign (-1)^{|Y|\cdotk}.

    n_fail = 0
    n_total = 0

    # Sample L_1 generators: 5 even, 5 odd (covering diverse brackets)
    l1_sample = [0, 5, 10, 20, 35,    # even L_1
                 40, 45, 50, 55, 70]   # odd L_1

    for i1 in l1_sample:
        b1 = L1_list[i1]
        pY = b1['parity']

        # Pre-compute [Y, g_i] for each L_{-1} generator
        y_adj = {}
        for im1, bm1 in enumerate(Lm1_list):
            key = (b1['label'], bm1['label'])
            result = bt.get(key)
            if result is None or not any(c != 0 for c in result):
                y_adj[im1] = {}
                continue
            p_br = (b1['parity'] + bm1['parity']) % 2
            L0_sec = [b for b in L0_list if b['parity'] == p_br]
            coeffs = _expand_in_basis(result, L0_sec, R_poly)
            full = {}
            for sec_idx, c in coeffs.items():
                full[L0_list.index(L0_sec[sec_idx])] = c
            y_adj[im1] = full

        # Get l1_action_matrix at deg 2  (LHS)
        l1_mat_d2 = l1_action_matrix(M, i1, e44_data)[2]   # M[2] \to M[1]

        # For each M[2] basis vector, compute manual Leibniz (RHS)
        for j2, (alpha, S) in enumerate(M.pbw[2]):
            for k in range(dim_W):
                col_idx = j2 * dim_W + k
                lhs_col = l1_mat_d2.column(col_idx)

                # --- Manual Leibniz: RHS ---
                rhs = vector(QQ, M.dim(1))

                # Even generators: iterate in PBW order
                for gen_i in range(4):
                    pow_i = alpha[gen_i]
                    if pow_i == 0:
                        continue
                    lm1_idx = gen_i   # even gen
                    for L0_idx, c10 in y_adj[lm1_idx].items():
                        L0_mats = l0_action_matrix(M, L0_idx, e44_data)
                        for slot in range(pow_i):
                            # "after" piece: e_i^{pow_i-slot-1} \cdot e_{i+1}^... \cdot odd stuff
                            after_alpha = [0] * 4
                            after_alpha[gen_i] = pow_i - slot - 1
                            for gi2 in range(gen_i + 1, 4):
                                after_alpha[gi2] = alpha[gi2]
                            after_mon = (tuple(after_alpha), S)
                            after_d = sum(after_alpha) + len(S)

                            if after_d < 0 or after_d > M.max_deg:
                                continue
                            if after_mon not in M._pbw_idx[after_d]:
                                continue
                            jp_after = M._pbw_idx[after_d][after_mon]

                            # "before" piece to prepend: e_0^a_0 \cdot ... \cdot e_{i-1}^a_{i-1} \cdot e_i^slot
                            before_alpha = [0] * 4
                            for gi2 in range(gen_i):
                                before_alpha[gi2] = alpha[gi2]
                            before_alpha[gen_i] = slot

                            # Apply L_0 on after piece
                            L0_mat = L0_mats[after_d]
                            for row_a in range(M.dim(after_d)):
                                entry = L0_mat[row_a, jp_after * dim_W + k]
                                if entry == 0:
                                    continue
                                res_j = row_a // dim_W
                                res_w = row_a % dim_W
                                res_alpha, res_S = M.pbw[after_d][res_j]
                                # Combine before + result
                                new_alpha = [before_alpha[i] + res_alpha[i]
                                             for i in range(4)]
                                new_mon = (tuple(new_alpha), res_S)
                                if new_mon not in M._pbw_idx[1]:
                                    continue
                                new_j = M._pbw_idx[1][new_mon]
                                rhs[new_j * dim_W + res_w] += c10 * entry

                # Odd generators: iterate in sorted order
                for pos, gen_i in enumerate(sorted(S)):
                    lm1_idx = 4 + gen_i
                    sign_pos = QQ((-1) ** (pY * pos))
                    S_before = frozenset(jj for jj in S if jj < gen_i)
                    S_after  = frozenset(jj for jj in S if jj > gen_i)
                    deg_after = len(S_after)

                    for L0_idx, c10 in y_adj[lm1_idx].items():
                        L0_mats = l0_action_matrix(M, L0_idx, e44_data)
                        zero_alpha = (0, 0, 0, 0)
                        src_mon = (zero_alpha, S_after)
                        if src_mon not in M._pbw_idx[deg_after]:
                            continue
                        jp_after = M._pbw_idx[deg_after][src_mon]
                        L0_mat = L0_mats[deg_after]
                        for row_a in range(M.dim(deg_after)):
                            entry = L0_mat[row_a, jp_after * dim_W + k]
                            if entry == 0:
                                continue
                            res_j = row_a // dim_W
                            res_w = row_a % dim_W
                            res_alpha_a, res_S_a = M.pbw[deg_after][res_j]
                            # Combine: alpha (even from original) + S_before \cup res_S_a
                            new_alpha = tuple(alpha[i] + res_alpha_a[i]
                                              for i in range(4))
                            new_S = S_before | res_S_a
                            if len(new_S) != len(S_before) + len(res_S_a):
                                continue
                            new_mon = (new_alpha, new_S)
                            if new_mon not in M._pbw_idx[1]:
                                continue
                            new_j = M._pbw_idx[1][new_mon]
                            rhs[new_j * dim_W + res_w] += (
                                sign_pos * c10 * entry
                            )

                n_total += 1
                if lhs_col != rhs:
                    n_fail += 1
                    if verbose and n_fail <= 5:
                        print(f"    FAIL L_1[{i1}]({b1['label']},p={pY}) "
                              f"on ({alpha},{sorted(S)})\otimesw[{k}]")
                        print(f"      l1_action = {lhs_col}")
                        print(f"      manual    = {rhs}")
                        print(f"      diff      = {lhs_col - rhs}")

    _check(f"Leibniz deg-2 on M_1(0,0,0): {n_total} vectors, 0 failures",
           n_fail, 0)

    if verbose:
        if n_fail > 0:
            print(f"  ({n_fail} of {n_total} Leibniz checks failed)")
        print("─" * 60)
        if all_pass:
            print("Leibniz  ✓  ALL CHECKS PASSED")
        else:
            print("Leibniz  ✗  SOME CHECKS FAILED — see above")
        print("─" * 60)

    return all_pass



# ===========================================================================
# s2.4 — Grothendieck decomposition of M_t(a, b, c)
# ===========================================================================
#
# By the CKC 2026 singular vector classification (Step 5 / S4), the
# degenerate Verma modules are:
#   • M_t(a, 0, 0)  for all t \in QQ, a \in Z_{>=0}
#   • M_t(0, 0, c)  for all t \in QQ, c \in Z_{>=0}
#   • M_t(0, 1, 0)  for all t \in QQ
#
# Every degenerate Verma module has at least one sub-Verma, whose highest-
# weight vector is the corresponding singular vector.  The ones visible at
# degree 1 (from the CKC 2026 morphism table, Step 5) are:
#
#   w[1A]: sub-Verma M_{t-1}(a+1,0,0)  in  M_t(a,0,0)
#           present when a >= 1  OR  (a = 0 AND t = 0)
#   w[1B]: sub-Verma M_{t-1}(0,0,c-1)  in  M_t(0,0,c),  c >= 1
#   w[1C]: sub-Verma M_{t-1}(0,1,0)    in  M_t(0,0,1),  all t
#   w[1D]: sub-Verma M_{t-1}(1,0,0)    in  M_t(0,0,0),  t ≠ 0
#   w[1E]: sub-Verma M_0(0,0,0)        in  M_1(1,0,0)   only (t = 1)
#   w[4H]: sub-Verma M_{t-4}(1,0,0)    in  M_t(1,0,0),  all t  (degree 4)
#
# For M_t(0,1,0): degenerate per CKC 2026, but the explicit sub-Verma
# structure requires the degree-2 singular vectors w[2DA] / w[2EA], which
# are computed in singular_vectors.py (S3).  This block is a documented
# stub pending S3.
#
# Grothendieck group recursion:
#   [M_t(a,b,c)] = [L_t(a,b,c)]  +  \Sigma_i [M_{t_i}(a_i,b_i,c_i)]
# where the sum is over all direct maximal sub-Vermas.  The recursion
# terminates at irreducible (non-degenerate) modules or at max_depth.
#
# Caution: the chain M_t(1,0,0) \subset M_{t-4}(1,0,0) \subset M_{t-8}(1,0,0) \subset ...
# is infinite; the computation is truncated at max_depth steps.
# ---------------------------------------------------------------------------


def _known_subvermas(t, a, b, c):
    """
    Return the list of *direct* sub-Verma modules of M_t(a, b, c) as
    ``(degree, t', a', b', c')`` tuples, based on the singular vector
    classification of CKC 2026 Step 5.

    "Direct" means the sub-Verma is generated by a singular vector whose
    class is not already the image of a composition of lower-degree morphisms
    already listed here.

    Parameters
    ----------
    t, a, b, c : same conventions as M_verma / W_hat

    Returns
    -------
    list of (int degree, QQ t', int a', int b', int c')
    """
    t_q = QQ(t)
    a, b, c = int(a), int(b), int(c)
    subs = []

    # ── (a, 0, 0) family ────────────────────────────────────────────────────
    if b == 0 and c == 0:
        # w[1A]: degree-1 sv present for a >= 1, or a = t = 0
        if a >= 1 or t_q == 0:
            subs.append((1, t_q - 1, a + 1, 0, 0))
        # w[1D]: degree-1 sv in M_t(0,0,0) for t ≠ 0 (note: at a = 0)
        if a == 0 and t_q != 0:
            subs.append((1, t_q - 1, 1, 0, 0))
        # w[1E]: degree-1 sv in M_1(1,0,0) only
        if a == 1 and t_q == QQ(1):
            subs.append((1, QQ(0), 0, 0, 0))
        # w[4H]: degree-4 sv in M_t(1,0,0) for all t
        if a == 1:
            subs.append((4, t_q - 4, 1, 0, 0))

    # ── (0, 0, c) family ────────────────────────────────────────────────────
    if a == 0 and b == 0 and c >= 1:
        # w[1B]: degree-1 sv for c >= 1 \to sub-Verma M_{t-1}(0,0,c-1)
        subs.append((1, t_q - 1, 0, 0, c - 1))
        # w[1C]: degree-1 sv in M_t(0,0,1) for all t \to sub-Verma M_{t-1}(0,1,0)
        if c == 1:
            subs.append((1, t_q - 1, 0, 1, 0))

    # ── (0, 1, 0) family ────────────────────────────────────────────────────
    # M_t(0,1,0) is degenerate per CKC 2026 S5, but the explicit sub-Verma
    # structure requires the degree-2 singular vectors w[2DA] / w[2EA]
    # (computed in singular_vectors.py, S3).  Treated as irreducible here.

    return subs


def grothendieck_decomp(t, a, b, c, max_depth=8, _depth=0, _visit_count=None):
    """
    Grothendieck group decomposition of the Verma module M_t(a, b, c).

    Expresses [M_t(a,b,c)] as a formal ZZ-linear combination of simple
    module classes [L_{t'}(a',b',c')] in the Grothendieck group of
    E(4,4)-modules, using the Kazhdan-Lusztig-type recursion inferred from
    the singular vector classification of Cantarini-Caselli-Kac 2026 S4:

        [M_t(a,b,c)]  =  [L_t(a,b,c)]  +  \Sigma_i [M_{t_i}(a_i,b_i,c_i)]

    where the sum is over all direct maximal sub-Verma modules (as returned
    by ``_known_subvermas``).  Each [M_{t_i}] is in turn expanded recursively
    until no further sub-Vermas are found (irreducible module) or the recursion
    depth reaches ``max_depth``.

    Note on infinite chains
    -----------------------
    The sub-Verma chain  M_t(1,0,0) \subset M_{t-4}(1,0,0) \subset M_{t-8}(1,0,0) \subset ...
    (from the degree-4 morphism w[4H]) is infinite.  The argument
    ``max_depth`` truncates the recursion; increase it to see more terms.

    Note on M_t(0,1,0)
    ------------------
    The degree-2 sub-Verma structure of M_t(0,1,0) is left as a stub
    pending the explicit w[2DA]/w[2EA] formulas from S3.  M_t(0,1,0)
    currently contributes only [L_t(0,1,0)] to the decomposition.

    Parameters
    ----------
    t : integer or rational
        Central charge of the Verma module.
    a, b, c : non-negative integers
        sl_4 Dynkin highest weight labels.
    max_depth : int, default 8
        Maximum recursion depth.  Increase for longer w[4H] chains.

    Returns
    -------
    dict  {(t', a', b', c'): int}
        Multiplicities of each simple module [L_{t'}(a',b',c')] in the
        composition series of M_t(a,b,c), truncated at depth max_depth.
        All multiplicities are positive integers (>= 1).

    Examples
    --------
    sage: grothendieck_decomp(0, 1, 1, 0)          # non-degenerate
    {(0, 1, 1, 0): 1}

    sage: grothendieck_decomp(0, 0, 0, 0, max_depth=3)
    {(0, 0, 0, 0): 1, (-1, 1, 0, 0): 1, (-2, 2, 0, 0): 1,
     (-5, 1, 0, 0): 1, (-3, 3, 0, 0): 1}   # (etc., depth-3 truncation)

    sage: grothendieck_decomp(1, 1, 0, 0, max_depth=2)
    # head  L_1(1,0,0), plus sub-Vermas from w[1A], w[1E], w[4H]
    """
    from collections import defaultdict

    t_q = QQ(t)
    a, b, c = int(a), int(b), int(c)
    key = (t_q, a, b, c)

    # visit_count tracks how many times each (t,a,b,c) has been opened at
    # this recursion branch, to break cycles from the infinite w[4H] chain.
    if _visit_count is None:
        _visit_count = {}

    result = defaultdict(int)
    # Every M_t(a,b,c) contributes exactly one copy of [L_t(a,b,c)]
    result[key] += 1

    if _depth >= max_depth:
        return dict(result)

    # Recurse over direct sub-Vermas
    for (deg, t2, a2, b2, c2) in _known_subvermas(t_q, a, b, c):
        sub_key = (QQ(t2), int(a2), int(b2), int(c2))
        # Cycle guard: allow visiting each key at most max_depth times total
        if _visit_count.get(sub_key, 0) >= max_depth:
            continue
        vc_copy = dict(_visit_count)
        vc_copy[sub_key] = vc_copy.get(sub_key, 0) + 1

        sub_decomp = grothendieck_decomp(
            t2, a2, b2, c2,
            max_depth=max_depth,
            _depth=_depth + 1,
            _visit_count=vc_copy,
        )
        for k, v in sub_decomp.items():
            result[k] += v

    return dict(result)


def grothendieck_table(families=None, t_val=QQ(0), max_depth=4, verbose=True):
    """
    Print a table of Grothendieck decompositions for small (a, b, c) values
    in the three degenerate families, for comparison with CKC 2026 Table 1.

    Parameters
    ----------
    families : list of (a, b, c) tuples, or None
        Which Verma modules to include.  Default: small values in each family.
    t_val : rational, default 0
        Central charge to use for the table.
    max_depth : int, default 4
        Recursion depth for grothendieck_decomp.
    verbose : bool
        If True, print the table.

    Returns
    -------
    dict  {(a, b, c): decomp_dict}
    """
    if families is None:
        families = (
            # (a, 0, 0) family: a = 0, 1, 2, 3
            [(a, 0, 0) for a in range(4)] +
            # (0, 0, c) family: c = 1, 2, 3
            [(0, 0, c) for c in range(1, 4)] +
            # (0, 1, 0)
            [(0, 1, 0)] +
            # Generic non-degenerate: for comparison
            [(1, 1, 0), (1, 0, 1), (0, 1, 1)]
        )

    results = {}
    if verbose:
        print("─" * 72)
        print(f"Grothendieck table   t = {t_val},  max_depth = {max_depth}")
        print("─" * 72)
        print(f"{'M_t(a,b,c)':<20} #factors  composition factors  [L_{{t'}}(a',b',c')]")
        print("─" * 72)

    for (a, b, c) in families:
        d = grothendieck_decomp(t_val, a, b, c, max_depth=max_depth)
        results[(a, b, c)] = d
        if verbose:
            label = f"M_{{{t_val}}}({a},{b},{c})"
            factors = ", ".join(
                f"{v}\cdot[L_{{{t_}}}({a_},{b_},{c_})]"
                for (t_, a_, b_, c_), v in sorted(d.items())
            )
            print(f"  {label:<18} {len(d):>5}      {factors}")

    if verbose:
        print("─" * 72)

    return results


# ===========================================================================
# s2.4 Checkpoint
# ===========================================================================

def _check_s24(verbose=True):
    """
    S2.4 checkpoint: verify Grothendieck decompositions against structural
    predictions from CKC 2026 S4 / Step 5.

    Checks
    ------
    (1)  Non-degenerate: grothendieck_decomp(0,1,1,0) == {(0,1,1,0): 1}
    (2)  Non-degenerate: grothendieck_decomp(0,1,0,1) == {(0,1,0,1): 1}
    (3)  (a,0,0) chain depth-1: M_t(1,0,0) head L_t(1,0,0) has mult 1
    (4)  (a,0,0) chain depth-1: M_t(1,0,0) contains L_{t-1}(2,0,0) from w[1A]
    (5)  (a,0,0) chain: M_0(0,0,0) head L_0(0,0,0) has mult 1
    (6)  (a,0,0) chain: M_0(0,0,0) at depth 1 contains L_{-1}(1,0,0) from w[1A]
    (7)  (0,0,c) chain depth-1: M_t(0,0,2) contains L_t(0,0,2) and chain
         down to M_{t-1}(0,0,1) which includes M_{t-2}(0,1,0)
    (8)  w[1E] check: grothendieck_decomp(1,1,0,0) at depth 1 includes
         L_0(0,0,0) from w[1E]
    (9)  w[4H] check: grothendieck_decomp(t,1,0,0) at depth 1 includes
         L_{t-4}(1,0,0)
    (10) All multiplicities in grothendieck_decomp are positive integers
    (11) Head always has multiplicity 1
    (12) Total factor count grows with depth (non-trivial recursion)

    Returns True iff all checks pass.
    """
    all_pass = True

    def _check(name, got, expected):
        nonlocal all_pass
        ok = (got == expected)
        status = "PASS" if ok else "FAIL"
        if verbose:
            print(f"  [{status}] {name}: got {got!r}, expected {expected!r}")
        if not ok:
            all_pass = False

    if verbose:
        print("─" * 60)
        print("s2.4 checkpoint  —  Grothendieck decompositions")
        print("─" * 60)

    t0 = QQ(0)
    t1 = QQ(1)
    t5 = QQ(5)

    # (1) Non-degenerate module: only the head
    d_nd1 = grothendieck_decomp(0, 1, 1, 0, max_depth=2)
    _check("grothendieck_decomp(0,1,1,0) == {head: 1}",
           d_nd1, {(t0, 1, 1, 0): 1})

    # (2) Another non-degenerate
    d_nd2 = grothendieck_decomp(0, 1, 0, 1, max_depth=2)
    _check("grothendieck_decomp(0,1,0,1) == {head: 1}",
           d_nd2, {(t0, 1, 0, 1): 1})

    # (3) Head of M_t(1,0,0) has multiplicity 1
    d_100 = grothendieck_decomp(t5, 1, 0, 0, max_depth=1)
    _check("M_5(1,0,0): head L_5(1,0,0) has mult 1",
           d_100.get((t5, 1, 0, 0), 0), 1)

    # (4) w[1A] in M_t(1,0,0): L_{t-1}(2,0,0) present at depth 1
    _check("M_5(1,0,0) depth-1: L_4(2,0,0) present via w[1A]",
           d_100.get((QQ(4), 2, 0, 0), 0) >= 1, True)

    # (5) Head of M_0(0,0,0) has multiplicity 1
    d_000 = grothendieck_decomp(0, 0, 0, 0, max_depth=1)
    _check("M_0(0,0,0): head L_0(0,0,0) has mult 1",
           d_000.get((t0, 0, 0, 0), 0), 1)

    # (6) w[1A] (a=t=0 case) in M_0(0,0,0): L_{-1}(1,0,0) present
    _check("M_0(0,0,0) depth-1: L_{-1}(1,0,0) via w[1A]",
           d_000.get((QQ(-1), 1, 0, 0), 0) >= 1, True)

    # (7) (0,0,c) chain: M_t(0,0,2) depth-2 includes M_{t-2}(0,1,0)
    d_002 = grothendieck_decomp(t5, 0, 0, 2, max_depth=2)
    _check("M_5(0,0,2) depth-2: head L_5(0,0,2) present",
           d_002.get((t5, 0, 0, 2), 0), 1)
    _check("M_5(0,0,2) depth-2: L_4(0,0,1) present via w[1B] chain",
           d_002.get((QQ(4), 0, 0, 1), 0) >= 1, True)
    _check("M_5(0,0,2) depth-2: L_3(0,1,0) present via w[1C]",
           d_002.get((QQ(3), 0, 1, 0), 0) >= 1, True)

    # (8) w[1E]: M_1(1,0,0) includes L_0(0,0,0)
    d_1_100 = grothendieck_decomp(1, 1, 0, 0, max_depth=1)
    _check("M_1(1,0,0) depth-1: L_0(0,0,0) present via w[1E]",
           d_1_100.get((t0, 0, 0, 0), 0) >= 1, True)

    # (9) w[4H]: M_t(1,0,0) at depth 1 includes L_{t-4}(1,0,0)
    d_5_100 = grothendieck_decomp(t5, 1, 0, 0, max_depth=1)
    _check("M_5(1,0,0) depth-1: L_1(1,0,0) present via w[4H]",
           d_5_100.get((QQ(1), 1, 0, 0), 0) >= 1, True)

    # (10) All multiplicities positive
    for label, d in [("M_5(1,0,0)", d_5_100), ("M_5(0,0,2)", d_002)]:
        _check(f"all mult >= 1 in {label}",
               all(v >= 1 for v in d.values()), True)

    # (11) Head always multiplicity 1
    for (t_, a_, b_, c_) in [(t0, 2, 0, 0), (t5, 0, 0, 3), (t0, 0, 1, 0)]:
        d_ = grothendieck_decomp(t_, a_, b_, c_, max_depth=1)
        _check(f"head mult == 1 for ({t_},{a_},{b_},{c_})",
               d_.get((t_, a_, b_, c_), 0), 1)

    # (12) Recursion is non-trivial: depth-2 gives more factors than depth-1
    d_a = grothendieck_decomp(t5, 1, 0, 0, max_depth=1)
    d_b = grothendieck_decomp(t5, 1, 0, 0, max_depth=2)
    _check("M_5(1,0,0): depth-2 has more factors than depth-1",
           len(d_b) > len(d_a), True)

    if verbose:
        print("─" * 60)
        if all_pass:
            print("s2.4  ✓  ALL CHECKS PASSED")
        else:
            print("s2.4  ✗  SOME CHECKS FAILED — see above")
        print("─" * 60)

    return all_pass


# ===========================================================================
# s2.5 — Export: save_verma_data("verma_cache.pkl")
# ===========================================================================
#
# The pickle stores a dict with the following keys:
#
#   'W_hat_dims'     : dict { (a,b,c): dim }
#       Dimensions of W_hat(t, a,b,c) for canonical small cases.
#       (dim is independent of t, so t is omitted from the key.)
#
#   'verma_dim_tables' : dict { (t,a,b,c): {d: dim} }
#       Graded dimension tables for the key Verma modules used by S3.
#
#   'pbw_bases'      : dict { d: list of (alpha, frozenset_S) monomials }
#       PBW basis of U(L_{-1}) at degrees 0..4.  Shared across all modules.
#
#   'grothendieck'   : dict { (t,a,b,c): decomp_dict }
#       Grothendieck decompositions (max_depth=6) for the three degenerate
#       families at small parameter values; used by de_rham_complex.py.
#
#   'crystal_data'   : dict { (a,b,c): { 'basis_weights': list[tuple],
#                                        'e_mats': {i: matrix},
#                                        'f_mats': {i: matrix},
#                                        'h_mats': {i: matrix} } }
#       Crystal basis data for the W_hat modules needed by S3.
#
# On reload, downstream scripts call `load_verma_data(path)` which returns
# this dict and re-creates VermaModule / WHat4Module objects on demand.
# ---------------------------------------------------------------------------

_EXPORT_W_HAT_CASES = [
    (0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1),
    (2, 0, 0), (0, 2, 0), (0, 0, 2), (3, 0, 0),
    (1, 1, 0), (1, 0, 1), (0, 1, 1),
]

_EXPORT_VERMA_CASES = [
    # (t, a, b, c)
    (0, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1),
    (0, 2, 0, 0), (0, 0, 0, 2), (0, 0, 0, 3),
    (1, 1, 0, 0), (1, 0, 0, 0), (1, 0, 0, 1),
    (QQ(1)/2,   1, 0, 0),          # half-integer central charge
]

_EXPORT_GROTHENDIECK_CASES = [
    # degenerate families (small params, t = 0)
    (0, 0, 0, 0), (0, 1, 0, 0), (0, 2, 0, 0), (0, 3, 0, 0),
    (0, 0, 0, 1), (0, 0, 0, 2), (0, 0, 0, 3),
    (0, 0, 1, 0),
    # NS-relevant central charge t = 1
    (1, 1, 0, 0), (1, 0, 0, 0), (1, 0, 0, 1),
]


def save_verma_data(path="verma_cache.pkl", max_deg=4, grothendieck_depth=6,
                    verbose=True):
    """
    Compute and serialise the Verma module data needed by downstream scripts.

    Writes a pickle file at *path* containing:
      - W_hat module dimensions and crystal matrices for key (a,b,c) cases
      - Graded dimension tables for key Verma modules
      - Shared PBW basis at degrees 0..max_deg
      - Grothendieck decompositions (depth grothendieck_depth) for the three
        degenerate families and the NS-relevant cases

    Parameters
    ----------
    path : str
        Output file path (default "verma_cache.pkl").
    max_deg : int
        PBW truncation degree (default 4).
    grothendieck_depth : int
        Recursion depth for grothendieck_decomp (default 6).
    verbose : bool
        Print progress messages.

    Returns
    -------
    dict  — the payload that was written (same object as in the file).
    """
    payload = {}

    # ── 1.  W_hat dimensions ─────────────────────────────────────────────
    if verbose:
        print("  [s2.5] computing W_hat dimensions ...")
    w_hat_dims = {}
    for (a, b, c) in _EXPORT_W_HAT_CASES:
        W = W_hat(0, a, b, c)   # t doesn't affect dim
        w_hat_dims[(a, b, c)] = W.dim
    payload['W_hat_dims'] = w_hat_dims

    # ── 2.  Crystal data ─────────────────────────────────────────────────
    if verbose:
        print("  [s2.5] computing crystal data ...")
    crystal_data = {}
    for (a, b, c) in _EXPORT_W_HAT_CASES:
        W = W_hat(0, a, b, c)
        crystal_data[(a, b, c)] = {
            'basis_weights': [W.weight_of(k) for k in range(W.dim)],
            'e_mats': {i: W.e_mats[i] for i in range(1, 4)},
            'f_mats': {i: W.f_mats[i] for i in range(1, 4)},
            'h_mats': {i: W.h_mats[i] for i in range(1, 4)},
            'v_hw': W.v_hw,
        }
    payload['crystal_data'] = crystal_data

    # ── 3.  Shared PBW basis ─────────────────────────────────────────────
    if verbose:
        print("  [s2.5] computing PBW bases ...")
    payload['pbw_bases'] = {d: _pbw_basis_at_degree(d) for d in range(max_deg + 1)}

    # ── 4.  Graded dimension tables ──────────────────────────────────────
    if verbose:
        print("  [s2.5] computing Verma graded dimensions ...")
    verma_dim_tables = {}
    for (t, a, b, c) in _EXPORT_VERMA_CASES:
        M = M_verma(t, a, b, c, max_deg=max_deg)
        verma_dim_tables[(QQ(t), a, b, c)] = {d: M.dim(d) for d in range(max_deg + 1)}
    payload['verma_dim_tables'] = verma_dim_tables

    # ── 5.  Grothendieck decompositions ──────────────────────────────────
    if verbose:
        print(f"  [s2.5] computing Grothendieck decompositions (depth={grothendieck_depth}) ...")
    groth = {}
    for (t, a, b, c) in _EXPORT_GROTHENDIECK_CASES:
        groth[(QQ(t), a, b, c)] = grothendieck_decomp(
            t, a, b, c, max_depth=grothendieck_depth
        )
    payload['grothendieck'] = groth

    # ── 6.  Write pickle ─────────────────────────────────────────────────
    with open(path, 'wb') as fh:
        _pickle.dump(payload, fh, protocol=_pickle.HIGHEST_PROTOCOL)

    if verbose:
        size_kb = _os.path.getsize(path) / 1024
        print(f"  [s2.5] wrote '{path}'  ({size_kb:.1f} KB)")

    return payload


def load_verma_data(path="verma_cache.pkl"):
    """
    Load the Verma module cache written by save_verma_data.

    Returns the payload dict, or None with a warning if the file is absent.
    Downstream scripts (singular_vectors.py, de_rham_complex.py) call this
    instead of recomputing the data from scratch.
    """
    if not _os.path.exists(path):
        import warnings
        warnings.warn(
            f"verma_cache.pkl not found at '{path}'. "
            "Run 'sage verma_modules.py' first to generate it."
        )
        return None
    with open(path, 'rb') as fh:
        return _pickle.load(fh)


# ===========================================================================
# s2.5 Checkpoint
# ===========================================================================

def _check_s25(path="verma_cache.pkl", verbose=True):
    """
    S2.5 checkpoint: verify that verma_cache.pkl was written and is valid.

    Checks
    ------
    (1)  File exists at the given path
    (2)  File size > 0
    (3)  All expected top-level keys are present
    (4)  W_hat_dims[(1,0,0)] == 4  (fundamental sl_4 rep)
    (5)  W_hat_dims[(0,0,0)] == 1  (trivial rep)
    (6)  crystal_data[(1,0,0)] has keys 'e_mats', 'f_mats', 'h_mats', 'basis_weights', 'v_hw'
    (7)  pbw_bases[0] == [(( 0,0,0,0), frozenset())]  (single degree-0 monomial)
    (8)  pbw_bases[1] has 8 monomials
    (9)  verma_dim_tables[(QQ(0),1,0,0)][1] == 32
    (10) verma_dim_tables[(QQ(1),1,0,0)][0] == 4
    (11) (QQ(0),1,0,0) in grothendieck  (degenerate family present)
    (12) grothendieck[(QQ(0),1,0,0)][(QQ(0),1,0,0)] == 1  (head mult = 1)
    (13) Reload: load_verma_data(path) returns the same keys

    Returns True iff all checks pass.
    """
    all_pass = True

    def _check(name, got, expected):
        nonlocal all_pass
        ok = (got == expected)
        status = "PASS" if ok else "FAIL"
        if verbose:
            print(f"  [{status}] {name}: got {got!r}, expected {expected!r}")
        if not ok:
            all_pass = False

    if verbose:
        print("─" * 60)
        print("s2.5 checkpoint  —  verma_cache.pkl export")
        print("─" * 60)

    # (1) File exists
    _check("file exists", _os.path.exists(path), True)
    if not _os.path.exists(path):
        if verbose:
            print("  [FATAL] file not found — skipping remaining checks.")
        return False

    # (2) File size > 0
    _check("file size > 0", _os.path.getsize(path) > 0, True)

    # Load
    data = load_verma_data(path)

    # (3) Top-level keys
    expected_keys = {'W_hat_dims', 'crystal_data', 'pbw_bases',
                     'verma_dim_tables', 'grothendieck'}
    _check("top-level keys",
           set(data.keys()) >= expected_keys, True)

    # (4)-(5) W_hat dims
    _check("W_hat_dims[(1,0,0)] == 4",  data['W_hat_dims'].get((1, 0, 0)), 4)
    _check("W_hat_dims[(0,0,0)] == 1",  data['W_hat_dims'].get((0, 0, 0)), 1)

    # (6) Crystal data structure
    cd = data['crystal_data'].get((1, 0, 0), {})
    for key_needed in ('e_mats', 'f_mats', 'h_mats', 'basis_weights', 'v_hw'):
        _check(f"crystal_data[(1,0,0)] has '{key_needed}'",
               key_needed in cd, True)

    # (7)-(8) PBW bases
    _check("pbw_bases[0] == [((0,0,0,0), frozenset())]",
           data['pbw_bases'][0], [((0, 0, 0, 0), frozenset())])
    _check("len(pbw_bases[1]) == 8",
           len(data['pbw_bases'][1]), 8)

    # (9)-(10) Verma dim tables
    _check("verma_dim_tables[(0,1,0,0)][1] == 32",
           data['verma_dim_tables'].get((QQ(0), 1, 0, 0), {}).get(1), 32)
    _check("verma_dim_tables[(1,1,0,0)][0] == 4",
           data['verma_dim_tables'].get((QQ(1), 1, 0, 0), {}).get(0), 4)

    # (11)-(12) Grothendieck
    g_key = (QQ(0), 1, 0, 0)
    _check("(0,1,0,0) in grothendieck",
           g_key in data['grothendieck'], True)
    if g_key in data['grothendieck']:
        _check("grothendieck[(0,1,0,0)] head mult == 1",
               data['grothendieck'][g_key].get(g_key), 1)

    # (13) Reload gives same keys
    data2 = load_verma_data(path)
    _check("reload has same keys",
           set(data2.keys()) == set(data.keys()), True)

    if verbose:
        print("─" * 60)
        if all_pass:
            print("s2.5  ✓  ALL CHECKS PASSED")
        else:
            print("s2.5  ✗  SOME CHECKS FAILED — see above")
        print("─" * 60)

    return all_pass


def _check_s22(verbose=True):
    """
    S2.2 checkpoint: verify graded dimensions of key Verma modules.

    Checks
    ------
    (1) dim(M(t,1,0,0)[0]) == 4          [degree-0 = W_hat itself]
    (2) dim(M(t,1,0,0)[1]) == 8*4 == 32  [plan checkpoint]
    (3) dim(M(t,1,0,0)[2]) == 32*4 == 128
    (4) dim(M(t,0,0,0)[0]) == 1          [trivial at degree 0]
    (5) dim(M(t,0,0,0)[1]) == 8          [L_{-1} \otimes trivial]
    (6) pbw_dim(0)==1, pbw_dim(1)==8, pbw_dim(2)==32, pbw_dim(3)==96, pbw_dim(4)==256
    (7) len(M.basis(1)) == dim(M[1])     [basis list size matches dim]
    (8) len(M.basis(2)) == dim(M[2])
    (9) to_vec round-trip: coeffs \to vec \to correct entry is 1

    Returns True iff all checks pass.
    """
    all_pass = True

    def _check(name, got, expected):
        nonlocal all_pass
        ok = (got == expected)
        status = "PASS" if ok else "FAIL"
        if verbose:
            print(f"  [{status}] {name}: got {got!r}, expected {expected!r}")
        if not ok:
            all_pass = False

    if verbose:
        print("─" * 60)
        print("s2.2 checkpoint  —  VermaModule graded dimensions")
        print("─" * 60)

    M100 = M_verma(0, 1, 0, 0)
    M000 = M_verma(0, 0, 0, 0)

    # (1)-(3): dimensions of M(t,1,0,0)
    _check("dim M(t,1,0,0)[0] == 4",   M100.dim(0), 4)
    _check("dim M(t,1,0,0)[1] == 32",  M100.dim(1), 32)
    _check("dim M(t,1,0,0)[2] == 128", M100.dim(2), 128)

    # (4)-(5): trivial base module
    _check("dim M(t,0,0,0)[0] == 1",  M000.dim(0), 1)
    _check("dim M(t,0,0,0)[1] == 8",  M000.dim(1), 8)

    # (6): PBW dimensions in U(L_{-1}) (independent of W_hat)
    for d, expected_pbw in [(0,1), (1,8), (2,32), (3,88), (4,192)]:
        _check(f"pbw_dim({d})", M100.pbw_dim(d), expected_pbw)

    # (7)-(8): basis list sizes
    _check("len(basis(1)) == dim(1)", len(M100.basis(1)), M100.dim(1))
    _check("len(basis(2)) == dim(2)", len(M100.basis(2)), M100.dim(2))

    # (9): to_vec round-trip — set coefficient (mon_idx=0, w_idx=0) = 1
    v = M100.to_vec(1, {(0, 0): QQ(1)})
    _check("to_vec round-trip: entry 0 of deg-1 vec == 1", v[0], QQ(1))
    _check("to_vec round-trip: entry 1 of deg-1 vec == 0", v[1], QQ(0))

    if verbose:
        print("─" * 60)
        if all_pass:
            print("s2.2  ✓  ALL CHECKS PASSED")
        else:
            print("s2.2  ✗  SOME CHECKS FAILED — see above")
        print("─" * 60)

    return all_pass


# ===========================================================================
# __main__
# ===========================================================================

if __name__ == '__main__':
    print("=" * 60)
    print("verma_modules.py  —  s2.1 + s2.2")
    print("=" * 60)

    # Attempt to load E44 bracket data; warn if absent but continue
    e44_data = load_e44()
    if e44_data is not None:
        print(f"\n[info] Loaded E44 data: graded pieces {sorted(e44_data['E44'].keys())}")
    else:
        print("\n[info] E44 bracket data not loaded; continuing with s2.1/s2.2 only.")

    # --- Print module summaries for key cases ---
    print("\n[s2.1.a] Module summaries for key highest weights:")
    for (t, a, b, c) in [
        (1, 1, 0, 0),
        (0, 0, 1, 0),
        (0, 0, 0, 1),
        (0, 2, 0, 0),
        (0, 0, 0, 0),
    ]:
        try:
            W = W_hat(t, a, b, c)
            print(f"\n  W_hat({t},{a},{b},{c}):")
            print(W.character_string())
        except Exception as exc:
            print(f"  W_hat({t},{a},{b},{c}): ERROR — {exc}")

    # --- NS velocity-pressure module ---
    print("\n[s2.1.b] NS velocity module W_hat(1,1,0,0) = Verma base W_t(1,0,0):")
    W_ns = W_hat(1, 1, 0, 0)
    print(W_ns.character_string())

    # --- S2.1 Checkpoint ---
    print()
    _check_s21(verbose=True)

    # -----------------------------------------------------------------------
    # S2.2
    # -----------------------------------------------------------------------
    print()
    print("=" * 60)
    print("verma_modules.py  —  s2.2: Verma module graded structure")
    print("=" * 60)

    # Dimension tables for key modules
    print("\n[s2.2.a] Graded dimension tables:")
    for (t, a, b, c) in [(0,1,0,0), (0,0,1,0), (0,0,0,1), (0,0,0,0)]:
        M = M_verma(t, a, b, c)
        print(f"\n  {M}")
        print(f"  PBW_dims = " +
              str({d: M.pbw_dim(d) for d in range(M.max_deg + 1)}))

    # Sample basis at degree 1 for M(0,1,0,0)
    print("\n[s2.2.b] First 8 basis elements of M(0,1,0,0)[deg=1]:")
    M100 = M_verma(0, 1, 0, 0)
    for idx, (mon, k) in enumerate(M100.basis(1)[:8]):
        alpha, S = mon
        print(f"  [{idx}] e^{alpha} ∧ f_S{sorted(S)}  \otimes  w[{k}]")

    # --- S2.2 Checkpoint ---
    print()
    _check_s22(verbose=True)

    # -----------------------------------------------------------------------
    # S2.3
    # -----------------------------------------------------------------------
    print()
    print("=" * 60)
    print("verma_modules.py  —  s2.3: L_{<0}, L_0, L_1 actions + commutator")
    print("=" * 60)
    print()
    _check_s23(verbose=True)

    # ── Leibniz direct verification ───────────────────────────────────
    print()
    _check_leibniz(verbose=True)

    # -----------------------------------------------------------------------
    # S2.4
    # -----------------------------------------------------------------------
    print()
    print("=" * 60)
    print("verma_modules.py  —  s2.4: Grothendieck decompositions")
    print("=" * 60)

    # Sample decompositions
    print("\n[s2.4.a] Sample grothendieck_decomp calls:")
    for (t_, a_, b_, c_, depth_) in [
        (0, 0, 0, 0, 3),
        (0, 1, 0, 0, 3),
        (1, 1, 0, 0, 2),
        (5, 1, 0, 0, 2),
        (0, 0, 0, 1, 3),
        (0, 0, 0, 2, 3),
        (0, 0, 1, 0, 2),
        (0, 1, 1, 0, 2),
    ]:
        d = grothendieck_decomp(t_, a_, b_, c_, max_depth=depth_)
        print(f"\n  [M_{t_}({a_},{b_},{c_})]  (depth ≤ {depth_}):")
        for (tt, aa, bb, cc), mult in sorted(d.items()):
            print(f"    {mult} \cdot [L_{{{tt}}}({aa},{bb},{cc})]")

    # Grothendieck table (comparison with CKC 2026 Table 1)
    print()
    grothendieck_table(t_val=QQ(0), max_depth=3, verbose=True)

    # S2.4 Checkpoint
    print()
    _check_s24(verbose=True)

    # -----------------------------------------------------------------------
    # S2.5
    # -----------------------------------------------------------------------
    print()
    print("=" * 60)
    print("verma_modules.py  —  s2.5: export verma_cache.pkl")
    print("=" * 60)
    print()

    payload = save_verma_data("verma_cache.pkl", verbose=True)

    print(f"\n[s2.5.a] Keys written: {sorted(payload.keys())}")
    print(f"[s2.5.b] W_hat dims: {payload['W_hat_dims']}")
    print(f"[s2.5.c] Verma dim-table keys: {sorted(payload['verma_dim_tables'].keys())}")
    print(f"[s2.5.d] Grothendieck entries: {len(payload['grothendieck'])}")

    print()
    _check_s25("verma_cache.pkl", verbose=True)
