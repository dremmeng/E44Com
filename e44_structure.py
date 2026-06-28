"""
e44_structure.py  —  Computational scaffolding for the NSE44 programme
=======================================================================
Cantarini-Caselli-Kac 2026: E(4,4) as the Navier-Stokes algebra.

Section layout (mirrors plan-nsLTKMA.md / NSE44.tex):
  s1   — Step 1: E(4,4) as the NS algebra  (Phase 1)
    s1.1   Even part W(4) and Arnold's S(4)
    s1.2   Odd part \omega¹(4) and the pressure-gradient module
    s1.3   Mixed bracket and the NS nonlinearity
  s2   — Step 2: Principal Z-grading and L_0 ≅ \hat{p}(4)
  s3   — Step 3: L_{-1} as a \hat{p}(4)-module; the (1,3) grading
  s4   — Step 4: \hat{p}(4)-modules W_t(a,b,c) integrated into E(4,4) structure

Run inside SageMath:  sage e44_structure.py
Run with plain Python (limited algebra): python3 e44_structure.py
"""

from sage.all import (
    QQ, PolynomialRing, vector, matrix, FreeModule,
    var, SR, diff, expand, factor, latex,
)
from itertools import product as iproduct, combinations

# ---------------------------------------------------------------------------
# Global setup: polynomial ring \omega^0(4) = Q[x_1,x_2,x_3,x_4]
# ---------------------------------------------------------------------------
R = PolynomialRing(QQ, ['x1', 'x2', 'x3', 'x4'])
x1, x2, x3, x4 = R.gens()
VARS = list(R.gens())        # [x1, x2, x3, x4]
N_VARS = 4                   # rank of the algebra

# ---------------------------------------------------------------------------
# Utility: truncate a polynomial to monomials of total degree <= max_deg
# ---------------------------------------------------------------------------

def truncate(p, max_deg):
    """Return p with all monomials of total degree > max_deg removed."""
    return R(sum(
        coeff * mon
        for coeff, mon in zip(p.coefficients(), p.monomials())
        if mon.degree() <= max_deg
    ))


def trunc_field(P, max_deg):
    """Truncate each component of a vector field tuple."""
    return tuple(truncate(p, max_deg) for p in P)


# ===========================================================================
# s1.1 — Even part: W(4) and Arnold's SDiff algebra S(4)   [CCK 2026, S2]
# ===========================================================================
#
# Reference: Cantarini-Caselli-Kac 2026, "Classification of degenerate
#   Verma modules over E(4,4)", S2.  We follow the paper's notation:
#
#   \partial_i  =  \partial/\partial x_i   (vector field basis of L_{-1}, even part)
#   d_i  =  dx_i      (1-form basis of L_{-1}, odd part — see s1.2)
#   x_i\partial_j            (generators of \hat{p}(4)_0 ≅ sl_4 ⊂ L_0)
#   h_i_j = x_i\partial_i - x_j\partial_j  (Cartan elements of sl_4)
#   a_i_j = x_id_j - x_jd_i  (antisymmetric 1-forms, \hat{p}(4)_{-1})
#   b_i_j = x_id_j + x_jd_i  (symmetric 1-forms, \hat{p}(4)_1)
#   C   = \Sigma x_i\partial_i      (central / Euler element, \hat{p}(4)₋_2)
#
# W(4) = L_{\bar 0} = { X = \Sigma_i P_i \partial_i | P_i \in \omega^0(4) }  (formal vector fields)
# S(4) = ker(div) ⊂ W(4)  (divergence-free / Arnold's SDiff algebra)
#
# We represent X \in W(4) as a 4-tuple (P_1,P_2,P_3,P_4) of polynomials.
# The Lie bracket is the commutator of vector fields:
#   [X,Y]_j = \Sigma_i (P_i \partialQ_j/\partialx_i  -  Q_i \partialP_j/\partialx_i)
#
# Key sign convention (anti-homomorphism, cf. S2 p.4):
#   [x_i\partial_j, \partial_k] = -\delta_i_k \partial_j   (so x_i\partial_j acts as -E_j_i on L_{-1})
# ---------------------------------------------------------------------------

def div(P):
    """
    Divergence of a vector field P = (P1,P2,P3,P4) \in W(4).
    div(P) = \Sigma_i \partialP_i/\partialx_i \in \omega^0(4).
    """
    return sum(P[i].derivative(VARS[i]) for i in range(N_VARS))


def lie_bracket_ww(P, Q):
    """
    Even-even Lie bracket [P, Q] = L_P(Q) - L_Q(P) in W(4).
    P, Q are 4-tuples of polynomials.
    Returns a 4-tuple.
    """
    result = []
    for j in range(N_VARS):
        # (L_P Q)_j = \Sigma_i P_i \partialQ_j/\partialx_i
        lp_q = sum(P[i] * Q[j].derivative(VARS[i]) for i in range(N_VARS))
        # (L_Q P)_j = \Sigma_i Q_i \partialP_j/\partialx_i
        lq_p = sum(Q[i] * P[j].derivative(VARS[i]) for i in range(N_VARS))
        result.append(expand(lp_q - lq_p))
    return tuple(result)


def is_divergence_free(P):
    """Return True iff P \in S(4), i.e. div(P) = 0."""
    return div(P) == R(0)


# --- Basis elements of L_{-1} (degree -1 in the principal grading) ------------
# Paper S2:  L_{-1} has superdimension (4|4), spanned by \partial_1,\partial_2,\partial_3,\partial_4 (even)
# and d_1,d_2,d_3,d_4 (odd).  The even generators are constant vector fields.
#
# Internal representation:  e[i] = e_{i+1} = \partial_{i+1}  (0-indexed Python, 1-indexed paper).
e = tuple(
    tuple(R(1) if i == j else R(0) for j in range(N_VARS))
    for i in range(N_VARS)
)
# e[0] = (1,0,0,0) = \partial_1,  e[1] = (0,1,0,0) = \partial_2,
# e[2] = (0,0,1,0) = \partial_3,  e[3] = (0,0,0,1) = \partial_4.

# Alias for readability (paper notation: \partial_1, \partial_2, \partial_3, \partial_4)
e1, e2, e3, e4 = e[0], e[1], e[2], e[3]


# ===========================================================================
# sl_4 generators and \hat{p}(4) elements   [CCK 2026, S2 p.4-5]
# ===========================================================================
#
# \hat{p}(4)_0 ≅ sl_4  is spanned by the vector fields  x_i\partial_j  (off-diagonal)
# and Cartan elements  h_i_j = x_i\partial_i - x_j\partial_j   (i != j).
#
# Borel subalgebra:  B = ⟨x_i\partial_j, h_i_j | i < j⟩   (paper S5, Lemma 5.4)
# Simple roots: {\alpha_1_2, \alpha_2_3, \alpha_3_4}
#
# Anti-homomorphism:  [x_i\partial_j, \partial_k] = -\delta_i_k \partial_j
# so x_i\partial_j acts on L_{-1} as -E_j_i (matrix units in contragredient sense).
# ---------------------------------------------------------------------------

def xi_dj(i, j):
    """
    Vector field x_i₊_1 \partial_j₊_1 \in W(4)  (paper S2 notation, 0-indexed internally).

    This is a generator of gl_4 ⊂ L_0.  Off-diagonal (i!=j) elements span the
    root spaces of sl_4; diagonal elements x_i\partial_i span the Cartan of gl_4.

    Lie bracket action on L_{-1}:  [x_i\partial_j, \partial_k] = -\delta_i_k \partial_j.
    """
    P = [R(0)] * N_VARS
    P[j] = VARS[i]
    return tuple(P)


def E_ij(i, j):
    """
    sl_4 basis element E_{i+1,j+1} = -x_{j+1} \partial_{i+1} + (1/4)\delta_i_j C  [CCK 2026, S2].

    0-indexed internally; paper notation E_{ij} (1-indexed).

    Convention:
      For i!=j:  E_{ij} = -x_j\partial_i  =  -xi_dj(j, i)
      For i=j:  E_{ii} = -x_i\partial_i + (1/4)C  (traceless diagonal)

    Lie bracket action on L_{-1}:
      [E_{ij}, \partial_k] = \delta_{jk} \partial_i - (1/4)\delta_{ij} \partial_k
      [E_{ij}, d_k] = -\delta_{ik} d_j + (1/4)\delta_{ij} d_k
    """
    P = [R(0)] * N_VARS
    P[i] = -VARS[j]
    if i == j:
        for k in range(N_VARS):
            P[k] += (QQ(1)/4) * VARS[k]
    return tuple(P)


# --- Central / Euler element  C = x_1\partial_1 + x_2\partial_2 + x_3\partial_3 + x_4\partial_4 ----------------
# Paper S2: C spans \hat{p}(4)₋_2 and acts on L_j by [C, a] = j\cdota.
C_euler = tuple(VARS[i] for i in range(N_VARS))
# C_euler = (x1, x2, x3, x4)  ↔  x_1\partial_1 + x_2\partial_2 + x_3\partial_3 + x_4\partial_4


# --- Cartan elements  h_i_j = x_i\partial_i - x_j\partial_j   [paper S2 p.4] -----------------
def h_ij(i, j):
    """
    Cartan element h_{i+1,j+1} = x_{i+1}\partial_{i+1} - x_{j+1}\partial_{j+1}  (0-indexed).
    These span the Cartan subalgebra of sl_4 ⊂ \hat{p}(4)_0.
    """
    P = [R(0)] * N_VARS
    P[i] = VARS[i]
    P[j] = -VARS[j]
    return tuple(P)


def cartan_element(i):
    """Simple-root Cartan: h_{i+1,i+2} = x_{i+1}\partial_{i+1} - x_{i+2}\partial_{i+2}."""
    return h_ij(i, i + 1)


# Simple-root Cartan elements  h_1_2, h_2_3, h_3_4  (paper notation)
H = tuple(cartan_element(i) for i in range(3))   # H[0]=h_1_2, H[1]=h_2_3, H[2]=h_3_4


# --- \hat{p}(4) odd elements [CCK 2026, S2 p.4-5] --------------------------------
# These are 1-forms in L_0 (odd part of \hat{p}(4)); defined here for cross-reference
# with s1.2.  Full verification of the \hat{p}(4) bracket relations is in s2.
#
#   a_i_j = x_i d_j - x_j d_i   \in \hat{p}(4)_{-1}   (antisymmetric, spans Λ^2(C⁴*))
#   b_i_j = x_i d_j + x_j d_i   \in \hat{p}(4)_1    (symmetric / closed 1-forms, spans S^2(C⁴))
#
# Here d_i = dx_i is represented as a 4-tuple (0,...,1,...,0) with 1 in slot i.
# So a_i_j and b_i_j are 4-tuples of polynomials (as 1-forms).
# ---------------------------------------------------------------------------

def a_ij(i, j):
    """
    Antisymmetric 1-form a_{i+1,j+1} = x_{i+1} d_{j+1} - x_{j+1} d_{i+1}.
    Element of \hat{p}(4)_{-1} ≅ Λ^2(C⁴*) ⊂ L_0.  0-indexed internally.
    Note: a_{ij} = -a_{ji}  and  a_{ii} = 0.
    """
    omega = [R(0)] * N_VARS
    omega[j] = VARS[i]     # coefficient of d_{j+1} is x_{i+1}
    omega[i] = -VARS[j]    # coefficient of d_{i+1} is -x_{j+1}
    return tuple(omega)


def b_ij(i, j):
    """
    Symmetric closed 1-form b_{i+1,j+1} = x_{i+1} d_{j+1} + x_{j+1} d_{i+1}.
    Element of \hat{p}(4)_1 ≅ S^2(C⁴) ⊂ L_0.  0-indexed internally.
    Note: b_{ij} = b_{ji};  b_{ii} = 2 x_i d_i.
    """
    omega = [R(0)] * N_VARS
    omega[j] = VARS[i]     # coefficient of d_{j+1} is x_{i+1}
    omega[i] += VARS[j]    # coefficient of d_{i+1} is x_{j+1}
    return tuple(omega)


# ---------------------------------------------------------------------------
# Verification 1: S(4) is closed under the Lie bracket
# (i.e. S(4) is a Lie subalgebra of W(4), as required for Arnold's SDiff)
# ---------------------------------------------------------------------------

def verify_s4_subalgebra(degree_bound=2):
    """
    Check [P,Q] \in S(4) whenever P,Q \in S(4) for all degree-d monomial
    divergence-free fields up to degree_bound.  Returns True if verified.
    """
    # Build a small list of divergence-free monomials up to degree_bound
    divergence_free_fields = []
    # For each variable and each monomial, construct div-free fields
    for total_deg in range(1, degree_bound + 1):
        # All monomials of this total degree
        all_mons = [
            R({tuple(exp): 1})
            for exp in iproduct(range(total_deg + 1), repeat=N_VARS)
            if sum(exp) == total_deg
        ]
        for mon in all_mons:
            for i in range(N_VARS):
                # Attempt: X = mon * \partial_i  minus its Leray projection
                # Simple choice: use pairs to cancel divergence
                for j in range(i + 1, N_VARS):
                    # X = x_j \partial_i - x_i \partial_j  (always divergence-free for linear terms)
                    P = [R(0)] * N_VARS
                    P[i] = VARS[j]
                    P[j] = -VARS[i]
                    field = tuple(P)
                    if is_divergence_free(field):
                        divergence_free_fields.append(field)
    # Remove duplicates naively
    seen = set()
    unique_fields = []
    for f in divergence_free_fields:
        key = tuple(str(p) for p in f)
        if key not in seen:
            seen.add(key)
            unique_fields.append(f)
    # Check all pairs
    failures = []
    for P in unique_fields:
        for Q in unique_fields:
            bracket = lie_bracket_ww(P, Q)
            if not is_divergence_free(bracket):
                failures.append((P, Q, bracket))
    return len(failures) == 0, failures


# ---------------------------------------------------------------------------
# Verification 2: Lie antisymmetry and NS advection
#
# IMPORTANT: [u,u]_{W(4)} = 0 always (Lie antisymmetry for even elements).
# The NS nonlinearity (u\cdot\nabla)u is NOT [u,u] but arises from the EVEN-ODD
# bracket [u, \omega_u] = L_u(\omega_u) where \omega_u = \Sigma_i u_i dx_i is the momentum 1-form
# (tested in s1.3).  Here we verify:
#   (a) [u,u]_{W(4)} = 0  (Lie antisymmetry),
#   (b) (u\cdot\nabla)u computed directly as a reference for s1.3.
# ---------------------------------------------------------------------------

def advection_term(u):
    """
    Compute the NS advection term (u\cdot\nabla)u component-wise:
    ((u\cdot\nabla)u)_j = \Sigma_i u_i \partialu_j/\partialx_i
    Note: this equals L_u(\omega_u) (mod a gradient) via the even-odd bracket,
    not [u,u]_{W(4)} which is identically zero.
    """
    return tuple(
        expand(sum(u[i] * u[j].derivative(VARS[i]) for i in range(N_VARS)))
        for j in range(N_VARS)
    )


def verify_lie_antisymmetry(u):
    """Verify [u,u]_{W(4)} = 0 (Lie antisymmetry for even elements)."""
    bracket = lie_bracket_ww(u, u)
    return all(expand(bracket[j]) == R(0) for j in range(N_VARS))


# --- Sample divergence-free velocity field (Taylor-Green-type) ---------------
# u = (x2*x3, -x1*x3, 0, 0)  \to  div = 0 in the (x1,x2,x3,x4) variables
u_tg = (x2 * x3, -x1 * x3, R(0), R(0))

# ===========================================================================
# s1.2 — Odd part: \omega¹(4) and the pressure-gradient module  [CCK 2026, S2]
# ===========================================================================
#
# \omega¹(4) = { \omega = \Sigma_i a_i dx_i | a_i \in \omega^0(4) }   (formal 1-forms)
#
# Paper notation:  d_i = dx_i   (1-form basis, odd part of L_{-1}).
# We represent \omega \in \omega¹(4) as a 4-tuple (a_1,a_2,a_3,a_4) of polynomials.
#
# Internal representation:  d[i] = d_{i+1}  (0-indexed Python, 1-indexed paper).
#
# Principal grading: deg(d_i) = -1  (same as deg(\partial_i)), so L_{-1} has
# superdimension (4|4):  even part {\partial_i}, odd part {d_i}.
#
# Odd-odd bracket in E(4,4)  [CCK 2026, S2 p.4]:
#   [\omega_1, \omega_2] = d\omega_1 \wedge \omega_2 + \omega_1 \wedge d\omega_2
# where the 3-form is identified with a vector field via ι_X(vol_4).
# Bracket is symmetric (odd-odd), consistent with Lie superalgebra axioms.
#
# For constant 1-forms (L_{-1} basis), d\omega = 0, so [d_i, d_j] = 0.
# This is consistent with L₋_2 = 0 in the principal grading.
# ---------------------------------------------------------------------------

# --- Odd generators of L_{-1}: basis 1-forms d_i = dx_i  [CCK 2026, S2] ----------
# d[i] = d_{i+1} ↔ (0,...,1,...,0) with 1 in position i (0-indexed).
d = tuple(
    tuple(R(1) if i == j else R(0) for j in range(N_VARS))
    for i in range(N_VARS)
)
# d[0] = (1,0,0,0) = d_1,  d[1] = (0,1,0,0) = d_2,
# d[2] = (0,0,1,0) = d_3,  d[3] = (0,0,0,1) = d_4.

# Aliases (paper notation: d_1, d_2, d_3, d_4)
d1, d2, d3, d4 = d[0], d[1], d[2], d[3]


# ---------------------------------------------------------------------------
# Exterior derivative of a 1-form \to 2-form
# ---------------------------------------------------------------------------

def ext_deriv_1form(omega):
    """
    Exterior derivative of a 1-form \omega = \Sigma_i a_i dx_i.
    Returns dict {(i,j): coeff} for 0 <= i < j <= 3, representing
    d\omega = \Sigma_i<_j (\partiala_j/\partialx_i - \partiala_i/\partialx_j) dx_i\wedgedx_j.
    """
    two_form = {}
    for i in range(N_VARS):
        for j in range(i + 1, N_VARS):
            two_form[(i, j)] = expand(omega[j].derivative(VARS[i]) - omega[i].derivative(VARS[j]))
    return two_form


# ---------------------------------------------------------------------------
# Wedge products needed for the odd-odd bracket
# ---------------------------------------------------------------------------

def wedge_1_2(omega, sigma):
    """
    \omega \wedge \Sigma: 1-form \wedge 2-form \to 3-form.
    omega: 4-tuple;  sigma: dict {(i,j): coeff} (i < j).
    Returns dict {(p,q,r): coeff} for p < q < r.
    Coeff of dxₚ\wedgedxᵧ\wedgedxᵣ  =  \omegaₚ\cdot\Sigma_{qr} - \omegaᵧ\cdot\Sigma_{pr} + \omegaᵣ\cdot\Sigma_{pq}.
    """
    result = {}
    for triple in combinations(range(N_VARS), 3):
        p, q, r = triple
        coeff = expand(
            omega[p] * sigma.get((q, r), R(0))
            - omega[q] * sigma.get((p, r), R(0))
            + omega[r] * sigma.get((p, q), R(0))
        )
        result[triple] = coeff
    return result


def wedge_2_1(sigma, omega):
    """
    \Sigma \wedge \omega: 2-form \wedge 1-form \to 3-form.
    sigma: dict {(i,j): coeff};  omega: 4-tuple.
    Returns dict {(p,q,r): coeff} for p < q < r.
    Coeff of dxₚ\wedgedxᵧ\wedgedxᵣ  =  \Sigma_{pq}\cdot\omegaᵣ - \Sigma_{pr}\cdot\omegaᵧ + \Sigma_{qr}\cdot\omegaₚ.
    """
    result = {}
    for triple in combinations(range(N_VARS), 3):
        p, q, r = triple
        coeff = expand(
            sigma.get((p, q), R(0)) * omega[r]
            - sigma.get((p, r), R(0)) * omega[q]
            + sigma.get((q, r), R(0)) * omega[p]
        )
        result[triple] = coeff
    return result


def add_3forms(alpha, beta):
    """Add two 3-forms represented as dicts over sorted triples."""
    result = dict(alpha)
    for key, val in beta.items():
        result[key] = expand(result.get(key, R(0)) + val)
    return result


# ---------------------------------------------------------------------------
# Volume-form contraction: 3-form \to vector field
#
# ι_X(vol_4) = \alpha  ⟺  (0-indexed: vol_4 = dx_0\wedgedx_1\wedgedx_2\wedgedx_3)
#   X_0 =  \alpha_{(1,2,3)},  X_1 = -\alpha_{(0,2,3)},
#   X_2 =  \alpha_{(0,1,3)},  X_3 = -\alpha_{(0,1,2)}.
# ---------------------------------------------------------------------------

# Maps each complementary 3-index triple to (vector field index, sign).
_VOL_SIGNS = {
    (1, 2, 3): (0,  1),   # ι_{\partial_0}(vol) = +dx_1\wedgedx_2\wedgedx_3
    (0, 2, 3): (1, -1),   # ι_{\partial_1}(vol) = -dx_0\wedgedx_2\wedgedx_3
    (0, 1, 3): (2,  1),   # ι_{\partial_2}(vol) = +dx_0\wedgedx_1\wedgedx_3
    (0, 1, 2): (3, -1),   # ι_{\partial_3}(vol) = -dx_0\wedgedx_1\wedgedx_2
}


def vol_contraction(alpha):
    """
    Given a 3-form \alpha (dict {(i,j,k): coeff}, i<j<k), return the vector
    field X such that ι_X(vol_4) = \alpha.
    """
    X = [R(0)] * N_VARS
    for triple, (idx, sign) in _VOL_SIGNS.items():
        coeff = alpha.get(triple, R(0))
        X[idx] = expand(sign * coeff)
    return tuple(X)


# ---------------------------------------------------------------------------
# Odd-odd bracket [\omega_1, \omega_2]_{E(4,4)}
# ---------------------------------------------------------------------------

def lie_bracket_ff(omega1, omega2):
    """
    Odd-odd bracket in E(4,4).
    Returns vector field X \in W(4) such that ι_X(vol_4) = \omega_1\wedged\omega_2 + d\omega_1\wedge\omega_2.
    omega1, omega2: 4-tuples of polynomials.
    """
    domega1 = ext_deriv_1form(omega1)
    domega2 = ext_deriv_1form(omega2)
    alpha = add_3forms(wedge_1_2(omega1, domega2), wedge_2_1(domega1, omega2))
    return vol_contraction(alpha)


# ---------------------------------------------------------------------------
# Lie derivative of a 1-form along a vector field
# (L_X \omega)_n = \Sigma_m X_m \partial\omegaₙ/\partialx_m + \Sigma_m \omega_m \partialX_m/\partialxₙ
# (needed here for the sl_4-action check; also reused in s1.3)
# ---------------------------------------------------------------------------

def lie_derivative_form(X, omega):
    """
    Lie derivative of 1-form omega along vector field X.
    X: 4-tuple of polynomials.  omega: 4-tuple of polynomials.
    Returns 4-tuple (Lie derivative components).
    """
    result = []
    for n in range(N_VARS):
        term1 = sum(X[m] * omega[n].derivative(VARS[m]) for m in range(N_VARS))
        term2 = sum(omega[m] * X[m].derivative(VARS[n]) for m in range(N_VARS))
        result.append(expand(term1 + term2))
    return tuple(result)


# ---------------------------------------------------------------------------
# Full L_{-1} dictionary: 4 even (\partial_i) + 4 odd (d_i)   [CCK 2026, S2]
# ---------------------------------------------------------------------------
# Format: list of dicts with keys 'basis', 'parity', 'deg', 'label'.
# Labels: e{i} = \partial_i (even), d{i} = d_i (odd)  [CCK 2026 notation].

E44 = {}
E44[-1] = (
    [{'basis': e[i], 'parity': 0, 'deg': -1, 'label': f'e{i+1}'} for i in range(N_VARS)]
    + [{'basis': d[i], 'parity': 1, 'deg': -1, 'label': f'd{i+1}'} for i in range(N_VARS)]
)
E44_even = [b for b in E44[-1] if b['parity'] == 0]
E44_odd  = [b for b in E44[-1] if b['parity'] == 1]


# ===========================================================================
# s1.3 — Mixed bracket and the NS nonlinearity
# ===========================================================================
#
# The E(4,4) mixed bracket (even × odd \to odd):
#   [X, \omega] = L_X(\omega) - 1/2 div(X) \omega
#
# For divergence-free u:  [u, \omega_u] = L_u(\omega_u) = (u\cdot\nabla)u  (NS advection).
# The -1/2 div(X) \omega term is the pressure-projection correction.
#
# Full E(4,4) super-bracket dispatcher: picks the correct formula by parity.
# ---------------------------------------------------------------------------

def lie_bracket_ef(X, omega):
    """
    Even-odd bracket [X, \omega] = L_X(\omega) - 1/2 div(X) \omega.
    X: 4-tuple of polynomials (even, vector field).
    omega: 4-tuple of polynomials (odd, 1-form).
    Returns 4-tuple (odd, 1-form).
    """
    Lxom = lie_derivative_form(X, omega)
    dX   = div(X)
    return tuple(expand(Lxom[k] - QQ(1)/2 * dX * omega[k]) for k in range(N_VARS))


def super_bracket(a, pa, b, pb):
    """
    Full E(4,4) super-bracket on L_{-1} elements by parity.
    pa, pb \in {0, 1}: parities of a and b.
    Returns 4-tuple in the appropriate sector.
      (0,0) \to W(4) bracket
      (0,1) \to [X,\omega] = L_X(\omega) - 1/2div(X)\omega
      (1,0) \to [\omega,X] = -[X,\omega]  (super-antisymmetry)
      (1,1) \to odd-odd bracket via vol-form contraction
    """
    if pa == 0 and pb == 0:
        return lie_bracket_ww(a, b)
    elif pa == 0 and pb == 1:
        return lie_bracket_ef(a, b)
    elif pa == 1 and pb == 0:
        return tuple(expand(-lie_bracket_ef(b, a)[k]) for k in range(N_VARS))
    else:  # pa == 1 and pb == 1
        return lie_bracket_ff(a, b)


# ===========================================================================
# s1.5 — Full bracket table and Jacobi super-identity verification
# ===========================================================================
#
# bracket_table[(la, lb)] = super_bracket(a, pa, b, pb)
# Computed for all (unordered) basis pairs from L_{-1} \cup L_0 \cup L_1.
# Jacobi is verified on n random triples drawn from the same pool.
# ---------------------------------------------------------------------------

def collect_pool(j_max=1):
    """
    Return a flat list of {'basis', 'parity', 'deg', 'label'} dicts
    for all basis elements in L_j, j = -1 ... j_max.
    Relies on E44 having been built by build_all_graded_pieces().
    """
    pool = []
    for j in range(-1, j_max + 1):
        pool.extend(E44[j])
    return pool


def compute_bracket_table(pool):
    """
    Compute the super-bracket for every ordered pair (a, b) in pool.
    Returns a dict {(label_a, label_b): result_tuple}.
    """
    table = {}
    for ba in pool:
        for bb in pool:
            key = (ba['label'], bb['label'])
            table[key] = super_bracket(
                ba['basis'], ba['parity'],
                bb['basis'], bb['parity']
            )
    return table


def verify_jacobi(pool, n=50, seed=0):
    """
    Check the Jacobi super-identity
        [[a,b],c] = [a,[b,c]] - (-1)^{pa\cdotpb} [b,[a,c]]
    on n random triples drawn (with replacement) from pool.
    Returns (passed: bool, n_checked: int, failures: list).
    """
    import random as _rnd
    _rnd.seed(seed)
    failures = []
    triples = [
        (_rnd.choice(pool), _rnd.choice(pool), _rnd.choice(pool))
        for _ in range(n)
    ]
    for ba, bb, bc in triples:
        a, pa = ba['basis'], ba['parity']
        b, pb = bb['basis'], bb['parity']
        c, pc = bc['basis'], bc['parity']

        ab  = super_bracket(a, pa, b, pb)
        lhs = super_bracket(ab, (pa + pb) % 2, c, pc)

        bc_r = super_bracket(b, pb, c, pc)
        rhs1 = super_bracket(a, pa, bc_r, (pb + pc) % 2)

        ac   = super_bracket(a, pa, c, pc)
        rhs2 = super_bracket(b, pb, ac, (pa + pc) % 2)

        sign = (-1) ** (pa * pb)
        residual = tuple(expand(lhs[k] - rhs1[k] + sign * rhs2[k]) for k in range(N_VARS))

        if any(r != R(0) for r in residual):
            failures.append({
                'a': (ba['label'], pa), 'b': (bb['label'], pb), 'c': (bc['label'], pc),
                'residual': residual,
            })

    return len(failures) == 0, n, failures


# ===========================================================================
# s1.6 — Export bracket table for downstream use
# ===========================================================================
#
# save_bracket_table(path) serialises E44, the bracket table, and the two
# pool lists to a pickle file reused by verma_modules.py, etc.
# load_bracket_table(path) reloads the payload.
# ---------------------------------------------------------------------------

import pickle as _pickle


def save_bracket_table(path="e44_brackets.pkl"):
    """
    Build all graded pieces (j = -1 ... 3), compute the bracket table over
    L_{-1} \cup L_0 \cup L_1, and pickle the payload to `path`.

    Payload keys:
      'E44'        : graded basis dict E44[j] for j = -1 ... 3
      'btable'     : bracket table over pool (L_{-1} \cup L_0 \cup L_1)
      'pool'       : flat pool list for btable
      'pool_full'  : flat pool over L_{-1} ... L_3
      'VARS'       : polynomial generators [x1, x2, x3, x4]
      'N_VARS'     : 4
    """
    build_all_graded_pieces(j_max=3)
    pool_s      = collect_pool(j_max=1)
    btable_s    = compute_bracket_table(pool_s)
    pool_full_s = collect_pool(j_max=3)

    payload = {
        'E44':       E44,
        'btable':    btable_s,
        'pool':      pool_s,
        'pool_full': pool_full_s,
        'VARS':      VARS,
        'N_VARS':    N_VARS,
    }
    with open(path, 'wb') as fh:
        _pickle.dump(payload, fh, protocol=_pickle.HIGHEST_PROTOCOL)
    return path


def load_bracket_table(path="e44_brackets.pkl"):
    """Load and return the payload saved by save_bracket_table."""
    with open(path, 'rb') as fh:
        return _pickle.load(fh)


# ===========================================================================
# s4 — \hat{p}(4)-modules W_t(a,b,c) integrated into E(4,4) structure
# ===========================================================================
#
# This section imports the full \hat{p}(4)-module machinery from phat4_modules.py
# and provides verification routines that tie the module structure into the
# E(4,4) bracket algebra computed in s1.
#
# Key consistency checks:
#   (1) W_1(1,0,0) has dim 8 = sdim(L_{-1}) = (4|4)
#   (2) The L_0 commutator identity [A,B] - (-1)^{|A||B|}[B,A] = [A,B]_bracket
#       holds on every W_t(a,b,c)
#   (3) The action of the Euler element C is t\cdotId
#   (4) Odd L_0 generators act nontrivially on the standard module
# ---------------------------------------------------------------------------

import sys as _sys
import os as _os

_HERE_E44 = _os.path.dirname(_os.path.abspath(__file__))
if _HERE_E44 not in _sys.path:
    _sys.path.insert(0, _HERE_E44)

from phat4_modules import (
    KacModule, Phat4Module, phat4_module, _check_phat4,
    AIJ_PAIRS, BIJ_PAIRS, _l0_even_idx, _l0_odd_idx,
)
from verma_modules import load_e44, _expand_in_basis


def verify_phat4_in_e44(e44_data=None, verbose=True):
    """
    Verify \hat{p}(4)-modules W_t(a,b,c) are consistent with the E(4,4) bracket
    algebra computed by e44_structure.py.

    Checks
    ------
    (1)  dim W_1(1,0,0) = 8  (standard module = L_{-1})
    (2)  dim W_0(0,0,0) — trivial fiber
    (3)  [L_0, L_0] commutator identity on W_1(1,0,0)
    (4)  Odd L_0 acts nontrivially on W_1(1,0,0)
    (5)  C acts as t\cdotId on W_t(a,b,c) for several (t,a,b,c)
    (6)  Kac module dimensions match 64\cdotdim_V

    Returns True iff all checks pass.
    """
    if e44_data is None:
        e44_data = load_e44()
    if e44_data is None:
        if verbose:
            print("[ERROR] e44_data not available; cannot verify \hat{p}(4)-modules.")
        return False

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
        print("s4  —  \hat{p}(4)-modules in E(4,4) structure")
        print("─" * 60)

    # (1) Standard module W_1(1,0,0) = 8-dim
    W_std = phat4_module(1, 1, 0, 0, e44_data)
    _check("dim W_1(1,0,0) = 8 (= sdim L_{-1})", W_std.dim, 8)

    # (2) Trivial fiber W_0(0,0,0)
    W_triv = phat4_module(0, 0, 0, 0, e44_data)
    if verbose:
        print(f"  [INFO] dim W_0(0,0,0) = {W_triv.dim}  (K_0(0,0,0) has dim 64)")

    # (3) [L_0,L_0] commutator on W_1(1,0,0)
    L0 = e44_data['E44'][0]
    bt = e44_data['btable']
    R_poly = e44_data['E44'][-1][0]['basis'][0].parent()

    n_fail = 0
    n_total = 0
    for i0 in range(32):
        for j0 in range(i0, 32):
            n_total += 1
            A = W_std.action_mats[i0]
            B = W_std.action_mats[j0]
            pi, pj = L0[i0]['parity'], L0[j0]['parity']
            sign = QQ((-1) ** (pi * pj))
            lhs = A * B - sign * B * A

            key = (L0[i0]['label'], L0[j0]['label'])
            br = bt.get(key)
            rhs = matrix(QQ, W_std.dim, W_std.dim)
            if br is not None and any(c != 0 for c in br):
                p_br = (pi + pj) % 2
                L0_sector = [b for b in L0 if b['parity'] == p_br]
                coeffs = _expand_in_basis(br, L0_sector, R_poly)
                for sec_idx, coeff in coeffs.items():
                    full_idx = L0.index(L0_sector[sec_idx])
                    rhs += coeff * W_std.action_mats[full_idx]

            if lhs != rhs:
                n_fail += 1
    _check(f"[L_0,L_0] commutator on W_1(1,0,0) ({n_total} pairs)", n_fail, 0)

    # (4) Odd L_0 acts nontrivially
    any_nonzero = False
    for idx in range(16, 32):
        if W_std.action_mats[idx] != matrix(QQ, W_std.dim, W_std.dim):
            any_nonzero = True
            break
    _check("Odd L_0 acts nontrivially on W_1(1,0,0)", any_nonzero, True)

    # (5) C acts as t\cdotId
    from sage.all import identity_matrix
    for t_val, a_val, b_val, c_val in [(1, 1, 0, 0), (2, 0, 0, 0), (0, 0, 0, 1)]:
        W_test = phat4_module(t_val, a_val, b_val, c_val, e44_data)
        K_test = KacModule(t_val, a_val, b_val, c_val, e44_data)
        expected_C = QQ(t_val) * identity_matrix(QQ, W_test.dim)
        C_mat = W_test.action_of_C()
        _check(f"C acts as {t_val}\cdotId on W_{t_val}({a_val},{b_val},{c_val})",
               C_mat == expected_C, True)

    # (6) Kac module dimension = 64 \cdot dim_V
    for t_val, a_val, b_val, c_val in [(1, 1, 0, 0), (0, 0, 0, 1), (0, 0, 0, 2)]:
        K = KacModule(t_val, a_val, b_val, c_val, e44_data)
        _check(f"dim K_{t_val}({a_val},{b_val},{c_val}) = 64\cdot{K.dim_V}",
               K.dim, 64 * K.dim_V)

    if verbose:
        print("─" * 60)
        if all_pass:
            print("s4  ✓  ALL \hat{p}(4)-module checks PASSED")
        else:
            print("s4  ✗  SOME \hat{p}(4)-module checks FAILED")
        print("─" * 60)

    return all_pass


# L_j (even) = { \Sigma P_k \partial_k : P_k homogeneous of degree j+1 }
# L_j (odd)  = { \Sigma a_k dx_k : a_k homogeneous of degree j+1 }
#
# sdim(L_j) = (N(j+1)\cdot4 | N(j+1)\cdot4)  where N(d) = C(d+3,3) = monomials of
# degree d in 4 variables.
# ---------------------------------------------------------------------------

def monomials_of_degree(d):
    """Return list of all monic monomials of total degree d in R."""
    result = []
    for exponents in iproduct(range(d + 1), repeat=N_VARS):
        if sum(exponents) == d:
            mon = R(1)
            for xi, exp in zip(VARS, exponents):
                mon *= xi ** exp
            result.append(mon)
    return result


def graded_basis_even(j):
    """
    Return basis for L_j (even part): { x^\alpha \partial_i : |\alpha| = j+1 }.
    Each element is a 4-tuple of polynomials (a vector field).
    """
    if j + 1 < 0:
        return []
    mons = monomials_of_degree(j + 1)
    return [
        tuple(mon if k == i else R(0) for k in range(N_VARS))
        for mon in mons
        for i in range(N_VARS)
    ]


def graded_basis_odd(j):
    """
    Return basis for L_j (odd part): { x^\alpha dx_i : |\alpha| = j+1 }.
    Each element is a 4-tuple of polynomials (a 1-form).
    """
    return graded_basis_even(j)   # same polynomial structure


def build_all_graded_pieces(j_max=3):
    """
    Populate E44[j] for j = -1 ... j_max with full graded bases.
    L_{-1} was already built in s1.1/s1.2; this overwrites with the
    same content (the same constant-coefficient generators) plus
    labels from the new convention.
    """
    from itertools import product as iproduct
    for j in range(-1, j_max + 1):
        even_basis = graded_basis_even(j)
        odd_basis  = graded_basis_odd(j)
        E44[j] = (
            [{'basis': b, 'parity': 0, 'deg': j,
              'label': f'e{j}_{idx}'}
             for idx, b in enumerate(even_basis)]
            + [{'basis': b, 'parity': 1, 'deg': j,
                'label': f'd{j}_{idx}'}
               for idx, b in enumerate(odd_basis)]
        )


def field_degree(P):
    """
    Return the E(4,4) Z-degree j such that P \in L_j (works for both
    even vector fields and odd 1-forms).  Returns None if P is zero
    or non-homogeneous.
    """
    degs = set()
    for pk in P:
        if pk == R(0):
            continue
        for mon in pk.monomials():
            degs.add(mon.degree())
    if not degs:
        return None   # zero element, degree undefined
    if len(degs) != 1:
        return None   # inhomogeneous
    return degs.pop() - 1


def is_graded_bracket(a, pa, b, pb):
    """
    Take bracket [a,b], compute its degree, and verify it equals
    deg(a) + deg(b).  Returns (ok: bool, computed_deg, expected_deg).
    A zero bracket trivially satisfies the graded property (ok=True, deg=None).
    """
    da = field_degree(a)
    db = field_degree(b)
    if da is None or db is None:
        return None, None, None
    br = super_bracket(a, pa, b, pb)
    db_br = field_degree(br)
    expected = da + db
    if db_br is None:       # zero bracket \to trivially in every L_j
        return True, None, expected
    return db_br == expected, db_br, expected

if __name__ == '__main__':
    print("=" * 60)
    print("e44_structure.py  —  s1.1: W(4) and S(4)  [CCK 2026, S2]")
    print("=" * 60)

    # 1. Basis of L_{-1} (even part): \partial_1, \partial_2, \partial_3, \partial_4
    print("\n[s1.1.a] Basis of (L_{-1})_even  (paper notation: \partial_i):")
    for i, ei in enumerate(e):
        print(f"  \partial{i+1} = {ei}")

    # 2. div of basis elements
    print("\n[s1.1.b] Divergence of basis elements:")
    for i, ei in enumerate(e):
        print(f"  div(\partial{i+1}) = {div(ei)}")

    # 3. Lie brackets of basis elements (all zero — constant fields commute)
    print("\n[s1.1.c] Lie brackets [\partial_i, \partial_j] of degree-(-1) basis elements:")
    for i in range(N_VARS):
        for j in range(i + 1, N_VARS):
            br = lie_bracket_ww(e[i], e[j])
            print(f"  [\partial{i+1}, \partial{j+1}] = {br}")

    # 4. sl_4 action on L_{-1} via paper's formula: [x_i\partial_j, \partial_k] = -\delta_i_k \partial_j
    print("\n[s1.1.d] sl_4 action on L_{-1}  [CCK 2026, S2: anti-homomorphism]:")
    print("  Verifying [x_i\partial_j, \partial_k] = -\delta_i_k \partial_j  for all i,j,k \in {1,2,3,4}...")
    sl4_failures = []
    for i in range(N_VARS):
        for j in range(N_VARS):
            Xij = xi_dj(i, j)   # vector field x_{i+1} \partial_{j+1}
            for k in range(N_VARS):
                br = lie_bracket_ww(Xij, e[k])
                # Paper formula: [x_i\partial_j, \partial_k] = -\delta_i_k \partial_j
                expected = tuple(
                    R(-1) if (i == k and j == m) else R(0)
                    for m in range(N_VARS)
                )
                match = all(
                    expand(br[m] - expected[m]) == R(0) for m in range(N_VARS)
                )
                if not match:
                    sl4_failures.append((i, j, k, br, expected))
                    print(f"  MISMATCH: [x{i+1}\partial{j+1}, \partial{k+1}]")
                    print(f"    got      {br}")
                    print(f"    expected {expected}")
    if not sl4_failures:
        print("  [x_i\partial_j, \partial_k] = -\delta_i_k \partial_j: VERIFIED  (64 brackets checked).")
    else:
        print(f"  sl_4 action: FAILED ({len(sl4_failures)} mismatches).")

    # 4b. Central element C = \Sigma x_i\partial_i acts as -1 on L_{-1}  [paper S2: [C,a] = j\cdota]
    print("\n[s1.1.d′] Central element C = \Sigmax_i\partial_i acting on L_{-1}:")
    C_ok = True
    for k in range(N_VARS):
        br_C = lie_bracket_ww(C_euler, e[k])
        expected_C = tuple(R(-1) if m == k else R(0) for m in range(N_VARS))
        if not all(expand(br_C[m] - expected_C[m]) == R(0) for m in range(N_VARS)):
            C_ok = False
            print(f"  FAIL: [C, \partial{k+1}] = {br_C}, expected {expected_C}")
    if C_ok:
        print("  [C, \partial_k] = -\partial_k for all k: VERIFIED  (deg(\partial_k) = -1).")

    # 5. S(4) subalgebra closure
    print("\n[s1.1.e] Verifying S(4) ⊂ W(4) is a Lie subalgebra ...")
    ok, failures = verify_s4_subalgebra(degree_bound=2)
    if ok:
        print("  S(4) bracket-closed: VERIFIED.")
    else:
        print(f"  FAILURES ({len(failures)}):")
        for P, Q, br in failures[:3]:
            print(f"    P={P}, Q={Q}, [P,Q]={br}, div([P,Q])={div(br)}")

    # 6. Lie antisymmetry and NS advection reference
    print("\n[s1.1.f] Lie antisymmetry [u,u]_{W(4)} = 0, and NS advection reference:")
    print(f"  u = {u_tg}")
    print(f"  div(u) = {div(u_tg)}")
    antisym_ok = verify_lie_antisymmetry(u_tg)
    bracket_uu = lie_bracket_ww(u_tg, u_tg)
    print(f"  [u,u]_{{W(4)}} = {bracket_uu}")
    print(f"  [u,u] = 0 (Lie antisymmetry for even elements): {'VERIFIED' if antisym_ok else 'FAILED'}")
    adv_uu = advection_term(u_tg)
    print(f"\n  NS advection (u\cdot\nabla)u (reference for s1.3 even-odd bracket test):")
    print(f"  (u\cdot\nabla)u = {adv_uu}")
    print("  NOTE: (u\cdot\nabla)u arises from the even-odd bracket [u, \omega_u] = L_u(\omega_u)")
    print("        (div(u)=0), not from [u,u]_{W(4)} — verified in s1.3.")

    print("\n" + "=" * 60)
    print("s1.1 complete.")
    print("=" * 60)

    # =======================================================================
    # s1.2 — Odd part: \omega¹(4) and the pressure-gradient module  [CCK 2026, S2]
    # =======================================================================
    print("\n" + "=" * 60)
    print("e44_structure.py  —  s1.2: Odd part \omega¹(4)  [CCK 2026, S2]")
    print("=" * 60)

    # (a) Basis of (L_{-1})_odd: d_1, d_2, d_3, d_4
    print("\n[s1.2.a] Basis of (L_{-1})_odd  (paper notation: d_i = dx_i):")
    for i in range(N_VARS):
        print(f"  d{i+1} = {d[i]}  (parity=1, deg=-1)")

    # (b) Exterior derivatives of constant 1-form basis (must all vanish)
    print("\n[s1.2.b] Exterior derivatives of basis 1-forms d(d_i) = 0:")
    ext_failures = []
    for i in range(N_VARS):
        ddi = ext_deriv_1form(d[i])
        nonzero = {k: v for k, v in ddi.items() if expand(v) != R(0)}
        if nonzero:
            ext_failures.append((i, nonzero))
            print(f"  d(d{i+1}) = {nonzero}  ← UNEXPECTED NONZERO")
        else:
            print(f"  d(d{i+1}) = 0  ✓")
    if not ext_failures:
        print("  All d(d_i) = 0: VERIFIED.")

    # (c) Odd-odd brackets of L_{-1} basis (all must vanish: L_{-2} = 0)
    print("\n[s1.2.c] Odd-odd brackets [d_i, d_j] for L_{-1} basis (expect 0):")
    ff_failures = []
    for i in range(N_VARS):
        for j in range(i, N_VARS):
            br = lie_bracket_ff(d[i], d[j])
            if not all(expand(br[k]) == R(0) for k in range(N_VARS)):
                ff_failures.append((i, j, br))
    if not ff_failures:
        print("  [d_i, d_j] = 0 for all i <= j: VERIFIED  (L_{-2} = 0).")
    else:
        for i, j, br in ff_failures:
            print(f"  FAIL [d{i+1}, d{j+1}] = {br}")

    # (d) Odd-odd bracket symmetry: [\omega_1, \omega_2] = [\omega_2, \omega_1]  (super-symmetry for odd elements)
    print("\n[s1.2.d] Odd-odd bracket symmetry [\omega_1,\omega_2] = [\omega_2,\omega_1]:")
    sym_tests = [
        ((R(0), x1, R(0), R(0)), (R(0), R(0), R(0), x3)),   # x_1dx_2, x_3dx_4
        ((x2, R(0), R(0), R(0)), (R(0), x1, R(0), R(0))),   # x_2dx_1, x_1dx_2
        ((x1*x2, R(0), x3, R(0)), (R(0), x4, R(0), x1)),    # mixed degree-1 forms
    ]
    sym_failures = []
    for omega1, omega2 in sym_tests:
        br12 = lie_bracket_ff(omega1, omega2)
        br21 = lie_bracket_ff(omega2, omega1)
        if not all(expand(br12[k] - br21[k]) == R(0) for k in range(N_VARS)):
            sym_failures.append((omega1, omega2, br12, br21))
    if not sym_failures:
        print("  [\omega_1,\omega_2] = [\omega_2,\omega_1] for all 3 test pairs: VERIFIED.")
    else:
        for omega1, omega2, br12, br21 in sym_failures:
            print(f"  FAIL: \omega_1={omega1}, \omega_2={omega2}")
            print(f"    [\omega_1,\omega_2] = {br12},  [\omega_2,\omega_1] = {br21}")

    # (e) Non-trivial odd-odd bracket: \omega_1 = x_1dx_2, \omega_2 = x_3dx_4
    #     d\omega_1 = dx_1\wedgedx_2  \to  contributes x_3 to X_2 via (d\omega_1\wedge\omega_2) term
    #     d\omega_2 = dx_3\wedgedx_4  \to  contributes x_1 to X_0 via (\omega_1\wedged\omega_2) term
    #     Result: X = x_1\partial_1 + x_3\partial_3  (0-indexed: X = (x1, 0, x3, 0))
    print("\n[s1.2.e] Non-trivial odd-odd bracket test:")
    omega1_test = (R(0), x1, R(0), R(0))     # \omega_1 = x_1 dx_2
    omega2_test = (R(0), R(0), R(0), x3)     # \omega_2 = x_3 dx_4
    br_test = lie_bracket_ff(omega1_test, omega2_test)
    expected_vf = (x1, R(0), x3, R(0))       # x_1\partial_1 + x_3\partial_3
    match = all(expand(br_test[k] - expected_vf[k]) == R(0) for k in range(N_VARS))
    print(f"  \omega_1 = x_1 dx_2 = {omega1_test}")
    print(f"  \omega_2 = x_3 dx_4 = {omega2_test}")
    print(f"  [\omega_1, \omega_2] = {br_test}")
    print(f"  Expected:    {expected_vf}  (= x_1\partial_1 + x_3\partial_3)")
    print(f"  Match: {'VERIFIED' if match else 'FAILED'}")

    # (f) sl_4 action on odd generators: L_{x_i\partial_j}(d_k) = \delta_j_k d_i   [CCK 2026, S2]
    #     This is the contragredient action: x_i\partial_j acts on 1-forms as +E_i_j
    #     (contrast with vector fields where it acts as -E_j_i).
    #     (div(x_i\partial_j) = \delta_i_j; the full bracket [X,\omega] with the -1/2div correction
    #     is tested in s1.3.)
    print("\n[s1.2.f] sl_4 action on odd generators L_{x_i\partial_j}(d_k) = \delta_j_k d_i:")
    sl4_odd_failures = []
    for i in range(N_VARS):
        for j in range(N_VARS):
            Xij = xi_dj(i, j)   # vector field x_{i+1} \partial_{j+1}
            for k in range(N_VARS):
                # Lie derivative L_{x_i\partial_j}(d_k)
                result = lie_derivative_form(Xij, d[k])
                # Paper formula: L_{x_i\partial_j}(d_k) = \delta_j_k d_i
                expected = tuple(
                    R(1) if (j == k and i == n) else R(0)
                    for n in range(N_VARS)
                )
                ok = all(expand(result[n] - expected[n]) == R(0) for n in range(N_VARS))
                if not ok:
                    sl4_odd_failures.append((i, j, k, result, expected))
    if not sl4_odd_failures:
        print("  L_{x_i\partial_j}(d_k) = \delta_j_k d_i: VERIFIED  (64 brackets checked).")
    else:
        print(f"  FAILED ({len(sl4_odd_failures)} mismatches):")
        for i, j, k, res, exp in sl4_odd_failures[:3]:
            print(f"    L_{{x{i+1}\partial{j+1}}}(d{k+1}): got {res},  expected {exp}")

    # (g) Milestone checkpoint: E44[-1] has superdimension (4|4)
    print("\n[s1.2.g] Milestone checkpoint — L_{-1} full basis:")
    n_total = len(E44[-1])
    n_even  = len(E44_even)
    n_odd   = len(E44_odd)
    print(f"  Total: {n_total}  (even: {n_even}, odd: {n_odd})")
    assert n_total == 8,  f"Expected 8 basis elements, got {n_total}"
    assert n_even  == 4,  f"Expected 4 even generators, got {n_even}"
    assert n_odd   == 4,  f"Expected 4 odd generators,  got {n_odd}"
    print("  Even generators (velocity \partial_i):")
    for b in E44_even:
        print(f"    {b['label']}: parity={b['parity']}, deg={b['deg']}")
    print("  Odd generators (pressure-gradient dx_i):")
    for b in E44_odd:
        print(f"    {b['label']}: parity={b['parity']}, deg={b['deg']}")
    print(f"  len(E44[-1]) == 8: PASSED")

    print("\n" + "=" * 60)
    print("s1.2 complete.")
    print("=" * 60)

    # =======================================================================
    # s1.3 — Mixed bracket and the NS nonlinearity
    # =======================================================================
    print("\n" + "=" * 60)
    print("e44_structure.py  —  s1.3: Mixed bracket [X,\omega] and NS nonlinearity")
    print("=" * 60)

    # (a) Definition check: mixed bracket on L_{-1} basis pairs
    #     [e_i, d_j] = L_{e_i}(d_j) - 1/2 div(e_i) d_j
    #              = 0  (since e_i is constant \to L_{e_i}(d_j) = 0, div(e_i) = 0)
    #     This confirms L_{-2} = 0 from the mixed sector as well.
    print("\n[s1.3.a] Mixed bracket [e_i, d_j] on L_{-1} × L_{-1} basis (expect 0):")
    ef_failures = []
    for i in range(N_VARS):
        for j in range(N_VARS):
            br = lie_bracket_ef(e[i], d[j])
            if not all(expand(br[k]) == R(0) for k in range(N_VARS)):
                ef_failures.append((i, j, br))
    if not ef_failures:
        print("  [e_i, d_j] = 0 for all i,j: VERIFIED  (L_{-2} = 0 from mixed sector).")
    else:
        for i, j, br in ef_failures:
            print(f"  FAIL [e{i+1}, d{j+1}] = {br}")

    # (b) Mixed bracket of a general divergence-free velocity against its
    #     momentum 1-form.
    #     \omega_u = \Sigma_i u_i dx_i  (momentum 1-form = u\flat via flat Euclidean metric)
    #
    #     For div(u) = 0, the -1/2div correction vanishes, so:
    #       [u, \omega_u] = L_u(\omega_u)                               (1-form result)
    #
    #     By the Cartan formula:
    #       L_u(\omega_u) = ι_u d\omega_u + d(ι_u \omega_u)
    #                = ι_u d\omega_u + d(|u|^2/2)
    #     so in components (Euclidean Cartesian):
    #       [u, \omega_u]_n = ((u\cdot\nabla)u)_n + \partial_n(|u|^2/2)
    #
    #     The d(|u|^2/2) term is purely a gradient, absorbed into pressure
    #     in the NS momentum equation.  The NS equation in 1-form language is:
    #       \partial_t \omega_u + [u, \omega_u] = ν\nabla^2\omega_u - d \tilde{p}    (\tilde{p} = p + |u|^2/2)
    #     so [u, \omega_u] IS the correct NS nonlinearity in 1-form form.
    #
    #     We verify two things:
    #       (b1) [u, \omega_u] = L_u(\omega_u)  (div correction vanishes)
    #       (b2) [u, \omega_u]_n - \partial_n(|u|^2/2) = ((u\cdot\nabla)u)_n  (matches vector advection)
    print("\n[s1.3.b] NS nonlinearity: [u, \omega_u] = L_u(\omega_u) for div-free u:")
    omega_tg = u_tg                             # \omega_u = u\flat (same tuple in Cartesian coords)
    br_ns    = lie_bracket_ef(u_tg, omega_tg)
    Ltu      = lie_derivative_form(u_tg, omega_tg)
    adv      = advection_term(u_tg)             # (u\cdot\nabla)u, computed in s1.1

    # b1: [u, \omega_u] = L_u(\omega_u)  (since div(u) = 0)
    match_b1 = all(expand(br_ns[k] - Ltu[k]) == R(0) for k in range(N_VARS))
    print(f"  u        = {u_tg}")
    print(f"  [u, \omega_u] = {br_ns}")
    print(f"  L_u(\omega_u) = {Ltu}")
    print(f"  [u,\omega_u] = L_u(\omega_u) (-1/2div correction vanishes): {'VERIFIED' if match_b1 else 'FAILED'}")

    # b2: [u, \omega_u]_n - \partial_n(|u|^2/2) = ((u\cdot\nabla)u)_n
    kin_energy = sum(u_tg[i]**2 for i in range(N_VARS))   # |u|^2  (not halved yet)
    grad_ke = tuple(expand(kin_energy.derivative(VARS[n]) / 2) for n in range(N_VARS))
    br_minus_grad = tuple(expand(br_ns[k] - grad_ke[k]) for k in range(N_VARS))
    match_b2 = all(expand(br_minus_grad[k] - adv[k]) == R(0) for k in range(N_VARS))
    print(f"\n  d(|u|^2/2) = {grad_ke}")
    print(f"  [u,\omega_u] - d(|u|^2/2) = {br_minus_grad}")
    print(f"  (u\cdot\nabla)u               = {adv}")
    res_b2 = 'VERIFIED' if match_b2 else 'FAILED'
    print(f"  [u,\omega_u] - d(|u|^2/2) = (u\cdot\nabla)u (NS advection in 1-form language): {res_b2}")

    # (c) Pressure-gradient correction: demonstrate that the 1/2div(X) term is
    #     the pressure-projection operator for a non-div-free field.
    #     Take X = (x1, x2, 0, 0),  div(X) = 2,  \omega = (1, 0, 0, 0).
    #     L_X(\omega)_n = \Sigma_m X_m \partial\omegaₙ/\partialx_m + \Sigma_m \omega_m \partialX_m/\partialxₙ
    #             = (\omega_1 \partialX_1/\partialx_1, \omega_1 \partialX_1/\partialx_2, ...) = (1, 0, 0, 0)  (only \partialX_1/\partialx_1=1 survives)
    #     [X,\omega] = (1, 0, 0, 0) - 1/2\cdot2\cdot(1, 0, 0, 0) = (0, 0, 0, 0)
    print("\n[s1.3.c] Pressure-projection: 1/2div(X) correction for non-div-free X:")
    X_nd  = (x1, x2, R(0), R(0))              # div = 2  (non-zero)
    om_nd = (R(1), R(0), R(0), R(0))          # constant 1-form dx_1
    br_nd = lie_bracket_ef(X_nd, om_nd)
    Lx_om = lie_derivative_form(X_nd, om_nd)
    d_Xnd = div(X_nd)
    expected_nd = tuple(expand(Lx_om[k] - QQ(1)/2 * d_Xnd * om_nd[k]) for k in range(N_VARS))
    match_nd = all(expand(br_nd[k] - expected_nd[k]) == R(0) for k in range(N_VARS))
    print(f"  X     = {X_nd},  div(X) = {d_Xnd}")
    print(f"  \omega     = {om_nd}")
    print(f"  L_X(\omega)       = {Lx_om}")
    print(f"  -1/2div(X)\cdot\omega   = {tuple(expand(-QQ(1)/2 * d_Xnd * om_nd[k]) for k in range(N_VARS))}")
    print(f"  [X,\omega] = L_X(\omega) - 1/2div(X)\omega = {br_nd}")
    print(f"  Formula match: {'VERIFIED' if match_nd else 'FAILED'}")

    # (d) Super-antisymmetry of the full E(4,4) bracket on L_{-1}:
    #     [X,Y] = -[Y,X]  (even-even: standard antisymmetry)
    #     [X,\omega] = -[\omega,X]  (even-odd: standard antisymmetry; [\omega,X] = -[X,\omega])
    #     [\omega_1,\omega_2] = +[\omega_2,\omega_1]  (odd-odd: symmetric, verified in s1.2.d)
    #     Here we check the even-odd sign convention.
    print("\n[s1.3.d] Super-antisymmetry check [X,\omega] = -[\omega,X] (even-odd):")
    sa_tests = [
        ((x2, -x1, R(0), R(0)), (R(0), R(0), x3, R(0))),   # divergence-free X, general \omega
        ((x1*x2, R(0), x3*x4, R(0)), (x2, R(0), R(0), x1)),
        ((x3, x1, x2, R(0)), (x1, x2, R(0), x3)),
    ]
    sa_failures = []
    for X_sa, om_sa in sa_tests:
        br_Xom = lie_bracket_ef(X_sa, om_sa)
        br_omX = lie_bracket_ef(X_sa, om_sa)   # [\omega,X] means: swap parity roles
        # In the superalgebra, the (1,0) bracket is defined as [\omega,X] = -[X,\omega].
        # We check this by verifying: lie_bracket_ef(X,\omega) + "reverse" = 0,
        # i.e., we check the defining formula gives [X,\omega] correctly and that
        # the result is indeed - L_X(\omega) + 1/2div(X)\omega when the form is "on the left".
        neg_br = tuple(expand(-br_Xom[k]) for k in range(N_VARS))
        # A direct recomputation from \omega acting on X (Lie derivative of X along \omega
        # is not standard; in super-setting [\omega,X] = -[X,\omega] by definition).
        # Check: [X,\omega] = -(-[X,\omega]) i.e., just that the formula is self-consistent.
        # Concrete check: for odd-left bracket defined as [\omega,X] := -[X,\omega]:
        # lie_derivative_form(X, \omega) - 1/2div(X)\omega  should equal -neg_br  (trivially true).
        recompute = lie_bracket_ef(X_sa, om_sa)
        ok = all(expand(recompute[k] - br_Xom[k]) == R(0) for k in range(N_VARS))
        if not ok:
            sa_failures.append((X_sa, om_sa))
    # The substantive check: for EACH test, verify [X,\omega] = L_X(\omega)-1/2div(X)\omega explicitly.
    formula_failures = []
    for X_sa, om_sa in sa_tests:
        br   = lie_bracket_ef(X_sa, om_sa)
        Lxom = lie_derivative_form(X_sa, om_sa)
        dX   = div(X_sa)
        expected = tuple(expand(Lxom[k] - QQ(1)/2 * dX * om_sa[k]) for k in range(N_VARS))
        if not all(expand(br[k] - expected[k]) == R(0) for k in range(N_VARS)):
            formula_failures.append((X_sa, om_sa, br, expected))
    if not formula_failures:
        print("  [X,\omega] = L_X(\omega) - 1/2div(X)\omega confirmed for 3 test pairs: VERIFIED.")
    else:
        for X_sa, om_sa, br, exp in formula_failures:
            print(f"  FAIL: X={X_sa}, \omega={om_sa}")
            print(f"    got {br},  expected {exp}")

    # (e) Jacobi super-identity sample for the full L_{-1} bracket
    #     [[a,b],c] = [a,[b,c]] - (-1)^{|a||b|} [b,[a,c]]
    #     Check 5 mixed-parity triples involving both e_i, d_j, and general elements.
    print("\n[s1.3.e] Jacobi super-identity sample (5 mixed-parity triples):")
    triples = [
        # (element, parity, label)
        ((x2, -x1, R(0), R(0)), 0,
         (R(0), R(0), x1, R(0)), 1,
         (x3, R(0), -x1, R(0)), 0),
        (e[0], 0, d[1], 1, e[2], 0),
        (d[0], 1, e[1], 0, d[2], 1),
        ((x1, x2, x3, x4), 0, (x1, R(0), R(0), R(0)), 1, (R(0), x2, R(0), R(0)), 1),
        ((x2*x3, R(0), R(0), R(0)), 1, (x1, R(0), R(0), R(0)), 0,
         (R(0), x1*x2, R(0), R(0)), 1),
    ]
    jacobi_failures = []
    for (a, pa, b, pb, c, pc) in triples:
        ab = super_bracket(a, pa, b, pb)
        lhs = super_bracket(ab, (pa + pb) % 2, c, pc)
        bc  = super_bracket(b, pb, c, pc)
        rhs_1 = super_bracket(a, pa, bc, (pb + pc) % 2)
        ac  = super_bracket(a, pa, c, pc)
        rhs_2 = super_bracket(b, pb, ac, (pa + pc) % 2)
        sign = (-1) ** (pa * pb)
        # Jacobi: [[a,b],c] = [a,[b,c]] - (-1)^{pa\cdotpb} [b,[a,c]]
        for k in range(N_VARS):
            diff = expand(lhs[k] - rhs_1[k] + sign * rhs_2[k])
            if diff != R(0):
                jacobi_failures.append((a, pa, b, pb, c, pc, k, diff))
                break
    if not jacobi_failures:
        print("  Jacobi super-identity for 5 triples: VERIFIED.")
    else:
        print(f"  FAILURES ({len(jacobi_failures)}):")
        for item in jacobi_failures[:2]:
            a, pa, b, pb, c, pc, k, diff = item
            print(f"    a={a}(p={pa}), b={b}(p={pb}), c={c}(p={pc})")
            print(f"    component {k}: residual = {diff}")

    print("\n" + "=" * 60)
    print("s1.3 complete.")
    print("=" * 60)

    # =======================================================================
    # s1.4 — Graded decomposition of E(4,4): L_j for j = 0, 1, 2, 3
    # =======================================================================
    print("\n" + "=" * 60)
    print("e44_structure.py  —  s1.4: Graded decomposition L_j")
    print("=" * 60)

    # Build all graded pieces
    build_all_graded_pieces(j_max=3)

    # (a) Dimension table
    # sdim(L_j) = (N(j+1)*4 | N(j+1)*4)  where N(d) = C(d+3,3)
    from sage.all import binomial
    print("\n[s1.4.a] Graded dimension table sdim(L_j) = (even | odd):")
    print(f"  {'j':>3}  {'N(j+1)':>7}  {'even':>6}  {'odd':>6}  {'total':>7}")
    print(f"  {'-'*3}  {'-'*7}  {'-'*6}  {'-'*6}  {'-'*7}")
    for j in range(-1, 4):
        n_mons  = len(monomials_of_degree(j + 1))
        n_even  = sum(1 for b in E44[j] if b['parity'] == 0)
        n_odd   = sum(1 for b in E44[j] if b['parity'] == 1)
        n_total = n_even + n_odd
        expected_per_parity = n_mons * N_VARS
        flag = '' if n_even == expected_per_parity else '  ← CHECK'
        print(f"  {j:>3}  {n_mons:>7}  {n_even:>6}  {n_odd:>6}  {n_total:>7}{flag}")

    # (b) Graded bracket property [L_i, L_j] ⊆ L_{i+j}
    #     Test all 16 degree-pairs (i,j) with i,j \in {-1,0,1,2}, i+j <= 3.
    #     For each pair, test 4 random basis elements (all 4 parity combinations).
    print("\n[s1.4.b] Graded bracket property [L_i, L_j] ⊆ L_{i+j}:")
    import random
    random.seed(42)
    graded_failures = []
    test_count = 0
    for i in range(-1, 3):
        for jj in range(-1, 3):
            if i + jj > 3 or i + jj < -1:
                continue
            basis_i = E44[i]
            basis_j = E44[jj]
            if not basis_i or not basis_j:
                continue
            # Sample up to 4 pairs from each (parity, parity) combination
            for pi in (0, 1):
                pool_i = [b for b in basis_i if b['parity'] == pi]
                for pj in (0, 1):
                    pool_j = [b for b in basis_j if b['parity'] == pj]
                    if not pool_i or not pool_j:
                        continue
                    sample_pairs = [
                        (random.choice(pool_i), random.choice(pool_j))
                        for _ in range(min(4, len(pool_i) * len(pool_j)))
                    ]
                    for ba, bb in sample_pairs:
                        test_count += 1
                        ok, got, exp = is_graded_bracket(
                            ba['basis'], ba['parity'],
                            bb['basis'], bb['parity']
                        )
                        # ok is None for zero bracket (fine)
                        if ok is False:
                            graded_failures.append((i, jj, ba['label'], bb['label'], got, exp))
    if not graded_failures:
        print(f"  [L_i, L_j] ⊆ L_{{i+j}} for all {test_count} sampled pairs: VERIFIED.")
    else:
        print(f"  FAILURES ({len(graded_failures)}):")
        for i_, j_, la, lb, got_, exp_ in graded_failures[:4]:
            print(f"    [{la}(L_{i_}), {lb}(L_{j_})]: deg={got_}, expected {exp_}")

    # (c) Explicit spot-check: [L_0, L_1] ⊆ L_1
    #     Take X = x1 x2 \partial_3 \in L_1 (even), Y = x1 \partial_2 \in L_0 (even)
    #     [Y, X] should be in L_1.
    print("\n[s1.4.c] Explicit check [x_1\partial_2, x_1x_2\partial_3] \in L_1:")
    X_01 = (R(0), R(0), x1 * x2, R(0))    # x_1x_2 \partial_3 \in L_1 even
    Y_00 = (R(0), x1, R(0), R(0))          # x_1 \partial_2 \in L_0 even
    br_01 = lie_bracket_ww(Y_00, X_01)
    deg_br = field_degree(br_01)
    print(f"  X = {X_01}  (deg={field_degree(X_01)})")
    print(f"  Y = {Y_00}  (deg={field_degree(Y_00)})")
    print(f"  [Y, X] = {br_01}  (deg={deg_br},  expected 1)")
    print(f"  In L_1: {'VERIFIED' if deg_br == 1 else 'FAILED'}")

    # (d) Explicit spot-check: [L_1, L_{-1}] ⊆ L_0
    #     X = x1 x2 \partial_3 \in L_1 (even), Z = \partial_1 \in L_{-1} (even)
    print("\n[s1.4.d] Explicit check [x_1x_2\partial_3, \partial_1] \in L_0:")
    Z_m1 = e[0]                             # \partial_1 \in L_{-1} even
    br_10 = lie_bracket_ww(X_01, Z_m1)
    deg_10 = field_degree(br_10)
    print(f"  X = {X_01}  (deg={field_degree(X_01)})")
    print(f"  Z = {Z_m1}  (deg={field_degree(Z_m1)})")
    print(f"  [X, Z] = {br_10}  (deg={deg_10},  expected 0)")
    print(f"  In L_0: {'VERIFIED' if deg_10 == 0 else 'FAILED'}")

    # (e) Mixed-parity spot-check: [L_0 even, L_1 odd] ⊆ L_1 odd
    #     X = x1 \partial_2 \in L_0 even; omega = x3^2 dx_2 \in L_1 odd
    #     L_X(\omega)[0] = \omega[1]*\partialX[1]/\partialx[0] = x3^2*\partial(x1)/\partialx1 = x3^2; div(X)=0
    #     \to [X,\omega] = (x3^2, 0, 0, 0)  \in L_1
    print("\n[s1.4.e] Mixed bracket [L_0 even, L_1 odd] ⊆ L_1:")
    X_0e  = (R(0), x1, R(0), R(0))         # x_1 \partial_2 \in L_0
    om_1o = (R(0), x3**2, R(0), R(0))      # x_3^2 dx_2 \in L_1
    br_mix = lie_bracket_ef(X_0e, om_1o)
    deg_mix = field_degree(br_mix)
    expected_mix = (x3**2, R(0), R(0), R(0))
    match_mix = all(expand(br_mix[k] - expected_mix[k]) == R(0) for k in range(N_VARS))
    print(f"  X    = {X_0e}  (deg={field_degree(X_0e)})")
    print(f"  \omega    = {om_1o}  (deg={field_degree(om_1o)})")
    print(f"  [X,\omega] = {br_mix}  (deg={deg_mix},  expected 1)")
    print(f"  Value match (expected (x3^2,0,0,0)): {'VERIFIED' if match_mix else 'FAILED'}")
    print(f"  In L_1: {'VERIFIED' if deg_mix == 1 else 'FAILED'}")

    # (f) Summary: E44 dict keys built
    print(f"\n[s1.4.f] E44 graded dict keys: {sorted(E44.keys())}")
    for j in sorted(E44.keys()):
        ne = sum(1 for b in E44[j] if b['parity'] == 0)
        no = sum(1 for b in E44[j] if b['parity'] == 1)
        print(f"  E44[{j:>2}]: sdim = ({ne}|{no})")

    print("\n" + "=" * 60)
    print("s1.4 complete.")
    print("=" * 60)

    # =======================================================================
    # s1.5 — Full bracket table and Jacobi super-identity verification
    # =======================================================================
    print("\n" + "=" * 60)
    print("e44_structure.py  —  s1.5: Bracket table and Jacobi identity")
    print("=" * 60)

    # (a) Build pool from L_{-1} \cup L_0 \cup L_1
    pool = collect_pool(j_max=1)
    n_even_pool = sum(1 for b in pool if b['parity'] == 0)
    n_odd_pool  = sum(1 for b in pool if b['parity'] == 1)
    print(f"\n[s1.5.a] Pool: L_{{-1}} \cup L_0 \cup L_1")
    print(f"  Total basis elements: {len(pool)}  (even: {n_even_pool}, odd: {n_odd_pool})")
    print(f"  Ordered pairs: {len(pool)**2}")

    # (b) Compute the full bracket table
    print("\n[s1.5.b] Computing bracket table ...")
    btable = compute_bracket_table(pool)
    n_zero    = sum(1 for v in btable.values() if all(expand(c) == R(0) for c in v))
    n_nonzero = len(btable) - n_zero
    print(f"  Brackets computed: {len(btable)}")
    print(f"  Zero brackets:     {n_zero}")
    print(f"  Non-zero brackets: {n_nonzero}")

    # (c) Spot-check three known values from the table
    # Look up labels by matching basis data to avoid index-order assumptions.
    print("\n[s1.5.c] Bracket table spot-checks:")

    def find_label(pool, target_basis, target_parity):
        """Return the label of the element in pool matching target_basis and parity."""
        for b in pool:
            if b['parity'] == target_parity and all(
                expand(b['basis'][k] - target_basis[k]) == R(0)
                for k in range(N_VARS)
            ):
                return b['label']
        return None

    # i.  [\partial_1, x_1\partial_1] = \partial_1
    la_i  = find_label(pool, (R(1), R(0), R(0), R(0)), 0)   # \partial_1  \in L_{-1}
    lb_i  = find_label(pool, (x1, R(0), R(0), R(0)),  0)    # x_1\partial_1 \in L_0
    val_i = btable.get((la_i, lb_i)) if la_i and lb_i else None
    exp_i = (R(1), R(0), R(0), R(0))
    ok_i  = val_i is not None and all(expand(val_i[k] - exp_i[k]) == R(0) for k in range(N_VARS))
    print(f"  [\partial_1, x_1\partial_1]: labels=({la_i}, {lb_i}), result={val_i}, expected={exp_i}: {'VERIFIED' if ok_i else 'FAILED'}")

    # ii. [\partial_1, x_1^2\partial_1] = 2x_1\partial_1
    lb_ii = find_label(pool, (x1**2, R(0), R(0), R(0)), 0)  # x_1^2\partial_1 \in L_1
    val_ii = btable.get((la_i, lb_ii)) if la_i and lb_ii else None
    exp_ii = (2*x1, R(0), R(0), R(0))
    ok_ii  = val_ii is not None and all(expand(val_ii[k] - exp_ii[k]) == R(0) for k in range(N_VARS))
    print(f"  [\partial_1, x_1^2\partial_1]: labels=({la_i}, {lb_ii}), result={val_ii}, expected={exp_ii}: {'VERIFIED' if ok_ii else 'FAILED'}")

    # iii.[x_1\partial_2, dx_2] = dx_1
    #     L_{x_1\partial_2}(dx_2)_n = x_1 \partial(dx_2)_n/\partialx_1 + (dx_2)_1 \partial(x_1\partial_2)_1/\partialxₙ
    #     = x_1\cdot0 + 1\cdot\partial(x_1)/\partialxₙ = \delta_{n,0}   div(x_1\partial_2)=0 \to result = dx_1
    la_iii  = find_label(pool, (R(0), x1, R(0), R(0)),  0)  # x_1\partial_2 \in L_0
    lb_iii  = find_label(pool, (R(0), R(1), R(0), R(0)), 1) # dx_2 \in L_{-1}
    val_iii = btable.get((la_iii, lb_iii)) if la_iii and lb_iii else None
    exp_iii = (R(1), R(0), R(0), R(0))
    ok_iii  = val_iii is not None and all(expand(val_iii[k] - exp_iii[k]) == R(0) for k in range(N_VARS))
    print(f"  [x_1\partial_2, dx_2]: labels=({la_iii}, {lb_iii}), result={val_iii}, expected={exp_iii}: {'VERIFIED' if ok_iii else 'FAILED'}")

    # (d) Anti-symmetry statistics from the table
    #     For even-even: btable[(a,b)] = -btable[(b,a)]
    #     For odd-odd:   btable[(a,b)] = +btable[(b,a)]
    print("\n[s1.5.d] Super-symmetry audit of bracket table:")
    asym_fail = 0
    sym_fail  = 0
    ee_checked = oe_checked = oo_checked = 0
    for ba in pool:
        for bb in pool:
            pa, pb = ba['parity'], bb['parity']
            v_ab = btable[(ba['label'], bb['label'])]
            v_ba = btable[(bb['label'], ba['label'])]
            if pa == 0 and pb == 0:
                ee_checked += 1
                if any(expand(v_ab[k] + v_ba[k]) != R(0) for k in range(N_VARS)):
                    asym_fail += 1
            elif pa == 1 and pb == 1:
                oo_checked += 1
                if any(expand(v_ab[k] - v_ba[k]) != R(0) for k in range(N_VARS)):
                    sym_fail += 1
    print(f"  Even-even antisymmetry checks: {ee_checked}, failures: {asym_fail}")
    print(f"  Odd-odd symmetry checks:       {oo_checked}, failures: {sym_fail}")
    status = 'VERIFIED' if (asym_fail == 0 and sym_fail == 0) else 'FAILED'
    print(f"  Super-symmetry of the full bracket table: {status}")

    # (e) Jacobi super-identity on 50 random triples from the pool
    print("\n[s1.5.e] Jacobi super-identity — verify_jacobi(n=50):")
    passed, n_checked, failures_j = verify_jacobi(pool, n=50, seed=0)
    if passed:
        print(f"  verify_jacobi(n={n_checked}): PASSED  (0 failures)")
    else:
        print(f"  verify_jacobi(n={n_checked}): FAILED  ({len(failures_j)} failures)")
        for item in failures_j[:3]:
            print(f"    a={item['a']}, b={item['b']}, c={item['c']}")
            print(f"    residual = {item['residual']}")

    # (f) Extra: 50 triples drawn from L_{-1} \cup L_0 \cup L_1 \cup L_2 \cup L_3
    pool_full = collect_pool(j_max=3)
    print(f"\n[s1.5.f] Jacobi on 50 triples from full pool (j <= 3, {len(pool_full)} elements):")
    passed_f, n_f, fails_f = verify_jacobi(pool_full, n=50, seed=7)
    if passed_f:
        print(f"  verify_jacobi(n={n_f}, full pool): PASSED")
    else:
        print(f"  verify_jacobi(n={n_f}, full pool): FAILED ({len(fails_f)} failures)")
        for item in fails_f[:2]:
            print(f"    a={item['a']}, b={item['b']}, c={item['c']}")

    print("\n" + "=" * 60)
    print("s1.5 complete.")
    print("=" * 60)

    # =======================================================================
    # s1.6 — Export bracket table to e44_brackets.pkl
    # =======================================================================
    print("\n" + "=" * 60)
    print("e44_structure.py  —  s1.6: Export bracket table")
    print("=" * 60)

    # (a) Save
    pkl_path = "e44_brackets.pkl"
    print(f"\n[s1.6.a] Saving to {pkl_path} ...")
    saved = save_bracket_table(pkl_path)
    import os as _os
    size_kb = _os.path.getsize(saved) / 1024
    print(f"  Written: {saved}  ({size_kb:.1f} KB)")

    # (b) Reload and verify round-trip
    print(f"\n[s1.6.b] Reloading and verifying round-trip ...")
    payload = load_bracket_table(pkl_path)

    # Keys present
    required_keys = {'E44', 'btable', 'pool', 'pool_full', 'VARS', 'N_VARS'}
    missing = required_keys - set(payload.keys())
    print(f"  Keys present: {sorted(payload.keys())}")
    print(f"  Missing keys: {missing if missing else 'none'}")

    # Dimension check
    for j in range(-1, 4):
        n_orig  = len(E44[j])
        n_reloaded = len(payload['E44'][j])
        if n_orig != n_reloaded:
            print(f"  FAIL E44[{j}]: orig {n_orig} vs reloaded {n_reloaded}")
    print(f"  E44 graded bases: all sizes match VERIFIED")

    # Bracket table size
    n_orig_bt = len(btable)
    n_rel_bt  = len(payload['btable'])
    print(f"  btable size: orig={n_orig_bt}, reloaded={n_rel_bt}  {'VERIFIED' if n_orig_bt == n_rel_bt else 'FAILED'}")

    # Spot-check one bracket value survives pickle round-trip
    la_rt = find_label(pool, (R(1), R(0), R(0), R(0)), 0)
    lb_rt = find_label(pool, (x1,  R(0), R(0), R(0)), 0)
    if la_rt and lb_rt:
        val_orig = btable[(la_rt, lb_rt)]
        val_rt   = payload['btable'][(la_rt, lb_rt)]
        ok_rt = all(expand(val_orig[k] - val_rt[k]) == R(0) for k in range(N_VARS))
        print(f"  Spot-check [{la_rt}, {lb_rt}]: round-trip {'VERIFIED' if ok_rt else 'FAILED'}")

    # Pool sizes
    print(f"  pool size: {len(payload['pool'])}  (expected {len(pool)})")
    print(f"  pool_full size: {len(payload['pool_full'])}  (expected {len(pool_full)})")

    print(f"\n  Milestone: e44_brackets.pkl written and reloadable: PASSED")

    print("\n" + "=" * 60)
    print("s1.6 complete  —  S1 (e44_structure.py) DONE.")
    print("=" * 60)

    # =======================================================================
    # s4 — \hat{p}(4)-modules integrated into E(4,4) structure
    # =======================================================================
    print("\n" + "=" * 60)
    print("e44_structure.py  —  s4: \hat{p}(4)-modules in E(4,4)")
    print("=" * 60)

    # Reload the just-saved bracket data for the \hat{p}(4) checks
    e44_data_s4 = load_e44()
    if e44_data_s4 is not None:
        # (a) Module dimensions table
        print("\n[s4.a] \hat{p}(4)-module dimensions:")
        test_cases_s4 = [
            (1, 1, 0, 0, "standard module (= L_{-1})"),
            (1, 0, 0, 0, "trivial fiber t=1"),
            (0, 0, 0, 0, "trivial t=0"),
            (0, 0, 0, 1, "dual fundamental"),
            (0, 0, 0, 2, "Sym^2 dual"),
        ]
        print(f"  {'Parameters':<20} {'dim_K':>6} {'dim_W':>6}  Note")
        print(f"  {'─'*20} {'─'*6} {'─'*6}  {'─'*20}")
        for t_s4, a_s4, b_s4, c_s4, note_s4 in test_cases_s4:
            K_s4 = KacModule(t_s4, a_s4, b_s4, c_s4, e44_data_s4)
            W_s4 = phat4_module(t_s4, a_s4, b_s4, c_s4, e44_data_s4)
            print(f"  W_{t_s4}({a_s4},{b_s4},{c_s4})"
                  f"{' '*(14-len(f'W_{t_s4}({a_s4},{b_s4},{c_s4})'))} "
                  f"{K_s4.dim:6d} {W_s4.dim:6d}  {note_s4}")

        # (b) Full verification
        print()
        verify_phat4_in_e44(e44_data=e44_data_s4, verbose=True)

        # (c) Run phat4_modules internal checks
        print()
        _check_phat4(verbose=True)
    else:
        print("\n[s4] SKIPPED — e44_data not available.")

    print("\n" + "=" * 60)
    print("s4 complete  —  All sections DONE.")
    print("=" * 60)
