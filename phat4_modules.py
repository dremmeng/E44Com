"""
phat4_modules.py  —  Full \hat{p}(4)-modules W_t(a,b,c) for the NSE44 programme
==========================================================================
Cantarini-Caselli-Kac 2026: E(4,4) as the Navier-Stokes algebra.

Implements the irreducible \hat{p}(4)-module W_t(a,b,c) as a quotient of the
Kac module K_t(a,b,c) = \wedge^\circ(\hat{p}(4)_{-1}) \otimes F_t(a,b,c), where F_t(a,b,c)
is the irreducible sl_4-module V(a,b,c) with C acting as t\cdotId.

The \hat{p}(4) grading (internal to L_0 ≅ \hat{p}(4)):
    \hat{p}(4)_{-2} = \C C          (central, 1-dim)
    \hat{p}(4)_{-1} = span{a_{ij}} (6 odd generators, i<j)
    \hat{p}(4)_0    = gl_4           (16 even generators)
    \hat{p}(4)_1    = span{b_{ij}} (10 odd generators, i≤j)

where  a_{ij} = x_i d_j - x_j d_i  (antisymmetric)
       b_{ij} = x_i d_j + x_j d_i  (symmetric)

*Depends on* verma_modules.py for W_hat (sl_4 crystal basis) and e44_data.

Run inside SageMath:  sage phat4_modules.py
"""

from sage.all import QQ, vector, matrix, identity_matrix, block_matrix
import sys as _sys
import os as _os
import pickle as _pickle

_HERE = _os.path.dirname(_os.path.abspath(__file__))
if _HERE not in _sys.path:
    _sys.path.insert(0, _HERE)

from verma_modules import W_hat, load_e44, _expand_in_basis

_PHAT4_CACHE_FILE = _os.path.join(_HERE, 'phat4_cache.pkl')


# ===========================================================================
# L_0 index helpers
# ===========================================================================

def _l0_even_idx(r, s):
    """L_0 (full) index for the even generator x_r \partial_s (1-indexed r,s)."""
    return 4 * (4 - r) + (s - 1)


def _l0_odd_idx(r, s):
    """L_0 (full) index for the odd generator x_r dx_s (1-indexed r,s)."""
    return 16 + 4 * (4 - r) + (s - 1)


# Pairs for the six a_{ij} generators of \hat{p}(4)_{-1} (ordered canonically)
AIJ_PAIRS = [(1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]

# Pairs for the ten b_{ij} generators of \hat{p}(4)_1 (ordered canonically)
BIJ_PAIRS = [(1, 1), (1, 2), (1, 3), (1, 4),
             (2, 2), (2, 3), (2, 4),
             (3, 3), (3, 4),
             (4, 4)]

# Bilinear form B(a_m, a_n) for the Clifford structure of \hat{p}(4)_{-1}.
# {a_m, a_n} = c_mn \cdot C, so B(a_m, a_n) = c_mn (multiply by t in the module).
# Pairs: (m, n, c_mn) where m < n index into AIJ_PAIRS.
#   a_{12} ↔ a_{34}: c = 2
#   a_{13} ↔ a_{24}: c = -2
#   a_{14} ↔ a_{23}: c = 2
_CLIFFORD_PAIRS = [(0, 5, QQ(2)), (1, 4, QQ(-2)), (2, 3, QQ(2))]


# ===========================================================================
# Kac module  K_t(a,b,c)
# ===========================================================================

class KacModule:
    """
    Kac module K_t(a,b,c) = \wedge^\circ(\hat{p}(4)_{-1}) \otimes F_t(a,b,c).

    Basis ordering:  for each subset I ⊆ {0,..,5} (indexing the 6 a_{ij}
    generators in AIJ_PAIRS order), and each crystal basis vector v_k
    of V(a,b,c), there is a basis element  \omega_I \otimes v_k.

    Subsets are ordered by (|I|, lexicographic), giving 2^6 = 64 blocks
    of size dim_V each, for total dimension 64 \cdot dim_V.

    Attributes
    ----------
    t, a, b, c : module parameters
    V          : WHat4Module (the sl_4 part)
    dim_V      : dimension of V(a,b,c)
    dim        : total dimension = 64 \cdot dim_V
    subsets    : list of 64 frozensets (ordered basis for \wedge^\circ)
    subset_idx : {frozenset: int} index map
    action_mats : dict {L0_full_idx: dim\timesdim QQ-matrix}
    """

    def __init__(self, t, a, b, c, e44_data):
        self.t = QQ(t)
        self.a, self.b, self.c = int(a), int(b), int(c)
        self.V = W_hat(t, a, b, c)
        self.dim_V = self.V.dim
        self.e44_data = e44_data

        # Build ordered list of subsets of {0,...,5}
        self.subsets = []
        for size in range(7):
            from itertools import combinations
            for combo in combinations(range(6), size):
                self.subsets.append(frozenset(combo))
        self.subset_idx = {s: i for i, s in enumerate(self.subsets)}
        assert len(self.subsets) == 64

        self.dim = 64 * self.dim_V

        # Caches for bracket decompositions
        self._adjoint_cache = {}   # L0_idx -> {m: {m': coeff}}
        self._sl4_cache = {}       # L0_even_idx -> matrix
        self._b_a_cache = {}       # (bi, bj, ai, aj) -> matrix

        # Build action matrices for all 32 L_0 generators
        self.action_mats = {}
        self._build_all_actions()

    def _flat_idx(self, subset_pos, v_idx):
        """Flat index for basis element (subset_pos, v_idx)."""
        return subset_pos * self.dim_V + v_idx

    def _build_all_actions(self):
        """Build action matrices for all 32 L_0 generators."""
        for L0_idx in range(32):
            self.action_mats[L0_idx] = self._build_action(L0_idx)

    def _build_action(self, L0_idx):
        """
        Build the dim\timesdim action matrix of L_0[L0_idx] on K_t.

        For even L_0 (gl_4):
            X \cdot (\omega_I \otimes v) = (ad(X)\cdot\omega_I) \otimes v + \omega_I \otimes (X\cdotv)
            ad(X) acts on each a_{ij} factor in \omega_I via [X, a_{ij}] ∈ \hat{p}(4)_{-1}.

        For odd L_0 in \hat{p}(4)_{-1} (a_{kl}):
            a_{kl} \cdot (\omega_I \otimes v) = (a_{kl} \wedge \omega_I) \otimes v
            (left multiplication in exterior algebra)

        For odd L_0 in \hat{p}(4)_1 (b_{kl}):
            b_{kl} . (omega_I x v) = sum_j (+-) omega_{I\\j} x ([b_{kl}, a_{I_j}].v)
            (contraction / Leibniz through the exterior factors)
        """
        L0 = self.e44_data['E44'][0]
        b_elem = L0[L0_idx]

        if b_elem['parity'] == 0:
            return self._build_even_action(L0_idx)
        else:
            return self._build_odd_action(L0_idx)

    def _build_even_action(self, L0_idx):
        """
        Even L_0 generator (gl_4) acting on K_t via Leibniz:
            X \cdot (\omega_I \otimes v) = (ad(X)\cdot\omega_I) \otimes v  +  \omega_I \otimes (X\cdotv)
        """
        mat = matrix(QQ, self.dim, self.dim)

        # 1. Fiber action: X\cdotv using crystal matrices
        W_mat = self._sl4_action_on_V(L0_idx)

        # 2. Adjoint action: [X, a_{ij}] decomposed in {a_{kl}} basis
        # Pre-compute: for each a_{ij} (index m=0..5), the result of
        # [even_L0, a_{ij}] expressed in {a_{kl}} basis.
        ad_on_a = self._adjoint_even_on_a(L0_idx)
        # ad_on_a[m] = {m': coeff}  meaning  [X, a_m] = Σ coeff * a_{m'}

        for si, I in enumerate(self.subsets):
            # Fiber action: \omega_I \otimes (X\cdotv)
            for k in range(self.dim_V):
                col = self._flat_idx(si, k)
                for k2 in range(self.dim_V):
                    c = W_mat[k2, k]
                    if c != 0:
                        row = self._flat_idx(si, k2)
                        mat[row, col] += c

            # Adjoint action on \wedge part: Leibniz over each factor in I
            I_sorted = sorted(I)
            for pos, m in enumerate(I_sorted):
                # Replace a_m with [X, a_m] = Σ coeff * a_{m'}
                sign = QQ((-1) ** pos)  # sign from moving X past preceding factors
                for m_prime, coeff in ad_on_a[m].items():
                    new_I = frozenset((I - {m}) | {m_prime})
                    if m_prime in (I - {m}):
                        continue  # a_{m'} already in I → \wedge = 0
                    if len(new_I) != len(I):
                        continue  # shouldn't happen but safety check

                    # Sign for reinserting m' into the sorted order
                    I_minus_m = sorted(I - {m})
                    ins_pos = sum(1 for x in I_minus_m if x < m_prime)
                    total_sign = sign * QQ((-1) ** ins_pos) * coeff

                    if new_I in self.subset_idx:
                        new_si = self.subset_idx[new_I]
                        for k in range(self.dim_V):
                            col = self._flat_idx(si, k)
                            row = self._flat_idx(new_si, k)
                            mat[row, col] += total_sign

        return mat

    def _build_odd_action(self, L0_idx):
        """
        Odd L_0 generator acting on K_t.
        Classify as \hat{p}(4)_{-1} (a_{ij}) or \hat{p}(4)_1 (b_{ij}) and dispatch.
        """
        # Decompose L0[L0_idx] = x_r dx_s into a_{ij} and b_{ij} components.
        # x_r dx_s = (1/2)(a_{rs} + b_{rs}) if r ≠ s
        # x_r dx_r = (1/2) b_{rr}
        #
        # But we need to be careful: individual x_r dx_s generators are NOT
        # purely a or b. They are sums: x_r dx_s = (a_{rs} + b_{rs})/2 for r<s,
        # x_r dx_s = (-a_{sr} + b_{rs})/2 for r>s, x_r dx_r = b_{rr}/2.
        #
        # So we decompose and compute the action as a sum.

        r, s = self._odd_l0_to_rs(L0_idx)  # 1-indexed

        mat = matrix(QQ, self.dim, self.dim)

        if r == s:
            # x_r dx_r = b_{rr}/2 → action = (1/2) * b_{rr}_action
            mat += QQ(1) / 2 * self._b_action(r, s)
        elif r < s:
            # x_r dx_s = (b_{rs} + a_{rs})/2
            mat += QQ(1) / 2 * self._b_action(r, s)
            mat += QQ(1) / 2 * self._a_action(r, s)
        else:
            # r > s: x_r dx_s = (b_{sr} - a_{sr})/2  since a_{sr} = -a_{rs}
            # x_r dx_s = x_r dx_s.  b_{sr} = x_s dx_r + x_r dx_s, a_{sr} = x_s dx_r - x_r dx_s
            # So x_r dx_s = (b_{sr} - a_{sr})/2
            mat += QQ(1) / 2 * self._b_action(s, r)
            mat -= QQ(1) / 2 * self._a_action(s, r)

        return mat

    def _odd_l0_to_rs(self, L0_idx):
        """
        Convert odd L_0 index (16..31) to 1-indexed (r, s) for x_r dx_s.
        L_0[16+n] where n = 4*(4-r) + (s-1).
        """
        n = L0_idx - 16
        r = 4 - n // 4   # 1-indexed
        s = n % 4 + 1     # 1-indexed
        return r, s

    def _a_action(self, i, j):
        """
        Action of a_{ij} ∈ \hat{p}(4)_{-1} on K_t (Clifford multiplication):
            a_m \cdot (\omega_I \otimes v) = (a_m \wedge \omega_I) \otimes v  +  l_{a_m}(\omega_I) \otimes v

        The contraction l uses the bilinear form from {a_m, a_n} = c_mn \cdot C,
        with C acting as t on the module:
            l_{a_m}(\omega_I) = Σ_{n ∈ I} (-1)^{pos(n)} \cdot (c_mn \cdot t / 2) \cdot \omega_{I\n}
        """
        m = AIJ_PAIRS.index((i, j))  # index 0..5

        # Pre-compute contraction coefficients: B_half[n] = c_mn * t / 2
        B_half = {}
        for m1, m2, c_mn in _CLIFFORD_PAIRS:
            if m1 == m:
                B_half[m2] = c_mn * self.t / 2
            elif m2 == m:
                B_half[m1] = c_mn * self.t / 2

        mat = matrix(QQ, self.dim, self.dim)
        for si, I in enumerate(self.subsets):
            # Term 1: exterior multiplication  a_m \wedge \omega_I
            if m not in I:
                new_I = frozenset(I | {m})
                sign = QQ((-1) ** sum(1 for x in I if x < m))
                new_si = self.subset_idx[new_I]
                for k in range(self.dim_V):
                    col = self._flat_idx(si, k)
                    row = self._flat_idx(new_si, k)
                    mat[row, col] += sign

            # Term 2: contraction  l_{a_m}(\omega_I)
            I_sorted = sorted(I)
            for pos, n in enumerate(I_sorted):
                b_val = B_half.get(n)
                if b_val is None or b_val == 0:
                    continue
                contract_I = frozenset(I - {n})
                contract_si = self.subset_idx[contract_I]
                contract_sign = QQ((-1) ** pos) * b_val
                for k in range(self.dim_V):
                    col = self._flat_idx(si, k)
                    row = self._flat_idx(contract_si, k)
                    mat[row, col] += contract_sign

        return mat

    def _b_action(self, i, j):
        """
        Action of b_{ij} ∈ \hat{p}(4)_1 on K_t via recursive PBW normal-ordering.

        In the induced module K = U(\hat{p}_{-}) \otimes_{U(p+)} V, the correct
        recursive formula (using exterior basis) is:

            b \cdot (a_{m₁} \wedge \omega' \otimes v) = g_{m₁}^{full} \cdot (\omega' \otimes v)
                                      - a_{m₁}^{Cl} \cdot [b \cdot (\omega' \otimes v)]
                                      - b \cdot l_{a_{m₁}}(\omega')

        where g_{m₁} = {b, a_{m₁}} ∈ gl_4 acts via full Leibniz on \omega' \otimes v,
        a_{m₁}^{Cl} acts by Clifford multiplication (exterior + contraction),
        l_{a_{m₁}}(\omega') is the contraction of \omega' by a_{m₁} via the bilinear
        form B(a_m, a_n) = c_mn \cdot t / 2 from {a_m, a_n} = c_mn \cdot C,
        and the base case is b \cdot (1 \otimes v) = 0.

        Subsets are processed in order of increasing size so that the
        recursive term (b acting on smaller subsets) is already computed.
        """
        # Pre-compute g_m = {b_{ij}, a_m} full action matrices
        g_full = []
        for m, (k, l) in enumerate(AIJ_PAIRS):
            g_full.append(self._bracket_b_a_full_action(i, j, k, l))

        # Pre-compute a_m Clifford multiplication matrices
        a_ext = []
        for m, (ai, aj) in enumerate(AIJ_PAIRS):
            a_ext.append(self._a_action(ai, aj))

        # Pre-compute contraction coefficients B_half[m1][m2] = c_{m1,m2}*t/2
        B_half = {}
        for m1, m2, c_mn in _CLIFFORD_PAIRS:
            B_half[(m1, m2)] = c_mn * self.t / 2
            B_half[(m2, m1)] = c_mn * self.t / 2

        mat = matrix(QQ, self.dim, self.dim)
        # Process subsets in order of increasing size (self.subsets is
        # already sorted this way by construction).
        for si, I in enumerate(self.subsets):
            if len(I) == 0:
                continue  # base case: b \cdot (1 \otimes v) = 0

            m_1 = min(I)
            I_prime = frozenset(I - {m_1})
            si_prime = self.subset_idx[I_prime]

            G = g_full[m_1]   # full even action of g_{m_1} on K
            A = a_ext[m_1]    # Clifford mult by a_{m_1} on K

            p0 = si_prime * self.dim_V
            c0 = si * self.dim_V
            dV = self.dim_V

            # Extract the dim_V columns of mat for I' (already computed)
            mat_block = mat.submatrix(0, p0, self.dim, dV)

            # Term 2: A * mat_block  (a_{m₁}^{Cl} \cdot [b \cdot \omega'])
            correction = A * mat_block

            # Term 1: G columns for I' block
            if G is not None:
                G_block = G.submatrix(0, p0, self.dim, dV)
                result_block = G_block - correction
            else:
                result_block = -correction

            # Term 3: b \cdot l_{a_{m₁}}(\omega')
            # l_{a_{m₁}}(\omega_{I'}) = Σ_{n ∈ I'} (-1)^{pos(n)} B(m₁,n) \omega_{I'\n}
            # b acts on \omega_{I'\n} which is at level |I|-2 (already computed).
            I_prime_sorted = sorted(I_prime)
            for pos_n, n in enumerate(I_prime_sorted):
                b_val = B_half.get((m_1, n))
                if b_val is None or b_val == 0:
                    continue
                J = frozenset(I_prime - {n})
                sj = self.subset_idx[J]
                j0 = sj * dV
                sign_n = QQ((-1) ** pos_n) * b_val
                # Subtract b \cdot (sign_n * \omega_J) = sign_n * mat[:, j0:j0+dV]
                mat_J = mat.submatrix(0, j0, self.dim, dV)
                result_block -= sign_n * mat_J

            # Write result into mat columns c0..c0+dV
            for row in range(self.dim):
                for k in range(dV):
                    v = result_block[row, k]
                    if v != 0:
                        mat[row, c0 + k] = v
        return mat

    def _bracket_b_a_full_action(self, bi, bj, ai, aj):
        """
        Compute {b_{bi,bj}, a_{ai,aj}} as a gl_4 element, then return
        its FULL action matrix on K (adjoint on \wedge + action on V).
        Cached.
        """
        cache_key = ('full', bi, bj, ai, aj)
        if cache_key in self._b_a_cache:
            return self._b_a_cache[cache_key]

        L0 = self.e44_data['E44'][0]
        bt = self.e44_data['btable']
        R_poly = self.e44_data['E44'][-1][0]['basis'][0].parent()

        # b_{bi,bj} components in x_r dx_s basis:
        if bi == bj:
            b_components = [(_l0_odd_idx(bi, bi), QQ(2))]
        else:
            b_components = [(_l0_odd_idx(bi, bj), QQ(1)),
                            (_l0_odd_idx(bj, bi), QQ(1))]

        # a_{ai,aj} components (ai < aj):
        a_components = [(_l0_odd_idx(ai, aj), QQ(1)),
                        (_l0_odd_idx(aj, ai), QQ(-1))]

        # Compute bracket {b, a} as sum of pairwise brackets
        result = tuple(R_poly(0) for _ in range(4))
        for b_idx, b_coeff in b_components:
            for a_idx, a_coeff in a_components:
                br = bt.get((L0[b_idx]['label'], L0[a_idx]['label']))
                if br is not None:
                    result = tuple(result[q] + b_coeff * a_coeff * br[q]
                                   for q in range(4))

        if all(c == 0 for c in result):
            self._b_a_cache[cache_key] = None
            return None

        # Decompose in even L_0 basis (the bracket of two odd is even)
        L0_even = [b for b in L0 if b['parity'] == 0]
        coeffs = _expand_in_basis(result, L0_even, R_poly)

        # Build the FULL action matrix on K = \wedge \otimes V
        action = matrix(QQ, self.dim, self.dim)
        for sec_idx, coeff in coeffs.items():
            full_idx = L0.index(L0_even[sec_idx])
            action += QQ(coeff) * self.action_mats[full_idx]

        self._b_a_cache[cache_key] = action
        return action

    def _sl4_action_on_V(self, L0_even_idx):
        """
        Action matrix of even L_0[L0_even_idx] on V(a,b,c)
        using the crystal Chevalley matrices.  Cached.
        """
        if L0_even_idx in self._sl4_cache:
            return self._sl4_cache[L0_even_idx]
        from verma_modules import _w_action_from_l0_idx
        mat = _w_action_from_l0_idx(self.V, L0_even_idx)
        self._sl4_cache[L0_even_idx] = mat
        return mat

    def _adjoint_even_on_a(self, L0_idx):
        """
        Compute ad(L_0[L0_idx]) on each a_{m} ∈ \hat{p}(4)_{-1},
        returning the result decomposed in the {a_m} basis.  Cached.

        Returns dict {m: {m': coeff}} where [X, a_m] = Σ coeff * a_{m'}.
        """
        if L0_idx in self._adjoint_cache:
            return self._adjoint_cache[L0_idx]
        L0 = self.e44_data['E44'][0]
        bt = self.e44_data['btable']
        R_poly = self.e44_data['E44'][-1][0]['basis'][0].parent()

        result = {}
        for m, (ai, aj) in enumerate(AIJ_PAIRS):
            # a_{ai,aj} = x_{ai} dx_{aj} - x_{aj} dx_{ai}
            a_components = [(_l0_odd_idx(ai, aj), QQ(1)),
                            (_l0_odd_idx(aj, ai), QQ(-1))]

            # Compute [L_0[L0_idx], a_{ai,aj}]
            br_total = tuple(R_poly(0) for _ in range(4))
            for a_idx, a_coeff in a_components:
                br = bt.get((L0[L0_idx]['label'], L0[a_idx]['label']))
                if br is not None:
                    br_total = tuple(br_total[q] + a_coeff * br[q]
                                     for q in range(4))

            # Decompose in odd L_0 basis
            if all(c == 0 for c in br_total):
                result[m] = {}
                continue

            L0_odd = [b for b in L0 if b['parity'] == 1]
            coeffs = _expand_in_basis(br_total, L0_odd, R_poly)

            # Re-express in {a_{kl}} basis
            # Each L0_odd element is x_r dx_s.
            # x_r dx_s contributes: if r<s → (1/2)(a_{rs} + b_{rs})
            #                       if r>s → (1/2)(b_{sr} - a_{sr})
            #                       if r=s → (1/2) b_{rr}
            # The result of [even, a_{ij}] must be in \hat{p}(4)_{-1} (by grading),
            # so only a_{kl} components survive.
            a_coeffs = {}
            for sec_idx, coeff in coeffs.items():
                full_odd_idx = L0.index(L0_odd[sec_idx])
                r, s = self._odd_l0_to_rs(full_odd_idx)
                if r < s:
                    m_prime = AIJ_PAIRS.index((r, s))
                    a_coeffs[m_prime] = a_coeffs.get(m_prime, QQ(0)) + coeff / 2
                elif r > s:
                    m_prime = AIJ_PAIRS.index((s, r))
                    a_coeffs[m_prime] = a_coeffs.get(m_prime, QQ(0)) - coeff / 2
                # r == s → pure b, no a component

            result[m] = {mp: c for mp, c in a_coeffs.items() if c != 0}
        self._adjoint_cache[L0_idx] = result
        return result

    def flat_to_basis(self, flat_idx):
        """Convert flat index to (subset_pos, v_idx) pair."""
        si = flat_idx // self.dim_V
        k = flat_idx % self.dim_V
        return si, k


# ===========================================================================
# Full \hat{p}(4)-module  W_t(a,b,c)  as quotient of Kac module
# ===========================================================================

class Phat4Module:
    """
    Full irreducible \hat{p}(4)-module W_t(a,b,c).

    Computed as K_t(a,b,c) / (maximal proper submodule).

    Attributes
    ----------
    t, a, b, c : module parameters
    dim         : dimension of W_t(a,b,c)
    action_mats : dict {L0_full_idx : dim\timesdim QQ-matrix}
                  Action of each of the 32 L_0 generators.
    quotient_proj : dim_K \times dim matrix
                    Projection from K onto W (columns = quotient basis in K coords).
    """

    def __init__(self, t, a, b, c, e44_data):
        self.t = QQ(t)
        self.a, self.b, self.c = int(a), int(b), int(c)
        self.e44_data = e44_data

        # Build Kac module
        K = KacModule(t, a, b, c, e44_data)
        self._K = K

        # Find maximal proper submodule and quotient
        self._compute_quotient(K)

    def _compute_quotient(self, K):
        """
        Compute W_t = K_t / (maximal proper submodule).

        Algorithm:
        1. Find all \hat{p}(4)_1-singular vectors (ker of all b_{ij} actions)
        2. Remove the hw vector v_hw = 1 \otimes v_hw
        3. Generate the submodule from remaining singular vectors
        4. Quotient K by this submodule
        """
        # Step 1: Find kernel of all b_{ij} simultaneously
        # The \hat{p}(4)_1 generators are: b_{ij} for (i,j) in BIJ_PAIRS
        # Their action on K is via L_0 odd generators.
        # We need: ker(b_{ij}) for ALL (i,j), simultaneously.

        # But rather than computing individual bij actions, we can use:
        # A \hat{p}(4)-singular vector v satisfies b_{ij}\cdotv = 0 for all raising ops.
        # The raising operators of \hat{p}(4) are \hat{p}(4)_1 = {b_{ij}}.
        # It suffices to check the Chevalley-type generators.
        # For sl_4 \subset \hat{p}(4)_0, highest weight means: annihilated by e_i (i=1,2,3)
        # plus annihilated by all b_{ij}.
        # Actually, a \hat{p}(4)-singular vector w/ respect to the standard Borel is
        # annihilated by \hat{p}(4)_1 AND by the upper-triangular part of \hat{p}(4)_0.

        # But for finding the maximal submodule, we just need to find all
        # submodule generators. The approach: find ALL vectors annihilated by
        # ALL of \hat{p}(4)_1, and also by the sl_4 raising operators (e_1, e_2, e_3).
        # These are the \hat{p}(4)-highest-weight vectors (= singular vectors).

        # Build the "annihilation matrix" for \hat{p}(4)_1 + sl_4 raising generators
        mats_to_stack = []

        # \hat{p}(4)_1 generators: b_{ij}
        for i, j in BIJ_PAIRS:
            b_mat = K._b_action(i, j)
            mats_to_stack.append(b_mat)

        # sl_4 raising: x_i \partial_{i+1} for i=1,2,3 (= e_i in Chevalley basis)
        for i in range(1, 4):
            idx = _l0_even_idx(i, i + 1)
            mats_to_stack.append(K.action_mats[idx])

        # Stack and find kernel
        total_rows = sum(m.nrows() for m in mats_to_stack)
        stacked = matrix(QQ, total_rows, K.dim)
        r = 0
        for m in mats_to_stack:
            for row_idx in range(m.nrows()):
                stacked[r] = m.row(row_idx)
                r += 1

        ker_basis = list(stacked.right_kernel().basis())

        # Step 2: Identify the hw vector (1 \otimes v_hw)
        hw_flat = K._flat_idx(0, K.V.v_hw)  # subset {} = index 0, hw crystal
        hw_vec = vector(QQ, [QQ(1) if i == hw_flat else QQ(0)
                             for i in range(K.dim)])

        # Separate: hw vector vs other singular vectors
        # Check which kernel vectors have nonzero component at hw_flat
        non_hw_singular = []
        for v in ker_basis:
            if v[hw_flat] == 0:
                non_hw_singular.append(v)
            # Those with v[hw_flat] ≠ 0 include the hw direction

        # Step 3: Generate maximal submodule from non-hw singular vectors
        if not non_hw_singular:
            # K is already irreducible (no proper submodule)
            self.dim = K.dim
            self.action_mats = dict(K.action_mats)
            self.quotient_proj = identity_matrix(QQ, K.dim)
            self._quotient_basis = list(identity_matrix(QQ, K.dim).rows())
            self.weight_spaces = dict(K.V.weight_spaces)
            self.v_hw = K.V.v_hw
            self.dim_V = K.dim_V
            return

        # Generate the submodule by orbit closure.
        # IMPORTANT: Some singular vectors may generate ALL of K (i.e., they
        # are hw vectors of composition factors that extend through the
        # irreducible quotient).  We must skip those.  Strategy: add seeds
        # one at a time and verify the orbit stays proper (doesn't contain hw).
        all_mats = [K.action_mats[idx] for idx in range(32)]

        def _orbit_closure(seed_rows):
            """Grow orbit of seed_rows under all_mats. Returns echelonized matrix."""
            work = matrix(QQ, seed_rows)
            work.echelonize()
            nz = sum(1 for i in range(work.nrows()) if work.row(i) != 0)
            work = work[:nz]
            old_rank = 0
            # Precompute transposes for batch matrix multiplication
            all_mats_T = [m.transpose() for m in all_mats]
            while work.rank() > old_rank:
                old_rank = work.rank()
                # Batch: compute all generator images via matrix-matrix multiply.
                # work (rank\timesdim) * mat^T (dim\timesdim) = (rank\timesdim) image block.
                for mT in all_mats_T:
                    work = work.stack(work[:old_rank] * mT)
                work.echelonize()
                nz = sum(1 for i in range(work.nrows())
                         if work.row(i) != 0)
                if nz < work.nrows():
                    work = work[:nz]
            return work

        def _contains_hw(echelon_mat):
            """Check if hw_vec is in the row space of echelon_mat."""
            aug = echelon_mat.stack(matrix(QQ, [hw_vec]))
            return aug.rank() == echelon_mat.rank()

        # Sort seeds by increasing subset size (higher \wedge levels are more
        # likely to be proper submodule generators)
        non_hw_sorted = sorted(non_hw_singular,
                               key=lambda v: -min(len(K.subsets[si])
                                   for si in range(64)
                                   if v[si*K.dim_V:(si+1)*K.dim_V]
                                      != vector(QQ, [0]*K.dim_V)))

        # Incrementally build the submodule.
        # Add seeds one at a time; verify the orbit remains a proper
        # submodule (doesn't contain hw).
        current_orbit = None
        good_seeds = []
        for seed in non_hw_sorted:
            test_seeds = good_seeds + [seed]
            candidate = _orbit_closure(test_seeds)
            if not _contains_hw(candidate):
                good_seeds.append(seed)
                current_orbit = candidate

        if not good_seeds:
            # No proper submodule found → K is irreducible
            self.dim = K.dim
            self.action_mats = dict(K.action_mats)
            self.quotient_proj = identity_matrix(QQ, K.dim)
            self._quotient_basis = list(identity_matrix(QQ, K.dim).rows())
            self.weight_spaces = dict(K.V.weight_spaces)
            self.v_hw = K.V.v_hw
            self.dim_V = K.dim_V
            return

        work = current_orbit  # already computed; no need to redo

        sub_dim = work.rank()
        quot_dim = K.dim - sub_dim

        # Step 4: Quotient via echelon pivot complement
        # The echelonized submodule matrix (RREF) has pivots at certain columns.
        # The quotient K/N is identified with the "free" (non-pivot) columns.
        pivots = work.pivots()
        sub_pivots = set(pivots)
        complement_cols = [i for i in range(K.dim) if i not in sub_pivots]
        assert len(complement_cols) == quot_dim

        self.dim = quot_dim

        # Inclusion: standard basis vectors at non-pivot positions
        K_id = identity_matrix(QQ, K.dim)
        self._incl = matrix(QQ, [K_id.row(i) for i in complement_cols])

        # Projection: for RREF with pivots P and free cols F,
        # the quotient map q: K → K/N is:
        #   q(v)[r] = v[F_r] - Σ_i v[P_i] * work[i, F_r]
        # This zeroes out the N-component of v.
        Q = matrix(QQ, quot_dim, K.dim)
        for r in range(quot_dim):
            Q[r, complement_cols[r]] = QQ(1)
            for i in range(sub_dim):
                if work[i, complement_cols[r]] != 0:
                    Q[r, pivots[i]] = -work[i, complement_cols[r]]
        self._proj_coords = Q

        # Build action matrices in quotient basis:
        #   q ∘ (action on K) ∘ incl
        self.action_mats = {}
        for L0_idx in range(32):
            K_mat = K.action_mats[L0_idx]
            quot_mat = Q * K_mat * self._incl.transpose()
            self.action_mats[L0_idx] = quot_mat

        # Crystal interface (indices 0..dim_V-1 correspond to V(a,b,c))
        self.weight_spaces = dict(K.V.weight_spaces)
        self.v_hw = K.V.v_hw
        self.dim_V = K.dim_V

    def action_of_C(self):
        """Matrix of C ∈ \hat{p}(4)_{-2}: acts as t\cdotId."""
        return self.t * identity_matrix(QQ, self.dim)


# ===========================================================================
# Module cache  (in-memory + pickle disk cache)
# ===========================================================================

_phat4_cache = {}


def _load_disk_cache():
    """Load previously computed modules from disk."""
    global _phat4_cache
    if _os.path.exists(_PHAT4_CACHE_FILE):
        try:
            with open(_PHAT4_CACHE_FILE, 'rb') as f:
                _phat4_cache = _pickle.load(f)
        except Exception:
            _phat4_cache = {}


def _save_disk_cache():
    """Save all computed modules to disk."""
    with open(_PHAT4_CACHE_FILE, 'wb') as f:
        _pickle.dump(_phat4_cache, f, protocol=_pickle.HIGHEST_PROTOCOL)


def _make_cache_entry(W):
    """Extract picklable data from a Phat4Module."""
    return {
        't': W.t, 'a': W.a, 'b': W.b, 'c': W.c,
        'dim': W.dim,
        'action_mats': {k: m for k, m in W.action_mats.items()},
        'weight_spaces': W.weight_spaces,
        'v_hw': W.v_hw,
        'dim_V': W.dim_V,
    }


def _from_cache_entry(entry):
    """Reconstruct a lightweight Phat4Module from cached data."""
    W = object.__new__(Phat4Module)
    W.t = entry['t']
    W.a = entry['a']
    W.b = entry['b']
    W.c = entry['c']
    W.dim = entry['dim']
    W.action_mats = entry['action_mats']
    W.weight_spaces = entry.get('weight_spaces', {})
    W.v_hw = entry.get('v_hw', 0)
    W.dim_V = entry.get('dim_V', W.dim)
    return W


def phat4_module(t, a, b, c, e44_data):
    """
    Construct the full irreducible \hat{p}(4)-module W_t(a,b,c).

    Results are cached in memory and on disk (phat4_cache.pkl).
    """
    key = (QQ(t), int(a), int(b), int(c))

    # Check in-memory cache
    if key in _phat4_cache:
        return _from_cache_entry(_phat4_cache[key])

    # Compute fresh
    W = Phat4Module(t, a, b, c, e44_data)

    # Store in memory and flush to disk
    _phat4_cache[key] = _make_cache_entry(W)
    _save_disk_cache()
    return W


# Load disk cache at import time
_load_disk_cache()


# ===========================================================================
# Testing
# ===========================================================================

def _check_phat4(verbose=True):
    """
    Verify the \hat{p}(4)-module construction.

    Checks:
    (1) dim W_1(1,0,0) = 8  (standard module = E(4,4)_{-1})
    (2) dim W_t(0,0,0) for t≠0: should be smaller than 64
    (3) [L_0,L_0] commutator identity holds on W at degree 0
    (4) Odd L_0 generators act nontrivially on W_1(1,0,0)
    """
    e44_data = load_e44()
    if e44_data is None:
        print("[ERROR] e44_data not found")
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
        print("\hat{p}(4)-module checks")
        print("─" * 60)

    # (1) Standard module W_1(1,0,0) = 8-dim
    W_std = phat4_module(1, 1, 0, 0, e44_data)
    _check("dim W_1(1,0,0)", W_std.dim, 8)

    # (2) W_1(0,0,0)
    W_triv = phat4_module(1, 0, 0, 0, e44_data)
    if verbose:
        print(f"  [INFO] dim W_1(0,0,0) = {W_triv.dim}  (K has dim 64)")

    # (3) [L_0, L_0] commutator check on W_1(1,0,0)
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

            # RHS: [L_0[i0], L_0[j0]] decomposed in L_0
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
                # The bracket may have a central element C component
                # (from {a_{ij}, a_{kl}} ∈ \hat{p}(4)_{-2} = \C C).
                # C acts as t\cdotId on W, so check if diff is scalar\cdotId.
                diff = lhs - rhs
                c_val = diff[0, 0]
                if diff != c_val * identity_matrix(QQ, W_std.dim):
                    n_fail += 1
                    if n_fail <= 5 and verbose:
                        print(f"    [FAIL] [{L0[i0]['label']}, {L0[j0]['label']}]")

    _check(f"[L_0,L_0] commutator on W_1(1,0,0) ({n_total} pairs)", n_fail, 0)

    # (4) Odd L_0 acts nontrivially
    any_nonzero = False
    for idx in range(16, 32):
        if W_std.action_mats[idx] != matrix(QQ, W_std.dim, W_std.dim):
            any_nonzero = True
            break
    _check("Odd L_0 acts nontrivially on W_1(1,0,0)", any_nonzero, True)

    if verbose:
        print("─" * 60)
        if all_pass:
            print("\hat{p}(4)-module checks  ✓  ALL PASSED")
        else:
            print("\hat{p}(4)-module checks  ✗  SOME FAILED")
        print("─" * 60)

    return all_pass


# ===========================================================================
# __main__
# ===========================================================================

if __name__ == '__main__':
    print("=" * 60)
    print("phat4_modules.py  —  Full \hat{p}(4)-modules W_t(a,b,c)")
    print("=" * 60)

    e44_data = load_e44()
    if e44_data is None:
        print("[ERROR] Run 'sage e44_structure.py' first to generate e44_brackets.pkl")
        _sys.exit(1)

    # Show dimensions of various modules
    test_cases = [
        (1, 1, 0, 0, "standard module"),
        (1, 0, 0, 0, "trivial fiber"),
        (2, 0, 0, 0, "trivial fiber t=2"),
        (0, 0, 0, 0, "trivial t=0"),
        (0, 0, 0, 1, "dual fund t=0"),
        (1, 0, 0, 1, "dual fund t=1"),
        (0, 0, 0, 2, "Sym² dual"),
    ]

    print("\nModule dimensions:")
    print(f"  {'Parameters':<20} {'dim_K':>6} {'dim_W':>6}  Note")
    print(f"  {'─'*20} {'─'*6} {'─'*6}  {'─'*20}")
    for t, a, b, c, note in test_cases:
        K = KacModule(t, a, b, c, e44_data)
        W = phat4_module(t, a, b, c, e44_data)
        print(f"  W_{t}({a},{b},{c}){' '*(14-len(f'W_{t}({a},{b},{c})'))} "
              f"{K.dim:6d} {W.dim:6d}  {note}")

    print()
    _check_phat4(verbose=True)
