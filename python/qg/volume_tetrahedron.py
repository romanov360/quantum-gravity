#!/usr/bin/env python3
"""G17. The volume prefactor fixed by the classical tetrahedron.

On the space of invariants of four spin-j legs the triple product Q = J1 . (J2 x J3) is a Hermitian matrix. In the coupled basis
|(12)K, (34)K; 0> it is  Q_{KK'} = i (x_K - x_K') <K| J2.J3 |K'>,  x_K = [K(K+1) - 2 j(j+1)]/2, with J2.J3 diagonal in the (23)L channel and the
two channels related by the orthogonal recoupling  U_{KL} = sqrt((2K+1)(2L+1)) {j j K; j j L}.  Here the 6j symbols are computed in high precision
(Racah's sum, mpmath), so spins up to j = 40 are reachable (sympy's exact route is limited to j ~ 3).

Classical limit. Face-area vectors A_i = 8 pi gamma l_P^2 J_i of a tetrahedron sum to zero. For the regular one the three unit normals have
Gram determinant 16/27, so |n1.(n2 x n3)| = 4/(3 sqrt 3) and  V = (sqrt 2 / 3) sqrt|A1.(A2 x A3)|  holds exactly (both sides equal
(2 sqrt 2/3) 3^{-3/4} A^{3/2}).  That fixes the constant in  V^ = (sqrt 2/3)(8 pi gamma l_P^2)^{3/2} sqrt|Q|:  any other factor fails the classical limit.

Result (computed here). With r(j) = max|Q| / [(4/(3 sqrt 3)) (j(j+1))^{3/2}]:
    r(j) = 1 - 3/(4j) + O(1/j^2),   r -> 1   (Richardson fit through j = 8..40: limit 0.9999994, coefficient -0.74994).
So the largest eigenvalue of Q reproduces the regular tetrahedron's triple product at large spin, with the stated prefactor.
Smallest nonzero volume (four spin-1/2 legs, Q = +-sqrt3/4, gamma = ln2/(pi sqrt3)):  V = 1.777 l_P^3.

Not claimed: that Q is the right quantisation of the volume (Rovelli-Smolin and Ashtekar-Lewandowski operators differ at finite spin),
only that the prefactor common to all of them is fixed by the classical limit. Run: python3 python/qg/volume_tetrahedron.py
"""
from __future__ import annotations

import math
from functools import lru_cache

import mpmath as mp
import numpy as np

mp.mp.dps = 40


@lru_cache(None)
def _lf(n: int):
    return mp.log(mp.factorial(n))


def _delta(a, b, c):
    return mp.mpf(1) / 2 * (_lf((a + b - c) // 2) + _lf((a - b + c) // 2) + _lf((-a + b + c) // 2) - _lf((a + b + c) // 2 + 1))


def sixj(j1, j2, j3, j4, j5, j6):
    """{j1 j2 j3; j4 j5 j6} with all arguments doubled (k = 2j), by Racah's formula."""
    tr = [(j1, j2, j3), (j1, j5, j6), (j4, j2, j6), (j4, j5, j3)]
    for a, b, c in tr:
        if abs(a - b) > c or c > a + b or (a + b + c) % 2:
            return mp.mpf(0)
    pre = sum(_delta(*t) for t in tr)
    a = [(j1 + j2 + j3) // 2, (j1 + j5 + j6) // 2, (j4 + j2 + j6) // 2, (j4 + j5 + j3) // 2]
    b = [(j1 + j2 + j4 + j5) // 2, (j2 + j3 + j5 + j6) // 2, (j3 + j1 + j6 + j4) // 2]
    s = mp.mpf(0)
    for t in range(max(a), min(b) + 1):
        s += (-1) ** t * mp.exp(_lf(t + 1) - sum(_lf(t - x) for x in a) - sum(_lf(y - t) for y in b))
    return mp.exp(pre) * s


def q_matrix(k: int):
    """Q on Inv(V_j^{x4}), j = k/2, in the coupled basis (dimension 2j+1 = k+1)."""
    j = k / 2
    Ks = list(range(0, 2 * k + 1, 2))
    x = [(K / 2 * (K / 2 + 1) - 2 * j * (j + 1)) / 2 for K in Ks]
    n = len(Ks)
    U = mp.matrix(n, n)
    for a, K in enumerate(Ks):
        for b, L in enumerate(Ks):
            U[a, b] = mp.sqrt((K + 1) * (L + 1)) * sixj(k, k, K, k, k, L)
    Y = U * mp.diag(x) * U.T
    Q = mp.matrix(n, n)
    for a in range(n):
        for b in range(n):
            Q[a, b] = 1j * (x[a] - x[b]) * Y[a, b]
    return Q


def q_spectrum(k: int) -> np.ndarray:
    return np.sort(np.array([float(v) for v in mp.eigh(q_matrix(k), eigvals_only=True)]))


def ratio(k: int) -> float:
    """max|Q| / [(4/(3 sqrt3)) (j(j+1))^{3/2}]"""
    j = k / 2
    return float(max(abs(q_spectrum(k)))) / ((4 / (3 * math.sqrt(3))) * (j * (j + 1)) ** 1.5)


def regular_tetrahedron_triple_product() -> float:
    """|n1.(n2 x n3)| for unit outward normals of a regular tetrahedron = sqrt(det Gram), Gram off-diagonals -1/3."""
    G = np.array([[1, -1 / 3, -1 / 3], [-1 / 3, 1, -1 / 3], [-1 / 3, -1 / 3, 1]])
    return math.sqrt(np.linalg.det(G))


def volume_prefactor_identity() -> tuple[float, float]:
    """(sqrt2/3 sqrt(|A^3 triple|), direct volume of the regular tetrahedron with face area A=1)"""
    A = 1.0
    lhs = math.sqrt(2) / 3 * math.sqrt(regular_tetrahedron_triple_product() * A**3)
    a = math.sqrt(4 * A / math.sqrt(3))          # edge from face area sqrt3 a^2/4
    return lhs, a**3 / (6 * math.sqrt(2))


def smallest_volume(gamma: float | None = None) -> float:
    """Volume (units of l_P^3) of the four-spin-1/2 intertwiner with |Q| = sqrt3/4."""
    g = math.log(2) / (math.pi * math.sqrt(3)) if gamma is None else gamma
    return math.sqrt(2) / 3 * (8 * math.pi * g) ** 1.5 * math.sqrt(math.sqrt(3) / 4)


def richardson(ks=(16, 24, 32, 40, 48, 60, 72, 80), deg=4):
    J = np.array([k / 2 for k in ks])
    y = np.array([ratio(k) for k in ks])
    A = np.vstack([J ** -i for i in range(deg + 1)]).T
    return np.linalg.lstsq(A, y, rcond=None)[0]


if __name__ == "__main__":
    print("regular tetrahedron triple product", regular_tetrahedron_triple_product(), 4 / (3 * math.sqrt(3)))
    print("volume identity", volume_prefactor_identity())
    print("smallest volume / l_P^3", smallest_volume())
    for k in (3, 4, 6, 16, 24, 40):
        print(k, ratio(k))
    c = richardson((16, 24, 32, 40))
    print("fit r = c0 + c1/j + c2/j^2 + ...", c[:3])
