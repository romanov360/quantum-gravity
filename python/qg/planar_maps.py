#!/usr/bin/env python3
"""G15. Planar maps, Fortuin-Kasteleyn weights and the SU(2)_k loop weight: counts that sit beside the spin-network Catalan numbers.

  * Rooted planar maps with n edges, by brute force over vertex rotations (genus 0, connected), against Tutte's
        M_n = 2 * 3^n (2n)! / (n! (n+2)!)  = 2, 9, 54, 378, ...
    and the same enumeration weighted by spanning trees (Kirchhoff) against Mullin's  Cat_n Cat_{n+1} = 2, 10, 70, 588, ...
  * The FK parametrisation q = 2 + 2 cos(pi gamma^2/2), kappa' = 16/gamma^2, c = 25 - 6 Q^2, Q = 2/gamma + gamma/2, and the
    KPZ string susceptibility gamma_str = 1 - 4/gamma^2 (exact rational values at q = 0, 1, 2, 3, 4).
  * Invariants of N spin-1/2 legs at level k of U_q(sl_2): closed walks of length N on the Dynkin diagram A_{k+1};
    they equal Cat_{N/2} for N <= 2k and grow as (2 cos(pi/(k+2)))^N = (sqrt q_k)^N, q_k = 4 cos^2(pi/(k+2)) = 1, 2, 3, 4 at k = 1, 2, 4, infinity.

Not claimed: any map from a spin network to a decorated planar map (the numbers agree because the same Temperley-Lieb algebra with loop
weight sqrt q acts on both), or any statement about the scaling limit.

Run: python3 python/qg/planar_maps.py
"""
from __future__ import annotations

import itertools
import math
from fractions import Fraction as Fr

import numpy as np
import sympy as sp


def catalan(n: int) -> int:
    return math.comb(2 * n, n) // (n + 1)


def tutte(n: int) -> int:
    return 2 * 3**n * math.factorial(2 * n) // (math.factorial(n) * math.factorial(n + 2))


def _cycles(perm):
    seen, out = [False] * len(perm), []
    for s in range(len(perm)):
        if not seen[s]:
            c, x = [], s
            while not seen[x]:
                seen[x] = True
                c.append(x)
                x = perm[x]
            out.append(c)
    return out


def rooted_counts(n: int) -> tuple[Fr, Fr]:
    """(maps, tree-weighted maps) with n edges: sum over vertex rotations sigma in S_{2n} for the fixed involution alpha = (0 1)(2 3)...;
    genus 0 (V - n + F = 2) and connected; rooted count = labelled count / (2^(n-1) (n-1)!)."""
    h = 2 * n
    alpha = [i ^ 1 for i in range(h)]
    count = trees = 0
    for sigma in itertools.permutations(range(h)):
        sc = _cycles(sigma)
        V = len(sc)
        F = len(_cycles([sigma[alpha[i]] for i in range(h)]))
        if V - n + F != 2:
            continue
        vert = [0] * h
        for vi, c in enumerate(sc):
            for x in c:
                vert[x] = vi
        parent = list(range(V))

        def find(u):
            while parent[u] != u:
                parent[u] = parent[parent[u]]
                u = parent[u]
            return u

        for i in range(0, h, 2):
            parent[find(vert[i])] = find(vert[i + 1])
        if len({find(u) for u in range(V)}) != 1:
            continue
        count += 1
        L = np.zeros((V, V))
        for i in range(0, h, 2):
            u, w = vert[i], vert[i + 1]
            if u != w:
                L[u, u] += 1; L[w, w] += 1; L[u, w] -= 1; L[w, u] -= 1
        trees += int(round(np.linalg.det(L[1:, 1:]))) if V > 1 else 1
    scale = Fr(1, 2 ** (n - 1) * math.factorial(n - 1))
    return count * scale, trees * scale


def mullin_shuffles(n: int) -> int:
    return sum(math.comb(2 * n, 2 * k) * catalan(k) * catalan(n - k) for k in range(n + 1))


# ---- FK dictionary (exact) ----------------------------------------------------------------------------------------
def fk_row(q: int):
    """(kappa', gamma^2, c, gamma_str) at q in {0,1,2,3,4}, from sqrt q = -2 cos(4 pi/kappa'), kappa' in [4, 8]."""
    kappa = {0: sp.Integer(8), 1: sp.Integer(6), 2: sp.Rational(16, 3), 3: sp.Rational(24, 5), 4: sp.Integer(4)}[q]
    assert sp.simplify(-2 * sp.cos(4 * sp.pi / kappa) - sp.sqrt(q)) == 0
    g2 = 16 / kappa
    c = (3 * kappa - 8) * (6 - kappa) / (2 * kappa)
    gs = sp.simplify((c - 1 - sp.sqrt((1 - c) * (25 - c))) / 12)
    return kappa, g2, c, gs


# ---- invariants at level k ------------------------------------------------------------------------------------------
def closed_walks(k: int, length: int) -> int:
    """closed walks from node 0 on the path graph with nodes 0..k (Dynkin diagram A_{k+1}); = dim Inv_k(V_{1/2}^{length})."""
    v = [0] * (k + 1)
    v[0] = 1
    for _ in range(length):
        v = [(v[i - 1] if i > 0 else 0) + (v[i + 1] if i < k else 0) for i in range(k + 1)]
    return v[0]


def q_of_level(k: int) -> float:
    return 4 * math.cos(math.pi / (k + 2)) ** 2


if __name__ == "__main__":
    for n in (1, 2, 3):
        print(n, rooted_counts(n), tutte(n), catalan(n) * catalan(n + 1))
    for q in range(5):
        print("q =", q, fk_row(q))
    print([closed_walks(6, 2 * m) for m in range(7)], [catalan(m) for m in range(7)])
