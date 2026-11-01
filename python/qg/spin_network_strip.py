#!/usr/bin/env python3
"""G16. SU(2) Kogut-Susskind Hamiltonian on an open strip of N plaquettes, in the spin-network basis (N = 3, 4, 5).

H = (g^2/2) sum_links C(j) + (1/g^2) sum_p (2 - W_p),  C(j) = j(j+1),  W_p = chi_{1/2}(U_p).

States are spin networks on a ladder of squares; labels in k = 2j: j_a, j_c (the three outer links of the end plaquettes), t_i, b_i
for the interior plaquettes i = 2..N-1, and the rungs v_1..v_{N-1}: 3N-3 labels.  W_p shifts every link on the boundary of p by +-1/2
with one signed 6j symbol per vertex on the loop (two for an end plaquette, four for an interior one), times sqrt(prod dims) and a phase
sigma_p.  The phase is the rule that was fixed for N = 3 against an independent Clebsch-Gordan construction, extended by a factor
(-1)^((t_q - b_q)/2) for each interior neighbour q.  That extension is a conjecture fitted by comparison, NOT a derivation.

Checks (python/tests): N = 3 ground state at g^2 = 2 is 2.753963; N = 4: 3.67192; N = 5: 4.5899 (these agree with an independent
Peter-Weyl calculation, wavefunctions of the N loop holonomies, with no 6j symbols and no phases, to 2e-6 at N = 4 and 1e-3 at N = 5;
that engine is in github.com/brentharts/openai_math_yang_mills, sympy_models/models/su2_peter_weyl.py).
Not tested here: the weak-coupling gap 2 sqrt(2 - cos(pi/(N+1))) for N = 4, 5 (the energy cut converges too slowly).

Run: python3 python/qg/spin_network_strip.py
"""
from __future__ import annotations

import itertools
from functools import lru_cache

import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as sla
from sympy import Rational
from sympy.physics.wigner import wigner_6j


def tri(a, b, c):
    return abs(a - b) <= c <= a + b and (a + b + c) % 2 == 0


@lru_cache(maxsize=None)
def _sixj(a, b, c, d, e, f):
    """{a/2 b/2 c/2; d/2 e/2 f/2} as a float (arguments are the doubled spins)."""
    return float(wigner_6j(*(Rational(x, 2) for x in (a, b, c, d, e, f))))


def vert(x, xp, y, yp, z):
    """signed 6j at a loop vertex: {z/2 yp/2 xp/2; 1/2 x/2 y/2} (the module su2_lattice_strip.half6j_exact(z, yp, xp, x, y))"""
    if not (tri(z, yp, xp) and tri(z, y, x) and tri(1, yp, y) and tri(1, x, xp)):
        return 0.0
    return _sixj(z, yp, xp, 1, x, y)


class Net:
    def __init__(self, N):
        self.N = N
    def unpack(self, s):
        """s = (ja, jc, t2..t_{N-1}, b2..b_{N-1}, v1..v_{N-1}) in k -> arrays t[1..N], b[1..N], v[0..N] (index = plaquette/rung)."""
        N = self.N
        ja, jc = s[0], s[1]
        mid = N - 2
        t = [None, ja] + list(s[2:2 + mid]) + [jc]
        b = [None, ja] + list(s[2 + mid:2 + 2 * mid]) + [jc]
        v = [ja] + list(s[2 + 2 * mid:2 + 2 * mid + N - 1]) + [jc]
        return t, b, v
    def pack(self, t, b, v):
        N = self.N
        return (t[1], t[N]) + tuple(t[2:N]) + tuple(b[2:N]) + tuple(v[1:N])
    def valid(self, s):
        t, b, v = self.unpack(s)
        N = self.N
        if min(s) < 0: return False
        for i in range(1, N):
            if not tri(t[i], t[i + 1], v[i]) or not tri(b[i], b[i + 1], v[i]):
                return False
        return True
    def e_el(self, s):
        t, b, v = self.unpack(s)
        N = self.N
        c = lambda k: k * (k + 2) / 4
        e = 3 * c(t[1]) + 3 * c(t[N])
        e += sum(c(t[i]) + c(b[i]) for i in range(2, N))
        e += sum(c(v[i]) for i in range(1, N))
        return e
    def states(self, lam):
        N = self.N
        kmax = int(np.sqrt(4 * lam)) + 2
        out = []
        n = 3 * N - 3
        def rec(prefix):
            # brute force with pruning on energy only; validity checked at end (small spaces)
            if len(prefix) == n:
                if self.valid(prefix) and self.e_el(prefix) <= lam + 1e-12:
                    out.append(prefix)
                return
            for k in range(kmax):
                new = prefix + (k,)
                # crude energy prune: ja,jc weighted 3, others 1
                w = [3, 3] + [1] * (n - 2)
                e = sum(wi * (ki * (ki + 2) / 4) for wi, ki in zip(w, new))
                if e <= lam + 1e-12:
                    rec(new)
        rec(())
        return out
    def amp(self, p, s, s2, rule):
        t, b, v = self.unpack(s); t2, b2, v2 = self.unpack(s2)
        N = self.N
        if p == 1 or p == N:
            if p == 1:
                x, xp, y, yp = t[1], t2[1], v[1], v2[1]
                zt, zb = t[2], b[2]
                ztn, zbn = t[2], b[2]
            else:
                x, xp, y, yp = t[N], t2[N], v[N - 1], v2[N - 1]
                zt, zb = t[N - 1], b[N - 1]
            val = vert(x, xp, y, yp, zt) * vert(x, xp, y, yp, zb) * np.sqrt((x + 1) * (xp + 1) * (y + 1) * (yp + 1))
            return val * rule('end', p, s, s2, self)
        val = (vert(t[p], t2[p], v[p - 1], v2[p - 1], t[p - 1]) * vert(t[p], t2[p], v[p], v2[p], t[p + 1])
               * vert(b[p], b2[p], v[p], v2[p], b[p + 1]) * vert(b[p], b2[p], v[p - 1], v2[p - 1], b[p - 1]))
        dims = (t[p] + 1) * (t2[p] + 1) * (b[p] + 1) * (b2[p] + 1) * (v[p - 1] + 1) * (v2[p - 1] + 1) * (v[p] + 1) * (v2[p] + 1)
        return val * np.sqrt(dims) * rule('mid', p, s, s2, self)
    def transitions(self, p, s):
        t, b, v = self.unpack(s)
        N = self.N
        out = []
        if p == 1 or p == N:
            idx = [0, 2 + 2 * (N - 2)] if p == 1 else [1, 2 + 2 * (N - 2) + N - 2]
            for d in itertools.product((1, -1), repeat=2):
                s2 = list(s)
                for i, dd in zip(idx, d): s2[i] += dd
                s2 = tuple(s2)
                if min(s2) >= 0 and self.valid(s2): out.append(s2)
        else:
            mid = N - 2
            idx = [2 + (p - 2), 2 + mid + (p - 2)]  # t_p, b_p
            idx += [2 + 2 * mid + (p - 2), 2 + 2 * mid + (p - 1)]   # v_{p-1}, v_p
            for d in itertools.product((1, -1), repeat=4):
                s2 = list(s)
                for i, dd in zip(idx, d): s2[i] += dd
                s2 = tuple(s2)
                if min(s2) >= 0 and self.valid(s2): out.append(s2)
        return out
    def hamiltonian(self, lam, g2, rule):
        S = self.states(lam)
        ix = {s: i for i, s in enumerate(S)}
        n = len(S)
        rows, cols, vals = [], [], []
        for s in S:
            for p in range(1, self.N + 1):
                for s2 in self.transitions(p, s):
                    if s2 in ix:
                        a = self.amp(p, s, s2, rule)
                        if a != 0:
                            rows.append(ix[s2]); cols.append(ix[s]); vals.append(a)
        W = sps.csr_matrix((vals, (rows, cols)), shape=(n, n))
        H = sps.diags([0.5 * g2 * self.e_el(s) + 2 * self.N / g2 for s in S]) - W / g2
        return H.tocsr(), W, S



def rule_general(kind, p, s, s2, net):
    """end plaquette: (-1)^((t-b)/2) of its interior neighbour; interior: same/opposite shift sign times the neighbours' factors."""
    t, b, v = net.unpack(s); t2, b2, v2 = net.unpack(s2)
    N = net.N
    if kind == 'end':
        nb = 2 if p == 1 else N - 1
        return -1.0 if ((t[nb] - b[nb]) // 2) % 2 else 1.0
    sg = 1.0 if (t2[p] - t[p]) * (b2[p] - b[p]) > 0 else -1.0
    for q in (p - 1, p + 1):
        if 1 < q < N:
            sg *= -1.0 if ((t[q] - b[q]) // 2) % 2 else 1.0
    return sg




def rule_general(kind, p, s, s2, net):
    t, b, v = net.unpack(s); t2, b2, v2 = net.unpack(s2)
    N = net.N
    if kind == 'end':
        nb = 2 if p == 1 else N - 1
        return -1.0 if ((t[nb] - b[nb]) // 2) % 2 else 1.0
    sg = 1.0 if (t2[p] - t[p]) * (b2[p] - b[p]) > 0 else -1.0
    for q in (p - 1, p + 1):
        if 1 < q < N:
            sg *= -1.0 if ((t[q] - b[q]) // 2) % 2 else 1.0
    return sg


def lowest(N, lam, g2, k=6):
    """lowest k levels of H_N at coupling g2 with electric-energy cut lam (in units of C(j))."""
    H, W, S = Net(N).hamiltonian(lam, g2, rule_general)
    if len(S) <= k + 2:
        return np.sort(np.linalg.eigvalsh(H.toarray()))[:k], len(S)
    return np.sort(sla.eigsh(H, k=k, sigma=0, which='LM', return_eigenvectors=False)), len(S)


if __name__ == "__main__":
    for N, lam in ((3, 20), (4, 16), (5, 12)):
        w, n = lowest(N, lam, 2.0, k=4)
        print(N, n, np.round(w, 5))
