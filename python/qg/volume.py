#!/usr/bin/env python3
"""G17. The quantum tetrahedron: prefactor, closed-form matrix elements, and the semiclassical limit.

Hartshorn (qgravity paper, Sec. 4) computes the triple product Q = eps_abc J1a J2b J3c on four-valent intertwiner
spaces and finds that lambda_max(Q) / ((4 / 3 sqrt 3) [j (j+1)]^(3/2)) rises from 0.794 (j = 3/2) to 0.921 (j = 8),
"a trend, not a proof of a semiclassical limit". He leaves the volume prefactor open.

1. Matrix elements. In the coupled basis |k> ((j1 j2) -> k), Q is purely imaginary and tridiagonal:
       <k-1| Q |k> = -i a_k,  a_k = sqrt(D(j1, j2, k) D(j3, j4, k)) / (4 sqrt(4 k^2 - 1)),
       D(a, b, k) = (a + b + k + 1)(-a + b + k)(a - b + k)(a + b - k + 1).
   This closed form is checked against brute-force tensor-product matrices (`brute_force_spectrum`) for
   (1/2)^4, 1^4, (3/2)^4, 2^4 and (1, 1, 1/2, 1/2), with the phase convention up to a diagonal unitary.
2. Prefactor. A classical tetrahedron with face area vectors A_i (sum zero) has V^2 = (2/9) |A1 . (A2 x A3)|. With
   A_i = 8 pi gamma l_P^2 J_i:   V = (sqrt 2 / 3) (8 pi gamma l_P^2)^(3/2) sqrt|Q|.
   A regular tetrahedron with face areas |A| has |A1 . (A2 x A3)| = (4 / (3 sqrt 3)) |A|^3, the denominator of
   Hartshorn's ratio.
3. The ratio tends to 1. Numerically (to j = 65536), j (1 - ratio) -> 3/4, i.e.
       lambda_max = (4 / (3 sqrt 3)) [j (j+1)]^(3/2) (1 - 3/(4 j) + O(j^-2))
   (j (1 - ratio) = 0.74894, 0.74974, 0.74993, 0.74998 at j = 1024, 4096, 16384, 65536). A conjecture from numerics;
   the Bohr-Sommerfeld analysis of Bianchi-Haggard should give it.

Run: python3 python/qg/volume.py
"""
from __future__ import annotations

import math

import numpy as np
from scipy import linalg


def D(a, b, k):
    return (a + b + k + 1) * (-a + b + k) * (a - b + k) * (a + b - k + 1)


def tridiagonal_a(j1, j2, j3, j4):
    kmin = max(abs(j1 - j2), abs(j3 - j4))
    kmax = min(j1 + j2, j3 + j4)
    ks = [kmin + n for n in range(int(round(kmax - kmin)) + 1)]
    a = [math.sqrt(max(D(j1, j2, k), 0) * max(D(j3, j4, k), 0)) / (4 * math.sqrt(4 * k * k - 1)) for k in ks[1:]]
    return ks, np.array(a)


def spectrum(j1, j2, j3, j4):
    """Eigenvalues of Q: those of the Hermitian tridiagonal matrix with zero diagonal and off-diagonal i a_k,
    equal to the eigenvalues of the real symmetric tridiagonal matrix with off-diagonal a_k (unitary diagonal
    conjugation by i^n)."""
    ks, a = tridiagonal_a(j1, j2, j3, j4)
    if len(a) == 0:
        return np.zeros(1)
    return linalg.eigh_tridiagonal(np.zeros(len(ks)), a, eigvals_only=True)


def brute_force_spectrum(js):
    def spin(j):
        d = int(round(2 * j + 1))
        m = np.array([j - k for k in range(d)])
        Jp = np.zeros((d, d), complex)
        for k in range(1, d):
            Jp[k - 1, k] = math.sqrt(j * (j + 1) - m[k] * (m[k] + 1))
        return [(Jp + Jp.conj().T) / 2, (Jp - Jp.conj().T) / 2j, np.diag(m).astype(complex)]

    dims = [int(round(2 * j + 1)) for j in js]

    def embed(o, i):
        M = np.array([[1]], complex)
        for k, d in enumerate(dims):
            M = np.kron(M, o if k == i else np.eye(d))
        return M

    J = [[embed(o, i) for o in spin(j)] for i, j in enumerate(js)]
    Jt = [sum(J[i][a] for i in range(4)) for a in range(3)]
    w, v = np.linalg.eigh(sum(x @ x for x in Jt))
    P = v[:, abs(w) < 1e-8]
    eps = [(0, 1, 2, 1), (1, 2, 0, 1), (2, 0, 1, 1), (0, 2, 1, -1), (2, 1, 0, -1), (1, 0, 2, -1)]
    Q = sum(s * J[0][a] @ J[1][b] @ J[2][c] for a, b, c, s in eps)
    return np.sort(np.linalg.eigvalsh(P.conj().T @ Q @ P))


def ratio_regular(j):
    lam = spectrum(j, j, j, j).max()
    return lam / (4 / (3 * math.sqrt(3)) * (j * (j + 1)) ** 1.5)


def main():
    for js in ((0.5,) * 4, (1,) * 4, (1.5,) * 4, (2,) * 4, (1, 1, 0.5, 0.5), (1.5, 1, 2, 0.5)):
        bf, cf = brute_force_spectrum(js), np.sort(spectrum(*js))
        print(f"{js}: closed form {np.round(cf, 6)}  max |diff| vs brute force {np.max(np.abs(bf - cf)):.1e}")
    print("Hartshorn's ratio lambda_max / ((4/3sqrt3)[j(j+1)]^(3/2)):")
    js = [1.5, 2, 3, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048]
    rs = [ratio_regular(j) for j in js]
    for j, r in zip(js, rs):
        print(f"  j = {j:<6}  ratio = {r:.6f}   1 - ratio = {1 - r:.3e}")
    big = np.array(js[-6:])
    slope = np.polyfit(np.log(big), np.log(1 - np.array(rs[-6:])), 1)[0]
    print(f"fitted 1 - ratio ~ j^({slope:.3f}) for j in [{big[0]}, {big[-1]}]")
    for j in (4096, 16384, 65536):
        print(f"  j (1 - ratio) at j = {j}: {j * (1 - ratio_regular(j)):.6f}   (conjecture: 3/4)")
    print("prefactor: V = (sqrt2/3) (8 pi gamma l_P^2)^(3/2) sqrt|Q|; smallest nonzero volume (j=1/2):",
          f"(sqrt2/3) sqrt(sqrt3/4) = {math.sqrt(2) / 3 * math.sqrt(math.sqrt(3) / 4):.6f} (8 pi gamma)^(3/2) l_P^3")


if __name__ == "__main__":
    main()
