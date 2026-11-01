#!/usr/bin/env python3
"""G18. Weak-coupling gaps of small SU(2) plaquette clusters from one matrix: K = d^T d.

Hartshorn's SU(2) Hamiltonian clusters (one plaquette, the two-plaquette ladder, the three-plaquette strip and the
L-shape) have weak-coupling (harmonic) gaps sqrt 8, sqrt 6, 2 sqrt(2 - cos(pi/4)), and for the L-shape the same
number as the strip, "as for the strip". This module explains all of them by one formula.

Abelianise at weak coupling: each link carries an su(2) angle A_l, each plaquette the flux B_p = (d A)_p
(signed sum around p), and H = (g^2/2) sum_l E_l^2 + (1/g^2) sum_p (2 - W_p) with 2 - W_p ~ |B_p|^2 / 4. In
plaquette variables the kinetic metric is K = d d^T restricted to plaquettes, i.e. K_pq = (number of links of p)
on the diagonal and +-1 for each link p and q share. The normal-mode frequencies are omega^2 = mu / 2, mu in spec K,
each three times over (three su(2) components). The lowest gauge-invariant excitation is two quanta of the softest
mode, so
      Delta_{g -> 0} = 2 sqrt(mu_min / 2) = sqrt(2 mu_min),
with mu_min the smallest *nonzero* eigenvalue of K. Zero modes of K are Bianchi identities: they occur only for
closed surfaces (on the boundary of a cube, the six fluxes sum to zero).

For planar clusters with plaquettes oriented alike, K = 4 I - A, with A the adjacency matrix of the graph of
plaquettes that share a link. So Delta_{g -> 0} = 2 sqrt(2 - lambda_max(A)/2):
  * it depends only on the adjacency graph: the L-shape and the 3-strip are both the path P3, hence equal gaps;
  * the N-strip (path P_N, lambda_max = 2 cos(pi/(N+1))) gives 2 sqrt(2 - cos(pi/(N+1))) -> 2, Hartshorn's
    "prediction" for general N, here a theorem about the harmonic limit (the eigenvectors are sines; see
    QuantumGravity/ClusterModes.lean);
  * an L x L square block has lambda_max(A) = 4 cos(pi/(L+1)) -> 4, so the gap 2 sqrt(2 - 2 cos(pi/(L+1))) closes like 2 pi / (L+1):
    the free massless gluons of weak coupling.

Not claimed: anything beyond the harmonic limit; finite-g corrections need the full Hamiltonian (Hartshorn's code).
Run: python3 python/qg/cluster_gap.py
"""
from __future__ import annotations

import itertools
import math

import numpy as np


def planar_cluster(cells):
    """cells: set of (i, j) unit squares. Returns the signed incidence matrix d (plaquettes x links),
    counterclockwise orientation."""
    links = {}
    rows = []
    for (i, j) in sorted(cells):
        edges = [((i, j), (i + 1, j), +1), ((i + 1, j), (i + 1, j + 1), +1),
                 ((i, j + 1), (i + 1, j + 1), -1), ((i, j), (i, j + 1), -1)]
        row = {}
        for a, b, sgn in edges:
            key = (a, b)
            if key not in links:
                links[key] = len(links)
            row[links[key]] = sgn
        rows.append(row)
    d = np.zeros((len(rows), len(links)))
    for p, row in enumerate(rows):
        for l, sgn in row.items():
            d[p, l] = sgn
    return d


def cube_surface():
    """The six faces of the unit cube, outward orientation."""
    verts = list(itertools.product((0, 1), repeat=3))
    links = {}

    def link(a, b):
        key = tuple(sorted((a, b)))
        if key not in links:
            links[key] = len(links)
        return links[key], (1 if key == (a, b) else -1)

    faces = []
    for axis in range(3):
        for side in (0, 1):
            others = [ax for ax in range(3) if ax != axis]
            def v(u, w):
                p = [0, 0, 0]
                p[axis] = side
                p[others[0]] = u
                p[others[1]] = w
                return tuple(p)
            loop = [v(0, 0), v(1, 0), v(1, 1), v(0, 1)]
            if (side == 0) ^ (axis == 1):
                loop = loop[::-1]
            faces.append([link(loop[k], loop[(k + 1) % 4]) for k in range(4)])
    d = np.zeros((6, len(links)))
    for p, f in enumerate(faces):
        for l, s in f:
            d[p, l] += s
    return d


def gap(d):
    K = d @ d.T
    mu = np.linalg.eigvalsh(K)
    nz = mu[mu > 1e-9]
    return math.sqrt(2 * nz.min()), mu, K


def main():
    cases = {
        "one plaquette": {(0, 0)},
        "ladder (2)": {(0, 0), (1, 0)},
        "strip (3)": {(0, 0), (1, 0), (2, 0)},
        "L-shape (3)": {(0, 0), (1, 0), (0, 1)},
        "T-shape (4)": {(0, 0), (1, 0), (2, 0), (1, 1)},
        "2x2 block": {(0, 0), (1, 0), (0, 1), (1, 1)},
        "strip (10)": {(i, 0) for i in range(10)},
        "8x8 block": {(i, j) for i in range(8) for j in range(8)},
    }
    print(f"{'cluster':<15} {'gap (harmonic)':<16} {'formula':<12} spectrum of K")
    for name, cells in cases.items():
        d = planar_cluster(cells)
        Delta, mu, K = gap(d)
        A = 4 * np.eye(len(cells)) - K
        lam = np.linalg.eigvalsh(A).max()
        pred = 2 * math.sqrt(2 - lam / 2)
        spec = np.round(mu, 4) if len(mu) <= 6 else f"{len(mu)} values, min {mu.min():.4f}"
        print(f"{name:<15} {Delta:<16.6f} {pred:<12.6f} {spec}")
    print(f"strip formula 2 sqrt(2 - cos(pi/4)) = {2 * math.sqrt(2 - math.cos(math.pi / 4)):.6f};"
          f" ladder sqrt 6 = {math.sqrt(6):.6f}; plaquette sqrt 8 = {math.sqrt(8):.6f}")
    d = cube_surface()
    Delta, mu, K = gap(d)
    print(f"cube surface (6 faces, closed): spectrum of K {np.round(mu, 6)}; one zero mode (Bianchi: fluxes sum to 0),"
          f" gap {Delta:.6f} = 2 sqrt 2")


if __name__ == "__main__":
    main()
