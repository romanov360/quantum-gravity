#!/usr/bin/env python3
"""G13. The Barbero-Immirzi parameter from three puncture countings, rederived.

Hartshorn's quantum gravity paper fixes gamma_I = ln 2 / (pi sqrt 3) ~ 0.1274 from spin-1/2 punctures only and
quotes the all-spin value 0.2375 (Domagala-Lewandowski) without rederiving it. Here all three standard countings are
solved. A horizon of area A = 8 pi gamma l_P^2 sum_i sqrt(j_i (j_i + 1)) with S = A / (4 l_P^2) gives, at leading order:

  (a) spin-1/2 only, 2 states per puncture (ABCK / Rovelli):  gamma = ln 2 / (pi sqrt 3)
  (b) all spins, 2 states per puncture, projection m = +-j only (Domagala-Lewandowski, Meissner 2004):
          1 = sum_{j = 1/2, 1, 3/2, ...} 2 exp(-2 pi gamma sqrt(j (j + 1)))
  (c) all spins, 2j + 1 states per puncture (Ghosh-Mitra):
          1 = sum_j (2 j + 1) exp(-2 pi gamma sqrt(j (j + 1)))

Each f(gamma) = sum - 1 is strictly decreasing in gamma, so the root is bracketed by evaluating f with an explicit
tail bound: for j >= J, sqrt(j (j+1)) >= j + 1/2 - 1/(8 j) >= j, so the terms are at most c_j e^(-2 pi gamma j)
and the tail is a geometric(-type) series summed in closed form. Arithmetic is 60-digit mpmath.

Not claimed: which counting is right. That is a question about the quantum horizon (and, in newer work, about
the SU(2) Chern-Simons treatment), not about this arithmetic.
Two independent implementations live here and are checked against each other in the tests:
  * `solve(kind)` / `f_bounds`: bisection with the explicit tail bound above, 60-digit arithmetic (Romanov);
  * `gamma_spin_half`, `gamma_domagala_lewandowski`, `gamma_all_spins`, `renewal_sum`: Newton root finding, and
    `growth_rate`: the direct count N(x) = 1 + sum_j (2j+1) N(x - a_j) of puncture sequences of area <= x, whose
    growth rate ln N(x)/x approaches the root lambda = 2 pi gamma from below (Hartshorn).
Run: python3 python/qg/immirzi.py
"""
from __future__ import annotations

import math

import mpmath as _mpm
from mpmath import mpf, findroot, sqrt, log, pi, exp

_mpm.mp.dps = 60
mp = _mpm


def _terms(gamma, weight, J):
    s = mp.mpf(0)
    for k in range(1, 2 * J + 1):  # j = k/2
        j = mp.mpf(k) / 2
        s += weight(j) * mp.e ** (-2 * mp.pi * gamma * mp.sqrt(j * (j + 1)))
    return s


def _tail_bound(gamma, weight_kind, J):
    """Upper bound on sum over j > J (j = k/2, k > 2J) of w(j) e^(-2 pi gamma j)."""
    q = mp.e ** (-mp.pi * gamma)  # e^(-2 pi gamma (k/2)) = q^k
    k0 = 2 * J + 1
    if weight_kind == "two":
        return 2 * q ** k0 / (1 - q)
    # weight 2j + 1 = k + 1: sum_{k >= k0} (k + 1) q^k
    return q ** k0 * (k0 + 1 - k0 * q) / (1 - q) ** 2


def f_bounds(gamma, weight_kind, J=400):
    w = (lambda j: 2) if weight_kind == "two" else (lambda j: 2 * j + 1)
    s = _terms(gamma, w, J)
    return s - 1, s - 1 + _tail_bound(gamma, weight_kind, J)


def solve(weight_kind, lo=0.05, hi=0.5, digits=40):
    lo, hi = mp.mpf(lo), mp.mpf(hi)
    for _ in range(200):
        mid = (lo + hi) / 2
        fl, fu = f_bounds(mid, weight_kind)
        if fl > 0:
            lo = mid
        elif fu < 0:
            hi = mid
        else:
            break
        if hi - lo < mp.mpf(10) ** (-digits):
            break
    # certify: f(lo) > 0 > f(hi) including the tail
    assert f_bounds(lo, weight_kind)[0] > 0 and f_bounds(hi, weight_kind)[1] < 0
    return lo, hi


def spin_half():
    return mp.log(2) / (mp.pi * mp.sqrt(3))


def a_j(two_j: int) -> mpf:
    j = mpf(two_j) / 2
    return sqrt(j * (j + 1))


def renewal_sum(lam, jmax_two: int = 4000):
    """F(lam) = sum_{j=1/2,1,...} (2j+1) exp(-lam sqrt(j(j+1))); terms decay like (2j) e^{-lam j}."""
    return sum((mpf(tj) + 1) * exp(-lam * a_j(tj)) for tj in range(1, jmax_two + 1))


def gamma_spin_half() -> mpf:
    return log(2) / (pi * sqrt(3))


def gamma_all_spins() -> mpf:
    lam = findroot(lambda x: renewal_sum(x) - 1, mpf("1.5"))
    return lam / (2 * pi)


def gamma_domagala_lewandowski() -> mpf:
    """Domagala-Lewandowski counting (only the extreme projections m = +-j per puncture, degeneracy 2): sum_j 2 exp(-lam a_j) = 1."""
    lam = findroot(lambda x: 2 * sum(exp(-x * a_j(tj)) for tj in range(1, 4000)) - 1, mpf("1.2"))
    return lam / (2 * pi)


def growth_rate(x_max: float = 40.0, h: float = 0.005):
    """Discretised renewal equation  N(x) = 1 + sum_j (2j+1) N(x - a_j)  (number of puncture sequences of total area <= x, including the empty one).
    Returns ln N(x)/x at x_max, which tends to lambda from below as x -> infinity (N ~ C e^{lambda x})."""
    import numpy as np
    n = int(x_max / h) + 1
    N = np.ones(n)
    steps = [(int(round(math.sqrt((tj / 2) * (tj / 2 + 1)) / h)), tj + 1) for tj in range(1, int(2 * x_max) + 2)]
    steps = [(s, d) for s, d in steps if 0 < s < n]
    for i in range(n):
        acc = 1.0
        for s, d in steps:
            if i - s < 0:
                break
            acc += d * N[i - s]
        N[i] = acc
    return math.log(N[-1]) / (n * h), N



def main():
    print(f"(a) spin-1/2 only:            gamma = ln2/(pi sqrt3) = {mp.nstr(gamma_spin_half(), 20)}")
    lo, hi = solve("two")
    print(f"(b) all spins, 2 states (DL): gamma in [{mp.nstr(lo, 25)}, {mp.nstr(hi, 25)}]  "
          f"(Newton: {mp.nstr(gamma_domagala_lewandowski(), 15)})")
    lo2, hi2 = solve("dim")
    print(f"(c) all spins, 2j+1 states:   gamma in [{mp.nstr(lo2, 25)}, {mp.nstr(hi2, 25)}]  "
          f"(Newton: {mp.nstr(gamma_all_spins(), 15)})")
    r, _ = growth_rate(25.0, 0.01)
    print(f"direct count (c): ln N(x)/x at x = 25 -> {r:.5f}, below lambda = {float(2 * pi * gamma_all_spins()):.5f}")


if __name__ == "__main__":
    main()
