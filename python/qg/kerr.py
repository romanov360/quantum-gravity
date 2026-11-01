#!/usr/bin/env python3
"""G14. Kerr horizon identities and the Penrose gap, in exact algebra (G = c = 1).

Parametrise Kerr by the Boyer-Lindquist radii r_+ > r_- > 0:  M = (r_+ + r_-)/2,  a^2 = r_+ r_-,  J = M a.  Then

    A_pm = 4 pi (r_pm^2 + a^2) = 8 pi M r_pm,      kappa_pm = (r_pm - r_mp) / (2 (r_pm^2 + a^2)),   Omega_pm = a/(r_pm^2 + a^2)
    Smarr         M = 2 Omega_pm J + kappa_pm A_pm / (4 pi)                       at both horizons
    first law     dM = kappa_pm dA_pm / (8 pi) + Omega_pm dJ                        as one-forms in (r_+, r_-)
    Penrose gap   M^2 - A_+/(16 pi) = M r_- / 2   (zero exactly for Schwarzschild)
    area product  A_+ A_- = (8 pi J)^2,    A_+ - 8 pi |J| = 8 pi M (r_+ - a) >= 0

These are the Kerr quantities used by the near-Kerr strong-cosmic-censorship papers (the inner surface gravity |kappa_-|
fixes the blueshift rate).  The last function composes the bound A >= 8 pi |J| with the spin-1/2 area spectrum: a heuristic,
labelled as such (nothing connects a classical Kerr area to a spin-network eigenvalue beyond the shared symbol A).

Run: python3 python/qg/kerr.py
"""
from __future__ import annotations

import math

import sympy as sp

rp, rm = sp.symbols("r_p r_m", positive=True)
M = (rp + rm) / 2
a = sp.sqrt(rp * rm)
J = M * a
A = {"+": 4 * sp.pi * (rp**2 + a**2), "-": 4 * sp.pi * (rm**2 + a**2)}
kappa = {"+": (rp - rm) / (2 * (rp**2 + a**2)), "-": (rm - rp) / (2 * (rm**2 + a**2))}
Omega = {"+": a / (rp**2 + a**2), "-": a / (rm**2 + a**2)}


def simplify0(e) -> bool:
    return sp.simplify(e) == 0


def horizon_identities() -> dict[str, bool]:
    out = {}
    for s in "+-":
        out[f"area {s}"] = simplify0(A[s] - 8 * sp.pi * M * (rp if s == "+" else rm))
        out[f"smarr {s}"] = simplify0(M - (2 * Omega[s] * J + kappa[s] * A[s] / (4 * sp.pi)))
        out[f"first law {s}"] = all(
            simplify0(sp.diff(M, x) - (kappa[s] / (8 * sp.pi) * sp.diff(A[s], x) + Omega[s] * sp.diff(J, x))) for x in (rp, rm))
    out["|kappa_-|/kappa_+ = r_+/r_-"] = simplify0(-kappa["-"] / kappa["+"] - rp / rm)
    return out


def penrose_identities() -> dict[str, bool]:
    Mirr2 = A["+"] / (16 * sp.pi)
    return {
        "Penrose gap M r_-/2": simplify0(M**2 - Mirr2 - M * rm / 2),
        "Christodoulou-Ruffini": simplify0(M**2 - (Mirr2 + J**2 / (4 * Mirr2))),
        "area product": simplify0(A["+"] * A["-"] - (8 * sp.pi * J) ** 2),
        "A_+ - 8 pi J": simplify0(A["+"] - 8 * sp.pi * J - 8 * sp.pi * M * (rp - a)),
    }


def limits() -> dict[str, object]:
    x = sp.Symbol("x", positive=True)
    # a -> 0 at M = 1: r_+ = 2 - x, r_- = x
    km = kappa["-"].subs({rp: 2 - x, rm: x}, simultaneous=True)
    a2 = (2 - x) * x
    return {
        "a^2 |kappa_-| -> M (a -> 0)": sp.limit(a2 * (-km), x, 0),
        "kappa_+ -> 0 (a -> M)": sp.limit(kappa["+"].subs({rp: 1 + x, rm: 1 - x}, simultaneous=True), x, 0),
    }


def punctures_needed(J_hbar: float, gamma: float | None = None) -> float:
    """Heuristic: N spin-1/2 punctures have area 4 sqrt(3) pi gamma N l_P^2; A >= 8 pi |J| requires N >= 2 |J| / (sqrt(3) gamma).
    With gamma = ln2/(pi sqrt 3) this is 2 pi |J| / ln 2."""
    g = math.log(2) / (math.pi * math.sqrt(3)) if gamma is None else gamma
    return 2 * abs(J_hbar) / (math.sqrt(3) * g)


if __name__ == "__main__":
    for d in (horizon_identities(), penrose_identities()):
        for k, v in d.items():
            print(f"{k:35s} {v}")
    print(limits())
    print("punctures for J = 1/2, 1:", punctures_needed(0.5), punctures_needed(1.0))
