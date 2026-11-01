# G17. The volume prefactor, fixed by the classical tetrahedron

Script: `python/qg/volume_tetrahedron.py`. Tests: `test_volume_classical`, `test_volume_semiclassical`.

On the intertwiner space of a four-valent vertex the triple product Q = J₁·(J₂×J₃) is Hermitian, and the volume is
V̂ = c·(8πγℓ_P²)^{3/2}·√|Q| for some constant c that different papers (Rovelli–Smolin, Ashtekar–Lewandowski) write
differently. The classical tetrahedron fixes it.

**Classical identity (exact).** For a regular tetrahedron the three unit outward normals have Gram determinant 16/27, so
|n₁·(n₂×n₃)| = 4/(3√3). With face area A, V = (√2/3)·√|A₁·(A₂×A₃)| equals the direct volume a³/(6√2), a² = 4A/√3,
to machine precision. Hence **c = √2/3**, and any other constant fails the classical limit.

**Quantum side (computed).** Q is built in the coupled basis from 6j recouplings (Racah's sum in mpmath, so spins up to
j = 40 are reachable; the earlier exact-radical route stops near j = 3). With
r(j) = max|Q| / [(4/(3√3))(j(j+1))^{3/2}]:

| j | 3/2 | 3 | 8 | 12 | 20 | 30 | 40 |
|---|---|---|---|---|---|---|---|
| r | 0.7937 | 0.8403 | 0.9212 | 0.9444 | 0.9651 | 0.9762 | 0.9819 |

A polynomial fit in 1/j through j = 8…40 gives r(j) = 1 − 3/(4j) + O(1/j²), limit 0.99999. So the largest eigenvalue of Q
reproduces the regular tetrahedron's triple product at large spin, and the prefactor √2/3 is the right one.

**Smallest volume.** Four spin-½ legs (Q = ±√3/4) with γ = ln 2/(π√3): V = 1.777 ℓ_P³.

**Not claimed.** That Q is *the* volume operator: the two standard quantisations differ at finite spin, and the
1/j correction above is specific to this one. Only the common prefactor is fixed. The value of γ used for the numerical
volume is the spin-½ count of [G13](G13-immirzi-counting.md); with the other values of G13 the volume changes by
(γ/γ₀)^{3/2}.
