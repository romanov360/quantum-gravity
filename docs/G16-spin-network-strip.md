# G16. The SU(2) strip of N plaquettes in the spin-network basis

Script: `python/qg/spin_network_strip.py`. Test: `test_strip_levels`.

The Kogut–Susskind Hamiltonian H = (g²/2)ΣC(j) + (1/g²)Σ_p(2 − W_p) on an open row of N plaquettes, in the basis of
spin networks on a ladder of squares: 3N − 3 labels, electric energy a sum of squared areas, W_p = χ_{1/2}(U_p) acting by
±½ shifts with one signed 6j symbol per vertex on the loop. This is the Hamiltonian counterpart of the 3D
Ponzano–Regge state sum of [G9](../PROGRESS.md): the same 6j symbols, but the weights are a gauge-theory Hamiltonian,
not BF theory.

Lowest levels at g² = 2 (cut chosen for convergence):

| N | ground state | next three |
|---|---|---|
| 3 | 2.75396 | 5.83491, 5.85682, 5.86795 |
| 4 | 3.67192 | 6.74827, 6.76384, 6.77939 |
| 5 | 4.58991 | 7.66634, 7.67567, 7.68949 |

**Status of the phases.** The sign rule of W_p was fitted at N = 3 against an independent Clebsch–Gordan construction. The
extension to N ≥ 4, a factor (−1)^((t−b)/2) for each interior neighbour, is a conjecture. It is tested by the independent
Peter–Weyl calculation (N loop holonomies, no 6j symbols, no phases; in `brentharts/openai_math_yang_mills`,
`sympy_models/models/su2_peter_weyl.py`), which agrees to 2·10⁻⁶ at N = 4 and 10⁻³ at N = 5 on the lowest four levels,
and the plaquette operators commute. The independent engine reproduces the two-plaquette ladder exactly.

**Weak coupling: consistent, not confirmed.** The harmonic limit predicts a gap 2√(2 − cos(π/(N+1))): 2.2741, 2.1827, 2.1297 for
N = 3, 4, 5. For N = 3 this was extrapolated directly (2.2749, from g² = 0.2–0.4). For N = 4 the gap at g² ≤ 0.5 is still
falling with the energy cut (g² = 0.4: 2.493, 2.298, 2.198 at cuts 32, 40, 48, ratio of successive drops ≈ 0.5), so a direct
extrapolation is not possible; N = 5 was too slow. What can be said at g² = 0.6–0.8, where the cut has converged:

| g² | N = 3 gap | N = 4 gap | ratio N = 4 / N = 3 |
|---|---|---|---|
| 0.8 | 2.06059 | ≈ 2.0295 (cuts 32, 40: 2.0325, 2.0302; asymptote estimated) | 0.9849 |
| 0.6 | 2.12735 | ≈ 2.0833 (cuts 32, 40, 48: 2.1282, 2.0914, 2.0848; asymptote estimated) | 0.9793 |

A straight line through the two ratios extrapolates to 0.9624 at g² = 0; the harmonic formula gives 0.9598. That is a 0.3%
agreement from two points with estimated asymptotes, so it is **consistent with** the formula, not a test of it. A real test
needs a faster state generator (the Python enumeration takes 140 s for 13,500 states at N = 4) and cuts of order 100.

**Relation to the guiding observation.** A finite-N lattice gap says nothing about a gap of a field theory
([G1](../PROGRESS.md)). The formula's N → ∞ limit is 2, so in the harmonic regime the gap would stay finite in lattice
units as the strip grows; this is a statement about the formula, not about the numerics.
