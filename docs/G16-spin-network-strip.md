# G16. The SU(2) strip of N plaquettes in the spin-network basis

Script: `python/qg/spin_network_strip.py` (gap scans: `python/qg/strip_gap_scan.py`, data: `docs/data/g16_gap_scan.jsonl`). Tests: `test_strip_*` in `python/tests/test_g13_g16.py`.

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

**Faster code.** The 140 s for 13,500 states at N = 4 was the state enumeration: it generated every label tuple under the
energy bound and checked the triangle conditions at the end (1.1 million recursive calls for 909 states at N = 4, cut 24). `states()`
now builds the labels one plaquette boundary at a time with both triangle conditions and the energy bound applied during
construction, in milliseconds (13,527 states in 0.04 s); the old routine is kept as `states_slow()` and the tests check that
the two return the same set of states. The second bottleneck was the eigensolver: `lowest()` used shift-invert at σ = 0,
whose sparse factorisation fills in (13.3 s at N = 4, cut 40), while plain Lanczos for the smallest eigenvalues takes
0.03 s and agrees to 2·10⁻¹⁴ (`method="shift_invert"` keeps the old path). What remains is the Python matrix assembly: 28 s
for 113,051 states (N = 4, cut 80), about 2 minutes for 475,000 states, about 2 GB. Nothing in the physics changed: the levels
in the table above and the G16 gaps at g² = 0.4 (2.4932, 2.2985, 2.1984 at cuts 32, 40, 48) are reproduced exactly.

**Converged gaps.** `python/qg/strip_gap_scan.py` raises the cut until the gap moves by less than 2·10⁻⁵ (or the state
count passes a limit) and Aitken-extrapolates the cut from the last three values. All numbers are in
`docs/data/g16_gap_scan.jsonl`. The gap E₁ − E₀:

| g² | N = 3 | N = 4 | N = 5 | N = 4 / N = 3 | N = 5 / N = 3 |
|---|---|---|---|---|---|
| 0.8 | 2.06059 | 2.03004 | 2.01442 | 0.98517 | 0.97759 |
| 0.6 | 2.12735 | 2.08381 | ≈ 2.0615 | 0.97953 | ≈ 0.969 |
| 0.5 | 2.16058 | 2.10969 | | 0.97645 | |
| 0.4 | 2.18908 | 2.13079 | | 0.97337 | |
| 0.3 | 2.21384 | ≈ 2.1480 | | ≈ 0.9703 | |
| 0.25, 0.2, 0.15 (N = 3 only) | 2.22515, 2.23587, 2.24609 | | | | |
| harmonic limit, g² → 0 | 2.27411 | 2.18264 | 2.12976 | 0.95978 | 0.93653 |

Convergence: every entry converged to better than 10⁻⁴ except N = 4 at g² = 0.3 (cuts 96, 104, 112 give 2.15047, 2.14888, 2.14832;
the Aitken limit 2.1480 is uncertain by a few 10⁻⁴) and N = 5 at g² = 0.6 (stopped at cut 52 with the last step 2·10⁻³; the
limit is uncertain by about 10⁻³). The hand estimates of the earlier version of this note (≈ 2.0295 and ≈ 2.0833 for N = 4 at
g² = 0.8 and 0.6) agree with these to 5·10⁻⁴.

**What this says about the harmonic formula.** The gap falls with g² (anharmonic corrections), so the value at g² = 0.4 is
not the harmonic one and the formula must be tested at g² → 0.

- *N = 3, direct.* Polynomial fits in g² to the eight converged points extrapolate to 2.2740–2.2751 when restricted to
  g² ≤ 0.5 (degree 3 to 5), against 2.27411. Fits that include g² = 0.8 are unstable (spread up to 0.009), so the quoted
  agreement is at the 10⁻³ level, not better.
- *N = 4, direct.* Inconclusive. Five points (g² = 0.3 to 0.8) give extrapolations scattered from 2.15 to 2.22 depending on
  the degree and the range.
- *N = 4 / N = 3.* The ratio is smooth. Fits of degree 1 to 3 over several ranges extrapolate to 0.9602–0.9619, against 0.95978:
  consistent at the 0.1–0.2% level. All the fits lie slightly above the harmonic ratio (by 0.0004 to 0.0021). That could be extrapolation error
  (the nearest point, g² = 0.3, is itself uncertain at the 10⁻³ level in the gap) or a small genuine deviation; the data here do not
  separate the two. It is a consistency check of the formula, not a confirmation, at roughly the same level as the two-point
  extrapolation it replaces (0.9624).
- *N = 5.* Values now exist at g² = 0.8 and 0.6 (they did not before). Two points do not support an extrapolation, so none is attempted.

The phase rule for N ≥ 4 is the same conjectured rule as before, so these runs do not test it independently; that remains
the Peter–Weyl comparison above. That the N = 4 ratio approaches the harmonic value is indirect support, since a wrong rule would
be unlikely to land within 0.2% of an independently derived number, but at that level it is not a proof.

To go lower in g² (0.2 needs cuts of order 130 and about a million states at N = 4) the Python assembly has to be replaced
(numba or C++), or the spin-network basis replaced by another.

**Relation to the guiding observation.** A finite-N lattice gap says nothing about a gap of a field theory
([G1](../PROGRESS.md)). The formula's N → ∞ limit is 2, so in the harmonic regime the gap would stay finite in lattice
units as the strip grows; this is a statement about the formula, not about the numerics.
