# G15. Planar maps, FK weights and the SU(2)_k loop weight

Script: `python/qg/planar_maps.py`. Tests: `test_tutte_mullin`, `test_fk_table`, `test_su2k_walks`.

This is the nearest the Catalan counts of spin networks come to Liouville quantum gravity, and it connects to
[G10](../PROGRESS.md) (YM₂ on the Brownian disk) through the same KPZ exponents.

1. **Counts.** Brute force over vertex rotations of 2n half-edges (genus 0, connected) gives the rooted planar maps
   2, 9, 54 for n = 1, 2, 3, equal to Tutte's 2·3ⁿ(2n)!/(n!(n+2)!); weighting each map by its spanning-tree count
   (Kirchhoff) gives 2, 10, 70 = Cat_n·Cat_{n+1} (Mullin). Mullin's number also equals the shuffle sum
   Σ_k C(2n,2k)Cat_k Cat_{n−k} (checked to n = 29): a tree-rooted map is a pair of Dyck words, each counted by the
   Catalan numbers that count the invariants of spin-½ legs.
2. **Dictionary.** q = 2 + 2cos(πγ²/2), κ′ = 16/γ², √q = −2cos(4π/κ′), c = (3κ′−8)(6−κ′)/(2κ′) = 25 − 6Q², and
   γ_str = 1 − 4/γ² (KPZ). At q = 0, 1, 2, 3, 4: c = −2, 0, ½, 4/5, 1 and γ_str = −1, −½, −⅓, −⅕, 0 (exact). Rooted
   counts grow like μⁿ n^{γ_str−2}, giving the exponents −3 (Mullin, q = 0) and −5/2 (Tutte, pure gravity).
3. **The SU(2)_k link.** The invariants of N spin-½ legs at level k of U_q(sl₂) are counted by closed walks of length N
   on the Dynkin diagram A_{k+1}. For N ≤ 2k this equals Cat_{N/2}, the undeformed SU(2) count; beyond that it is
   smaller, and it grows as (2cos(π/(k+2)))^N = (√q_k)^N with q_k = 4cos²(π/(k+2)) = 1, 2, 3, 4 at k = 1, 2, 4, ∞
   (percolation, Ising, three-state Potts, the critical end). The loop weight of the spin-½ line is the FK loop weight √q.

**Not claimed.** No map is constructed from a spin network to a decorated planar map: the numbers agree because the
same Temperley–Lieb algebra with loop weight √q acts on both. Nothing here tests a theorem about scaling limits
(convergence to Liouville quantum surfaces). The growth rate ½ln q_k is the exponential rate of a count, not an
Immirzi-type constant.
