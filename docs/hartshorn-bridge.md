# Bridge to Hartshorn, *Spin Networks on the SU(2) Ladder* (openai_math_yang_mills, QUANTUM_GRAVITY v3)

Brent's paper reads the exactly solved SU(2) Hamiltonian clusters (one plaquette, ladder, strip, L-shape) as
loop-quantum-gravity kinematics, and builds the knowledge graph around it. This note records:
- what this repository adds to his open items;
- where his paper and ours connect;
- one small correction to a remark in his paper.

Code and Lean live here; the numbers are written so that they can be dropped into his `models/` as new rows.

## His open items, answered

| His open item (§12, QUANTUM_GRAVITY.md) | Result here | Where |
|---|---|---|
| "The all-spin Immirzi value (about 0.2375) is quoted from Domagala–Lewandowski, not rederived" | Rederived and certified by bisection with an explicit tail bound: **γ_DL = 0.23753 29579 65920 66435 8377** (all spins, 2 states per puncture), and the Ghosh–Mitra counting (2j+1 states) **γ_GM = 0.27406 68582 32829 99621 6763**, beside his ln 2/(π√3) = 0.12738 40231 4094 | `python/qg/immirzi.py` (G13) |
| "a volume prefactor for the triple product" | V = (√2/3)(8πγℓ_P²)^{3/2}√\|Q\|, from the classical tetrahedron V² = (2/9)\|A₁·(A₂×A₃)\| with A_i = 8πγℓ_P²J_i. The smallest volume is (√2/3)(√3/4)^{1/2}(8πγ)^{3/2}ℓ_P³ = 0.3102 (8πγ)^{3/2}ℓ_P³ | `python/qg/volume.py` (G17) |
| The largest-volume ratio "rises from 0.794 at j = 3/2 to 0.921 at j = 8 … a trend, not a proof of a semiclassical limit" | Closed-form tridiagonal matrix elements a_k = √(Δ(j₁,j₂,k)Δ(j₃,j₄,k))/(4√(4k²−1)), matched against brute force. The ratio is taken to j = 65536: it tends to 1, with **j(1 − ratio) → 3/4** (0.74894, 0.74974, 0.74993, 0.74998 at j = 1024, 4096, 16384, 65536). Conjecture: λ_max = (4/3√3)[j(j+1)]^{3/2}(1 − 3/(4j) + O(j⁻²)). All his printed ratios are reproduced | `python/qg/volume.py` (G17) |
| "the strip formula for general N is a prediction, tested at three values" | **Proved in the harmonic limit for every N**: sine eigenvectors of the tridiagonal kinetic matrix, the lowest mode is k = 1, gap 2√(2 − cos(π/(N+1))) | `QuantumGravity/ClusterModes.lean` (G18) |
| The L-shape gap equals the strip's "(as for the strip)" | **Explained.** The harmonic kinetic matrix is K = ∂∂ᵀ on plaquettes; for planar clusters K = 4I − A, with A the adjacency of the plaquette graph. The L-shape and the 3-strip both have plaquette graph P₃, so their gaps agree; the corner vertex never enters | `python/qg/cluster_gap.py` (G18) |

## New predictions for his code (harmonic limit, g → 0)

- 2×2 block: Δ = **2** exactly (K spectrum 2, 4, 4, 6). Its centre is a four-valent vertex, like the L corner.
- T-shape (four plaquettes): Δ = 2√(2 − √3/2) = **2.12977**.
- Closed cube surface (six faces): K has one zero mode, the Bianchi identity (4Y-M `cube_bianchi`); the gap of the
  remaining modes is 2√2.
- An L×L block: Δ = 2√(2 − 2cos(π/(L+1))) → 0, the massless gluons of weak coupling.

## Where the two papers meet

1. **The self-dual SU(2).** His §8: the SU(2) of a spin-network link is SU(2)₊ ⊂ Spin(4) = SU(2)₊ × SU(2)₋. Our
   `GravityGroups.lean`: the forest-gauge range for SU(2)₊ is β_W < 4/29, the same as for SO(4) in its vector
   representation. The strong-coupling exclusion of G1/G8 therefore applies to a Wilson-type lattice action for the
   self-dual (Ashtekar) connection as well. This is not the Plebanski action, which also carries simplicity
   constraints.
2. **The area quantum in a Wilson loop.** G10 (`ym2_brownian.py`): on a Brownian disk the YM₂ Wilson loop is
   d_F(1 + gL√C_F)e^{−gL√C_F}. Its decay rate is g√(j(j+1)), the LQG area eigenvalue √(j(j+1)) of his §3. This is
   an arrow, not an equals sign: the √C comes from averaging an area law over the Brownian area distribution.
3. **Planar maps and Catalan counts.** His §(fifth round): Catalan counts of spin-½ intertwiners and Tutte counts
   of planar maps. Our G10 uses the continuum limit of those maps (the Brownian disk), so both papers use the same
   random surface at different levels.
4. **Dimension four, twice.** His observation that (d−1)(d−4) = 0 decides both the existence of the forest with one
   link per plaquette and the splitting Λ² = Λ⁺ ⊕ Λ⁻ is now a Lean lemma (`GravityGroups.dim_four_twice`, G19).
5. **Gauge fixing.** The ladder's spanning-tree gauge and the 4Y-M forest gauge remove the same redundancy. In our
   G2 count (`MMPartners.lean`), the 4Y-M forest leaves 30 partners per link for the MacDowell–Mansouri plaquette
   pairs, against 12 for the Wilson action.

## One correction

His §10 says that −½ log 0.9653 = 0.017658 "is printed as 0.0177" and that the printed value "does not follow".
0.017658 rounded to three significant figures is 0.0177, so the printed number is a rounding of the stated
quantity, not an error. (His QUANTUM_GRAVITY.md phrases it as "0.00004 above", which is accurate.)

## After PR #1

Brent's PR added G13–G17 (`docs/G13`–`G17`, `kerr.py`, `planar_maps.py`, `spin_network_strip.py`,
`volume_tetrahedron.py`). Two topics overlapped and are merged:
- **G13, Immirzi:** both implementations are kept in `immirzi.py` and agree. His uses Newton and a direct renewal
  count; ours uses bisection with a certified tail bound.
- **G17, volume:** his `volume_tetrahedron.py` (prefactor, to j = 40) and our `volume.py` (closed-form matrix
  elements, to j = 65536) find the same 1 − 3/(4j).

His G16 tests the strip formula at N = 4 with the full Hamiltonian. Our G18 proves it in the harmonic limit for
every N.
