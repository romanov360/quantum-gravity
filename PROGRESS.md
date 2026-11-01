# Progress on G0–G19

Ideas are as numbered in [docs/03-ideas.md](docs/03-ideas.md). Status as of 2026-11-01. G13–G19 come from, or answer, Brent Hartshorn's quantum gravity paper (`docs/hartshorn-bridge.md`); G13–G17 were numbered in his PR (#1) and are kept.

| Idea | What | Status | Where |
|---|---|---|---|
| G0 | Four open problems of lattice QG | **written** | `OPEN.md` |
| G1 | A physical gap excludes a massless mode | **proved in Lean** | `QuantumGravity/NoGraviton.lean` |
| G2 | Pregeometric phase of MacDowell–Mansouri lattice gravity | group constants and **pair-action partner count 30 (42 unfixed) proved in Lean**; Lipschitz constant of a pair term and the S⁴ compensator **to do** | `GravityGroups.lean`, `MMPartners.lean`, `docs/G2-pregeometric-phase.md` |
| G3 | CDT through the transfer matrix | 2D: **spectrum proved**, λ_n = x₋^{2n+2}/g (Möbius identity and fixed-point algebra in Lean; triangular/completeness argument in the doc); 3+1 programme written | `CDT2D.lean`, `python/qg/cdt2d.py`, `docs/G3-cdt.md` |
| G4 | YM on any geometry; induced gravity | count lemma, ℤ³/ℤ⁴ reproduced, **600-cell optimum 9 proved**, general lower bound; induced-gravity cluster expansion **to do** | `ForestRowSum.lean`, `python/qg/forest_partners.py`, `docs/G4-geometry-uniform.md` |
| G5 | Stress tensor of continuum YM | programme written (4 steps), all open | `docs/G5-stress-tensor.md` |
| G6 | Gravity and asymptotic freedom | **proved in Lean** | `QuantumGravity/RunningWithGravity.lean` |
| G7 | Cosmological billiards | AE₃/E₁₀ **in Lean**; BKL = continued fractions **checked** | `Billiards.lean`, `python/qg/bkl_billiard.py` |
| G8 | Confinement of the curvature-squared sector | corollary of 4Y-M + G1 + group constants; character coefficients computed; rigorous area-law constants **to do** | `paper/`, `python/qg/lorentz_strong.py` |
| G9 | 3D: Ponzano–Regge ↔ forest gauge | exact 6j identities and PR asymptotics **checked**; **quantum 6j of SU(2)_k** (orthogonality, Biedenharn–Elliott, total dimension, q → 1 limit) **checked**; Turaev–Viro invariants and the interpolating family **to do** | `python/qg/sixj.py`, `python/qg/quantum6j.py` |
| G10 | YM₂ on Liouville quantum gravity | Brownian-disk Wilson loop **derived and checked**; sphere sectors shift μ | `python/qg/ym2_brownian.py`, `paper/` |
| G11 | Decoherence–diffusion trade-off for QGL | saturating model **built and checked**; port into QGL **to do** (Brent's repo) | `python/qg/cq_tradeoff.py` |
| G12 | Small QGL tests | Gauss–Kuzmin for BKL eras **checked** | `python/qg/bkl_billiard.py` |
| G13 | Immirzi parameter from counting punctures | three conventions, **two independent implementations agree**: 0.12738… (spin ½), **0.2375329579659206643583771** (Domagala–Lewandowski, 2 states), **0.2740668582328299962167637** (Ghosh–Mitra, 2j+1), certified by bisection with a tail bound; the direct renewal count approaches the root from below | `python/qg/immirzi.py`, `docs/G13-immirzi-counting.md` |
| G14 | Kerr horizons and the Penrose gap (Hartshorn) | identities (Smarr, first law at both horizons, Penrose gap M r₋/2, A₊A₋ = (8πJ)²) **checked exactly**; the heuristic composition with spin-½ area is labelled as such | `python/qg/kerr.py`, `docs/G14-kerr-horizons.md` |
| G15 | Planar maps, FK weights, SU(2)_k loop weight (Hartshorn) | Tutte and Mullin counts by brute force, FK dictionary, closed walks on A_{k+1} **checked**; no map from spin networks to maps constructed | `python/qg/planar_maps.py`, `docs/G15-planar-maps-su2k.md` |
| G16 | SU(2) strip of N = 3, 4, 5 plaquettes in the spin-network basis (Hartshorn) | N = 4, 5 levels **agree** with an independent Peter–Weyl calculation (2e-6, 1e-3) under a fitted phase rule; N = 4 weak-coupling gap ratio to N = 3 **consistent** with the formula (0.9624 vs 0.9598, two-point extrapolation); N = 5 not done. The harmonic limit of the formula is proved for every N (G18) | `python/qg/spin_network_strip.py`, `docs/G16-spin-network-strip.md` |
| G17 | Volume operator: prefactor and semiclassical limit | prefactor c = √2/3 **fixed** by the regular tetrahedron (exact); smallest volume 1.777 ℓ_P³ (Hartshorn, to j = 40). Closed-form tridiagonal matrix elements **checked against brute force**; j(1 − ratio) → **3/4** to j = 65536 (Romanov); proof to do | `python/qg/volume_tetrahedron.py`, `python/qg/volume.py`, `docs/G17-volume-prefactor.md` |
| G18 | Harmonic gaps of plaquette clusters | **one formula** Δ = √(2 μ_min(∂∂ᵀ)); strip formula for all N **proved in Lean**; L-shape = strip explained; predictions 2×2 → 2, T → 2.12977 | `ClusterModes.lean`, `python/qg/cluster_gap.py` |
| G19 | Dimension four, twice (Hartshorn's observation) | **Lean** | `GravityGroups.dim_four_twice` |

## Next steps, in order

0. G16: a faster state generator to reach g² < 0.4 at N = 4, 5 and test 2√(2 − cos(π/(N+1))) beyond the harmonic
   limit (G18); SU(2)_k deformation of the strip with quantum 6j (the link to G9).
1. G8: rigorous area law for SO(4) vector loops at β_W < 4/29, with explicit constants (polymer expansion), to add to
   the paper.
2. G17: prove j(1 − ratio) → 3/4 (Bohr–Sommerfeld for the tridiagonal a_k, or a WKB analysis of the Jacobi matrix).
3. G18: run the 2×2 and T-shape predictions through Brent's full-Hamiltonian code (`spin_network_strip.py` generalised to planar clusters).
4. G2: the Lipschitz constant of the pair term and Hessian bound of the MacDowell–Mansouri plaquette-pair action on the 4Y-M forest, in Lean
   like `GaugeFix.lean`.
5. G4: the induced-gravity expansion on degree-capped complexes.
6. G11: open a PR to `brentharts/Quantum_Gravity_Lab` with a `qgl_cq.py` module, if Brent wants it.
