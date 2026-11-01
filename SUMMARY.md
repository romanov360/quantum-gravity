# Summary: who did what

Status 2026-11-01. Details are in [PROGRESS.md](PROGRESS.md); the draft paper is [paper/main.pdf](paper/main.pdf)
(*A gap is not a graviton*, Romanov & Hartshorn).

**The main idea.** In Yang–Mills a mass gap is the goal; in gravity it rules a phase out. The graviton is massless, so
the strong-coupling gap theorems of 4Y-M, applied to lattice gravity written as an SO(4) or SO(5) gauge theory, prove
that those phases are not gravitational.

## Tristan Romanov

**Proved in Lean** (no `sorry`, standard axioms only):
- **G1:** a lattice theory whose gap stays positive in physical units has no continuum limit with a massless mode,
  so no graviton.
- **G2/G8:** the group constants that make the 4Y-M theorems apply to lattice gravity: SO(4) β_W < 4/29, SO(5)
  15/58, Sp(2) 12/29. In the SO(4) range there is a unique state, a gap, and no ghost states. Under the
  MacDowell–Mansouri pair action each link has 30 partners, against 12 for the Wilson action.
- **G3:** the 2D causal-triangulation transfer matrix has an exactly geometric spectrum,
  λₙ = x₋^{2n+2}/g, so its lattice gap is of order the spacing. The key identity and the algebra are in Lean; the
  last two steps of the proof are in `docs/G3-cdt.md`.
- **G4:** on the 600-cell, 9 partners per link is optimal for any gauge-fixing forest (β_W < 4/9).
- **G6:** if gravity antiscreens, Yang–Mills stays asymptotically free; if it screens, the coupling goes to a fixed
  point instead.
- **G7:** the AE₃ and E₁₀ billiard algebras of near-singularity gravity have Lorentzian signature.
- **G18:** one formula, Δ = √(2μ_min(∂∂ᵀ)), gives the weak-coupling gap of any plaquette cluster. It proves the
  strip formula for every N, explains why the L-shape matches the strip, and predicts 2 for a 2×2 block.
- **G19:** Brent's "dimension four, twice" observation, as a lemma.

**Computed and checked:**
- **G10:** Yang–Mills on a Brownian disk turns the area law into d_F(1 + gL√C)e^{−gL√C}.
- **G13:** the Immirzi parameter, certified to 25 digits: 0.23753… (2 states per puncture), 0.27407… (2j+1 states).
- **G17:** exact closed-form matrix elements of the volume operator. Its largest eigenvalue approaches the classical
  value as 1 − 3/(4j); this rate is numerical (to j = 65536), not proved.
- **G9:** ordinary and quantum 6j symbols pass their standard identities.
- **G11:** a toy classical–quantum model shows Oppenheim's decoherence–diffusion trade-off.
- **G12:** the BKL era lengths are continued-fraction digits, with Gauss–Kuzmin statistics.

**Written up:** a survey of quantum gravity as of October 2026, four open problems of lattice quantum gravity
(G0, `OPEN.md`), and plans for the remaining items. These include G5, the Yang–Mills stress tensor, which is how
Yang–Mills couples to gravity.

## Brent Hartshorn

**His paper**, *Spin Networks on the SU(2) Ladder* (github.com/brentharts/openai_math_yang_mills), reads the exactly
solved small SU(2) lattice systems as loop-quantum-gravity spin networks. It covers the area spectrum, the Immirzi
parameter from spin-½ counting (0.1274), the quantum-tetrahedron volume spectra, Einstein 4-manifolds and
Penrose-inequality checks, all built into his knowledge graph.

**His PR (#1, merged):**
- **G13:** an independent Immirzi solver and a direct count of puncture states.
- **G14:** exact Kerr-horizon identities, including the Penrose gap M r₋/2.
- **G15:** planar-map counts by brute force and the SU(2)_k loop weights.
- **G16:** the N = 3, 4, 5 plaquette strips with the full Hamiltonian, cross-checked against an independent
  Peter–Weyl calculation.
- **G17:** the volume prefactor √2/3 and a smallest volume of 1.777 ℓ_P³. He found the same 1 − 3/(4j) rate
  independently (to j = 40).

## Where they meet

Tristan's work answers several of Brent's open items:
- the all-spin Immirzi value he had only quoted;
- the volume prefactor and its large-spin limit;
- the strip formula, proved for every N.

Brent's full-Hamiltonian strips test that formula beyond the weak-coupling limit, the only regime the proof covers.
See [docs/hartshorn-bridge.md](docs/hartshorn-bridge.md).
