# G13. The Immirzi parameter from counting punctures, rederived

Script: `python/qg/immirzi.py`. Test: `python/tests/test_g13_g16.py::test_immirzi_*`.

A puncture of spin j contributes area 8πγ√(j(j+1)) ℓ_P² to the area spectrum of loop quantum gravity. Requiring the
horizon entropy to be S = A/4 fixes γ, but *which* γ depends on what is counted. Three conventions are in the
literature, and all three are computed here, not quoted:

| Counting | Condition | γ |
|---|---|---|
| spin-½ punctures only (the original count of this type) | 2^N states at area 4√3πγN | ln 2/(π√3) = **0.127384** |
| all spins, all 2j+1 projections (Ghosh–Mitra, gr-qc/0411035 eq. 7) | Σ_j (2j+1) e^{−λ√(j(j+1))} = 1, γ = λ/(2π) | λ = 1.722013, γ = **0.274067** |
| all spins, only m = ±j (Domagala–Lewandowski) | Σ_j 2 e^{−λ√(j(j+1))} = 1 | γ = **0.237533** |

The renewal equation behind the second and third rows: the number N(x) of sequences of punctures with total
Σ√(j(j+1)) ≤ x obeys N(x) = 1 + Σ_j d_j N(x − √(j(j+1))), so N(x) ~ C e^{λx} with λ the root of Σ_j d_j e^{−λ a_j} = 1.
`growth_rate` iterates this recursion on a grid and shows ln N(x)/x increasing to λ from below (1.683 at x = 25,
against λ = 1.722). Then S = A/4 requires λ/(8πγ) = 1/4, i.e. γ = λ/(2π).

**Correction to a number used elsewhere.** The value "about 0.2375" is the Domagala–Lewandowski count (degeneracy 2),
not the all-projection count. With all 2j+1 projections the root is 0.2741. Papers and notes that quote 0.2375 should
say which count they mean.

**Not claimed.** These are leading-order counts for sequences of distinguishable punctures. They ignore the
total-spin constraint, the logarithmic correction (−½ ln A), and the debate about which states are physical. The
three values show that γ is a property of the counting, not a prediction of the area spectrum.
