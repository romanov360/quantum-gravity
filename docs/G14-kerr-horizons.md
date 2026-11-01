# G14. Kerr horizons, the Penrose gap, and the area of spin-½ punctures

Script: `python/qg/kerr.py`. Test: `test_kerr_*`.

Parametrise Kerr by r₊ > r₋ > 0 (G = c = 1): M = (r₊ + r₋)/2, a² = r₊r₋, J = Ma. All of the following are checked
symbolically (sympy) as exact identities.

- A± = 4π(r±² + a²) = 8πM r±; κ± = (r± − r∓)/(2(r±² + a²)); Ω± = a/(r±² + a²).
- Smarr M = 2Ω±J + κ±A±/4π and the first law dM = κ±dA±/8π + Ω±dJ at **both** horizons (the latter as an identity of
  one-forms in (r₊, r₋), with the signed κ₋).
- |κ₋|/κ₊ = r₊/r₋ > 1.
- **Penrose gap:** M² − A₊/16π = M r₋/2. It vanishes only when the inner horizon does (Schwarzschild), so the Penrose
  inequality is strict for every rotating Kerr.
- Christodoulou–Ruffini M² = M_irr² + J²/(4M_irr²); A₊A₋ = (8πJ)²; A₊ − 8π|J| = 8πM(r₊ − a) ≥ 0, zero only at
  extremality.
- Limits: a²|κ₋| → M as a → 0 and κ± → 0 as a → M, so constants in near-Kerr stability and censorship arguments that are
  expressed through κ₋ cannot be uniform toward zero rotation or toward extremality.

**Heuristic composition with G13.** N spin-½ punctures have area 4√3πγN ℓ_P²; requiring A ≥ 8π|J| (ħ = 1) with
γ = ln 2/(π√3) gives N ≥ 2π|J|/ln 2: 4.53 (so five punctures) for J = ½, 9.06 (so ten) for J = 1.

**Not claimed.** The Kerr bound is a statement about a classical solution; a spin-network area is an eigenvalue of an
operator, and nothing here connects the two beyond the shared symbol A. The inner horizon is not a surface of the
spin-network states. No theorem about strong cosmic censorship is tested; only the Kerr background quantities that
such arguments use.
