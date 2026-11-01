import Mathlib

/-!
# G18. The harmonic modes of a strip of `N` plaquettes

At weak coupling the plaquette fluxes of a planar SU(2) cluster are oscillators with kinetic matrix `K = 4 I − A`,
where `A` is the adjacency matrix of the graph of plaquettes sharing a link (`python/qg/cluster_gap.py`). For a
strip of `N` plaquettes that graph is the path `P_N`, and `K` is tridiagonal with `4` on the diagonal and `−1`
next to it. This file proves:

* `strip_mode`: `v_m = sin(m θ)` satisfies `4 v_{m+1} − v_m − v_{m+2} = (4 − 2 cos θ) v_{m+1}` for every `m`;
* `strip_mode_boundary`: for `θ_k = k π/(N+1)`, `v_0 = 0` and `v_{N+1} = 0`, so `(v_1, …, v_N)` is an eigenvector of
  `K` with eigenvalue `4 − 2 cos θ_k` (`k = 1, …, N`);
* `strip_lowest`: `k = 1` gives the smallest of these eigenvalues;
* `strip_gap_formula`: the harmonic gap `2 √(μ/2)` at `μ = 4 − 2 cos (π/(N+1))` is `2 √(2 − cos (π/(N+1)))`.

This is Hartshorn's strip formula for general `N` (his "prediction, tested at three values"), as a statement
about the harmonic limit. The `N` eigenvalues are distinct, so they are the whole spectrum of `K`.
-/

open Real

namespace QuantumGravity.Cluster

theorem strip_mode (θ : ℝ) (m : ℕ) :
    4 * sin ((m + 1) * θ) - sin (m * θ) - sin ((m + 2) * θ) = (4 - 2 * cos θ) * sin ((m + 1) * θ) := by
  have h1 : sin ((m + 2) * θ) = sin ((m + 1) * θ + θ) := by ring_nf
  have h2 : sin (m * θ) = sin ((m + 1) * θ - θ) := by ring_nf
  rw [h1, h2, sin_add, sin_sub]
  ring

theorem strip_mode_boundary (N k : ℕ) :
    sin ((0 : ℕ) * (k * π / (N + 1))) = 0 ∧ sin (((N + 1 : ℕ) : ℝ) * (k * π / (N + 1))) = 0 := by
  refine ⟨by simp, ?_⟩
  have hN : ((N : ℝ) + 1) ≠ 0 := by positivity
  have : ((N + 1 : ℕ) : ℝ) * (k * π / (N + 1)) = (k : ℝ) * π := by
    push_cast
    field_simp
  rw [this]
  exact Real.sin_nat_mul_pi k

/-- The lowest mode is `k = 1`. -/
theorem strip_lowest (N k : ℕ) (hk1 : 1 ≤ k) (hkN : k ≤ N) :
    4 - 2 * cos (π / (N + 1)) ≤ 4 - 2 * cos (k * π / (N + 1)) := by
  have hpos : (0 : ℝ) < N + 1 := by positivity
  have hk : (1 : ℝ) ≤ k := by exact_mod_cast hk1
  have hkN' : (k : ℝ) ≤ N := by exact_mod_cast hkN
  have h0 : 0 ≤ π / (N + 1) := by positivity
  have hle : π / (N + 1) ≤ k * π / (N + 1) := by
    rw [div_le_div_iff_of_pos_right hpos]
    nlinarith [Real.pi_pos]
  have hpi : k * π / (N + 1) ≤ π := by
    rw [div_le_iff₀ hpos]
    nlinarith [Real.pi_pos]
  have := Real.cos_le_cos_of_nonneg_of_le_pi h0 hpi hle
  linarith

theorem strip_gap_formula (N : ℕ) :
    2 * Real.sqrt ((4 - 2 * cos (π / (N + 1))) / 2) = 2 * Real.sqrt (2 - cos (π / (N + 1))) := by
  congr 2
  ring

end QuantumGravity.Cluster
