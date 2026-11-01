"""Fast checks of G13-G16 (run: pytest python/tests)."""
import math

import pytest
from mpmath import mpf, pi, log, sqrt

from qg import immirzi, kerr, planar_maps, spin_network_strip


# G13 Immirzi parameter from three counting conventions
def test_immirzi_values():
    assert float(immirzi.gamma_spin_half()) == pytest.approx(math.log(2) / (math.pi * math.sqrt(3)), rel=1e-12)
    assert float(immirzi.gamma_domagala_lewandowski()) == pytest.approx(0.2375330, abs=1e-6)
    g = immirzi.gamma_all_spins()
    assert float(g) == pytest.approx(0.2740669, abs=1e-6)
    assert float(immirzi.renewal_sum(2 * pi * g)) == pytest.approx(1.0, abs=1e-12)
    assert float(immirzi.gamma_spin_half()) < float(immirzi.gamma_domagala_lewandowski()) < float(g)


def test_immirzi_growth_rate_approaches_root():
    r, _ = immirzi.growth_rate(20.0, 0.01)
    lam = float(2 * pi * immirzi.gamma_all_spins())
    assert 0.9 * lam < r < lam          # ln N(x)/x increases to lambda from below


# G14 Kerr
def test_kerr_identities():
    for d in (kerr.horizon_identities(), kerr.penrose_identities()):
        assert all(d.values()), d
    lim = kerr.limits()
    assert all(v in (0, 1) for v in lim.values())


def test_kerr_punctures():
    assert kerr.punctures_needed(0.5) == pytest.approx(math.pi / math.log(2), rel=1e-12)
    assert kerr.punctures_needed(1.0) == pytest.approx(2 * math.pi / math.log(2), rel=1e-12)


# G15 planar maps
def test_tutte_mullin():
    for n in (1, 2, 3):
        m, t = planar_maps.rooted_counts(n)
        assert m == planar_maps.tutte(n) and t == planar_maps.catalan(n) * planar_maps.catalan(n + 1)
    assert all(planar_maps.mullin_shuffles(n) == planar_maps.catalan(n) * planar_maps.catalan(n + 1) for n in range(30))


def test_fk_table():
    import sympy as sp
    rows = [planar_maps.fk_row(q) for q in range(5)]
    assert [r[2] for r in rows] == [-2, 0, sp.Rational(1, 2), sp.Rational(4, 5), 1]
    assert [r[3] for r in rows] == [-1, sp.Rational(-1, 2), sp.Rational(-1, 3), sp.Rational(-1, 5), 0]


def test_su2k_walks():
    assert all(planar_maps.closed_walks(6, 2 * m) == planar_maps.catalan(m) for m in range(7))
    assert planar_maps.closed_walks(2, 6) == 4 < planar_maps.catalan(3)
    for k, q in ((1, 1), (2, 2), (4, 3)):
        assert planar_maps.q_of_level(k) == pytest.approx(q, abs=1e-12)
        assert planar_maps.closed_walks(k, 122) / planar_maps.closed_walks(k, 120) == pytest.approx(q, rel=1e-6)


# G16 strip, spin-network basis; reference values from the independent Peter-Weyl construction
def test_strip_levels():
    w3, _ = spin_network_strip.lowest(3, 20, 2.0, k=2)
    assert w3[0] == pytest.approx(2.753963, abs=2e-5)
    w4, _ = spin_network_strip.lowest(4, 16, 2.0, k=2)
    assert w4[0] == pytest.approx(3.67192, abs=2e-5)


# G17 volume prefactor
def test_volume_classical():
    from qg import volume_tetrahedron as vt
    assert vt.regular_tetrahedron_triple_product() == pytest.approx(4 / (3 * math.sqrt(3)), rel=1e-12)
    lhs, rhs = vt.volume_prefactor_identity()
    assert lhs == pytest.approx(rhs, rel=1e-12)
    assert vt.smallest_volume() == pytest.approx(1.7769, abs=2e-3)


def test_volume_semiclassical():
    from qg import volume_tetrahedron as vt
    s = vt.q_spectrum(1)                     # four spin-1/2: Q = +-sqrt3/4
    assert max(abs(s)) == pytest.approx(math.sqrt(3) / 4, rel=1e-9)
    assert vt.ratio(3) == pytest.approx(0.79373, abs=2e-5)       # four spin-3/2, as in the exact radicals
    r = [vt.ratio(k) for k in (8, 16, 24, 32)]
    assert all(a < b < 1 for a, b in zip(r, r[1:]))
    j = 16.0
    assert r[-1] == pytest.approx(1 - 3 / (4 * j), abs=0.005)
