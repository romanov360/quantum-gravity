"""Fast checks of every numbered idea that has code (run: pytest python/tests)."""
import math

import numpy as np
from fractions import Fraction

import pytest

from qg import bkl_billiard, cdt2d, cq_tradeoff, forest_partners, lorentz_strong, sixj, ym2_brownian
from qg import gravity_group_checks as ggc


# G2 / G8 group constants
def test_group_ranges():
    for name, N, roots, wts, Hs in ggc.groups():
        c = ggc.killing_over_trace(roots, wts, Hs[0])
        bw = 2 * N * (c / 4) / 29
        expect = {"SU(2)": Fraction(4, 29), "SO(4)": Fraction(4, 29), "SO(5)": Fraction(15, 58), "Sp(2)=Spin(5)": Fraction(12, 29)}
        assert bw == expect[name.split()[0]]


# G7 Cartan matrices
def test_cartan():
    assert ggc.det([[2, -2, 0], [-2, 2, -1], [0, -1, 2]]) == -2
    e10 = ggc.cartan_from_edges(10, [(i, i + 1, 1) for i in range(8)] + [(6, 9, 1)])
    assert ggc.det(e10) == -1 and ggc.inertia(e10) == (9, 1)


# G4 forests
def test_forests():
    assert forest_partners.z4_4ym_pattern()[0] == 12
    assert forest_partners.z3_romanov_pattern()[0] == 7
    V, E, T = forest_partners.cell600()
    assert (len(V), len(E), len(T)) == (120, 720, 1200)
    assert forest_partners.triangle_lower_bound(120, 720, 5) == 9
    eidx = {e: i for i, e in enumerate(E)}
    plaqs = [[eidx[(a, b)], eidx[(a, c)], eidx[(b, c)]] for a, b, c in T]
    (m, _), F = forest_partners.greedy_cover_forest(V, E, plaqs, seeds=5)
    assert forest_partners.is_forest(len(V), E, sorted(F)) and m == 9


# G10 YM2 on the Brownian disk
@pytest.mark.parametrize("s", [0.3, 2.0])
def test_disk_laplace(s):
    assert ym2_brownian.laplace_closed(s, 1.1) == pytest.approx(ym2_brownian.laplace_numeric(s, 1.1), rel=1e-9)


def test_disk_wilson_mc():
    m, e = ym2_brownian.wilson_brownian_disk_mc(0.5, 1.0, 1.0, n=400_000, seed=3)
    assert abs(m - ym2_brownian.wilson_brownian_disk(0.5, 1.0, 1.0)) < 5 * e


# G3 2D CDT
def test_cdt():
    assert cdt2d.series_check() < 1e-12
    a = 0.1
    g = 0.5 * math.exp(-a * a)
    w = cdt2d.spectrum(a, lmax=400, k=4)
    s = math.sqrt(1 - 4 * g * g)
    assert w[1] / w[0] == pytest.approx((1 - s) / (1 + s), rel=1e-8)
    assert -math.log(w[1] / w[0]) / a == pytest.approx(2 * math.sqrt(2), rel=2e-3)


# G11 trade-off
def test_cq_product():
    _, varp = cq_tradeoff.simulate(0.2, F=1.0, T=1.0, dt=2e-4, n_traj=6000, seed=5, z0=1.0)
    assert 2 * 0.2 * varp / 2 == pytest.approx(0.25, rel=0.08)


# G12 BKL
def test_bkl():
    p = bkl_billiard.kasner(Fraction(5, 2))
    assert sum(p) == 1 and sum(q * q for q in p) == 1
    import mpmath
    with mpmath.workdps(300):
        u0 = mpmath.e + 2
        assert bkl_billiard.era_lengths(u0, 60, dps=300) == bkl_billiard.cf_digits(u0, 60, dps=300)


# G9 6j
def test_sixj():
    assert sixj.orthogonality(1, 1, 1, 1)[0] == 0
    assert sixj.biedenharn_elliott(2, 1, 1, 2, 1, 2, 1, 2, 2)[0] == 0
    assert sixj.racah_sum_rule(1, 1, 1, 1, 1, 2) == 0
    js = (48, 48, 48, 48, 48, 48)
    exact = float(sixj.wigner_6j(js[0], js[1], js[3], js[5], js[4], js[2]))
    assert exact == pytest.approx(sixj.ponzano_regge(js), rel=2e-2)


# G8 character coefficients
def test_strong_coupling():
    assert lorentz_strong.su2_aF_numeric(0.2) == pytest.approx(lorentz_strong.su2_aF(0.2), rel=1e-9)
    assert lorentz_strong.so4_aF(0.05) == pytest.approx(0.05 / 16, rel=1e-4)


# G13 Immirzi
def test_immirzi():
    from qg import immirzi
    lo, hi = immirzi.solve("two", digits=20)
    assert abs(float(lo) - 0.23753295796592066) < 1e-15
    lo, hi = immirzi.solve("dim", digits=20)
    assert abs(float(lo) - 0.27406685823282999) < 1e-15


# G17 volume
def test_volume():
    from qg import volume
    import numpy as np
    for js in ((0.5,) * 4, (1.5,) * 4, (1, 1, 0.5, 0.5)):
        assert np.max(np.abs(volume.brute_force_spectrum(js) - np.sort(volume.spectrum(*js)))) < 1e-12
    assert volume.ratio_regular(8) == pytest.approx(0.921229, abs=1e-6)
    assert 4096 * (1 - volume.ratio_regular(4096)) == pytest.approx(0.75, abs=1e-3)


# G18 clusters
def test_clusters():
    from qg import cluster_gap as cg
    g = lambda cells: cg.gap(cg.planar_cluster(cells))[0]
    assert g({(0, 0)}) == pytest.approx(math.sqrt(8))
    assert g({(0, 0), (1, 0)}) == pytest.approx(math.sqrt(6))
    assert g({(0, 0), (1, 0), (2, 0)}) == pytest.approx(g({(0, 0), (1, 0), (0, 1)}))
    assert g({(0, 0), (1, 0), (0, 1), (1, 1)}) == pytest.approx(2.0)
    N = 7
    assert g({(i, 0) for i in range(N)}) == pytest.approx(2 * math.sqrt(2 - math.cos(math.pi / (N + 1))))
    mu = cg.gap(cg.cube_surface())[1]
    assert abs(mu.min()) < 1e-9 and sorted(np.round(mu, 9))[1] == pytest.approx(4)


# G9 quantum 6j
def test_quantum6j():
    from qg import quantum6j as q6
    assert q6.classical_limit_check() < 1e-5
    from fractions import Fraction as Fr
    for k in (3, 5):
        Q = q6.SU2k(k)
        assert Q.orthogonality(1, 1, 1, 1)[0] < 1e-12
        assert Q.orthogonality(Fr(1, 2), 1, Fr(1, 2), 1)[0] < 1e-12
        assert abs(Q.biedenharn_elliott(1, 1, 1, 1, 1, 1, 1, 1, 1)[0]) < 1e-12
        assert Q.total_dimension_sq() == pytest.approx((k + 2) / (2 * math.sin(math.pi / (k + 2)) ** 2))
