import random

import pytest
from flint import fmpz_poly

from noninterlacing.bruteforce import brute_force_interlacers, integer_box
from noninterlacing.interlace import Interlacing, siegel_ceiling
from noninterlacing.poly import compose, is_irreducible, poly, real_cyclotomic, to_str
from noninterlacing.realroots import interlaces
from noninterlacing.search import totally_real_census

F21 = real_cyclotomic(21)


def test_paper_base_counts():
    I = Interlacing(F21)
    assert I.count(1) == 0
    assert I.count(2) == 141
    assert I.count(3) == 1324
    assert I.min_level() == 2
    S = I.levels()
    assert S.generators == [2, 3] and S.gaps == [1]
    _, box = integer_box(F21, 1)
    assert box == 22680


@pytest.mark.parametrize("level", [1, 2])
def test_engine_matches_bruteforce_on_census(level):
    polys = list(totally_real_census(4, 4.2)) + list(totally_real_census(5, 4.0))[:25]
    checked = 0
    for f in polys:
        _, size = integer_box(f, level)
        if size > 200_000:
            continue
        bf, _ = brute_force_interlacers(f, level)
        en = Interlacing(f).interlacers(level)
        assert [to_str(g) for g in bf] == [to_str(g) for g in en], to_str(f)
        checked += 1
    assert checked >= 15


def test_engine_matches_bruteforce_level2_F21():
    bf, _ = brute_force_interlacers(F21, 2)
    assert bf == Interlacing(F21).interlacers(2)


def test_listed_interlacers_are_exact():
    for f in (real_cyclotomic(28), real_cyclotomic(13), poly("x^4 - 4*x^2 + 2")):
        I = Interlacing(f)
        for j in (1, 2):
            gs = I.interlacers(j)
            assert len(gs) == I.count(j)
            assert all(interlaces(f, g) and int(g.leading_coefficient()) == j for g in gs)


def test_reducible_and_shifted_inputs():
    # reducible separable f: strict interlacing still makes sense
    f = poly("(x^2 - 2)*(x^2 - 3)")
    I = Interlacing(f)
    bf, _ = brute_force_interlacers(f, 1)
    assert I.interlacers(1) == bf
    # translation invariance of the counts
    g = compose(F21, fmpz_poly([5, 1]))
    assert Interlacing(g).count(2) == 141


def test_derivative_is_level_d():
    for f in (F21, real_cyclotomic(11), poly("x^3 - 3*x - 1")):
        I = Interlacing(f)
        assert I.is_interlaced(f.degree())
        assert I.min_level() <= siegel_ceiling(f.degree())


def test_milp_margin_agrees_with_dfs():
    polys = list(totally_real_census(6, 4.0))
    for f in polys:
        I = Interlacing(f)
        mm = I.margin(1)
        dfs = I.find(1, method="dfs") is not None
        if mm.exists is None:
            # only on the boundary (possible for reducible f: g(r) = 0 at a common root)
            assert abs(mm.value) < 1e-6 and not is_irreducible(f), to_str(f)
        else:
            assert mm.exists is dfs, to_str(f)
        assert I.find(1) is not None or not dfs


@pytest.mark.parametrize("n,expected", [(21, 2), (28, 1), (27, 1), (33, 2), (84, 2), (63, 2)])
def test_cyclotomic_levels(n, expected):
    assert Interlacing(real_cyclotomic(n)).min_level() == expected


def test_high_degree_milp():
    I = Interlacing(real_cyclotomic(99))  # degree 30
    assert I.find(1) is None
    g = I.find(2)
    assert g is not None and interlaces(I.f, g)


def test_f_N_family():
    for N in (2, 3, 7):
        f = compose(F21, fmpz_poly([-N, -1, 1]))
        I = Interlacing(f)
        assert I.count(1) == 0
        assert I.find(2) is not None


def test_feasibility_and_certificates():
    I = Interlacing(real_cyclotomic(21))
    assert I.lattice.feasible(1).exists is False
    fr = I.lattice.feasible(2)
    assert fr.exists and interlaces(I.f, I.coords_to_poly(fr.point))
    cert = I.level_certificate()
    assert cert["m"] == 2 and cert["empty_levels"] and interlaces(I.f, poly(cert["witness"]))
    assert Interlacing(real_cyclotomic(28)).level_certificate()["m"] == 1
