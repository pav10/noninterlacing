import random
from fractions import Fraction

import pytest

from noninterlacing.cyclotomic import (chebyshev_composition, is_admissible, odd_degree_conductors,
                                       reachable_degrees)
from noninterlacing.families import check_main_theorem
from noninterlacing.poly import is_irreducible, poly, real_cyclotomic
from noninterlacing.search import totally_real_census
from noninterlacing.siegel import dim_M, siegel_coefficients
from noninterlacing.symmetric import cauchy_interlacer, find_symmetric_realization, random_symmetric, charpoly

from .conftest import needs_pari


def test_siegel_coefficients():
    assert siegel_coefficients(4) == (Fraction(1, 240),)
    assert siegel_coefficients(6) == (Fraction(-1, 504),)
    assert [dim_M(2 * d) for d in (2, 3, 4, 5, 6, 7, 12, 13)] == [1, 1, 1, 1, 2, 1, 3, 2]
    b = siegel_coefficients(12)
    assert b[0] > 0 and b[1] > 0  # contradicts Remark 5.9's "b_1 b_2 < 0 in degree 6"


@needs_pari
def test_siegel_formula_matches_pari():
    from noninterlacing.siegel import siegel_zeta

    for f in (real_cyclotomic(21), poly("x^2 - x - 1"), real_cyclotomic(7)):
        assert siegel_zeta(f)["agree"]


def test_main_theorem_rows():
    r = check_main_theorem(5)
    assert r["totally_real"] and r["irreducible"] and r["span_ok"] and r["disc_formula_ok"]
    assert r["mahler_formula"] == 5**6 + 5**5 - 6 * 5**4 - 6 * 5**3 + 8 * 25 + 40 + 1
    assert r["mahler_formula_ok"] and r["noninterlacing"]
    assert not check_main_theorem(1)["totally_real"]


def test_cyclotomic_helpers():
    assert [n for n in range(21, 60) if is_admissible(n)] == [21, 33, 35, 39, 42, 51, 55, 57]
    assert [k for k in range(2, 12) if is_irreducible(chebyshev_composition(21, k))] == [3, 7, 9]
    D = reachable_degrees(250)
    assert all(x % 2 == 0 for x in D) and not {202, 206, 214, 218} & D and 12 in D and 8 not in D
    assert {27, 49, 81, 121} <= set(odd_degree_conductors(130))


def _brute_census(d, S, coeff_bound):
    import itertools

    from noninterlacing.poly import canonical_form
    from noninterlacing.realroots import RootData, is_totally_real
    from flint import fmpz_poly

    out = set()
    for top in range(-(d - 1), 1):
        for rest in itertools.product(range(-coeff_bound, coeff_bound + 1), repeat=d - 1):
            f = fmpz_poly(list(rest) + [top, 1])
            if is_totally_real(f) and RootData(f).span <= S:
                out.add(tuple(int(c) for c in canonical_form(f).coeffs()))
    return out


def test_census_matches_bruteforce():
    for d, S, B in ((3, 4.0, 8), (4, 3.6, 8)):
        cen = {tuple(int(c) for c in f.coeffs()) for f in totally_real_census(d, S)}
        assert cen == _brute_census(d, S, B)


def test_cauchy_lemma_random():
    rng = random.Random(7)
    n_ok = 0
    for _ in range(80):
        A = random_symmetric(rng.randint(2, 5), 2, rng)
        try:
            r = cauchy_interlacer(A)
        except ValueError:
            continue
        assert r["interlaces"]
        n_ok += 1
    assert n_ok > 30


def test_realization_search():
    M = find_symmetric_realization("x^3 - 3*x - 1")
    assert M is not None and charpoly(M) == poly("x^3 - 3*x - 1")
    assert find_symmetric_realization(real_cyclotomic(21), max_nodes=50) is None


def test_cli_smoke(capsys):
    from noninterlacing.cli import main

    main(["level", "x^6 - x^5 - 6*x^4 + 6*x^3 + 8*x^2 - 8*x + 1"])
    out = capsys.readouterr().out
    assert "min level m(f): 2" in out and "<2, 3>" in out
    main(["check", "x^2 - 2", "x"])
    assert "strictly interlaces: True" in capsys.readouterr().out
