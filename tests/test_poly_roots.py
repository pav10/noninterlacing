import random

from flint import fmpq_poly, fmpz_poly

from noninterlacing.poly import (canonical_form, chebyshev_C, compose, poly, real_cyclotomic,
                                 reflect, to_str, translate, z_equivalent)
from noninterlacing.realroots import (RootData, cauchy_index, interlaces, is_totally_real,
                                      sturm_count, tarski_query)

F21 = "x^6 - x^5 - 6*x^4 + 6*x^3 + 8*x^2 - 8*x + 1"


def test_parse_print_roundtrip():
    f = poly(F21)
    assert to_str(f) == F21
    assert poly([1, -1, -6, 6, 8, -8, 1]) == f
    assert poly("(x^2 - x - 3)^2 + 2x - 1") == poly("x^4 - 2*x^3 - 5*x^2 + 8*x + 8")
    assert poly("x**3 - 2 x + 1") == poly([1, 0, -2, 1])


def test_real_cyclotomic_and_chebyshev():
    assert real_cyclotomic(21) == poly(F21)
    assert chebyshev_C(5) == poly("x^5 - 5*x^3 + 5*x")
    # C_k(beta_n) = beta_{n/gcd}: m_beta_21(C_3) = m_beta_63
    assert compose(real_cyclotomic(21), chebyshev_C(3)) == real_cyclotomic(63)


def test_translation_classes():
    f = poly(F21)
    assert canonical_form(translate(f, 7)) == canonical_form(f)
    assert z_equivalent(f, translate(f, -3))
    assert canonical_form(reflect(f), reflections=True) == canonical_form(f, reflections=True)


def test_sturm_and_interlacing_exact():
    f = poly(F21)
    assert sturm_count(f) == 6 and is_totally_real(f)
    assert not is_totally_real(poly("x^2 + 1")) and not is_totally_real(poly("(x-1)^2*(x+2)"))
    fp = f.derivative()
    assert interlaces(f, fp) and cauchy_index(fp, f) == 6
    assert tarski_query(fmpz_poly([1]), f) == 6
    assert not interlaces(f, poly("x^5"))
    # residue signs agree with the exact test on random candidates
    rd = RootData(f)
    rng = random.Random(1)
    for _ in range(200):
        g = fmpz_poly([rng.randint(-15, 15) for _ in range(5)] + [rng.randint(1, 3)])
        assert rd.is_positive_interlacer(g) == (interlaces(f, g) and int(g.leading_coefficient()) > 0)


def test_rootdata_signs_with_common_roots():
    f = poly("(x^2 - 2)*(x^2 - 3)")
    rd = RootData(f)
    s = rd.signs(poly("x^2 - 2"))
    assert sorted(s) == [-1, -1, 0, 0] or s.count(0) == 2
