from noninterlacing.bruteforce import lp_search_interlacers
from noninterlacing.cyclotomic import (compositum_splits, compositum_witness, factorint,
                                       psi_interlacing_predicted)
from noninterlacing.interlace import Interlacing
from noninterlacing.kronecker import (amgm_lower_bound, check_disc_polynomial, components, convolve_counts,
                                      cyclotomic_products, factor_as_psi, linked, ms_all_pairs,
                                      ms_monic_interlacers, product_poly)
from noninterlacing.poly import compose, poly, real_cyclotomic
from noninterlacing.realroots import interlaces

from .conftest import needs_pari


def test_level1_equals_mckee_smyth():
    pairs = ms_all_pairs(12)
    for d in range(1, 7):
        for t in cyclotomic_products(d):
            I = Interlacing(product_poly(t), check=False)
            mine = sorted(tuple(sorted(factor_as_psi(g))) for g in I.interlacers(1))
            assert mine == ms_monic_interlacers(t, pairs), t


def test_psi_interlacing_formula():
    members = {m for pr in ms_all_pairs(80) for m in pr}

    def pred(n):
        odd = [p for p in factorint(n) if p != 2]
        r = 1
        for p in odd:
            r *= p
        return r in (1, 15) or len(odd) == 1

    assert all(pred(n) == psi_interlacing_predicted(n) for n in range(3, 1000))
    for n in range(3, 160):
        if real_cyclotomic(n).degree() <= 40:
            assert (frozenset([n]) in members) == pred(n), n
    for n in (21, 28, 33, 35, 45, 60, 75):
        assert (Interlacing(real_cyclotomic(n)).min_level() == 1) == pred(n), n


def test_compositum_witness():
    # m(Psi_n) <= 2 m(Psi_n1) m(Psi_n2) for coprime n = n1 n2 (n1, n2 > 2)
    I21 = Interlacing(real_cyclotomic(21))
    g = compositum_witness(21, 3, 7)
    assert g.coeffs()[-1] == 2 and g in I21.interlacers(2)
    assert compositum_witness(84, 7, 12).coeffs()[-1] == 2
    assert compositum_witness(115, 5, 23).degree() == 43
    # levels multiply: Psi_4 (level 1) and Psi_21 (level 2) give level 2 * 1 * 2 for Psi_84
    g = compositum_witness(84, 4, 21, None, I21.find(2))
    assert g.coeffs()[-1] == 4 and interlaces(real_cyclotomic(84), g)
    assert compositum_splits(231) == [] and compositum_splits(165) == [(11, 15)]
    assert compositum_splits(396) == [(9, 44), (11, 36)]


def test_unit_resultant_additivity():
    assert linked(3, 6) and linked(21, 63) and not linked(21, 33) and not linked(3, 4)
    for t in [(5, 7), (3, 4, 5), (1, 6, 10), (3, 6, 5, 7)]:
        comps = components(t)
        per = []
        for c in comps:
            I = Interlacing(product_poly(c), check=False)
            per.append([0] + [I.count(j) for j in range(1, 4)])
        If = Interlacing(product_poly(t), check=False)
        assert convolve_counts(*per)[:4] == [0] + [If.count(j) for j in range(1, 4)]
        assert If.min_level() == sum(Interlacing(product_poly(c), check=False).min_level() for c in comps)
    # not unit resultant: strictly smaller level than the sum
    assert Interlacing(product_poly((21, 42)), check=False).min_level() == 2


def test_amgm_bound_and_unbounded_levels():
    for t in [(3, 4, 5, 7), (5, 7, 11), (21,), (3, 4, 6)]:
        f = product_poly(t)
        assert amgm_lower_bound(f)[0] <= Interlacing(f, check=False).min_level()
    S = Interlacing(product_poly((3, 4, 5, 7, 11)), check=False).levels()
    assert S.multiplicity == 5 and S.gaps == [1, 2, 3, 4]


def test_level_j_disc_polynomial():
    for t in [(28,), (3, 4, 6), (1, 5), (21,)]:
        I = Interlacing(product_poly(t), check=False)
        for j in (1, 2, 3):
            for g in I.interlacers(j, limit=5):
                r = check_disc_polynomial(t, g)
                assert r["lead"] == j + 1 and r["max|root|"] < 1 and r["two_reciprocal"] and r["F_recovered"]


def test_noninterlacing_octic():
    f = poly("x^8 - 4*x^7 - x^6 + 17*x^5 - 5*x^4 - 23*x^3 + 6*x^2 + 9*x - 1")
    assert compose("x^4 - 7*x^3 + 15*x^2 - 9*x - 1", "x^2 - x") == f
    I = Interlacing(f)
    assert I.count(1) == 0 and I.margin(1).exists is False
    lp1, _ = lp_search_interlacers(f, 1)
    assert lp1 == []
    assert I.count(2) == 406 and I.levels().generators == [2, 3]


@needs_pari
def test_noninterlacing_octic_siegel():
    from noninterlacing.orders import Order
    from noninterlacing.siegel import siegel_zeta

    f = poly("x^8 - 4*x^7 - x^6 + 17*x^5 - 5*x^4 - 23*x^3 + 6*x^2 + 9*x - 1")
    assert Order.maximal(f).min_trace() == 2
    z = siegel_zeta(f)
    assert z["s_K"][0] == 0 and z["agree"]
