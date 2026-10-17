import json
import os

from flint import fmpq_poly

from noninterlacing.fibers import CyclotomicField, tower_levels
from noninterlacing.interlace import Interlacing
from noninterlacing.poly import poly, real_cyclotomic
from noninterlacing.realroots import interlaces

from .conftest import needs_pari

RESULTS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results")


def _chain(K, *ms):
    chain, fix = [K.G], set(K.G)
    for m in ms:
        fix &= set(K.fixer(m))
        chain.append(sorted(fix))
    return chain + [[1]]


def test_galois_action_and_subfields():
    K = CyclotomicField(21)
    a = fmpq_poly([0, 1])
    # sigma_2 sigma_4 = sigma_8 and sigma_a(alpha) = C_a(alpha)
    assert K.act(2, K.act(4, a)) == K.act(8, a) == K.cheb(8)
    H7 = K.fixer(7)
    beta7 = K.cheb(3)  # 2cos(2 pi 3/21) generates F_7
    assert all(K.act(h, beta7) == beta7 for h in H7)
    assert len(K.codifferent(H7)) == 3 and len(K.codifferent([1])) == 6


@needs_pari
def test_fibre_towers_match_engine():
    for n, ms, levels in ((21, (7,), (1, 2, 3)), (35, (7, 5), (1, 2)), (33, (11,), (2,))):
        K = CyclotomicField(n)
        I = Interlacing(real_cyclotomic(n))
        for j in levels:
            assert tower_levels(K, _chain(K, *ms), j)[-1].count == I.count(j), (n, j)


@needs_pari
def test_fibre_over_non_positive_target_is_empty():
    K = CyclotomicField(21)
    H7 = K.fixer(7)
    nu = fmpq_poly([1]) - K.cheb(3)  # 1 - beta_7 lies in F_7 and is not totally positive
    assert K.fibre([1], H7, nu).points == []


def test_psi231_psi273_psi399_witnesses():
    # m = 3 for n = 231, 273, 399: level 2 is empty (fibres; too slow for the test suite), level 3 attained
    for n in (231, 273, 399):
        w = json.load(open(os.path.join(RESULTS, f"psi{n}.json")))
        f = real_cyclotomic(n)
        for j, s in w["interlacers"].items():
            g = poly(s)
            assert int(g.leading_coefficient()) == int(j) and interlaces(f, g), (n, j)


@needs_pari
def test_central_search_finds_witness():
    from noninterlacing.fibers import central_search

    K = CyclotomicField(21)
    g = central_search(K, K.fixer(7), 2, 0.3)
    assert g is not None and int(g.leading_coefficient()) == 2 and interlaces(K.f, g)


def test_values_survive_cancellation():
    # 10^40 Psi_n(x) + 1 is 1 at every embedding but cancels catastrophically at fixed precision
    K = CyclotomicField(105)
    x = fmpq_poly(K.f) * 10**40 + 1
    for v in K.values(x, K.G, 128):
        assert abs(float(v.mid()) - 1) < 1e-30 and float(v.rad()) < 2.0 ** -100


def test_parse_rational_poly_round_trip():
    from noninterlacing.poly import parse_rational_poly

    K = CyclotomicField(21)
    for x in K.codifferent(K.fixer(7)) + K.codifferent([1]):
        assert parse_rational_poly(str(x)) == x


@needs_pari
def test_lp_and_fincke_pohst_fibres_agree():
    K = CyclotomicField(35)
    H7 = K.fixer(7)
    H57 = sorted(set(H7) & set(K.fixer(5)))
    for nu in K.fibre(H7, K.G, fmpq_poly([3])).points[:6]:
        a = sorted(map(str, K.fibre(H57, H7, nu, method="fp").points))
        b = sorted(map(str, K.fibre(H57, H7, nu, method="lp").points))
        assert a == b


@needs_pari
def test_subfield_fixers_and_towers_through_them():
    K = CyclotomicField(39)
    assert [len(K.G) // len(K.subfield_fixer(13, e)) for e in (1, 2, 3, 6)] == [1, 2, 3, 6]
    chain = [K.G, K.subfield_fixer(13, 2), K.fixer(13), [1]]
    assert [s.count for s in tower_levels(K, chain, 2)] == [1, 7, 156, 2214]



@needs_pari
def test_subfield_chains_find_interlacing_subfields():
    # the three index-2 subfields of F_65 (degree 12) are interlacing (abelian table), so m(Psi_65) <= 2 by lifting
    from noninterlacing.fibers import subfield_chain, subgroups_of_prime_order

    K = CyclotomicField(65)
    subs = subgroups_of_prime_order(K, 2)
    assert len(subs) == 3
    for H in subs:
        chain = subfield_chain(K, H)
        assert chain[0] == K.G and chain[-1] == sorted(H)
        assert tower_levels(K, chain, 1)[-1].count > 0
