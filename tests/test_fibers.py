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


def test_psi231_witnesses():
    # m(Psi_231) = 3: level 2 is empty (fibres; too slow for the test suite), levels 3 and 5 are attained
    w = json.load(open(os.path.join(RESULTS, "psi231.json")))
    f = real_cyclotomic(231)
    for j, s in w["interlacers"].items():
        g = poly(s)
        assert int(g.leading_coefficient()) == int(j) and interlaces(f, g)
