import random

from flint import fmpq_poly

from noninterlacing.fibers import CyclotomicField, interlacer_from_element
from noninterlacing.parity import (block_array, corner_parity, level_two_census, parity_terms,
                                   power_sums, three_q_family, trace_array)

from .conftest import needs_pari


@needs_pari
def test_parity_identity_on_random_codifferent_elements():
    # the arrays are even with zero line sums, and X + Y + Z is even (notes/parity.tex, Lemma 3.1)
    K = CyclotomicField(105)
    s, B = power_sums(K), K.codifferent([1])
    rng = random.Random(1)
    for _ in range(3):
        x = sum((rng.randint(-3, 3) * b for b in B), fmpq_poly([0]))
        c = trace_array(K, x, s)
        for blocks in ((3, 5, 7), (7, 3, 5), (5, 7, 3)):
            C = block_array(c, blocks)
            for key in C:
                assert C[key] == C[tuple(-a % q for a, q in zip(key, blocks))]
                for i, q in enumerate(blocks):
                    if key[i] == 0:
                        assert sum(C[key[:i] + (t,) + key[i + 1:]] for t in range(q)) == 0
            assert sum(parity_terms(C, blocks)) % 2 == 0


@needs_pari
def test_products_have_odd_corner_parity():
    # level-1 elements a = (2 - alpha_5)/5, b = (2 - alpha_7)/7: ab has axes 2a, 2b and tau = 1
    K = CyclotomicField(35)
    al = fmpq_poly([0, 1])
    ab = ((2 - K.embed(al, 5)) / 5 * (2 - K.embed(al, 7)) / 7) % K.F
    C = block_array(trace_array(K, ab), (5, 7))
    assert C[(0, 0)] == 4 and all(C[(x, 0)] % 2 == 0 for x in range(5))
    assert corner_parity(C, 5, 7) == 1


@needs_pari
def test_level_two_census_small_fields():
    assert level_two_census(3, 5) == {"field": "F_15", "degree": 4, "level2": 50, "product": 2,
                                      "tau=0": 12, "tau=1": 12, "other": 24}
    assert level_two_census(3, 7) == {"field": "F_21", "degree": 6, "level2": 141, "product": 3,
                                      "tau=0": 18, "tau=1": 18, "other": 102}


def test_three_q_family_positive_exactly_for_q_at_most_7():
    for q, positive in ((5, True), (7, True), (11, False), (13, False)):
        K = CyclotomicField(3 * q)
        for g in three_q_family(q):
            assert trace_array(K, g)[0] == 4  # in the codifferent, trace 2
            try:
                interlacer_from_element(K, g)
                ok = True
            except ArithmeticError:
                ok = False
            assert ok == positive, q
