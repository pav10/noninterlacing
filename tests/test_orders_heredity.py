import pytest
from flint import fmpq_poly

from noninterlacing.heredity import composition_bound, intersect_with_subfield, relative_trace
from noninterlacing.interlace import Interlacing
from noninterlacing.orders import Order, integer_kernel, power_sums
from noninterlacing.poly import compose, poly, real_cyclotomic
from noninterlacing.semigroup import NumericalSemigroup

from .conftest import needs_pari

F21 = real_cyclotomic(21)


def test_semigroups():
    S = NumericalSemigroup.from_generators([3, 5])
    assert S.gaps == [1, 2, 4, 7] and S.frobenius == 7 and S.generators == [3, 5]
    assert NumericalSemigroup.from_generators([1]).gaps == []
    T = NumericalSemigroup.from_generators([4, 6, 9])
    assert T.generators == [4, 6, 9] and T.gaps == [1, 2, 3, 5, 7, 11] and 13 in T


def test_power_sums():
    assert power_sums(poly("x^2 - x - 1"), 5) == [2, 1, 3, 4, 7]


def test_integer_kernel_is_saturated():
    from flint import fmpz_mat

    # kernel of (2, 4, 6) on Z^3; with p = (1, 0, 0) (t.p = gcd = 2) it must complete to a basis
    K = integer_kernel([[2, 4, 6]], 3)
    assert len(K) == 2 and all(2 * a + 4 * b + 6 * c == 0 for a, b, c in K)
    assert abs(int(fmpz_mat([K[0], K[1], [1, 0, 0]]).det())) == 1
    # a rank-2 system in Z^4: every small kernel vector is an integer combination of the basis
    import itertools
    from fractions import Fraction

    K2 = integer_kernel([[3, 6, 9, 12], [1, 1, 1, 1]], 4)
    assert len(K2) == 2
    a, b = K2
    i, j = next((i, j) for i, j in itertools.combinations(range(4), 2) if a[i] * b[j] != a[j] * b[i])
    det = a[i] * b[j] - a[j] * b[i]
    for v in itertools.product(range(-4, 5), repeat=4):
        if 3 * v[0] + 6 * v[1] + 9 * v[2] + 12 * v[3] == 0 and sum(v) == 0:
            x = Fraction(v[i] * b[j] - v[j] * b[i], det)
            y = Fraction(a[i] * v[j] - a[j] * v[i], det)
            assert x.denominator == 1 and y.denominator == 1
            assert all(x * a[k] + y * b[k] == v[k] for k in range(4))


def test_monogenic_order_matches_interlacing():
    for f in (F21, real_cyclotomic(28), poly("x^4 - 4*x^2 + 2")):
        O = Order.monogenic(f)
        I = Interlacing(f)
        assert O.min_trace() == I.min_level()
        assert O.count(2) == I.count(2)
        assert O.codifferent_traces.count(1) >= 1


def test_suborder_monotonicity():
    f = poly("x^3 - 3*x - 1")
    O = Order.monogenic(f)
    sub = Order.generated_by(f, [fmpq_poly([0, 2])])  # Z[2 alpha]
    assert O.contains_order(sub) and sub.index_in(O) == 8
    assert sub.min_trace() <= O.min_trace()


def test_heredity_composition_and_trace_down():
    b = composition_bound(F21, "x^2 - x - 4")
    assert b.m == 2 and (b.maximal in (True, None))
    f = compose(F21, poly("x^2 - x - 4"))
    O = Order.monogenic(f)
    OF = Order.monogenic(F21)
    for mu in O.totally_positive_codifferent(2, limit=5):
        nu = relative_trace(O, mu, fmpq_poly(poly("x^2 - x - 4")), F21)
        assert OF.trace(nu) == 2
        assert OF.roots.is_totally_positive(nu)
        assert all(OF.trace((nu * b) % fmpq_poly(F21)).q == 1 for b in OF.basis)


@needs_pari
def test_maximal_orders_and_subfields():
    from noninterlacing.heredity import heredity_report

    OK = Order.maximal(F21)
    assert OK.discriminant == 453789 and OK.min_trace() == 2
    rep = heredity_report(real_cyclotomic(84))
    assert rep.lower_bound == 2
    nonmono = Order.maximal(poly("x^3 - x^2 - 10*x + 8"))  # disc(poly) != disc(K)
    assert Order.monogenic(poly("x^3 - x^2 - 10*x + 8")).index_in(nonmono) == 2
