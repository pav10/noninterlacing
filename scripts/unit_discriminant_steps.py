"""Check totally real reciprocal lifts with identical polynomial root discriminant."""

import argparse
from math import exp, log

from noninterlacing.pari import pari, to_pari
from noninterlacing.poly import X, is_irreducible
from noninterlacing.realroots import sturm_count
from noncyclotomic_certificates import homogeneous_value


def least_second_cost(g, units):
    """Certify the minimum for the a=4 example using its explicit upper bound."""
    bound = 1801
    best = None
    for x, y in sorted(units):
        radius = 1 + max(abs(x), abs(y))
        entries = range(-radius, radius + 1)
        b, d = next((b, d) for b in entries for d in entries if x*d-y*b == 1)
        h = homogeneous_value(g, x*X+b, y*X+d)
        h *= h.leading_coefficient()  # unit leading coefficient
        def cost(s):
            return abs(int(h(s-2) * h(s+2)))

        left = right = 0
        while sturm_count(h, None, left+2) or cost(left) <= bound:
            left -= 1
        while sturm_count(h, right-2, None) or cost(right) <= bound:
            right += 1
        for s in range(left, right+1):
            value = cost(s)
            if not 0 < value <= bound or sturm_count(h, s-2, s+2):
                continue
            lifted = homogeneous_value(h, X**2+s*X+1, X)
            if is_irreducible(lifted):
                assert lifted.discriminant() == g.discriminant()**2 * value
                best = value if best is None else min(best, value)
    if best != bound:
        raise ArithmeticError(f"unexpected minimum {best}")
    return best


def check(a):
    f = (X**2 - 4) * (X - a) + 1
    g = (X**2 - 1)**2 * (X**2 - a*X + 1) + X**3
    assert f(2) == f(-2) == 1 and sturm_count(g) == 6
    assert is_irreducible(g) and g.discriminant() == f.discriminant()**2
    p = pari()
    thue = p.thueinit(to_pari(g), 1)
    units = {tuple(map(int, v)) for sign in (-1, 1) for v in p.thue(thue, sign)}
    # The leading coefficient and the two endpoint evaluations must be units.
    # Write the -2 endpoint as (b,d) so the +2 endpoint is (b+4*x,d+4*y).
    triples = [(x, y, b, d) for x, y in units for b, d in units
               if abs(x*d-y*b) == 1 and (b+4*x, d+4*y) in units]
    print(f"a={a}  cubic discriminant={f.discriminant()}  "
          f"unit vectors={len(units)}  candidate transformations={len(triples)}")
    if a == 4:
        print("least second-step norm", least_second_cost(g, units))
        for x, b, y, d, s in [(1, -1, -1, 2, -1), (0, -1, 1, -2, -3),
                              (0, -1, 1, -2, 6)]:
            h = homogeneous_value(g, x*X+b, y*X+d)
            assert abs(h.leading_coefficient()) == 1
            h *= h.leading_coefficient()
            cost = abs(int(h(s-2)*h(s+2)))
            lifted = homogeneous_value(h, X**2+s*X+1, X)
            assert cost and not sturm_count(h, s-2, s+2) and is_irreducible(lifted)
            disc = int(lifted.discriminant())
            assert disc == g.discriminant()**2 * cost
            print("degree", lifted.degree(), "cost", cost, "root discriminant",
                  exp(log(disc)/lifted.degree()), "below 200", disc < 200**lifted.degree())
            g = lifted
    return triples


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("parameters", nargs="*", type=int, default=[4, 5, 6, 7, 8, 9, 10, 11, 12, 16, 24])
    args = parser.parse_args()
    if min(args.parameters) < 4:
        parser.error("parameters must be at least 4")
    for parameter in args.parameters:
        check(parameter)
