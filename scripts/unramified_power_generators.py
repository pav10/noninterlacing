"""Exact discriminant certificates for generators of ordinary Hilbert class fields."""

from math import exp, isqrt, log, prod

from noninterlacing.pari import pari


def collision_valuation(n, residue_size):
    total = 0
    classes = residue_size
    while classes < n:
        a, r = divmod(n, classes)
        total += classes * a * (a - 1) + 2 * r * a
        classes *= residue_size
    return total


def quadratic_obstruction(discriminant, prime):
    p = pari()
    b = p.bnfinit(p.quadpoly(discriminant, "y"))
    assert p.bnfcertify(b) and int(p.nfdisc(p.quadpoly(discriminant, "y"))) == discriminant
    assert p.isprime(prime) and p.kronecker(discriminant, prime) == -1
    n = 2 * prod(map(int, b.bnf_get_cyc()))
    v = collision_valuation(n, prime**2)
    lower_bound = discriminant**(n//2) * prime**v
    assert lower_bound > 200**n
    print(discriminant, n, prime, v, exp(log(lower_bound) / n))


def cubic_generators():
    p = pari()
    base = p("y^3 + y^2 - 17*y + 16")
    b = p.bnfinit(base)
    assert p.bnfcertify(b) and list(b.bnf_get_cyc()) == [4]
    assert p.nfdisc(base) == 8069
    for relative_degree in (2, 4):
        f = p.polredbest(p.bnrclassfield(b, relative_degree, 2))
        n = 3 * relative_degree
        field_disc = 8069**relative_degree
        disc = int(p.poldisc(f))
        index = isqrt(disc // field_disc)
        assert p.poldegree(f) == n and p.pollead(f) == 1 and p.polisirreducible(f)
        assert p.polsturm(f) == n and p.nfdisc(f) == field_disc
        assert field_disc * index**2 == disc < 200**n
        print(n, index, disc**(1/n), f)


if __name__ == "__main__":
    for D, prime in ((7057, 5), (32401, 7), (36865, 7)):
        quadratic_obstruction(D, prime)
    cubic_generators()
