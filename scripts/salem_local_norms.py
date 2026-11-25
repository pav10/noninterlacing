"""Check the local norm estimate on Salem traces and give finite splitting benchmarks."""

from math import ceil, log

from flint import arb, ctx
from noninterlacing.pari import pari, to_pari
from noninterlacing.poly import X
from noninterlacing.realroots import RootData


def lower_bound(n, prime, f, places, k):
    mass = f*places*arb(prime).log()
    return 3*(mass - (n-2)*arb(k).log()/(k-1))/((k+1)*prime**f)


def check(trace):
    gp = pari()
    t = to_pari(trace)
    d = trace.degree()
    assert gp.polisirreducible(t) and gp.polsturm(t, [-2, 2]) == d-1 and trace(2) < 0
    nf = gp.nfinit(t)
    r = RootData(trace).balls(128)[-1]
    L = ((r+(r*r-4).sqrt())/2).log()
    z = gp.Mod(gp("x"), t)
    for p in (2, 3, 5, 7):
        degrees = [int(pr[3]) for pr in gp.idealprimedec(nf, p)]
        for f in set(degrees):
            q = p**f
            if q > 500:
                continue
            # C_q(x) = 2*T_q(x/2); binary doubling in the trace field.
            a, b = z**0*2, z
            for bit in bin(q)[2:]:
                ab = a*b-z
                a, b = (ab, b*b-2) if bit == '1' else (a*a-2, ab)
            norm = abs(int(gp.norm(a-z)))
            count = degrees.count(f)
            assert norm and norm % (p**(f*count)) == 0
            for k in (2, 8, 32, 128):
                rhs = (2*d-2)*arb(k).log()/(k-1) + (k+1)*q*L/3
                assert arb(norm).log() < rhs
    print("checked trace degree", d, "Salem measure", float(L.exp()))


if __name__ == "__main__":
    with ctx.workprec(160):
        check(X**5+X**4-5*X**3-5*X**2+4*X+3)
        t = (X-2)*(X**3+X**2-2*X-1)
        check(t*t-t-1)
        gp = pari()
        b = gp.bnfinit(gp.quadpoly(14401, "y"), 1)
        assert gp.bnfcertify(b) and list(b.bnf_get_cyc()) == [43]
        assert gp.kronecker(14401, 11) == -1
        bound = lower_bound(172, 11, 2, 43, 2)
        print("degree-86 Hilbert trace field, Salem measure >=", bound.exp())
        for n in (4096, 1048576):
            places = int(n**0.5)
            k = 1+ceil(2*n*log(n)/(places*log(2)))
            bound = lower_bound(n, 2, 1, places, k)
            assert bound > 0
            print("conditional n", n, "norm-2 places", places, "k", k, "log-measure >=", bound)
