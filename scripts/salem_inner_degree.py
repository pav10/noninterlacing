"""Certify inner degrees excluded by the power-trace bound for lambda <= 5/4."""

from flint import arb, ctx
from noninterlacing.poly import X, is_irreducible, real_cyclotomic
from noninterlacing.realroots import RootData, sturm_count


def excluded(n, b, cap):
    m = n - 2
    r = 1 - 1/arb(n)
    energy = (m*arb(n).log() - (m*(m-1)//2)*r.log()
              + 2*(m+1)*cap.log() + 2*m*arb(4).log())
    for k in range(1, b):
        upper = cap**k + 1 + (k*energy/r**k).sqrt()
        if not arb(k*n)/(2*b) > upper:
            return False
    return True


def quadratic_minimum():
    f = (X**2-3*X+1)(X**2-X)
    assert is_irreducible(f) and sturm_count(f, -2, 2) == 3 and f(2) < 0
    r = (1 + (7 + 2*arb(5).sqrt()).sqrt())/2
    assert r < arb(11)/5
    for m in (7, 8, 9):
        assert RootData(real_cyclotomic(m)).balls(96)[0] < -arb(5)/4
    print("quadratic inner minimum", (r + (r*r-4).sqrt())/2)


if __name__ == "__main__":
    with ctx.workprec(128):
        quadratic_minimum()
        for n in (128, 1024, 4096, 65536, 1048576, 10**8):
            cap = arb(5)/4
            b = 2
            while excluded(n, b, cap):
                b += 1
            result = f"2 through {b-1}" if b > 2 else "none"
            print("Salem degree", n, "inner degrees excluded by the energy test", result)
