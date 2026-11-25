"""Screen C_d-C_a-C_b+1 for irreducible Salem traces of measure below 5/4."""

import argparse
from math import gcd

from flint import fmpq, fmpq_poly

from noninterlacing.pari import pari, to_pari
from noninterlacing.poly import chebyshev_C
from noninterlacing.realroots import RootData, precision


def scan(degrees):
    chebyshev = [chebyshev_C(i) for i in range(max(degrees) + 1)]
    seen = set()
    hits = []
    for d in degrees:
        for a in range(1, d):
            for b in range(a):
                if gcd(gcd(a, b), d) != 1:
                    continue
                f = chebyshev[d] - chebyshev[a] - chebyshev[b] + 1
                for h, _ in f.factor()[1]:
                    if h.degree() < 2 or h(2) >= 0:
                        continue
                    key = tuple(map(int, h.coeffs()))
                    if key in seen:
                        continue
                    seen.add(key)
                    # 5/4 + 4/5 = 41/20. Both endpoint signs are necessary.
                    if fmpq_poly(h)(fmpq(41, 20)) <= 0:
                        continue
                    if pari().polsturm(to_pari(h), [-2, 2]) != h.degree() - 1:
                        continue
                    with precision(96):
                        r = RootData(h, check=False).balls(96)[-1]
                        measure = float((r + (r*r - 4).sqrt()) / 2)
                    hits.append((measure, h.degree(), d, a, b))
    print("distinct irreducible factors of degree >= 2 with negative value at 2", len(seen))
    for measure, degree, d, a, b in sorted(hits):
        print(f"measure {measure:.12f}, trace degree {degree}, first parameters {(d, a, b)}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("degrees", nargs="*", type=int,
                        default=[12, 16, 20, 24, 30, 32, 36, 40, 48, 60, 64, 80, 96])
    args = parser.parse_args()
    if min(args.degrees) < 2:
        parser.error("degrees must be at least 2")
    scan(args.degrees)
