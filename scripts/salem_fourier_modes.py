"""Certify the best single Fourier mode in the Salem correlation estimate."""

from flint import arb, ctx, fmpz_poly

from noninterlacing.pari import pari, to_pari
from noninterlacing.poly import X, chebyshev_C


def unique_peak(values):
    i = max(range(len(values)), key=lambda j: float(values[j]))
    assert all(v < values[i] for j, v in enumerate(values) if i != j)
    return i+1, values[i]


def peaks(p):
    n = p.degree()
    m = n-2
    d = n//2
    assert p.coeffs() == p.coeffs()[::-1] and pari().polisirreducible(to_pari(p))
    trace = sum((int(p[d+j])*chebyshev_C(j) for j in range(1, d+1) if p[d+j]),
                fmpz_poly([p[d]]))
    assert trace(2) < 0 and pari().polsturm(to_pari(trace), [-2, 2]) == d-1
    roots = [z for z, _ in p.complex_roots() if z.imag > 0]
    assert len(roots) == m//2
    r = 1-1/arb(n)
    penalty = arb(n).log()-(m-1)*r.log()/2+2*arb(4).log()
    scale = 2*(1+1/arb(m))
    raw, adjusted = [], []
    powers = roots.copy()
    best_raw = best_adjusted = arb(0)
    j = 1
    while True:
        moment = 2*sum(z.real for z in powers)/m
        square = moment**2
        raw.append(square/j)
        adjusted.append((r*m*square-penalty)/(scale*j))
        if raw[-1] > best_raw:
            best_raw = raw[-1]
        if adjusted[-1] > best_adjusted:
            best_adjusted = adjusted[-1]
        adjusted_done = (r*m <= penalty or
                         (r*m-penalty)/(scale*(j+1)) < best_adjusted)
        if arb(1)/(j+1) < best_raw and adjusted_done:
            break
        # Restart at regular intervals to avoid interval wrapping under rotations.
        if j % 16 == 0:
            powers = [z**(j+1) for z in roots]
        else:
            powers = [v*z for v, z in zip(powers, roots)]
        j += 1
    adjusted_peak = unique_peak(adjusted) if r*m > penalty else (None, arb(0))
    return unique_peak(raw), adjusted_peak, j


def shifted_bound(n, support, shifts):
    exponents = sorted({i+j for i in support for j in shifts})
    k = len(exponents)
    h = len(shifts)
    # Sum of v-u over unordered pairs u<v, computed in sorted order.
    avg = 2*arb(sum((2*i-k+1)*v for i, v in enumerate(exponents)))/(k*(k-1))
    m = n-2
    r = 1-1/arb(n)
    penalty = arb(n).log()-(m-1)*r.log()/2+2*arb(4).log()
    return (r*m*h/((k-h)*(k-1))-penalty)/(2*(1+1/arb(m))*avg)


def main():
    seed = X**3-X-1
    reverse = fmpz_poly(seed.coeffs()[::-1])
    for exponent in (8, 31, 63, 127, 255, 511):
        p, = [f for f, _ in (X**exponent*seed-reverse).factor()[1] if f(1) < 0]
        (j, u), (k, lower), stop = peaks(p)
        support = {0, 2, 3, exponent, exponent+1, exponent+3}
        bound = shifted_bound(p.degree(), support, range(4))
        assert bound < lower
        print("seed exponent", exponent, "degree", p.degree(), "raw mode", j,
              "J^2 B / n >=", float((1/(p.degree()*u)*10**8).floor())/10**8, "adjusted mode", k,
              "log measure >=", float((lower*10**8).floor())/10**8, "last tested", stop,
              "shift bound", float(bound), flush=True)
    theta, = [z.real for z, _ in seed.complex_roots() if z.imag.is_zero()]
    mean = -1-reverse(theta)/(theta*seed.derivative()(theta))
    assert mean < -arb(45488)/100000
    print("limiting ordinary moment at the seed exponent", mean)


if __name__ == "__main__":
    with ctx.workprec(192):
        main()
