"""Interlacing levels of all squarefree cyclotomic polynomials (roots in [-2, 2]).

Usage: python scripts/kronecker_census.py OUT.json DMAX

For every product f of distinct Psi_n of degree <= DMAX (one of f(x),
(-1)^d f(-x)) this records m(f), the level semigroup S(f), the number of
monic interlacers, the McKee--Smyth prediction for them, the AM--GM lower
bound for m(f) and the level-1 margin.
"""
import json
import os
import sys
import time

from noninterlacing.interlace import Interlacing
from noninterlacing.kronecker import (amgm_lower_bound, cyclotomic_products, factor_as_psi, label,
                                      ms_all_pairs, ms_monic_interlacers, product_poly)


def save(rows, out):
    json.dump(rows, open(out + ".tmp", "w"))
    os.replace(out + ".tmp", out)


def main(out, dmax):
    rows = json.load(open(out)) if os.path.exists(out) else []
    done = {tuple(r["indices"]) for r in rows}
    pairs = ms_all_pairs(2 * dmax)
    t0 = last = time.time()
    for d in range(1, dmax + 1):
        for t in cyclotomic_products(d):
            if t in done:
                continue
            f = product_poly(t)
            I = Interlacing(f, check=False)
            S = I.levels()
            mine = sorted(tuple(sorted(factor_as_psi(g))) for g in I.interlacers(1))
            ms = ms_monic_interlacers(t, pairs)
            lb, s = amgm_lower_bound(f)
            rows.append({"indices": list(t), "factors": label(t), "degree": d, "irreducible": len(t) == 1,
                         "m": S.multiplicity, "semigroup": str(S), "gaps": S.gaps,
                         "interval": S.gaps == list(range(1, S.multiplicity)),
                         "I_1": len(mine), "ms_agrees": mine == ms, "amgm_bound": lb, "amgm_sum": s,
                         "margin1": I.margin(1).value})
            if time.time() - last > 120:
                save(rows, out)
                last = time.time()
        save(rows, out)
        print(d, len(rows), round(time.time() - t0, 1), flush=True)


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]))
