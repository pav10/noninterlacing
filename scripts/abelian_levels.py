"""m(O_K) for real abelian fields K (real subfields of Q(zeta_n)), searching for m >= 3.

Usage: python scripts/abelian_levels.py OUT.json NMAX MINDEG MAXDEG

Fields are deduplicated by polredabs.  For each field: degree, conductor
(first n where it appears), discriminant, m(O_K) with certificate, and the
Siegel ceiling dim M_2d.
"""
import json
import os
import sys
import time

from flint import fmpz_poly

from noninterlacing.cyclotomic import euler_phi
from noninterlacing.orders import Order
from noninterlacing.pari import from_pari, nfdisc, pari, polredabs
from noninterlacing.poly import to_str
from noninterlacing.realroots import is_totally_real
from noninterlacing.siegel import dim_M


def main(out, nmax, mindeg, maxdeg):
    P = pari()
    rows = json.load(open(out)) if os.path.exists(out) else []
    seen = {r["poly"] for r in rows}
    done_n = max([r["conductor"] for r in rows], default=2)
    for n in range(max(3, done_n), nmax + 1):
        if n % 4 == 2:
            continue
        phi = euler_phi(n)
        for e in range(mindeg, min(maxdeg, phi // 2) + 1):
            if phi % e:
                continue
            res = P.polsubcyclo(n, e)
            polys = list(res) if str(P.type(res)) == "t_VEC" else [res]
            for gp in polys:
                g = from_pari(gp)
                g = fmpz_poly([int(c) for c in g.numer().coeffs()])
                if not is_totally_real(g):
                    continue
                g = polredabs(g)
                key = to_str(g)
                if key in seen:
                    continue
                seen.add(key)
                t = time.time()
                O = Order.maximal(g)
                cert = O.level_certificate(time_limit=200)
                row = {"poly": key, "degree": e, "conductor": n, "disc": O.discriminant,
                       "m": cert["m"], "empty_levels": cert["empty_levels"],
                       "witness_method": cert.get("witness_method"), "ceiling": dim_M(2 * e),
                       "seconds": round(time.time() - t, 2)}
                rows.append(row)
                print(n, e, row["m"], row["witness_method"], row["seconds"], key[:60], flush=True)
                json.dump(rows, open(out + ".tmp", "w"), indent=1)
                os.replace(out + ".tmp", out)


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]))
