"""m(O_{F_n}) for real cyclotomic fields F_n = Q(zeta_n)^+ (one conductor per field).

Usage: python scripts/cyclotomic_levels.py OUT.json NMAX [MAXDEG]
"""
import json
import sys
import time

from noninterlacing.cyclotomic import (admissible_divisors, is_admissible,
                                       real_cyclotomic_degree)
from noninterlacing.interlace import Interlacing
from noninterlacing.poly import real_cyclotomic


def main(out, nmax, maxdeg):
    import os

    rows = json.load(open(out)) if os.path.exists(out) else []  # resume
    done = {r["n"] for r in rows}
    for n in range(3, nmax + 1):
        if n % 4 == 2 or n in done:
            continue
        e = real_cyclotomic_degree(n)
        if e < 2 or e > maxdeg:
            continue
        t = time.time()
        I = Interlacing(real_cyclotomic(n))
        cert = I.level_certificate(time_limit=200)
        mm = I.margin(1, time_limit=60)
        row = {"n": n, "degree": e, "m": cert["m"], "empty_levels": cert["empty_levels"],
               "witness_method": cert.get("witness_method"), "witness": cert["witness"],
               "margin1": mm.value if mm.status == "optimal" else None,
               "admissible": is_admissible(n), "admissible_divisors": admissible_divisors(n),
               "seconds": round(time.time() - t, 2)}
        print(n, e, cert["m"], row["admissible"], row["admissible_divisors"][:3], row["margin1"],
              cert.get("witness_method"), row["seconds"], flush=True)
        rows.append(row)
        rows.sort(key=lambda r: r["n"])
        json.dump(rows, open(out, "w"), indent=1)


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]) if len(sys.argv) > 3 else 60)
