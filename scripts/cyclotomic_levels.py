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
        margins = {}
        j = 1
        while True:
            mm = I.margin(j)
            margins[j] = mm.value
            ex = mm.exists
            if ex is None:  # MILP undecided: fall back to the complete search
                ex = I.find(j, method="dfs") is not None
            if ex:
                break
            j += 1
        row = {"n": n, "degree": e, "m": j, "margins": margins,
               "admissible": is_admissible(n), "admissible_divisors": admissible_divisors(n),
               "seconds": round(time.time() - t, 2)}
        rows.append(row)
        rows.sort(key=lambda r: r["n"])
        print(n, e, j, row["admissible"], row["admissible_divisors"][:3], round(margins[1], 6), row["seconds"], flush=True)
        json.dump(rows, open(out, "w"), indent=1)


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]) if len(sys.argv) > 3 else 60)
