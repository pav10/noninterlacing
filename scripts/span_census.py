"""Census of irreducible totally real polynomials of degree d and span <= S, with levels.

Usage: python scripts/span_census.py DEGREE MAX_SPAN OUT.json

Every monic irreducible totally real f in Z[x] of degree d with span <= S is
listed once up to x -> +-x + n, with m(f) (complete search), the level-1
margin and, when PARI is available, m(O_K).
"""
import json
import sys
import time

from noninterlacing.interlace import Interlacing
from noninterlacing.poly import to_str
from noninterlacing.search import totally_real_census


def main(d, S, out):
    t0 = time.time()
    rows = []
    for f in totally_real_census(d, S, irreducible_only=True):
        I = Interlacing(f, check=False)
        row = {"f": to_str(f), "span": round(I.roots.span, 6), "m": I.min_level(),
               "margin1": round(I.margin(1).value, 8)}
        try:
            from noninterlacing.orders import Order

            row["m_OK"] = Order.maximal(f).min_trace()
        except Exception:
            pass
        rows.append(row)
        print(json.dumps(row), flush=True)
        json.dump({"degree": d, "max_span": S, "classes": rows, "seconds": round(time.time() - t0, 1)},
                  open(out, "w"), indent=1)
    json.dump({"degree": d, "max_span": S, "classes": rows, "seconds": round(time.time() - t0, 1)},
              open(out, "w"), indent=1)


if __name__ == "__main__":
    main(int(sys.argv[1]), float(sys.argv[2]), sys.argv[3])
