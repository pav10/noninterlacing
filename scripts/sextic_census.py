"""Interlacing census of the totally real sextic fields computed by sails.

Usage: python -I scripts/sextic_census.py OUT.json FILE.jsonl[.xz] ...

For every field record (sails JSON format) this computes m(O_K) (Yatsyna's
minimal trace of a totally positive codifferent element), the level-1
margin, the trace semigroup of O_K, m(Z[a]) for the record's defining
polynomial, and compares m(O_K) with the bound min{d_F : 1 in F} read off
the sail facets through the vertex 1.
"""
import json
import sys
import time

from flint import fmpz_poly

from noninterlacing.orders import Order
from noninterlacing.sailsdata import facet_bound_at_one, parse_element, read_records


def main(out, paths):
    rows = []
    t0 = time.time()
    for path in paths:
        for rec in read_records(path):
            f = fmpz_poly(rec["coeffs"])
            OK = Order(f, [parse_element(t) for t in rec["integral_basis"]], name=rec["label"])
            Za = Order.monogenic(f)
            mK = OK.min_trace()
            row = {
                "label": rec["label"],
                "disc": rec["disc"],
                "coeffs": rec["coeffs"],
                "index": int(Za.index_in(OK)),
                "m_OK": mK,
                "m_Za": Za.min_trace(),
                "margin_OK": OK.margin(1).value,
                "trace_semigroup_OK": str(OK.trace_semigroup()),
                "trace1_count_OK": OK.count(1),
                "sails_bound": facet_bound_at_one(rec),
                "num_facets": next(o for o in rec["orthants"] if not any(o["signature"]))["num_facets"],
            }
            b = row["sails_bound"]
            row["consistent"] = b is None or (mK <= b and (b != 1 or mK == 1))
            rows.append(row)
            print(row["label"], row["m_OK"], row["m_Za"], b, round(row["margin_OK"], 5), flush=True)
    json.dump({"fields": rows, "seconds": round(time.time() - t0, 1)}, open(out, "w"), indent=1)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2:])
