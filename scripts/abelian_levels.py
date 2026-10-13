"""m(O_K) for real abelian fields K (real subfields of Q(zeta_n)), searching for m >= 3.

Usage: python scripts/abelian_levels.py OUT.json NMAX MINDEG MAXDEG

Fields are deduplicated by polredabs.  For each field: degree, conductor
(first n where it appears), discriminant, m(O_K) with certificate, and the
Siegel ceiling dim M_2d.  Each field runs in a child process with a 5 GB
memory limit and a 30-minute wall-clock limit.  If that fails and Psi_n
(n = conductor) is interlacing, the trace-down of a witness of Q(zeta_n)^+
gives m = 1 (heredity.trace_down_witness); otherwise the row is recorded with
m = None and the reason.

Trace-down gives m(O_K) <= m(Psi_n) for every K of conductor n, so m >= 3 can
only occur for conductors with m(Psi_n) >= 3; there is none below n = 231
(see results/higher_interlacing_cyclotomic.md, Section 4).
"""
import json
import os
import sys
import time

from flint import fmpz_poly

from noninterlacing.cyclotomic import euler_phi, psi_interlacing_predicted
from noninterlacing.heredity import trace_down_witness
from noninterlacing.jobs import run_limited
from noninterlacing.orders import Order
from noninterlacing.pari import from_pari, nfdisc, pari, polredabs
from noninterlacing.poly import poly, real_cyclotomic, to_str
from noninterlacing.realroots import is_totally_real
from noninterlacing.siegel import dim_M


def field_row(key, e, n):
    t = time.time()
    O = Order.maximal(poly(key))
    cert = O.level_certificate(time_limit=200)
    return {"poly": key, "degree": e, "conductor": n, "disc": O.discriminant,
            "m": cert["m"], "empty_levels": cert["empty_levels"],
            "witness_method": cert.get("witness_method"), "ceiling": dim_M(2 * e),
            "seconds": round(time.time() - t, 2)}


def trace_down(key, n):
    w = trace_down_witness(poly(key), real_cyclotomic(n))
    return {"m": w["m_upper"], "empty_levels": [], "witness_method": f"trace-down from Q(zeta_{n})^+",
            "witness": w["nu"]}


def fill(row):
    """Apply the trace-down fallback to an undecided row (in place); True if it decided it."""
    if row["m"] is not None or not psi_interlacing_predicted(row["conductor"]):
        return False
    t = time.time()
    res, err = run_limited(trace_down, (row["poly"], row["conductor"]), mem_bytes=5 << 30, wall=1800)
    if res is None:
        return False
    if row.get("error"):
        row["direct_error"] = row.pop("error")
    row.update(res, seconds=round(time.time() - t, 2))
    return True


def save(rows, out):
    json.dump(rows, open(out + ".tmp", "w"), indent=1)
    os.replace(out + ".tmp", out)


def main(out, nmax, mindeg, maxdeg):
    P = pari()
    rows = json.load(open(out)) if os.path.exists(out) else []
    for r in rows:
        if fill(r):
            print(r["conductor"], r["degree"], r["m"], r["witness_method"], flush=True)
    save(rows, out)
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
                row, err = run_limited(field_row, (key, e, n), mem_bytes=5 << 30, wall=1800)
                if row is None:
                    row = {"poly": key, "degree": e, "conductor": n, "disc": int(nfdisc(g)), "m": None,
                           "empty_levels": [], "witness_method": None, "ceiling": dim_M(2 * e),
                           "seconds": None, "error": err}
                fill(row)
                rows.append(row)
                print(n, e, row["m"], row["witness_method"], row["seconds"], key[:60], flush=True)
                save(rows, out)


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]))
