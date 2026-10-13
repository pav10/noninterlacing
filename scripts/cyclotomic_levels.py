"""m(O_{F_n}) for real cyclotomic fields F_n = Q(zeta_n)^+ (one conductor per field).

Usage: python scripts/cyclotomic_levels.py OUT.json NMAX [MAXDEG]

Resumes from OUT.json.  Each field runs in a child process with a 5 GB memory
limit and a wall-clock limit.  If that fails and Psi_n is non-interlacing by
McKee--Smyth, a level-2 witness is built from a coprime split n = n1 n2 with
Psi_{n1}, Psi_{n2} interlacing (cyclotomic.compositum_witness), giving m = 2;
otherwise the row is recorded with m = None and the reason.
"""
import json
import os
import sys
import time

from noninterlacing.cyclotomic import (admissible_divisors, compositum_splits, compositum_witness,
                                       is_admissible, psi_interlacing_predicted, real_cyclotomic_degree)
from noninterlacing.interlace import Interlacing
from noninterlacing.jobs import run_limited
from noninterlacing.poly import real_cyclotomic


def _empty(n, e, reason):
    return {"n": n, "degree": e, "m": None, "empty_levels": [], "witness_method": None,
            "witness": None, "margin1": None, "error": reason}


def _field_row(n, e):
    t = time.time()
    I = Interlacing(real_cyclotomic(n))
    cert = I.level_certificate(time_limit=200)
    mm = I.margin(1, time_limit=60)
    return {"n": n, "degree": e, "m": cert["m"], "empty_levels": cert["empty_levels"],
            "witness_method": cert.get("witness_method"), "witness": cert["witness"],
            "margin1": mm.value if mm.status == "optimal" else None,
            "seconds": round(time.time() - t, 2)}


def compositum_row(n, e, previous=None):
    """m = 2 from McKee--Smyth (level 1 empty) and a compositum witness at level 2."""
    if psi_interlacing_predicted(n) or not compositum_splits(n):
        return None
    t = time.time()
    n1, n2 = min(compositum_splits(n), key=lambda s: s[1])
    g = compositum_witness(n, n1, n2)
    row = {"n": n, "degree": e, "m": 2, "empty_levels": [[1, "McKee-Smyth"]],
           "witness_method": f"compositum Psi_{n1} x Psi_{n2}", "witness": str(g), "margin1": None,
           "seconds": round(time.time() - t, 2)}
    if previous and previous.get("error"):
        row["direct_error"] = previous["error"]
    return row


def compute_row(n, e, wall=1800):
    row, err = run_limited(_field_row, (n, e), mem_bytes=5 << 30, wall=wall)
    return row if row is not None else _empty(n, e, err)


def save(rows, out):
    json.dump(rows, open(out + ".tmp", "w"), indent=1)
    os.replace(out + ".tmp", out)


def main(out, nmax, maxdeg):
    rows = json.load(open(out)) if os.path.exists(out) else []
    for i, r in enumerate(rows):
        if r["m"] is None and (c := compositum_row(r["n"], r["degree"], r)):
            c.update(admissible=r["admissible"], admissible_divisors=r["admissible_divisors"])
            rows[i] = c
            print(r["n"], r["degree"], "m = 2 by", c["witness_method"], flush=True)
    save(rows, out)
    done = {r["n"] for r in rows}
    for n in range(3, nmax + 1):
        if n % 4 == 2 or n in done:
            continue
        e = real_cyclotomic_degree(n)
        if e < 2 or e > maxdeg:
            continue
        row = compute_row(n, e)
        if row["m"] is None:
            row = compositum_row(n, e, row) or row
        row.update(admissible=is_admissible(n), admissible_divisors=admissible_divisors(n))
        print(n, e, row["m"], row["admissible"], row["admissible_divisors"][:3], row["margin1"],
              row.get("witness_method"), row.get("seconds"), row.get("error", ""), flush=True)
        rows.append(row)
        rows.sort(key=lambda r: r["n"])
        save(rows, out)


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]) if len(sys.argv) > 3 else 60)
