"""m(O_{F_n}) for real cyclotomic fields F_n = Q(zeta_n)^+ (one conductor per field).

Usage: python scripts/cyclotomic_levels.py OUT.json NMAX [MAXDEG]

Resumes from OUT.json.  Each field runs in a child process with a 5 GB memory
limit and a wall-clock limit; failures are recorded with m = None and the reason.
"""
import json
import os
import sys
import time

from noninterlacing.cyclotomic import admissible_divisors, is_admissible, real_cyclotomic_degree
from noninterlacing.interlace import Interlacing
from noninterlacing.poly import real_cyclotomic


def _empty(n, e, reason):
    return {"n": n, "degree": e, "m": None, "empty_levels": [], "witness_method": None,
            "witness": None, "margin1": None, "error": reason}


def _field_row(n, e, q):
    import resource

    resource.setrlimit(resource.RLIMIT_AS, (5 << 30, 5 << 30))
    try:
        t = time.time()
        I = Interlacing(real_cyclotomic(n))
        cert = I.level_certificate(time_limit=200)
        mm = I.margin(1, time_limit=60)
        q.put({"n": n, "degree": e, "m": cert["m"], "empty_levels": cert["empty_levels"],
               "witness_method": cert.get("witness_method"), "witness": cert["witness"],
               "margin1": mm.value if mm.status == "optimal" else None,
               "seconds": round(time.time() - t, 2)})
    except MemoryError:
        q.put(_empty(n, e, "memory limit"))


def compute_row(n, e, wall=1800):
    import multiprocessing as mp

    q = mp.Queue()
    p = mp.Process(target=_field_row, args=(n, e, q))
    p.start()
    p.join(wall)
    if p.is_alive():
        p.kill()
        p.join()
        return _empty(n, e, "wall-clock limit")
    if q.empty():
        return _empty(n, e, f"child exit code {p.exitcode}")
    return q.get()


def main(out, nmax, maxdeg):
    rows = json.load(open(out)) if os.path.exists(out) else []
    done = {r["n"] for r in rows}
    for n in range(3, nmax + 1):
        if n % 4 == 2 or n in done:
            continue
        e = real_cyclotomic_degree(n)
        if e < 2 or e > maxdeg:
            continue
        row = compute_row(n, e)
        row.update(admissible=is_admissible(n), admissible_divisors=admissible_divisors(n))
        print(n, e, row["m"], row["admissible"], row["admissible_divisors"][:3], row["margin1"],
              row.get("witness_method"), row.get("seconds"), row.get("error", ""), flush=True)
        rows.append(row)
        rows.sort(key=lambda r: r["n"])
        json.dump(rows, open(out + ".tmp", "w"), indent=1)
        os.replace(out + ".tmp", out)


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]) if len(sys.argv) > 3 else 60)
