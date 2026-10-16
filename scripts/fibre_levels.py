"""Level-j elements of Q(zeta_n)^+ by fibres of the relative trace along a tower.

Usage: python scripts/fibre_levels.py N J M1,M2,... [OUT.json]

The tower is Q subset F_{M1} subset F_{M1} F_{M2} subset ... subset F_n, where the
k-th field is the compositum of F_{M1}, ..., F_{Mk} (F_m = Q(zeta_m)^+, m | N),
followed by F_N itself.  Example: 231 2 11,7,77 uses
Q subset F_11 subset F_7 F_11 subset F_77 subset F_231.
Every stage is complete (see noninterlacing/fibers.py); the counts and one
element per Galois orbit are written to OUT.json at the end.  Each finished
stage is also saved to OUT.json.partial, from which a rerun resumes.
"""
import json
import os
import sys
import time

from flint import fmpq, fmpq_poly

from noninterlacing.fibers import CyclotomicField, Stage, tower_levels


def _enc(x):
    return [f"{c.p}/{c.q}" for c in (fmpq(x[i]) for i in range(x.degree() + 1))]


def _dec(cs):
    return fmpq_poly([fmpq(*map(int, c.split("/"))) for c in cs])


def _save_partial(path, stages):
    rows = [{"H": s.H, "degree": s.degree, "count": s.count, "candidates": s.candidates,
             "seconds": s.seconds, "reps": [_enc(x) for x in s.reps]} for s in stages]
    json.dump(rows, open(path + ".tmp", "w"))
    os.replace(path + ".tmp", path)


def _load_partial(path):
    rows = json.load(open(path))
    return [Stage(r["H"], r["degree"], [_dec(c) for c in r["reps"]], r["count"], r["candidates"], r["seconds"])
            for r in rows]


def main(n, j, ms, out):
    K = CyclotomicField(n)
    chain, fix = [K.G], set(K.G)
    for m in ms:
        fix &= set(K.fixer(m))
        chain.append(sorted(fix))
    if chain[-1] != [1]:
        chain.append([1])
    t0 = time.time()
    log = []

    def note(msg):
        print(msg, flush=True)
        log.append(msg)

    part = (out or f"fibre_{n}_{j}.json") + ".partial"
    start = _load_partial(part) if os.path.exists(part) else None
    if start:
        note(f"resuming after degree {start[-1].degree}")
    stages = tower_levels(K, chain, j, log=note, start=start, on_stage=lambda st: _save_partial(part, st))
    rows = [{"degree": s.degree, "subgroup": s.H, "count": s.count, "orbits": len(s.reps),
             "candidates": s.candidates, "seconds": s.seconds,
             "reps": [str(x) for x in s.reps] if s.degree < K.d or s.count else []}
            for s in stages]
    res = {"n": n, "level": j, "tower": ms, "stages": rows, "seconds": round(time.time() - t0, 1)}
    if out:
        json.dump(res, open(out, "w"), indent=1)
    if os.path.exists(part):
        os.remove(part)
    return res


if __name__ == "__main__":
    out = sys.argv[4] if len(sys.argv) > 4 else None
    main(int(sys.argv[1]), int(sys.argv[2]), [int(x) for x in sys.argv[3].split(",")], out)
