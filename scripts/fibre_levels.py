"""Level-j elements of Q(zeta_n)^+ by fibres of the relative trace along a tower.

Usage: python scripts/fibre_levels.py N J M1,M2,... [OUT.json]

An entry P:E stands for the subfield of degree E of F_P (P prime).

The tower is Q subset F_{M1} subset F_{M1} F_{M2} subset ... subset F_n, where the
k-th field is the compositum of F_{M1}, ..., F_{Mk} (F_m = Q(zeta_m)^+, m | N),
followed by F_N itself.  Example: 231 2 11,7,77 uses
Q subset F_11 subset F_7 F_11 subset F_77 subset F_231.
Every stage is complete (see noninterlacing/fibers.py); the counts and one
element per Galois orbit are written to OUT.json at the end.  Each finished
stage is also saved to OUT.json.partial, from which a rerun resumes.
"""
import json
import math
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


def _chain(K, ms):
    chain, fix = [K.G], set(K.G)
    for m in ms:
        fix &= set(K.tower_fixer(str(m)))
        chain.append(sorted(fix))
    if chain[-1] != [1]:
        chain.append([1])
    return chain


def main(n, j, ms, out):
    """The stages below ``F_L`` (``L = lcm(ms)``) are computed inside ``F_L``, which is much
    cheaper than inside ``F_n``; only the last step works in ``F_n``."""
    L = math.lcm(*(int(str(m).split(":")[0]) for m in ms))
    K = CyclotomicField(n)
    small = L != n and CyclotomicField(L)
    K0 = small or K
    chain = _chain(K0, ms)
    t0 = time.time()
    log = []

    def note(msg):
        print(msg, flush=True)
        log.append(msg)

    part = (out or f"fibre_{n}_{j}.json") + ".partial"
    start = _load_partial(part) if os.path.exists(part) else None
    if start:
        note(f"resuming after degree {start[-1].degree}")
    stages = tower_levels(K0, chain, j, log=note, start=start, on_stage=lambda st: _save_partial(part, st))
    if small:  # last step F_L -> F_n, in the coordinates of F_n
        H = K.fixer(L)
        reps = [K.embed(x, L) for x in stages[-1].reps]
        t1, found, cand = time.time(), {}, 0
        for nu in reps:
            fr = K.fibre([1], H, nu)
            cand += fr.candidates
            for mu in fr.points:
                found.setdefault(K.orbit_key(mu, [1]), mu)
        top = list(found.values())
        st = Stage([1], K.d, top, sum(K.orbit_size(mu, [1]) for mu in top), cand, round(time.time() - t1, 2))
        stages = stages + [st]
        note(f"degree {st.degree}: {st.count} elements of trace {j} in {len(top)} orbits "
             f"({cand} candidates, {st.seconds} s)")
    amb = [L if small and i < len(stages) - 1 else n for i in range(len(stages))]
    rows = [{"degree": s.degree, "ambient": a, "subgroup": s.H, "count": s.count, "orbits": len(s.reps),
             "candidates": s.candidates, "seconds": s.seconds,
             "reps": [str(x) for x in s.reps] if s.degree < K.d or s.count else []}
            for s, a in zip(stages, amb)]
    res = {"n": n, "level": j, "tower": ms, "stages": rows, "seconds": round(time.time() - t0, 1)}
    if out:
        json.dump(res, open(out, "w"), indent=1)
    if os.path.exists(part):
        os.remove(part)
    return res


if __name__ == "__main__":
    out = sys.argv[4] if len(sys.argv) > 4 else None
    main(int(sys.argv[1]), int(sys.argv[2]), [x if ":" in x else int(x) for x in sys.argv[3].split(",")], out)
