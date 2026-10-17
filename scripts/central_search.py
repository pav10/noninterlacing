"""Look for a level-j interlacer of Psi_N in fibres over central elements of F_M (M | N).

Usage: python scripts/central_search.py N J M [SCALE] [OUT.json]

See noninterlacing.fibers.central_search.  A result is a checked interlacer; a
failure proves nothing.  F_M may also be a compositum: M1xM2 means F_M1 F_M2.
"""
import json
import sys
import time

from noninterlacing.fibers import CyclotomicField, central_search


def main(n, j, spec, scale, out):
    K = CyclotomicField(n)
    H = set(K.G)
    for m in spec.split("x"):
        H &= set(K.fixer(int(m)))
    t = time.time()
    g = central_search(K, sorted(H), j, scale, log=lambda s: print(s, flush=True))
    row = {"n": n, "level": j, "subfield": spec, "scale": scale, "interlacer": str(g) if g is not None else None,
           "seconds": round(time.time() - t, 1)}
    print("found" if g is not None else "none found", row["seconds"], "s", flush=True)
    if out:
        json.dump(row, open(out, "w"), indent=1)


if __name__ == "__main__":
    main(int(sys.argv[1]), int(sys.argv[2]), sys.argv[3], float(sys.argv[4]) if len(sys.argv) > 4 else 0.007,
         sys.argv[5] if len(sys.argv) > 5 else None)
