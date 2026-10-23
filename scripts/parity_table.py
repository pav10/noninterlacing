"""Level-2 elements of F_pq by trace type and corner parity tau (table in notes/parity.tex).

usage: python scripts/parity_table.py 3,5 3,7 3,11 5,7 3,13 [--json OUT]
"""

import json
import sys
import time

from noninterlacing.parity import level_two_census


def main(argv):
    out = None
    if "--json" in argv:
        i = argv.index("--json")
        out, argv = argv[i + 1], argv[:i] + argv[i + 2:]
    rows = []
    print("| field | degree | level-2 | products | tau = 0 | tau = 1 | other | seconds |")
    print("|---|---|---|---|---|---|---|---|")
    for spec in argv:
        p, q = map(int, spec.split(","))
        t0 = time.time()
        r = level_two_census(p, q)
        r["seconds"] = round(time.time() - t0, 1)
        rows.append(r)
        print(f"| {r['field']} | {r['degree']} | {r['level2']} | {r['product']} | {r['tau=0']} "
              f"| {r['tau=1']} | {r['other']} | {r['seconds']} |", flush=True)
    if out:
        with open(out, "w") as fh:
            json.dump(rows, fh, indent=1)


if __name__ == "__main__":
    main(sys.argv[1:])
