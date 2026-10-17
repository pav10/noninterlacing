"""Level-j elements of the subfields of index r (r prime) of Q(zeta_n)^+, by fibre towers.

Usage: python scripts/subfield_levels.py N R J [OUT.json]

A level-1 element of a subfield F of index r lifts to a level-r element of F_n,
so an interlacing subfield of index r gives m(Psi_n) <= r.  Each subfield K^H is
reached by the tower K^H cap F_{q_1}, K^H cap F_{q_1 q_2}, ... (fibers.subfield_chain);
every stage is complete.  For each subgroup H the output lists how H acts on each
prime-power factor (the residues of its elements), the tower degrees and the count.
"""
import json
import sys
import time

from noninterlacing.cyclotomic import factorint
from noninterlacing.fibers import (CyclotomicField, interlacer_from_element, subfield_chain,
                                   subgroups_of_prime_order, tower_levels)


def main(n, r, j, out):
    K = CyclotomicField(n)
    rows = []
    for H in subgroups_of_prime_order(K, r):
        qs = [p ** a for p, a in factorint(n).items()]
        acts = {q: sorted({min(h % q, q - h % q) for h in H}) for q in qs}
        chain = subfield_chain(K, H)
        t = time.time()
        st = tower_levels(K, chain, j)
        row = {"H": H, "action": {str(q): v for q, v in acts.items()},
               "tower_degrees": [K.d // len(c) for c in chain], "counts": [s.count for s in st],
               "level_j_elements": st[-1].count, "seconds": round(time.time() - t, 1)}
        if j == 1 and st[-1].count:
            g = interlacer_from_element(K, st[-1].reps[0])  # lifted: level r in F_n
            row["lifted_interlacer"] = str(g)
            row["lifted_level"] = int(g.leading_coefficient())
        print(json.dumps({k: v for k, v in row.items() if k != "lifted_interlacer"}), flush=True)
        rows.append(row)
    if out:
        json.dump({"n": n, "index": r, "level": j, "subfields": rows}, open(out, "w"), indent=1)


if __name__ == "__main__":
    main(int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4] if len(sys.argv) > 4 else None)
