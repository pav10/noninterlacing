"""Command line interface: ``python -m noninterlacing COMMAND ...``."""

from __future__ import annotations

import argparse
import json
import sys


def _range(text: str) -> list:
    """``"2..30"`` or ``"2,5,9"`` or ``"7"`` -> list of ints."""
    out = []
    for part in text.split(","):
        if ".." in part:
            a, b = part.split("..")
            out.extend(range(int(a), int(b) + 1))
        elif part:
            out.append(int(part))
    return out


def _emit(obj, as_json: bool):
    if as_json:
        print(json.dumps(obj, indent=1, default=str))
    else:
        for k, v in obj.items():
            if isinstance(v, list) and v and isinstance(v[0], str) and len(v) > 3:
                print(f"{k}:")
                for x in v:
                    print(f"  {x}")
            else:
                print(f"{k}: {v}")


def cmd_interlacers(a):
    from .interlace import Interlacing
    from .poly import to_str

    I = Interlacing(a.poly)
    if a.count:
        res = I.enumerate(a.level, mode="count")
        _emit({"f": to_str(I.f), "level": a.level, "count": res.count, "search_nodes": res.nodes,
               "exact_checks": res.exact_checks, "seconds": round(res.seconds, 3)}, a.json)
        return
    gs = I.interlacers(a.level, limit=a.limit)
    _emit({"f": to_str(I.f), "level": a.level, "number": len(gs),
           "interlacers": [to_str(g) for g in gs]}, a.json)


def cmd_level(a):
    from .interlace import Interlacing, siegel_ceiling
    from .poly import is_irreducible, to_str

    I = Interlacing(a.poly)
    S = I.levels(method=a.method)
    m = S.multiplicity
    out = {"f": to_str(I.f), "degree": I.d, "irreducible": is_irreducible(I.f),
           "span": round(I.roots.span, 6), "min level m(f)": m,
           "interlacing": m == 1, "level semigroup": str(S),
           "non-interlacing levels (gaps)": S.gaps, "Frobenius number": S.frobenius,
           "Siegel ceiling dim M_2d": siegel_ceiling(I.d)}
    for j in range(1, m + 1):
        out[f"margin eps_{j}"] = round(I.margin(j).value, 8)
    g = I.find(m, method=a.method)
    out[f"level-{m} interlacer"] = to_str(g)
    _emit(out, a.json)


def cmd_check(a):
    from .poly import poly, to_str
    from .realroots import interlaces

    f, g = poly(a.poly), poly(a.g)
    ok = interlaces(f, g)
    _emit({"f": to_str(f), "g": to_str(g), "strictly interlaces": ok,
           "level": int(g.leading_coefficient()) if ok else None}, a.json)


def cmd_order(a):
    from .orders import Order
    from .poly import to_str

    O = Order.maximal(a.poly) if a.maximal else Order.monogenic(a.poly)
    S = O.trace_semigroup()
    out = {"order": repr(O), "m(O)": S.multiplicity, "trace semigroup": str(S), "gaps": S.gaps,
           "margin(1)": round(O.margin(1).value, 8)}
    if a.trace:
        els = O.totally_positive_codifferent(a.trace, limit=a.limit)
        out[f"totally positive codifferent elements of trace {a.trace}"] = [to_str(e, "a") for e in els]
    _emit(out, a.json)


def cmd_heredity(a):
    from .heredity import heredity_report

    _emit(heredity_report(a.poly).as_dict(), a.json)


def cmd_cyclotomic(a):
    from .cyclotomic import cyclotomic_level

    rows = []
    for n in _range(a.n):
        r = cyclotomic_level(n).as_dict()
        rows.append(r)
        if not a.json:
            print(f"n={r['n']:4d}  degree={r['degree']:3d}  m={r['m']}  margin={r['margin']:+.6f}  "
                  f"admissible={r['admissible']}  admissible divisors={r['admissible_divisors']}")
    if a.json:
        print(json.dumps(rows, indent=1))


def cmd_family(a):
    from .families import check_main_theorem

    for N in _range(a.N):
        r = check_main_theorem(N, a.base)
        if a.json:
            print(json.dumps(r, default=str))
        else:
            r.pop("f", None)
            print(r)


def cmd_siegel(a):
    from .siegel import coefficient_table, noninterlacing_sign_test, siegel_zeta

    if a.zeta:
        _emit(siegel_zeta(a.zeta), a.json)
        return
    for row in coefficient_table(_range(a.degrees)):
        print(f"d={row['d']:3d}  dim M_2d={row['dim M_2d']}  signs {row['signs']:6s}  b = {', '.join(row['b'])}"
              f"   [{noninterlacing_sign_test(row['d'])}]")


def cmd_census(a):
    from .search import levels_census

    rows = levels_census(a.degree, a.span, irreducible_only=not a.all)
    non = [r for r in rows if r.level > 1]
    if a.json:
        print(json.dumps([r.as_dict() for r in rows], indent=1))
        return
    print(f"{len(rows)} classes (degree {a.degree}, span <= {a.span}); non-interlacing: {len(non)}")
    for r in (rows if a.verbose else non):
        print(r.as_dict())


def cmd_random(a):
    from .search import random_search

    _emit(random_search(a.degree, a.samples, seed=a.seed, cluster=a.cluster), a.json)


def cmd_paper(a):
    from .paper import run, write_report

    res = run(full=a.full, only=a.only.split(",") if a.only else None)
    path = write_report(res, a.out)
    print(f"report written to {path}")


def main(argv=None):
    p = argparse.ArgumentParser(prog="noninterlacing", description=__doc__)
    p.add_argument("--json", action="store_true", help="JSON output")
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("interlacers", help="list or count the level-j interlacers of f")
    s.add_argument("poly")
    s.add_argument("--level", "-j", type=int, default=1)
    s.add_argument("--count", action="store_true")
    s.add_argument("--limit", type=int)
    s.set_defaults(fn=cmd_interlacers)

    s = sub.add_parser("level", help="interlacing level m(f), semigroup of levels, margins")
    s.add_argument("poly")
    s.add_argument("--method", choices=["auto", "dfs", "milp"], default="auto")
    s.set_defaults(fn=cmd_level)

    s = sub.add_parser("check", help="exact test: does g strictly interlace f?")
    s.add_argument("poly")
    s.add_argument("g")
    s.set_defaults(fn=cmd_check)

    s = sub.add_parser("order", help="m(O) and trace semigroup for Z[a] or O_K")
    s.add_argument("poly")
    s.add_argument("--maximal", action="store_true", help="use the maximal order (PARI)")
    s.add_argument("--trace", type=int)
    s.add_argument("--limit", type=int)
    s.set_defaults(fn=cmd_order)

    s = sub.add_parser("heredity", help="subfield lower bounds m(Z[a]) >= m(Z[a] cap F) (PARI)")
    s.add_argument("poly")
    s.set_defaults(fn=cmd_heredity)

    s = sub.add_parser("cyclotomic", help="levels of real cyclotomic fields F_n")
    s.add_argument("n", help="e.g. 21 or 21,84 or 20..60")
    s.set_defaults(fn=cmd_cyclotomic)

    s = sub.add_parser("family", help="check Theorem 3.1 for f_N = m_beta(x^2 - x - N)")
    s.add_argument("--N", default="2..10")
    s.add_argument("--base", default=None, help="base polynomial (default m_beta for n = 21)")
    s.set_defaults(fn=cmd_family)

    s = sub.add_parser("siegel", help="Siegel coefficients b_k(2d); --zeta f compares with PARI")
    s.add_argument("--degrees", default="2..20")
    s.add_argument("--zeta", default=None)
    s.set_defaults(fn=cmd_siegel)

    s = sub.add_parser("census", help="all totally real polys of degree d and span <= S, with levels")
    s.add_argument("--degree", "-d", type=int, required=True)
    s.add_argument("--span", "-S", type=float, required=True)
    s.add_argument("--all", action="store_true", help="include reducible polynomials")
    s.add_argument("--verbose", "-v", action="store_true")
    s.set_defaults(fn=cmd_census)

    s = sub.add_parser("random", help="random search for non-interlacing polynomials")
    s.add_argument("--degree", "-d", type=int, required=True)
    s.add_argument("--samples", "-n", type=int, default=100)
    s.add_argument("--seed", type=int, default=0)
    s.add_argument("--cluster", type=float, default=None)
    s.set_defaults(fn=cmd_random)

    s = sub.add_parser("paper", help="reproduce all computational claims of the paper")
    s.add_argument("--full", action="store_true")
    s.add_argument("--only", default=None)
    s.add_argument("--out", default="results")
    s.set_defaults(fn=cmd_paper)

    a = p.parse_args(argv)
    if not hasattr(a, "json"):
        a.json = False
    a.fn(a)


if __name__ == "__main__":  # pragma: no cover
    main()
