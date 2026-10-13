"""Bridge to the data of the ``sails`` project (https://github.com/pav10/sails).

``sails`` computes the sail (boundary of the convex hull of the totally
positive integers) of totally real orders.  Its facets ``F`` have normals
``delta_F in O^vee`` (totally positive, primitive) and integer distances
``d_F = min { Tr(delta_F x) : x in O^+ }``.  For a facet through the vertex
``1`` one has ``d_F = Tr(delta_F)``, so

    m(O) <= min { d_F : F a facet of the totally positive sail, 1 in F },

with equality ``= 1`` as soon as one such facet has distance 1 (then
``delta_F`` is a totally positive codifferent element of trace one, i.e. an
interlacer of the defining polynomial when ``O = Z[alpha]``).  The converse
fails: a trace-one ``delta`` need not be a facet normal (it can support the
sail along a lower-dimensional face through 1).

This module reads ``sails`` summary lines and JSON records and compares
these bounds with :mod:`noninterlacing.orders`.  It never modifies the data.
"""

from __future__ import annotations

import ast
import json
import lzma
from typing import Iterator, Optional

from flint import fmpq_poly, fmpz_poly

from .poly import _Parser


def _open(path: str):
    return lzma.open(path, "rt") if path.endswith(".xz") else open(path)


def parse_element(text: str, var: str = "a") -> fmpq_poly:
    """Parse an element printed by sails/Sage (``-2*a^5 + 15*a^3 + ...``)."""
    return _Parser(text, var).parse()


def read_summary(path: str) -> Iterator[dict]:
    """Rows of a ``summary.txt.*`` file (header line skipped)."""
    with _open(path) as fh:
        for line in fh:
            if line.startswith("label|"):
                continue
            p = line.rstrip("\n").split("|")
            if len(p) < 4:
                continue
            coeffs = [int(c) for c in ast.literal_eval(p[2])]
            basis_txt = p[3].strip()[1:-1]
            yield {
                "label": p[0],
                "disc": int(p[1]),
                "coeffs": coeffs,  # ascending
                "basis": [parse_element(t) for t in _split_top(basis_txt)],
                "min_norm_indecomposables": p[13] if len(p) > 13 else None,
            }


def _split_top(s: str) -> list:
    out, depth, cur = [], 0, ""
    for ch in s:
        if ch in "([":
            depth += 1
        elif ch in ")]":
            depth -= 1
        if ch == "," and depth == 0:
            out.append(cur.strip())
            cur = ""
        else:
            cur += ch
    if cur.strip():
        out.append(cur.strip())
    return out


def read_records(path: str) -> Iterator[dict]:
    with _open(path) as fh:
        for line in fh:
            line = line.strip()
            if line:
                yield json.loads(line)


def facet_bound_at_one(record: dict) -> Optional[int]:
    """``min d_F`` over facets of the totally positive orthant containing ``1``."""
    basis = [parse_element(t) for t in record["integral_basis"]]
    d = record["degree"]
    one = None
    for i, b in enumerate(basis):
        if b == fmpq_poly([1]):
            one = [int(i == j) for j in range(d)]
    if one is None:
        return None
    orth = next(o for o in record["orthants"] if not any(o["signature"]))
    ds = [f["d"] for f in orth["facets"] if any(v == one for v in f["vertices"])]
    return min(ds) if ds else None


def order_from_summary(row: dict):
    from .orders import Order

    f = fmpz_poly(row["coeffs"])
    return Order(f, row["basis"], name=row["label"])
