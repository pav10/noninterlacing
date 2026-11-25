"""Test periodic reciprocal towers for finite critical orbits, using integer pairs."""

from itertools import product
from functools import lru_cache


def matrices(bound):
    """Represent the integral unimodular matrices up to overall sign."""
    for matrix in product(range(-bound, bound + 1), repeat=4):
        a, b, c, d = matrix
        if abs(a*d - b*c) == 1 and next(v for v in matrix if v) > 0:
            yield matrix


def image(matrix, point):
    a, b, c, d = matrix
    x, y = point
    numerator, denominator = x*x + y*y, x*y
    x, y = a*numerator + b*denominator, c*numerator + d*denominator
    # Both maps preserve primitivity, including at infinity.
    if y < 0 or (y == 0 and x < 0):
        return -x, -y
    return x, y


def postcritically_finite(word):
    escape = max(max(abs(d)+abs(b), abs(c)+abs(a)) for a, b, c, d in word)
    for j in range(len(word)):
        for point in ((1, 1), (-1, 1)):
            for matrix in reversed(word[:j+1]):
                point = image(matrix, point)
                if max(map(abs, point)) > escape:
                    return False
            seen = set()
            while point not in seen:
                seen.add(point)
                for matrix in reversed(word):
                    point = image(matrix, point)
                    if max(map(abs, point)) > escape:
                        return False
    return True


def safe_frontiers(bound):
    """Read maps from inner to outer, tracking all critical values at once."""
    alphabet = list(matrices(bound))
    cutoff = max(max(abs(d)+abs(b), abs(c)+abs(a)) for a, b, c, d in alphabet)
    critical = {(1, 1), (-1, 1)}

    @lru_cache(None)
    def successors(point):
        choices = {}
        for i, matrix in enumerate(alphabet):
            result = image(matrix, point)
            if max(map(abs, result)) <= cutoff:
                choices[i] = result
        return choices

    frontier = frozenset([frozenset()])
    visited = set()
    counts = []
    while frontier and frontier not in visited:
        visited.add(frontier)
        next_frontier = set()
        for state in frontier:
            choices = sorted((successors(p) for p in state | critical), key=len)
            for i in choices[0]:
                if all(i in options for options in choices[1:]):
                    next_frontier.add(frozenset(options[i] for options in choices))
        frontier = frozenset(next_frontier)
        counts.append(len(frontier))
    return cutoff, counts, bool(frontier)


if __name__ == "__main__":
    for bound, length in ((16, 1), (3, 2), (1, 3), (1, 4)):
        alphabet = list(matrices(bound))
        count = sum(postcritically_finite(w) for w in product(alphabet, repeat=length))
        print(f"period {length}, matrix entries <= {bound} in absolute value, "
              f"{len(alphabet)**length} words, {count} postcritically finite")
    for bound in (1, 2, 3, 4, 8, 16):
        cutoff, counts, cycle = safe_frontiers(bound)
        print(f"all words, entry bound {bound}, escape above {cutoff}, "
              f"successive frontier sizes {counts}, safe cycle {cycle}")
