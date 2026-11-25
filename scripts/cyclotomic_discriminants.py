"""Compare real cyclotomic discriminants with d**d.

Usage: python scripts/cyclotomic_discriminants.py 105 1155 15015 45045
"""

import argparse
import math

from noninterlacing.cyclotomic import (
    real_cyclotomic_degree,
    real_cyclotomic_discriminant_factors,
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("conductors", nargs="+", type=int)
    args = parser.parse_args()
    if any(n < 1 for n in args.conductors):
        parser.error("conductors must be positive")

    print("| conductor | degree | root discriminant | root disc / degree | disc < d^d |")
    print("|---:|---:|---:|---:|:---:|")
    for n in args.conductors:
        factors = real_cyclotomic_discriminant_factors(n)
        d = real_cyclotomic_degree(n)
        log_rd = math.fsum((e / d) * math.log(p) for p, e in factors.items())
        rd = math.exp(log_rd)
        disc = math.prod(p**e for p, e in factors.items())
        print(f"| {n} | {d} | {rd:.8g} | {rd / d:.8g} | {'yes' if disc < d**d else 'no'} |")


if __name__ == "__main__":
    main()
