"""Interlacing polynomials of every level and non-interlacing families.

See :mod:`noninterlacing.interlace` for definitions and the codifferent
dictionary, :mod:`noninterlacing.simplex` for the enumeration algorithm.
"""

from .poly import poly, to_str, real_cyclotomic, chebyshev_C, compose
from .realroots import interlaces, is_totally_real, RootData
from .interlace import (Interlacing, count_interlacers, find_interlacer, interlacers,
                        is_interlacing, level_semigroup, min_level, siegel_ceiling)
from .orders import Order
from .semigroup import NumericalSemigroup

__version__ = "0.1.0"

__all__ = [
    "poly", "to_str", "real_cyclotomic", "chebyshev_C", "compose",
    "interlaces", "is_totally_real", "RootData",
    "Interlacing", "count_interlacers", "find_interlacer", "interlacers",
    "is_interlacing", "level_semigroup", "min_level", "siegel_ceiling",
    "Order", "NumericalSemigroup",
]
