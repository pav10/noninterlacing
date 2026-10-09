import pytest

try:
    import cypari2  # noqa: F401

    HAVE_PARI = True
except ImportError:  # pragma: no cover
    HAVE_PARI = False

needs_pari = pytest.mark.skipif(not HAVE_PARI, reason="cypari2 (PARI/GP) not installed")
