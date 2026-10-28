# Make sure Python treats the test directory as a package.

import importlib.util
import unittest


def requires(*modules: str):
    """Skip a test if one of its optional dependencies is not installed."""
    missing = [m for m in modules if importlib.util.find_spec(m) is None]
    return unittest.skipIf(
        bool(missing), f"optional dependency not installed: {missing}"
    )
