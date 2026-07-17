#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Unit tests for the exceptions submodule."""

import builtins
import unittest

import numpy as np

import pypoman
from pypoman import (
    ConversionError,
    EmptyPolytopeError,
    SolverError,
    compute_chebyshev_center,
)


class TestExceptions(unittest.TestCase):
    """Test fixture for pypoman exceptions."""

    def test_base_class(self):
        """All pypoman exceptions derive from the same base class."""
        self.assertTrue(issubclass(pypoman.Exception, builtins.Exception))
        for cls in (ConversionError, EmptyPolytopeError, SolverError):
            self.assertTrue(issubclass(cls, pypoman.Exception))

    def test_flat_hierarchy(self):
        """Concrete exceptions are plain pypoman errors, not ``ValueError``."""
        for cls in (ConversionError, EmptyPolytopeError, SolverError):
            self.assertFalse(issubclass(cls, builtins.ValueError))

    def test_caught_as_pypoman_exception(self):
        """An empty polytope raises a catchable pypoman exception."""
        G = np.array([[1.0, 0.0], [-1.0, 0.0]])
        h = np.array([0.0, -1.0])
        with self.assertRaises(EmptyPolytopeError):
            compute_chebyshev_center(G, h)
        with self.assertRaises(pypoman.Exception):
            compute_chebyshev_center(G, h)


if __name__ == "__main__":
    unittest.main()
