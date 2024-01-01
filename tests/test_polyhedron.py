#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright 2018 Stéphane Caron

"""Unit tests for polyhedron submodule."""

import unittest

import numpy as np

from pypoman import compute_chebyshev_center


class TestPolyhedron(unittest.TestCase):
    """Test fixture for polyhedron submodule."""

    def test_compute_chebyshev_center(self):
        G = np.array([[1.0, 0.0], [0.0, 1.0], [-1.0, 0.0], [0.0, -2.0]])
        h = np.zeros(4)
        x = compute_chebyshev_center(G, h)
        self.assertEqual(x.shape, (2,))

    def test_chebyshev_center_unbounded(self):
        G = np.array([[1.0, 0.0], [0.0, 1.0]])
        h = np.zeros(2)
        with self.assertRaises(ValueError):
            compute_chebyshev_center(G, h)

    def test_chebyshev_center_empty(self):
        G = np.array([[1.0, 0.0], [-1.0, 0.0]])
        h = np.array([0.0, -1.0])
        with self.assertRaises(ValueError):
            compute_chebyshev_center(G, h)
