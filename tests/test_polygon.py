#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2018-2020 Stephane Caron <stephane.caron@normalesup.org>
#
# This file is part of pypoman <https://github.com/stephane-caron/pypoman>.
#
# pypoman is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# pypoman is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# pypoman. If not, see <http://www.gnu.org/licenses/>.

"""Unit tests for polygon submodule."""

import unittest

import numpy as np

from pypoman import compute_polygon_hull, plot_polygon


class TestPolygon(unittest.TestCase):
    """Test fixture for polygon submodule."""

    def test_compute_polygon_hull(self):
        B = np.array([[1.0, 0.0], [0.0, 1.0], [-1.0, 0.0], [0.0, -2.0]])
        c = np.ones(4)
        x = compute_polygon_hull(B, c)
        self.assertGreater(len(x), 3)
        self.assertLess(len(x), 30)

    def test_compute_polygon_hull_offset(self):
        B = np.array([[1.0, 0.0], [0.0, 1.0], [-1.0, 0.0], [0.0, -2.0]])
        c = np.array([1.0, 1.0, -0.1, 2.0])
        x = compute_polygon_hull(B, c)
        self.assertGreater(len(x), 3)
        self.assertLess(len(x), 30)

    def test_compute_polygon_hull_empty(self):
        B = np.array([[1.0, 0.0], [0.0, 1.0], [-1.0, 0.0], [0.0, -2.0]])
        c = np.array([1.0, 1.0, -1.0, 2.0])
        with self.assertRaises(ValueError):
            compute_polygon_hull(B, c)

    def test_plot_polygon(self):
        vertices = [
            np.array([0.5, 0.5]),
            np.array([1.5, 0.5]),
            np.array([1.5, 1.5]),
            np.array([0.5, 1.5]),
        ]
        plot_polygon(vertices)
