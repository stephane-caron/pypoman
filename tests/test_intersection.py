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

"""Unit tests for intersection functions."""

import unittest

import numpy as np

from pypoman import intersect_line_polygon


class TestIntersection(unittest.TestCase):
    """Test fixture for intersection algorithms."""

    def test_intersect_line_polygon(self):
        vertices = [
            (np.cos(theta), np.sin(theta))
            for theta in np.arange(0, 2 * np.pi, np.pi / 6)
        ]
        line = (np.array([2.1, 1.9]), np.array([0.0, 0.0]))
        inter = intersect_line_polygon(line, vertices, apply_hull=False)
        self.assertEqual(len(inter), 1)

    def test_intersect_line_polygon_with_hull(self):
        vertices = [
            (np.cos(theta), np.sin(theta))
            for theta in np.arange(0, 2 * np.pi, np.pi / 6)
        ]
        line = (np.array([2.1, 1.9]), np.array([0.0, 0.0]))
        inter = intersect_line_polygon(line, vertices, apply_hull=True)
        self.assertEqual(len(inter), 1)
