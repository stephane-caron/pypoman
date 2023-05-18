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

"""Unit tests for halfspace enumeration."""

import unittest

import numpy as np

from pypoman import compute_polytope_halfspaces


class TestHalfspaceEnumeration(unittest.TestCase):
    """Test fixture for halfspace enumeration problems."""

    def test_halfspace_enumeration(self):
        vertices = map(
            np.array,
            [[1, 0, 0], [0, 1, 0], [1, 1, 0], [0, 0, 1], [0, 1, 1]],
        )
        A, b = compute_polytope_halfspaces(vertices)
        self.assertEqual(A.shape[0], b.shape[0])
        self.assertGreater(len(b), 4)
        self.assertLess(len(b), 10)
