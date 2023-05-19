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

"""Unit tests for linear programming."""

import unittest

import numpy as np

from pypoman.lp import solve_lp


class TestLP(unittest.TestCase):
    """Test fixture for linear programming."""

    def test_solve_lp(self):
        c = np.array([-1.0, 0.0])
        G = np.array([[1.0, 0.0], [0.0, 1.0]])
        h = np.zeros(2)
        x = solve_lp(c, G, h)
        self.assertEqual(x.shape, (2,))

    def test_unfeasible_lp(self):
        c = np.array([1.0, 0.0])
        G = np.array([[1.0, 0.0], [0.0, 1.0]])
        h = np.zeros(2)
        with self.assertRaises(ValueError):
            solve_lp(c, G, h)
