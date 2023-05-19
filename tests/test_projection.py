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

"""Unit tests for polytope projection."""

import unittest

import numpy as np

from pypoman import (
    compute_polytope_halfspaces,
    project_point_to_polytope,
    project_polytope,
)


class TestProjection(unittest.TestCase):
    """Test fixture for projection algorithms."""

    def test_project_polytope(self, n: int = 10, p: int = 2):
        r"""Test polytope projection.

        The original polytope is defined by:

        - Inequality constraints: \forall i, |x_i| <= 1
        - Equality constraint: sum_i x_i = 0

        Parameters
        ----------
        n :
            Dimension of the original polytope
        p :
            Dimension of the projected polytope
        """
        A = np.vstack([+np.eye(n), -np.eye(n)])
        b = np.ones(2 * n)
        C = np.ones(n).reshape((1, n))
        d = np.array([0])
        ineq = (A, b)  # A * x <= b
        eq = (C, d)  # C * x == d

        # Projection is proj(x) = [x_0 x_1]
        E = np.zeros((p, n))
        E[0, 0] = 1.0
        E[1, 1] = 1.0
        f = np.zeros(p)
        proj = (E, f)  # proj(x) = E * x + f

        vertices = project_polytope(proj, ineq, eq, method="bretl")
        self.assertGreater(len(vertices), 3)
        self.assertLess(len(vertices), 10)

    def test_project_point_to_polytope(self):
        vertices = [
            (np.cos(theta), np.sin(theta))
            for theta in np.arange(0, 2 * np.pi, np.pi / 6)
        ]
        A, b = compute_polytope_halfspaces(vertices)
        point = np.array([2.1, 1.9])
        proj = project_point_to_polytope(point, (A, b), qpsolver="cvxopt")
        self.assertEqual(proj.shape, (2,))
