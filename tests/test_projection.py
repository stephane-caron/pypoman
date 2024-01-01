#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright 2018 Stéphane Caron

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

        vertices_bretl = project_polytope(proj, ineq, eq, method="bretl")
        vertices_cdd = project_polytope(proj, ineq, eq, method="cdd")
        self.assertGreater(len(vertices_bretl), 3)
        self.assertGreater(len(vertices_cdd), 3)
        self.assertLess(len(vertices_bretl), 20)
        self.assertLess(len(vertices_cdd), 1000)

    def test_project_point_to_polytope(self):
        vertices = [
            (np.cos(theta), np.sin(theta))
            for theta in np.arange(0, 2 * np.pi, np.pi / 6)
        ]
        A, b = compute_polytope_halfspaces(vertices)
        point = np.array([2.1, 1.9])
        proj = project_point_to_polytope(point, (A, b), qpsolver="cvxopt")
        self.assertEqual(proj.shape, (2,))
