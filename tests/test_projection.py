# SPDX-License-Identifier: GPL-3.0-or-later

"""Unit tests for polytope projection."""

import unittest

import numpy as np

from pypoman import (
    compute_polytope_halfspaces,
    project_point_to_polytope,
    project_polytope,
)

from . import requires


class TestProjection(unittest.TestCase):
    """Test fixture for projection algorithms."""

    @staticmethod
    def get_problem(n: int = 10, p: int = 2):
        r"""Get a polytope projection problem.

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
        return proj, ineq, eq

    def test_project_polytope_bretl(self):
        proj, ineq, eq = self.get_problem()
        vertices = project_polytope(proj, ineq, eq, method="bretl")
        self.assertGreater(len(vertices), 3)
        self.assertLess(len(vertices), 20)

    @requires("cdd")
    def test_project_polytope_cdd(self):
        proj, ineq, eq = self.get_problem()
        vertices = project_polytope(proj, ineq, eq, method="cdd")
        self.assertGreater(len(vertices), 3)
        self.assertLess(len(vertices), 1000)

    @requires("cdd", "qpsolvers")
    def test_project_point_to_polytope(self):
        vertices = [
            (np.cos(theta), np.sin(theta))
            for theta in np.arange(0, 2 * np.pi, np.pi / 6)
        ]
        A, b = compute_polytope_halfspaces(vertices)
        point = np.array([2.1, 1.9])
        proj = project_point_to_polytope(point, (A, b), qpsolver="cvxopt")
        self.assertEqual(proj.shape, (2,))
