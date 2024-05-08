#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright 2018 Stéphane Caron

"""Unit tests for duality conversions."""

import unittest

import numpy as np

from pypoman import (
    compute_cone_face_matrix,
    compute_polytope_halfspaces,
    compute_polytope_vertices,
    convex_hull,
)


class TestDuality(unittest.TestCase):
    """Test fixture for duality conversions."""

    def test_halfspace_enumeration(self):
        vertices = [
            np.array(a)
            for a in [[1, 0, 0], [0, 1, 0], [1, 1, 0], [0, 0, 1], [0, 1, 1]]
        ]
        A, b = compute_polytope_halfspaces(vertices)
        self.assertEqual(A.shape[0], b.shape[0])
        self.assertGreater(len(b), 4)
        self.assertLess(len(b), 10)

    def test_vertex_enumeration(self):
        A = np.array(
            [
                [-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0],
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1],
                [1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1],
                [1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0],
                [0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0],
                [0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1],
            ]
        )
        b = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 2, 2, 1, 2, 3])
        vertices = compute_polytope_vertices(A, b)
        self.assertGreater(len(vertices), 50)
        self.assertLess(len(vertices), 500)

    def test_compute_cone_face_matrix(self):
        S = np.array([[1.0, 0.0], [0.0, 1.0]])
        F = compute_cone_face_matrix(S)
        self.assertEqual(len(F), 2)

    def test_convex_hull(self):
        polygon = [
            np.array([0.0, 0.0]),
            np.array([1.0, 0.0]),
            np.array([1.0, 1.0]),
            np.array([0.0, 1.0]),
            np.array([0.5, 0.5]),
        ]
        hull = convex_hull(polygon)
        self.assertEqual(len(hull), 4)
