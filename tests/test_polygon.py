#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright 2018 Stéphane Caron

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
        plot_polygon(vertices, resize=True)
