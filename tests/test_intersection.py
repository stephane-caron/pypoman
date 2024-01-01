#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright 2018 Stéphane Caron

"""Unit tests for intersection functions."""

import unittest

import numpy as np

from pypoman import (
    intersect_line_cylinder,
    intersect_line_polygon,
    intersect_polygons,
)


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

    def test_intersect_line_cylinder(self):
        vertices = [
            (np.cos(theta), np.sin(theta))
            for theta in np.arange(0, 2 * np.pi, np.pi / 6)
        ]
        line = (np.array([2.1, 1.9, -1.1]), np.array([0.0, 0.0, 0.0]))
        inter = intersect_line_cylinder(line, vertices)
        self.assertEqual(len(inter), 1)

    def test_intersect_polygons(self):
        polygon1 = [
            np.array([0.0, 0.0]),
            np.array([1.0, 0.0]),
            np.array([1.0, 1.0]),
            np.array([0.0, 1.0]),
        ]
        polygon2 = [
            np.array([0.5, 0.5]),
            np.array([1.5, 0.5]),
            np.array([1.5, 1.5]),
            np.array([0.5, 1.5]),
        ]
        inter = intersect_polygons(polygon1, polygon2)
        self.assertEqual(len(inter), 4)

    def test_empty_polygon_intersection(self):
        polygon1 = [
            np.array([0.0, 0.0]),
            np.array([1.0, 0.0]),
            np.array([1.0, 1.0]),
            np.array([0.0, 1.0]),
        ]
        polygon2 = [
            np.array([10.5, 10.5]),
            np.array([11.5, 10.5]),
            np.array([11.5, 11.5]),
            np.array([10.5, 11.5]),
        ]
        inter = intersect_polygons(polygon1, polygon2)
        self.assertEqual(len(inter), 0)
