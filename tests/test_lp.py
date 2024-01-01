#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright 2018 Stéphane Caron

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
