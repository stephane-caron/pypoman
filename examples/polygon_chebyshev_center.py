#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright 2025 Inria

import numpy as np

from pypoman import compute_chebyshev_center

A = np.array([[1.0, 0.0], [-1.0, 0.0], [0.0, 1.0], [0.0, -1]])
b = np.array([1.0, 0.0, 1.0, 0.0])

if __name__ == "__main__":
    x = compute_chebyshev_center(A, b)
    print("The polygon is a square: 0 <= x <= 1 and 0 <= y <= 1")
    print(f"Chebyshev center of the polygon: {x=}")
