#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright 2018 Stéphane Caron

"""General polyhedron-related functions."""

from __future__ import division

from warnings import warn

import numpy as np

from .lp import solve_lp
from .misc import norm

try:
    import cdd
except ImportError:
    warn("Could not import cdd, some functions will not be available")
    cdd = None


try:
    import cvxopt
except ImportError:
    warn("Could not import CVXOPT, some functions will not be available")
    cvxopt = None


def compute_chebyshev_center(A: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Compute the Chebyshev center of a polyhedron.

    The Chebyshev center is the point furthest away from all inequalities.

    Parameters
    ----------
    A :
        Matrix of halfspace representation.
    b :
        Vector of halfspace representation.

    Returns
    -------
    :
        Point further away from all inequalities.

    Notes
    -----
    The Chebyshev center is discussed in [Boyd04]_, Section 4.3.1, p. 148.
    """
    cost = np.zeros(A.shape[1] + 1)
    cost[-1] = -1.0
    a_cheby = np.array([norm(A[i, :]) for i in range(A.shape[0])])
    A_cheby = np.hstack([A, a_cheby.reshape((A.shape[0], 1))])
    z = solve_lp(cost, A_cheby, b)
    if z[-1] < -1e-1:  # last coordinate is distance to boundaries
        raise ValueError("Polytope is empty (margin violation %.2f)" % z[-1])
    return z[:-1]


__all__ = [
    "cdd",
    "compute_chebyshev_center",
    "cvxopt",
]
