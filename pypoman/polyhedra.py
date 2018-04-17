#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2018 Stephane Caron <stephane.caron@lirmm.fr>
#
# This file is part of pypoman <https://github.com/stephane-caron/pypoman>.
#
# pypoman is free software: you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# pypoman is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with pypoman. If not, see <http://www.gnu.org/licenses/>.

from __future__ import division

from numpy import array, hstack, zeros

from .lp import solve_lp
from .misc import norm, warn


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


def compute_chebyshev_center(A, b):
    """
    Compute the Chebyshev center of a polyhedron, that is, the point furthest
    away from all inequalities.

    Parameters
    ----------
    A : array, shape=(m, k)
        Matrix of halfspace representation.
    b : array, shape=(m,)
        Vector of halfspace representation.

    Returns
    -------
    z : array, shape=(k,)
        Point further away from all inequalities.

    Notes
    -----
    The Chebyshev center is discussed in [Boyd04]_, Section 4.3.1, p. 148.
    """
    cost = zeros(A.shape[1] + 1)
    cost[-1] = -1.
    a_cheby = array([norm(A[i, :]) for i in range(A.shape[0])])
    A_cheby = hstack([A, a_cheby.reshape((A.shape[0], 1))])
    z = solve_lp(cost, A_cheby, b)
    if z[-1] < -1e-1:  # last coordinate is distance to boundaries
        raise Exception("Polytope is empty (margin violation %.2f)" % z[-1])
    return z[:-1]
