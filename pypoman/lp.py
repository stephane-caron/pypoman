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

"""Functions for linear programming."""

from typing import Optional, Union

import cvxopt
import cvxopt.solvers
import numpy as np
from cvxopt.solvers import lp

from .misc import warn

cvxopt.solvers.options["show_progress"] = False  # disable cvxopt output

GLPK_IF_AVAILABLE: Optional[str] = None
try:
    import cvxopt.glpk

    GLPK_IF_AVAILABLE = "glpk"
    # GLPK is the fastest LP solver I could find so far:
    # <https://scaron.info/blog/linear-programming-in-python-with-cvxopt.html>
    # ... however, it's verbose by default, so tell it to STFU:
    cvxopt.solvers.options["glpk"] = {"msg_lev": "GLP_MSG_OFF"}  # cvxopt 1.1.8
    cvxopt.solvers.options["msg_lev"] = "GLP_MSG_OFF"  # cvxopt 1.1.7
    cvxopt.solvers.options["LPX_K_MSGLEV"] = 0  # previous versions
except ImportError:
    # issue a warning as GLPK is the best LP solver in practice
    warn("GLPK solver not found")


def cvxmat(M: Union[np.ndarray, cvxopt.matrix]) -> cvxopt.matrix:
    """Convert a NumPy array to a CVXOPT matrix."""
    if isinstance(M, cvxopt.matrix):
        return M
    return cvxopt.matrix(M)


def solve_lp(
    c: np.ndarray,
    G: np.ndarray,
    h: np.ndarray,
    A: Optional[np.ndarray] = None,
    b: Optional[np.ndarray] = None,
    solver: Optional[str] = GLPK_IF_AVAILABLE,
) -> np.ndarray:
    r"""Solve a linear program (LP).

    The linear program is defined by:

    .. math::

        \mathrm{minimize} \ & c^T x \\
        \mathrm{subject\ to} \ & G x \leq h \\
            & A x = b

    using the `CVXOPT
    <http://cvxopt.org/userguide/coneprog.html#linear-programming>`_ interface
    to LP solvers.

    Parameters
    ----------
    c :
        Linear-cost vector.
    G :
        Linear inequality constraint matrix.
    h :
        Linear inequality constraint vector.
    A :
        Linear equality constraint matrix.
    b :
        Linear equality constraint vector.
    solver :
        Solver to use, default is GLPK if available

    Returns
    -------
    :
        Optimal solution to the LP.

    Raises
    ------
    ValueError
        If the LP is not feasible.
    """
    args = [cvxmat(c), cvxmat(G), cvxmat(h)]
    if A is not None:
        args.extend([cvxmat(A), cvxmat(b)])
    sol = lp(*args, solver=solver)
    if "optimal" not in sol["status"]:
        raise ValueError("LP optimum not found: %s" % sol["status"])
    return np.array(sol["x"]).reshape((np.array(c).shape[0],))
