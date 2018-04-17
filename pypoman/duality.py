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

import cdd

from numpy import array, hstack, ones, vstack, zeros

from .misc import norm


def compute_cone_face_matrix(S):
    """
    Compute the face matrix of a polyhedral convex cone from its span matrix.

    Parameters
    ----------
    S : array, shape=(n, m)
        Span matrix defining the cone as :math:`x = S \\lambda` with
        :math:`\\lambda \\geq 0`.

    Returns
    -------
    F : array, shape=(k, n)
        Face matrix defining the cone equivalently by :math:`F x \\leq 0`.
    """
    V = vstack([
        hstack([zeros((S.shape[1], 1)), S.T]),
        hstack([1, zeros(S.shape[0])])])
    # V-representation: first column is 0 for rays
    mat = cdd.Matrix(V, number_type='float')
    mat.rep_type = cdd.RepType.GENERATOR
    P = cdd.Polyhedron(mat)
    ineq = P.get_inequalities()
    H = array(ineq)
    if H.shape == (0,):  # H == []
        return H
    A = []
    for i in range(H.shape[0]):
        # H matrix is [b, -A] for A * x <= b
        if norm(H[i, 1:]) < 1e-10:
            continue
        elif abs(H[i, 0]) > 1e-10:  # b should be zero for a cone
            raise Exception("Polyhedron is not a cone")
        elif i not in ineq.lin_set:
            A.append(-H[i, 1:])
    return array(A)


def compute_polytope_halfspaces(vertices):
    """
    Compute the halfspace representation (H-rep) of a polytope defined as
    convex hull of a set of vertices:

    .. math::

        A x \\leq b
        \\quad \\Leftrightarrow \\quad
        x \\in \\mathrm{conv}(\\mathrm{vertices})

    Parameters
    ----------
    vertices : list of arrays
        List of polytope vertices.

    Returns
    -------
    A : array, shape=(m, k)
        Matrix of halfspace representation.
    b : array, shape=(m,)
        Vector of halfspace representation.
    """
    V = vstack(vertices)
    t = ones((V.shape[0], 1))  # first column is 1 for vertices
    tV = hstack([t, V])
    mat = cdd.Matrix(tV, number_type='float')
    mat.rep_type = cdd.RepType.GENERATOR
    P = cdd.Polyhedron(mat)
    bA = array(P.get_inequalities())
    if bA.shape == (0,):  # bA == []
        return bA
    # the polyhedron is given by b + A x >= 0 where bA = [b|A]
    b, A = array(bA[:, 0]), -array(bA[:, 1:])
    return (A, b)


def compute_polytope_vertices(A, b):
    """
    Compute the vertices of a polytope given in halfspace representation by
    :math:`A x \\leq b`.

    Parameters
    ----------
    A : array, shape=(m, k)
        Matrix of halfspace representation.
    b : array, shape=(m,)
        Vector of halfspace representation.

    Returns
    -------
    vertices : list of arrays
        List of polytope vertices.
    """
    b = b.reshape((b.shape[0], 1))
    mat = cdd.Matrix(hstack([b, -A]), number_type='float')
    mat.rep_type = cdd.RepType.INEQUALITY
    P = cdd.Polyhedron(mat)
    g = P.get_generators()
    V = array(g)
    vertices = []
    for i in range(V.shape[0]):
        if V[i, 0] != 1:  # 1 = vertex, 0 = ray
            raise Exception("Polyhedron is not a polytope")
        elif i not in g.lin_set:
            vertices.append(V[i, 1:])
    return vertices
