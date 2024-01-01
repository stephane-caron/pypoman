#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright 2018 Stéphane Caron

"""Functions to switch between halfspace and vertex representations."""

from typing import List, Tuple, Union

import cdd
import numpy as np
from scipy.spatial import ConvexHull

from .misc import norm


def compute_cone_face_matrix(S: np.ndarray) -> np.ndarray:
    r"""Compute the face matrix of a polyhedral cone from its span matrix.

    Parameters
    ----------
    S :
        Span matrix defining the cone as :math:`x = S \lambda` with
        :math:`\lambda \geq 0`.

    Returns
    -------
    :
        Face matrix defining the cone equivalently by :math:`F x \leq 0`.
    """
    V = np.vstack(
        [
            np.hstack([np.zeros((S.shape[1], 1)), S.T]),
            np.hstack([1, np.zeros(S.shape[0])]),
        ]
    )
    # V-representation: first column is 0 for rays
    mat = cdd.Matrix(V, number_type="float")
    mat.rep_type = cdd.RepType.GENERATOR
    P = cdd.Polyhedron(mat)
    ineq = P.get_inequalities()
    H = np.array(ineq)
    if H.shape == (0,):  # H == []
        return H
    A = []
    for i in range(H.shape[0]):
        # H matrix is [b, -A] for A * x <= b
        if norm(H[i, 1:]) < 1e-10:
            continue
        elif abs(H[i, 0]) > 1e-10:  # b should be zero for a cone
            raise ValueError("Polyhedron is not a cone")
        elif i not in ineq.lin_set:
            A.append(-H[i, 1:])
    return np.array(A)


def compute_polytope_halfspaces(
    vertices: List[np.ndarray],
) -> Union[np.ndarray, Tuple[np.ndarray, np.ndarray]]:
    r"""Compute the halfspace representation (H-rep) of a polytope.

    The polytope is defined as convex hull of a set of vertices:

    .. math::

        A x \leq b
        \quad \Leftrightarrow \quad
        x \in \mathrm{conv}(\mathrm{vertices})

    Parameters
    ----------
    vertices :
        List of polytope vertices.

    Returns
    -------
    :
        Tuple ``(A, b)`` of the halfspace representation, or empty array if it
        is empty.
    """
    V = np.vstack(vertices)
    t = np.ones((V.shape[0], 1))  # first column is 1 for vertices
    tV = np.hstack([t, V])
    mat = cdd.Matrix(tV, number_type="float")
    mat.rep_type = cdd.RepType.GENERATOR
    P = cdd.Polyhedron(mat)
    bA = np.array(P.get_inequalities())
    if bA.shape == (0,):  # bA == []
        return bA
    # the polyhedron is given by b + A x >= 0 where bA = [b|A]
    b, A = np.array(bA[:, 0]), -np.array(bA[:, 1:])
    return (A, b)


def compute_polytope_vertices(
    A: np.ndarray, b: np.ndarray
) -> List[np.ndarray]:
    r"""Compute the vertices of a polytope.

    The polytope is given in halfspace representation by :math:`A x \leq b`.

    Parameters
    ----------
    A :
        Matrix of halfspace representation.
    b :
        Vector of halfspace representation.

    Returns
    -------
    :
        List of polytope vertices.

    Notes
    -----
    This method won't work well if your halfspace representation includes
    equality constraints :math:`A x = b` written as :math:`(A x \leq b \wedge
    -A x \leq -b)`. If this is your use case, consider using directly the
    linear set ``lin_set`` of `equality-constraint generatorsin pycddlib
    <https://pycddlib.readthedocs.io/en/latest/matrix.html>`_.
    """
    b = b.reshape((b.shape[0], 1))
    mat = cdd.Matrix(np.hstack([b, -A]), number_type="float")
    mat.rep_type = cdd.RepType.INEQUALITY
    P = cdd.Polyhedron(mat)
    g = P.get_generators()
    V = np.array(g)
    vertices = []
    for i in range(V.shape[0]):
        if V[i, 0] != 1:  # 1 = vertex, 0 = ray
            raise ValueError("Polyhedron is not a polytope")
        elif i not in g.lin_set:
            vertices.append(V[i, 1:])
    return vertices


def convex_hull(points: List[np.ndarray]) -> List[np.ndarray]:
    """Compute the convex hull of a set of points.

    Parameters
    ----------
    points :
        Set of points.

    Returns
    -------
    :
        List of polytope vertices.
    """
    hull = ConvexHull(points)
    return [points[i] for i in hull.vertices]
