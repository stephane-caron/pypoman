#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2016 Quang-Cuong Pham
# Copyright (C) 2017 Stéphane Caron

"""Iterative projection algorithm by [Bretl08]_."""

from typing import Any, List, Optional, Tuple, Union

import numpy as np
from numpy.random import random
from scipy.linalg import norm

from .lp import GLPK_IF_AVAILABLE, solve_lp


class Vertex:
    """Vertex of the projected polygon, with a pointer to its successor."""

    expanded: bool
    next: Optional[Any]
    x: float
    y: float

    def __init__(self, p: Union[List[float], np.ndarray]):
        """
        Initialize vertex from point coordinates.

        Parameters
        ----------
        p :
            2D coordinates of the vertex.
        """
        self.x = p[0]
        self.y = p[1]
        self.next = None
        self.expanded = False

    def expand(
        self,
        lp: Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray],
    ):
        """
        Expand the edge from the vertex to its successor.

        Parameters
        ----------
        lp :
            Tuple `(q, G, h, A, b)` defining the linear program. See
            :func:`pypoman.lp.solve_lp` for details.
        """
        v1 = self
        v2 = self.next
        if v2 is None:
            raise ValueError("cannot expand vertex as it has no successor")
        v = np.array(
            [v2.y - v1.y, v1.x - v2.x]
        )  # orthogonal direction to edge
        v /= norm(v)
        try:
            z = optimize_direction(v, lp)
        except ValueError:
            self.expanded = True
            return None
        xopt: float = z[0]
        yopt: float = z[1]
        if (
            abs(
                np.cross(
                    [xopt - v1.x, yopt - v1.y], [v1.x - v2.x, v1.y - v2.y]
                )
            )
            < 1e-4
        ):
            self.expanded = True
            return None
        vnew = Vertex([xopt, yopt])
        vnew.next = self.next
        self.next = vnew
        self.expanded = False
        return vnew


class Polygon:
    """Polygon, that is, 2D polyhedron."""

    vertices: List[Vertex]

    def __init__(self, v1: Vertex, v2: Vertex, v3: Vertex):
        """
        Initialize polygon from inscribed triangle.

        Parameters
        ----------
        v1 :
            First vertex of the initial triangle.
        v2 :
            Second vertex of the initial triangle.
        v3 :
            Third vertex of the initial triangle.
        """
        v1.next = v2
        v2.next = v3
        v3.next = v1
        self.vertices = [v1, v2, v3]

    def all_expanded(self) -> bool:
        """
        Check for unexpanded vertices.

        Returns
        -------
        :
            True if and only if all vertices have been expanded.
        """
        for v in self.vertices:
            if not v.expanded:
                return False
        return True

    def iter_expand(
        self,
        lp: Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray],
        max_iter: int,
    ) -> None:
        """Extend polygon until the termination condition is reached.

        Termination condition: there is no other vertex to expand, or the
        maximum number of iterations has been reached.

        Parameters
        ----------
        lp :
            Tuple `(q, G, h, A, b)` defining the linear program. See
            :func:`pypoman.lp.solve_lp` for details.
        max_iter :
            Maximum number of iterations.
        """
        nb_iter = 0
        v = self.vertices[0]
        while not self.all_expanded() and nb_iter < max_iter:
            if v.expanded:
                if v.next is None:
                    raise ValueError(
                        "Invalid expanded vertex with no successor"
                    )
                v = v.next
                continue
            vnew = v.expand(lp)
            if vnew is None:
                continue
            self.vertices.append(vnew)
            nb_iter += 1

    def sort_vertices(self):
        """Export vertices starting from the leftmost one and going clockwise.

        Note
        ----
        Assumes all vertices are on the positive halfplane.
        """
        minsd = 1e10
        ibottom = 0
        for i in range(len(self.vertices)):
            vertex = self.vertices[i]
            next_vertex = vertex.next
            if next_vertex is None:
                raise ValueError("Invalid expanded vertex with no successor")
            if (vertex.y + next_vertex.y) < minsd:
                ibottom = i
                minsd = vertex.y + next_vertex.y
        for vertex in self.vertices:
            vertex.checked = False
        vcur = self.vertices[ibottom]
        newvertices = []
        while not vcur.checked:
            vcur.checked = True
            newvertices.append(vcur)
            vcur = vcur.next
        newvertices.reverse()
        vfirst = newvertices.pop(-1)
        newvertices.insert(0, vfirst)
        self.vertices = newvertices

    def export_vertices(self, min_dist: float = 1e-2) -> List[Vertex]:
        """Get list of vertices.

        Parameters
        ----------
        min_dist :
            Minimum distance between two consecutive vertices.

        Returns
        -------
        :
            List of vertices.
        """
        vertices: List[Vertex] = [self.vertices[0]]
        for i in range(1, len(self.vertices) - 1):
            vcur = self.vertices[i]
            vlast = vertices[-1]
            if norm([vcur.x - vlast.x, vcur.y - vlast.y]) > min_dist:
                vertices.append(vcur)
        vertices.append(self.vertices[-1])  # always add last vertex
        return vertices


def optimize_direction(
    vdir: np.ndarray,
    lp: Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray],
    solver: Optional[str] = GLPK_IF_AVAILABLE,
) -> np.ndarray:
    """Optimize in one direction.

    Parameters
    ----------
    vdir :
        Direction (3D vector) in which the optimization is performed.
    lp :
        Tuple `(q, G, h, A, b)` defining the LP. See
        :func:`pypoman.lp..solve_lp` for details.
    solver :
        Backend LP solver to call.

    Returns
    -------
    :
        Vector ``z`` representing the maximum vertex of the polygon in the
        direction `vdir`.
    """
    lp_q, lp_Gextended, lp_hextended, lp_A, lp_b = lp
    lp_q[-2] = -vdir[0]
    lp_q[-1] = -vdir[1]
    x = solve_lp(lp_q, lp_Gextended, lp_hextended, lp_A, lp_b, solver=solver)
    return x[-2:]


def optimize_angle(
    theta: float,
    lp: Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray],
    solver: Optional[str] = GLPK_IF_AVAILABLE,
) -> np.ndarray:
    """Optimize in one direction.

    Parameters
    ----------
    theta :
        Angle of the direction in which the optimization is performed.
    lp :
        Tuple `(q, G, h, A, b)` defining the LP. See
        :func:`pypoman.lp..solve_lp` for details.
    solver :
        Backend LP solver to call.

    Returns
    -------
    :
        Vector ``z`` representing the maximum vertex of the polygon in the
        direction `vdir`.
    """
    d = np.array([np.cos(theta), np.sin(theta)])
    z = optimize_direction(d, lp, solver=solver)
    return z


def compute_polygon(
    lp: Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray],
    max_iter: int = 1000,
    solver: Optional[str] = GLPK_IF_AVAILABLE,
    init_angle: Optional[float] = None,
) -> Polygon:
    """Expand a polygon iteratively.

    Parameters
    ----------
    lp :
        Tuple `(q, G, h, A, b)` defining the linear program. See
        :func:`pypoman.lp.solve_lp` for details.
    max_iter :
        Maximum number of calls to the LP solver.
    solver :
        Name of backend LP solver.
    init_angle :
        Angle in [rad] giving the direction of the initial ray cast.

    Returns
    -------
    :
        Output polygon.
    """
    theta = init_angle if init_angle is not None else np.pi * random()
    init_vertices = [optimize_angle(theta, lp, solver)]
    step = 2.0 * np.pi / 3.0
    while len(init_vertices) < 3 and max_iter >= 0:
        theta += step
        if theta >= 2.0 * np.pi:
            step *= 0.25 + 0.5 * random()
            theta += step - 2.0 * np.pi
        z = optimize_angle(theta, lp, solver)
        if all([norm(z - z0) > 1e-5 for z0 in init_vertices]):
            init_vertices.append(z)
        max_iter -= 1
    if len(init_vertices) < 3:
        raise ValueError("problem is not linearly feasible")
    v0 = Vertex(init_vertices[0])
    v1 = Vertex(init_vertices[1])
    v2 = Vertex(init_vertices[2])
    polygon = Polygon(v0, v1, v2)
    polygon.iter_expand(lp, max_iter)
    return polygon
