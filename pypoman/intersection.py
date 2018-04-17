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

from numpy import array
from scipy.spatial import ConvexHull

from .misc import norm


PREC_TOL = 1e-10  # numerical tolerance


def intersect_line_polygon(line, vertices, apply_hull):
    """
    Intersect a line segment with a polygon.

    Parameters
    ----------
    line : couple of arrays
        End points of the line segment (2D or 3D).
    vertices : list of arrays
        Vertices of the polygon.
    apply_hull : bool
        Set to `True` to apply a convex hull algorithm to `vertices`.
        Otherwise, the function assumes that vertices are already sorted in
        clockwise or counterclockwise order.

    Returns
    -------
    inter_points : list of array
        List of intersection points between the line segment and the polygon.

    Notes
    -----
    This code is adapted from
    <https://stackoverflow.com/questions/20677795/how-do-i-compute-the-intersection-point-of-two-lines-in-python/20679579#20679579>.
    On the same setting with `apply_hull=False`, it %timeits to 6 us.
    """
    def line_coordinates(p1, p2):
        A = (p1[1] - p2[1])
        B = (p2[0] - p1[0])
        C = (p1[0] * p2[1] - p2[0] * p1[1])
        return A, B, -C

    def intersection(L1, L2):
        D = L1[0] * L2[1] - L1[1] * L2[0]
        Dx = L1[2] * L2[1] - L1[1] * L2[2]
        Dy = L1[0] * L2[2] - L1[2] * L2[0]
        if abs(D) < 1e-5:
            return None
        x = Dx / D
        y = Dy / D
        return x, y

    if apply_hull:
        points = vertices
        hull = ConvexHull(points)
        vertices = [points[i] for i in hull.vertices]

    n = len(vertices)
    p1, p2 = line
    L1 = line_coordinates(p1, p2)
    x_min, x_max = min(p1[0], p2[0]), max(p1[0], p2[0])
    y_min, y_max = min(p1[1], p2[1]), max(p1[1], p2[1])
    inter_points = []
    for i, v1 in enumerate(vertices):
        v2 = vertices[(i + 1) % n]
        L2 = line_coordinates(v1, v2)
        p = intersection(L1, L2)
        if p is not None:
            if not (x_min <= p[0] <= x_max and y_min <= p[1] <= y_max):
                continue
            vx_min, vx_max = min(v1[0], v2[0]), max(v1[0], v2[0])
            vy_min, vy_max = min(v1[1], v2[1]), max(v1[1], v2[1])
            if not (vx_min - PREC_TOL <= p[0] <= vx_max + PREC_TOL and
                    vy_min - PREC_TOL <= p[1] <= vy_max + PREC_TOL):
                continue
            inter_points.append(array(p))
    return inter_points


def intersect_line_cylinder(line, vertices):
    """
    Intersect the line segment [p1, p2] with a vertical cylinder of polygonal
    cross-section. If the intersection has two points, returns the one closest
    to p1.

    Parameters
    ----------
    line : couple of (3,) arrays
        End points of the 3D line segment.
    vertices : list of (3,) arrays
        Vertices of the polygon.

    Returns
    -------
    inter_points : list of (3,) arrays
        List of intersection points between the line segment and the cylinder.
    """
    inter_points = []
    inter_2d = intersect_line_polygon(line, vertices, apply_hull=True)
    for p in inter_2d:
        p1, p2 = array(line[0]), array(line[1])
        alpha = norm(p - p1[:2]) / norm(p2[:2] - p1[:2])
        z = p1[2] + alpha * (p2[2] - p1[2])
        inter_points.append(array([p[0], p[1], z]))
    return inter_points


def intersect_polygons(polygon1, polygon2):
    """
    Intersect two polygons.

    Parameters
    ----------
    polygon1 : list of arrays
        Vertices of the first polygon in counterclockwise order.
    polygon1 : list of arrays
        Vertices of the second polygon in counterclockwise order.

    Returns
    -------
    intersection : list of arrays
        Vertices of the intersection in counterclockwise order.
    """
    from pyclipper import Pyclipper, PT_CLIP, PT_SUBJECT, CT_INTERSECTION
    from pyclipper import scale_to_clipper, scale_from_clipper
    # could be accelerated by removing the scale_to/from_clipper()
    subj, clip = (polygon1,), polygon2
    pc = Pyclipper()
    pc.AddPath(scale_to_clipper(clip), PT_CLIP)
    pc.AddPaths(scale_to_clipper(subj), PT_SUBJECT)
    solution = pc.Execute(CT_INTERSECTION)
    if not solution:
        return []
    return scale_from_clipper(solution)[0]
