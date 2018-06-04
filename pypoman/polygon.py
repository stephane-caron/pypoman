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

from matplotlib.patches import Polygon
from numpy import array, dot, hstack
from pylab import axis, gca
from scipy.spatial import ConvexHull

from .polyhedron import compute_chebyshev_center


def __compute_polygon_hull(B, c):
    """
    Compute the vertex representation of a polygon defined by:

    .. math::

        B x \\leq c

    where `x` is a 2D vector.

    Parameters
    ----------
    B : array, shape=(2, K)
        Linear inequality matrix.
    c : array, shape=(K,)
        Linear inequality vector with positive coordinates.

    Returns
    -------
    vertices : list of arrays
        List of 2D vertices in counterclowise order.

    Notes
    -----
    The origin [0, 0] should lie inside the polygon (:math:`c \\geq 0`) in
    order to build the polar form. If you don't have this guarantee, call
    ``compute_polar_polygon()`` instead.

    Checking that :math:`c > 0` is not optional. The rest of the algorithm can
    be executed when some coordinates :math:`c_i < 0`, but the result would be
    wrong.
    """
    assert B.shape[1] == 2, \
        "Input (B, c) is not a polygon: B.shape = %s" % str(B.shape)
    assert all(c > 0), \
        "Polygon should contain the origin, but min(c) = %.2f" % min(c)

    B_polar = hstack([
        (B[:, column] * 1. / c).reshape((B.shape[0], 1))
        for column in range(2)])

    def axis_intersection(i, j):
        ai, bi = c[i], B[i]
        aj, bj = c[j], B[j]
        x = (ai * bj[1] - aj * bi[1]) * 1. / (bi[0] * bj[1] - bj[0] * bi[1])
        y = (bi[0] * aj - bj[0] * ai) * 1. / (bi[0] * bj[1] - bj[0] * bi[1])
        return array([x, y])

    # QHULL OPTIONS:
    #
    # - ``Pp`` -- do not report precision problems
    # - ``Q0`` -- no merging with C-0 and Qx
    #
    # ``Q0`` avoids [this bug](https://github.com/scipy/scipy/issues/6484).
    # It slightly diminishes computation times (0.9 -> 0.8 ms on my machine)
    # but raises QhullError at the first sight of precision errors.
    #
    hull = ConvexHull([row for row in B_polar], qhull_options='Pp Q0')
    #
    # contrary to hull.simplices (which was not in practice), hull.vertices is
    # guaranteed to be in counterclockwise order for 2-D (see scipy doc)
    #
    simplices = [(hull.vertices[i], hull.vertices[i + 1])
                 for i in range(len(hull.vertices) - 1)]
    simplices.append((hull.vertices[-1], hull.vertices[0]))
    vertices = [axis_intersection(i, j) for (i, j) in simplices]
    return vertices


def compute_polygon_hull(B, c):
    """
    Compute the vertex representation of a polygon defined by:

    .. math::

        B x \\leq c

    where `x` is a 2D vector.

    Parameters
    ----------
    B : array, shape=(2, K)
        Linear inequality matrix.
    c : array, shape=(K,)
        Linear inequality vector.

    Returns
    -------
    vertices : list of arrays
        List of 2D vertices in counterclockwise order.
    """
    x = None
    if not all(c > 0):
        x = compute_chebyshev_center(B, c)
        c = c - dot(B, x)
    if not all(c > 0):
        raise Exception("Polygon is empty (min. dist. to edge %.2f)" % min(c))
    vertices = __compute_polygon_hull(B, c)
    if x is not None:
        vertices = [v + x for v in vertices]
    return vertices


def plot_polygon(points, alpha=.4, color='g', linestyle='solid', fill=True,
                 linewidth=None):
    """
    Plot a polygon in matplotlib.

    Parameters
    ----------
    points : list of arrays
        List of poitns.
    alpha : scalar, optional
        Transparency value.
    color : string, optional
        Color in matplotlib format.
    linestyle : scalar, optional
        Line style in matplotlib format.
    fill : bool, optional
        When ``True``, fills the area inside the polygon.
    linewidth : scalar, optional
        Line width in matplotlib format.
    """
    if type(points) is list:
        points = array(points)
    ax = gca()
    hull = ConvexHull(points)
    points = points[hull.vertices, :]
    xmin1, xmax1, ymin1, ymax1 = axis()
    xmin2, ymin2 = 1.5 * points.min(axis=0)
    xmax2, ymax2 = 1.5 * points.max(axis=0)
    axis((min(xmin1, xmin2), max(xmax1, xmax2),
          min(ymin1, ymin2), max(ymax1, ymax2)))
    patch = Polygon(
        points, alpha=alpha, color=color, linestyle=linestyle, fill=fill,
        linewidth=linewidth)
    ax.add_patch(patch)
