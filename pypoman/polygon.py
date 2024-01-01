#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright 2018 Stéphane Caron

"""Functions on polygons, that is, 2D polyhedra."""

from typing import List, Optional

import numpy as np
from matplotlib.patches import Polygon
from pylab import axis, gca
from scipy.spatial import ConvexHull

from .polyhedron import compute_chebyshev_center


def __compute_polygon_hull(B: np.ndarray, c: np.ndarray):
    r"""Compute the vertex representation of a polygon.

    The polygon is defined by:

    .. math::

        B x \leq c

    where :math:`x` is a 2D vector.

    Parameters
    ----------
    B :
        Linear inequality matrix of size :math:`2 \times K`.
    c :
        Linear inequality vector of with positive coordinates.

    Returns
    -------
    vertices : list of arrays
        List of 2D vertices in counterclowise order.

    Notes
    -----
    The origin [0, 0] should lie inside the polygon (:math:`c \geq 0`) in
    order to build the polar form. If you don't have this guarantee, call
    ``compute_polar_polygon()`` instead.

    Checking that :math:`c > 0` is not optional. The rest of the algorithm can
    be executed when some coordinates :math:`c_i < 0`, but the result would be
    wrong.
    """
    assert (
        B.shape[1] == 2
    ), "Input (B, c) is not a polygon: B.shape = %s" % str(B.shape)
    assert all(
        c > 0
    ), "Polygon should contain the origin, but min(c) = %.2f" % min(c)

    B_polar = np.hstack(
        [
            (B[:, column] * 1.0 / c).reshape((B.shape[0], 1))
            for column in range(2)
        ]
    )

    def axis_intersection(i, j):
        ai, bi = c[i], B[i]
        aj, bj = c[j], B[j]
        x = (ai * bj[1] - aj * bi[1]) * 1.0 / (bi[0] * bj[1] - bj[0] * bi[1])
        y = (bi[0] * aj - bj[0] * ai) * 1.0 / (bi[0] * bj[1] - bj[0] * bi[1])
        return np.array([x, y])

    # QHULL OPTIONS:
    #
    # - ``Pp`` -- do not report precision problems
    # - ``Q0`` -- no merging with C-0 and Qx
    #
    # ``Q0`` avoids [this bug](https://github.com/scipy/scipy/issues/6484).
    # It slightly diminishes computation times (0.9 -> 0.8 ms on my machine)
    # but raises QhullError at the first sight of precision errors.
    #
    hull = ConvexHull([row for row in B_polar], qhull_options="Pp Q0")
    #
    # contrary to hull.simplices (which was not in practice), hull.vertices is
    # guaranteed to be in counterclockwise order for 2-D (see scipy doc)
    #
    simplices = [
        (hull.vertices[i], hull.vertices[i + 1])
        for i in range(len(hull.vertices) - 1)
    ]
    simplices.append((hull.vertices[-1], hull.vertices[0]))
    vertices = [axis_intersection(i, j) for (i, j) in simplices]
    return vertices


def compute_polygon_hull(B: np.ndarray, c: np.ndarray) -> List[np.ndarray]:
    r"""Compute the vertex representation of a polygon.

    The polygon is defined by:

    .. math::

        B x \leq c

    where :math:`x` is a 2D vector.

    Parameters
    ----------
    B :
        Linear inequality matrix of size :math:`2 \times K`.
    c :
        Linear inequality vector.

    Returns
    -------
    :
        List of 2D vertices in counterclockwise order.
    """
    x = None
    if not all(c > 0):
        x = compute_chebyshev_center(B, c)
        c = c - np.dot(B, x)
    if not all(c > 0):
        raise ValueError("Polygon is empty (min. dist. to edge %.2f)" % min(c))
    vertices = __compute_polygon_hull(B, c)
    if x is not None:
        vertices = [v + x for v in vertices]
    return vertices


def plot_polygon(
    points: np.ndarray,
    alpha: float = 0.4,
    color: str = "g",
    linestyle: str = "solid",
    fill: bool = True,
    linewidth: Optional[float] = None,
    resize: bool = False,
) -> None:
    """
    Plot a polygon in matplotlib.

    Parameters
    ----------
    points :
        Array or list of points.
    alpha :
        Transparency value.
    color :
        Color in matplotlib format.
    linestyle :
        Line style in matplotlib format.
    fill :
        When ``True``, fills the area inside the polygon.
    linewidth :
        Line width in matplotlib format.
    resize :
        When ``True``, resets axis limits to center on the polygon.
    """
    if isinstance(points, list):
        points = np.array(points)
    ax = gca()
    hull = ConvexHull(points)
    points = points[hull.vertices, :]
    if resize:
        xmin1, xmax1, ymin1, ymax1 = axis()
        xmin2, ymin2 = 1.5 * points.min(axis=0)
        xmax2, ymax2 = 1.5 * points.max(axis=0)
        axis(
            (
                min(xmin1, xmin2),
                max(xmax1, xmax2),
                min(ymin1, ymin2),
                max(ymax1, ymax2),
            )
        )
    patch = Polygon(
        points,
        alpha=alpha,
        color=color,
        linestyle=linestyle,
        fill=fill,
        linewidth=linewidth,
    )
    ax.add_patch(patch)
