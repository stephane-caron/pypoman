#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright 2021 Stéphane Caron

import pylab
from numpy import arange, array, cos, pi, sin

import pypoman

vertices = [(cos(theta), sin(theta)) for theta in arange(0, 2 * pi, pi / 6)]
A, b = pypoman.compute_polytope_halfspaces(vertices)

point = array([2.1, 1.9])
proj = pypoman.project_point_to_polytope(point, (A, b), qpsolver="cvxopt")


if __name__ == "__main__":  # plot projected polytope
    pylab.ion()
    pylab.figure()
    pylab.gca().set_aspect("equal")
    pypoman.plot_polygon(vertices)
    pylab.plot([point[0]], [point[1]], marker="o", markersize=3, color="r")
    pylab.plot([proj[0]], [proj[1]], marker="o", markersize=3, color="b")
    pylab.plot([point[0], proj[0]], [point[1], proj[1]], "b--")
    pylab.show(block=True)
