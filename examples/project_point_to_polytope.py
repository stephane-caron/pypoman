#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2021 Stephane Caron <stephane.caron@normalesup.org>
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

import IPython
import pylab

from numpy import arange, array, cos, pi, sin

import pypoman

vertices = [(cos(theta), sin(theta)) for theta in arange(0, 2 * pi, pi / 6)]
A, b = pypoman.compute_polytope_halfspaces(vertices)

point = array([2.1, 1.9])
proj = pypoman.project_point_to_polytope(point, (A, b))


if __name__ == "__main__":  # plot projected polytope
    pylab.ion()
    pylab.figure()
    pylab.gca().set_aspect("equal")
    pypoman.plot_polygon(vertices)
    pylab.plot([point[0]], [point[1]], marker='o', markersize=3, color='r')
    pylab.plot([proj[0]], [proj[1]], marker='o', markersize=3, color='b')
    pylab.plot([point[0], proj[0]], [point[1], proj[1]], 'b--')
    if IPython.get_ipython() is None:  # give the user a prompt
        IPython.embed()
