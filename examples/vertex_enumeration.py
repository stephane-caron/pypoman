#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2016 Quang-Cuong Pham <cuong.pham@normalesup.org>
# Copyright (C) 2017-2018 Stephane Caron <stephane.caron@lirmm.fr>
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

import numpy
import pypoman

A = numpy.array([
    [-1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
    [0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
    [0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0],
    [0,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0],
    [0,  0,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0],
    [0,  0,  0,  0,  0, -1,  0,  0,  0,  0,  0,  0],
    [0,  0,  0,  0,  0,  0, -1,  0,  0,  0,  0,  0],
    [0,  0,  0,  0,  0,  0,  0, -1,  0,  0,  0,  0],
    [0,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0,  0],
    [0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0],
    [0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  0],
    [0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1],
    [1,  1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0],
    [0,  0,  0,  1,  1,  1,  0,  0,  0,  0,  0,  0],
    [0,  0,  0,  0,  0,  0,  1,  1,  1,  0,  0,  0],
    [0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  1,  1],
    [1,  0,  0,  1,  0,  0,  1,  0,  0,  1,  0,  0],
    [0,  1,  0,  0,  1,  0,  0,  1,  0,  0,  1,  0],
    [0,  0,  1,  0,  0,  1,  0,  0,  1,  0,  0,  1]])
b = numpy.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 2, 2, 1, 2, 3])
vertices = pypoman.compute_polytope_vertices(A, b)

if __name__ == "__main__":
    print("Polytope has %d vertices:" % len(vertices))
    for i, vertex in enumerate(vertices):
        print("Vertex %3d:" % i, vertex)
    if IPython.get_ipython() is None:  # give the user a prompt
        IPython.embed()
