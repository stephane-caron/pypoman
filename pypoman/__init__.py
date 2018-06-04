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

from .duality import compute_cone_face_matrix
from .duality import compute_polytope_halfspaces
from .duality import compute_polytope_vertices
from .intersection import intersect_line_cylinder
from .intersection import intersect_line_polygon
from .intersection import intersect_polygons
from .lp import solve_lp
from .polygon import compute_polygon_hull
from .polygon import plot_polygon
from .projection import project_polytope
from .projection import project_polytope_bretl

__all__ = [
    'compute_cone_face_matrix',
    'compute_polygon_hull',
    'compute_polytope_halfspaces',
    'compute_polytope_vertices',
    'intersect_line_cylinder',
    'intersect_line_polygon',
    'intersect_polygons',
    'plot_polygon',
    'project_polytope',
    'project_polytope_bretl',
    'solve_lp',
]
