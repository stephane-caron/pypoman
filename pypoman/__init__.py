#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2018-2020 Stephane Caron <stephane.caron@normalesup.org>
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

"""Python module for polyhedral geometry."""

from .duality import (
    compute_cone_face_matrix,
    compute_polytope_halfspaces,
    compute_polytope_vertices,
    convex_hull,
)
from .intersection import (
    intersect_line_cylinder,
    intersect_line_polygon,
    intersect_polygons,
)
from .lp import solve_lp
from .polygon import compute_polygon_hull, plot_polygon
from .polyhedron import compute_chebyshev_center
from .projection import (
    project_point_to_polytope,
    project_polytope,
    project_polytope_bretl,
)

__version__ = "0.6.0rc0"

__all__ = [
    "compute_chebyshev_center",
    "compute_cone_face_matrix",
    "compute_polygon_hull",
    "compute_polytope_halfspaces",
    "compute_polytope_vertices",
    "convex_hull",
    "intersect_line_cylinder",
    "intersect_line_polygon",
    "intersect_polygons",
    "plot_polygon",
    "project_polytope",
    "project_polytope_bretl",
    "project_point_to_polytope",
    "solve_lp",
]
