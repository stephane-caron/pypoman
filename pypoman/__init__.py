# SPDX-License-Identifier: GPL-3.0-or-later

"""Python module for polyhedral geometry."""

from .duality import (
    compute_cone_face_matrix,
    compute_polytope_halfspaces,
    compute_polytope_vertices,
    convex_hull,
)
from .exceptions import (  # pylint: disable=redefined-builtin
    ConversionError,
    EmptyPolytopeError,
    Exception,
    SolverError,
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

__version__ = "2.0.0"

__all__ = [
    "ConversionError",
    "EmptyPolytopeError",
    "Exception",
    "SolverError",
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
