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

import cdd
import cvxopt

from numpy import array, dot, hstack, zeros

from .bretl import compute_polygon as bretl_compute_polygon


def project_polyhedron(proj, ineq, eq=None, canonicalize=True):
    """
    Apply the affine projection :math:`y = E x + f` to the polyhedron defined
    by:

    .. math::

        A x & \\leq b \\\\
        C x & = d

    Parameters
    ----------
    proj : pair of arrays
        Pair (`E`, `f`) describing the affine projection.
    ineq : pair of arrays
        Pair (`A`, `b`) describing the inequality constraint.
    eq : pair of arrays, optional
        Pair (`C`, `d`) describing the equality constraint.
    canonicalize : bool, optional
        Apply equality constraints from `eq` to reduce the dimension of the
        input polyhedron. May be a blessing or a curse, see notes below.

    Returns
    -------
    vertices : list of arrays
        List of vertices of the projection.
    rays : list of arrays
        List of rays of the projection.

    Notes
    -----
    When the equality set `eq` of the input polytope is not empty, it is
    usually faster to use these equality constraints to reduce the dimension of
    the input polytope (cdd function: `canonicalize()`) before enumerating
    vertices (cdd function: `get_generators()`). Yet, on some descriptions this
    operation may be problematic: if it fails, or if you get empty outputs when
    the output is supposed to be non-empty, you can try setting
    `canonicalize=False`.

    See also
    --------
    This webpage: https://scaron.info/teaching/projecting-polytopes.html
    """
    # the input [b, -A] to cdd.Matrix represents (b - A * x >= 0)
    # see ftp://ftp.ifor.math.ethz.ch/pub/fukuda/cdd/cddlibman/node3.html
    (A, b) = ineq
    b = b.reshape((b.shape[0], 1))
    linsys = cdd.Matrix(hstack([b, -A]), number_type='float')
    linsys.rep_type = cdd.RepType.INEQUALITY

    # the input [d, -C] to cdd.Matrix.extend represents (d - C * x == 0)
    # see ftp://ftp.ifor.math.ethz.ch/pub/fukuda/cdd/cddlibman/node3.html
    if eq is not None:
        (C, d) = eq
        d = d.reshape((d.shape[0], 1))
        linsys.extend(hstack([d, -C]), linear=True)
        if canonicalize:
            linsys.canonicalize()

    # Convert from H- to V-representation
    P = cdd.Polyhedron(linsys)
    generators = P.get_generators()
    if generators.lin_set:
        print("Generators have linear set: {}".format(generators.lin_set))
    V = array(generators)

    # Project output wrenches to 2D set
    (E, f) = proj
    vertices, rays = [], []
    free_coordinates = []
    for i in range(V.shape[0]):
        if generators.lin_set and i in generators.lin_set:
            free_coordinates.append(list(V[i, 1:]).index(1.))
        elif V[i, 0] == 1:  # vertex
            vertices.append(dot(E, V[i, 1:]) + f)
        else:  # ray
            rays.append(dot(E, V[i, 1:]))
    return vertices, rays


def project_polytope(proj, ineq, eq=None, method='cdd', **kwargs):
    """
    Apply the affine projection :math:`y = E x + f` to the polytope defined by:

    .. math::

        A x & \\leq b \\\\
        C x & = d

    Parameters
    ----------
    proj : pair of arrays
        Pair (`E`, `f`) describing the affine projection.
    ineq : pair of arrays
        Pair (`A`, `b`) describing the inequality constraint.
    eq : pair of arrays, optional
        Pair (`C`, `d`) describing the equality constraint.
    method : string, optional
        Choice between 'bretl' and 'cdd'.

    Returns
    -------
    vertices : list of arrays
        List of vertices of the projection.

    Note
    ----
    Additional keyword arguments can be provided when the method is 'bretl'.
    They are passed directly to the corresponding function
    :func:`pypoman.projection.project_polytope_bretl`.

    Notes
    -----
    The number of columns of all matrices `A`, `C` and `E` corresponds to the
    dimension of the input space, while the number of lines of `E` corresponds
    to the dimension of the output space.
    """
    if method == 'bretl':
        assert eq is not None, "Bretl method requires = constraints for now"
        return project_polytope_bretl(proj, ineq, eq, **kwargs)
    vertices, rays = project_polyhedron(proj, ineq, eq)
    assert not rays, "Projection is not a polytope"
    return vertices


def project_polytope_bretl(proj, ineq, eq, max_radius=1e5, max_iter=1000,
                           init_angle=None):
    """
    Project a polytope into a 2D polygon using the incremental projection
    algorithm from [Bretl08]_. The 2D affine projection :math:`y = E x + f` is
    applied to the polyhedron defined by:

    .. math::

        A x & \\leq b \\\\
        C x & = d

    Parameters
    ----------
    proj : pair of arrays
        Pair (`E`, `f`) describing the affine projection.
    ineq : pair of arrays
        Pair (`A`, `b`) describing the inequality constraint.
    eq : pair of arrays, optional
        Pair (`C`, `d`) describing the equality constraint.
    max_radius : scalar, optional
        Maximum distance from origin (in [m]) used to make sure the output
        is bounded.
    max_iter : integer, optional
        Maximum number of calls to the LP solver.
    init_angle : scalar, optional
        Angle in [rad] giving the direction of the initial ray cast.

    Returns
    -------
    vertices : list of arrays
        List of vertices of the projected polygon.
    """
    (E, f), (A, b), (C, d) = proj, ineq, eq
    assert E.shape[0] == f.shape[0] == 2

    # Inequality constraints: A_ext * [ x  u  v ] <= b_ext iff
    # (1) A * x <= b and (2) |u|, |v| <= max_radius
    A_ext = zeros((A.shape[0] + 4, A.shape[1] + 2))
    A_ext[:-4, :-2] = A
    A_ext[-4, -2] = 1
    A_ext[-3, -2] = -1
    A_ext[-2, -1] = 1
    A_ext[-1, -1] = -1
    A_ext = cvxopt.matrix(A_ext)

    b_ext = zeros(b.shape[0] + 4)
    b_ext[:-4] = b
    b_ext[-4:] = array([max_radius] * 4)
    b_ext = cvxopt.matrix(b_ext)

    # Equality constraints: C_ext * [ x  u  v ] == d_ext iff
    # (1) C * x == d and (2) [ u  v ] == E * x + f
    C_ext = zeros((C.shape[0] + 2, C.shape[1] + 2))
    C_ext[:-2, :-2] = C
    C_ext[-2:, :-2] = E[:2]
    C_ext[-2:, -2:] = array([[-1, 0], [0, -1]])
    C_ext = cvxopt.matrix(C_ext)

    d_ext = zeros(d.shape[0] + 2)
    d_ext[:-2] = d
    d_ext[-2:] = -f[:2]
    d_ext = cvxopt.matrix(d_ext)

    lp_obj = cvxopt.matrix(zeros(A.shape[1] + 2))
    lp = lp_obj, A_ext, b_ext, C_ext, d_ext
    polygon = bretl_compute_polygon(lp, max_iter=max_iter,
                                    init_angle=init_angle)
    polygon.sort_vertices()
    vertices_list = polygon.export_vertices()
    vertices = [array([v.x, v.y]) for v in vertices_list]
    return vertices
