#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright 2018 Stéphane Caron

"""Polytope projection functions."""

from typing import List, Optional, Tuple

import cdd
import cvxopt
import numpy as np

from .bretl import compute_polygon as bretl_compute_polygon


def project_polyhedron(
    proj: Tuple[np.ndarray, np.ndarray],
    ineq: Tuple[np.ndarray, np.ndarray],
    eq: Optional[Tuple[np.ndarray, np.ndarray]] = None,
    canonicalize: bool = True,
) -> Tuple[List[np.ndarray], List[np.ndarray]]:
    r"""Apply the affine projection :math:`y = E x + f` to a polyhedron.

    The polyhedron is defined by:

    .. math::

        \begin{split}\begin{array}{ll}
            A x & \leq b \\
            C x & = d
        \end{array}\end{split}

    Parameters
    ----------
    proj :
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
    This webpage: https://scaron.info/blog/projecting-polytopes.html
    """
    # the input [b, -A] to cdd.Matrix represents (b - A * x >= 0)
    # see ftp://ftp.ifor.math.ethz.ch/pub/fukuda/cdd/cddlibman/node3.html
    (A, b) = ineq
    b = b.reshape((b.shape[0], 1))
    linsys = cdd.matrix_from_array(np.hstack([b, -A]))  # type: ignore
    linsys.rep_type = cdd.RepType.INEQUALITY

    # the input [d, -C] to the cdd function represents (d - C * x == 0)
    # see ftp://ftp.ifor.math.ethz.ch/pub/fukuda/cdd/cddlibman/node3.html
    if eq is not None:
        (C, d) = eq
        d = d.reshape((d.shape[0], 1))
        new_matrix = cdd.matrix_from_array(np.hstack([d, -C]))  # type: ignore
        cdd.matrix_append_to(linsys, new_matrix)
        linsys.rep_type = cdd.RepType.INEQUALITY
        if canonicalize:
            cdd.matrix_canonicalize(linsys)

    # Convert from H- to V-representation
    P = cdd.polyhedron_from_matrix(linsys)
    generators = cdd.copy_generators(P)
    if generators.lin_set:
        print("Generators have linear set: {}".format(generators.lin_set))
    V = np.array(generators.array)

    # Project output wrenches to 2D set
    (E, f) = proj
    vertices, rays = [], []
    free_coordinates = []
    for i in range(V.shape[0]):
        if generators.lin_set and i in generators.lin_set:
            free_coordinates.append(list(V[i, 1:]).index(1.0))
        elif V[i, 0] == 1:  # vertex
            vertices.append(np.dot(E, V[i, 1:]) + f)
        else:  # ray
            rays.append(np.dot(E, V[i, 1:]))
    return vertices, rays


def project_polytope(proj, ineq, eq=None, method="cdd", **kwargs):
    r"""Apply the affine projection :math:`y = E x + f` to a polytope.

    The polytope is defined by:

    .. math::

        A x & \leq b \\
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
    if method == "bretl":
        assert eq is not None, "Bretl method requires = constraints for now"
        return project_polytope_bretl(proj, ineq, eq, **kwargs)
    vertices, rays = project_polyhedron(proj, ineq, eq)
    assert not rays, "Projection is not a polytope"
    return vertices


def project_polytope_bretl(
    proj: Tuple[np.ndarray, np.ndarray],
    ineq: Tuple[np.ndarray, np.ndarray],
    eq: Tuple[np.ndarray, np.ndarray],
    max_radius: float = 1e5,
    max_iter: int = 1000,
    init_angle: Optional[float] = None,
) -> List[np.ndarray]:
    r"""Project a polytope into a 2D polygon using the IP algorithm.

    The incremental projection algorithm is detailed in [Bretl08]_. The 2D
    affine projection :math:`y = E x + f` is applied to the polyhedron defined
    by:

    .. math::

        A x & \leq b \\
        C x & = d

    Parameters
    ----------
    proj :
        Pair (`E`, `f`) describing the affine projection.
    ineq :
        Pair (`A`, `b`) describing the inequality constraint.
    eq :
        Pair (`C`, `d`) describing the equality constraint.
    max_radius :
        Maximum distance from origin (in [m]) used to make sure the output
        is bounded.
    max_iter :
        Maximum number of calls to the LP solver.
    init_angle :
        Angle in [rad] giving the direction of the initial ray cast.

    Returns
    -------
    :
        List of vertices of the projected polygon.
    """
    (E, f), (A, b), (C, d) = proj, ineq, eq
    assert E.shape[0] == f.shape[0] == 2

    # Inequality constraints: A_ext * [ x  u  v ] <= b_ext iff
    # (1) A * x <= b and (2) |u|, |v| <= max_radius
    A_ext = np.zeros((A.shape[0] + 4, A.shape[1] + 2))
    A_ext[:-4, :-2] = A
    A_ext[-4, -2] = 1
    A_ext[-3, -2] = -1
    A_ext[-2, -1] = 1
    A_ext[-1, -1] = -1
    A_ext = cvxopt.matrix(A_ext)

    b_ext = np.zeros(b.shape[0] + 4)
    b_ext[:-4] = b
    b_ext[-4:] = np.array([max_radius] * 4)
    b_ext = cvxopt.matrix(b_ext)

    # Equality constraints: C_ext * [ x  u  v ] == d_ext iff
    # (1) C * x == d and (2) [ u  v ] == E * x + f
    C_ext = np.zeros((C.shape[0] + 2, C.shape[1] + 2))
    C_ext[:-2, :-2] = C
    C_ext[-2:, :-2] = E[:2]
    C_ext[-2:, -2:] = np.array([[-1, 0], [0, -1]])
    C_ext = cvxopt.matrix(C_ext)

    d_ext = np.zeros(d.shape[0] + 2)
    d_ext[:-2] = d
    d_ext[-2:] = -f[:2]
    d_ext = cvxopt.matrix(d_ext)

    lp_obj = cvxopt.matrix(np.zeros(A.shape[1] + 2))
    lp = lp_obj, A_ext, b_ext, C_ext, d_ext
    polygon = bretl_compute_polygon(
        lp, max_iter=max_iter, init_angle=init_angle
    )
    polygon.sort_vertices()
    vertices_list = polygon.export_vertices()
    vertices = [np.array([v.x, v.y]) for v in vertices_list]
    return vertices


def project_point_to_polytope(
    point: np.ndarray,
    ineq: Tuple[np.ndarray, np.ndarray],
    qpsolver: str,
    **kwargs,
) -> np.ndarray:
    """
    Projet a point onto a polytope in H-representation.

    Parameters
    ----------
    point :
        Point to project.
    ineq :
        Pair (`A`, `b`) describing the inequality constraint.
    qpsolver :
        Name of the backend quadratic programming solver to use, to be picked
        in ``qpsolvers.available_solvers``.

    Returns
    -------
    :
        Projected point.

    Note
    ----
    This function requires `qpsolvers <https://pypi.org/project/qpsolvers/>`_.
    """
    try:
        from qpsolvers import solve_ls
    except ImportError as e:
        raise ImportError(
            "This function requires qpsolvers: pip install qpsolvers"
        ) from e

    P = np.eye(len(point))
    return solve_ls(P, point, G=ineq[0], h=ineq[1], solver=qpsolver, **kwargs)
