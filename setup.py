#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2018 Stephane Caron <stephane.caron@lirmm.fr>
#
# This file is part of pypoman.
#
# pypoman is free software: you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# pypoman is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with pypoman. If not, see <http://www.gnu.org/licenses/>.

from distutils.core import setup

classifiers = """\
Development Status :: 5 - Production/Stable
License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)
Intended Audience :: Developers
Intended Audience :: Science/Research
Topic :: Scientific/Engineering :: Mathematics
Programming Language :: Python
Programming Language :: Python :: 2
Programming Language :: Python :: 3
Operating System :: OS Independent"""

links = [
    'https://en.wikipedia.org/wiki/Convex_polyhedron',
    'https://scaron.info/doc/pypoman/index.html#module-pypoman.projection',
    'https://scaron.info/doc/pypoman/index.html#module-pypoman.duality',
    'https://scaron.info/doc/pypoman/index.html#chebyshev-center',
    'https://scaron.info/doc/pypoman/']

packages = ["cython libglpk-dev python python-dev python-pip python-scipy"]

long_description = """\
This library implements common operations over `convex polyhedra <%s>`_:

- `Polytope projection <%s>`_
- `Duality <%s>`_: conversion between halfspace and vertex representations
- `Chebyshev center <%s>`_ calculation

See the complete `API documentation <%s>`_ for details.

Installation
------------

First, install all dependencies:

.. code-block::

    sudo apt-get install %s
    CVXOPT_BUILD_GLPK=1 pip install cvxopt --user
    pip install pycddlib --user

Then, install the module itself:

.. code-block::

    pip install pypoman --user

You can remove the ``--user`` options to install modules system-wide.

Examples
--------

Vertex enumeration
~~~~~~~~~~~~~~~~~~

We can compute the list of vertices of a polytope described in halfspace
representation by ``A * x <= b``:

.. code:: python

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

Polytope projection
~~~~~~~~~~~~~~~~~~~

Let us project an n-dimensional polytope over ``x = [x_1 ... x_n]`` onto its
first two coordinates ``proj(x) = [x_1 x_2]``:

.. code:: python

    import pypoman
    from numpy import array, eye, ones, vstack, zeros

    n = 10  # dimension of the original polytope
    p = 2   # dimension of the projected polytope

    # Original polytope:
    # - inequality constraints: \\forall i, |x_i| <= 1
    # - equality constraint: sum_i x_i = 0
    A = vstack([+eye(n), -eye(n)])
    b = ones(2 * n)
    C = ones(n).reshape((1, n))
    d = array([0])
    ineq = (A, b)  # A * x <= b
    eq = (C, d)    # C * x == d

    # Projection is proj(x) = [x_0 x_1]
    E = zeros((p, n))
    E[0, 0] = 1.
    E[1, 1] = 1.
    f = zeros(p)
    proj = (E, f)  # proj(x) = E * x + f

    vertices = pypoman.project_polytope(proj, ineq, eq, method='bretl')

    if __name__ == "__main__":   # plot projected polytope
        import pylab
        pylab.ion()
        pylab.figure()
        pypoman.plot_polygon(vertices)
""" % tuple(links + packages)

setup(
    name='pypoman',
    version='0.5.0',
    description="Python Polyhedron Manipulation",
    long_description=long_description,
    url="https://github.com/stephane-caron/pypoman",
    author="Stéphane Caron",
    author_email="stephane.caron@lirmm.fr",
    license="LGPL",
    keywords="convex, polyhedron, polytope, projection, duality",
    classifiers=classifiers.split('\n'),
    packages=['pypoman']
)
