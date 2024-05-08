# Polyhedron manipulation in Python

[![Build](https://img.shields.io/github/actions/workflow/status/stephane-caron/pypoman/ci.yml?branch=main)](https://github.com/stephane-caron/pypoman/actions)
[![Coverage](https://coveralls.io/repos/github/stephane-caron/pypoman/badge.svg?branch=main)](https://coveralls.io/github/stephane-caron/pypoman?branch=main)
[![Documentation](https://img.shields.io/badge/docs-online-brightgreen?logo=read-the-docs&style=flat)](https://scaron.info/doc/pypoman/)
[![PyPI version](https://img.shields.io/pypi/v/pypoman)](https://pypi.org/project/pypoman/)
[![PyPI downloads](https://pepy.tech/badge/pypoman/month)](https://pepy.tech/project/pypoman)

This library allows common operations over [convex polyhedra](https://en.wikipedia.org/wiki/Convex_polyhedron) such as [polytope projection](https://scaron.info/doc/pypoman/index.html#module-pypoman.projection) and [vertex enumeration](https://scaron.info/doc/pypoman/index.html#module-pypoman.duality). Check out the [API documentation](https://scaron.info/doc/pypoman/) for details.

## Installation

Install system packages for Python and GLPK, for instance for Debian-based Linux distributions:

```console
$ sudo apt-get install cython libglpk-dev python python-dev python-pip
```

Then, install the library by:

```console
$ pip install pypoman
```

Some functions, such as point-polytope projection and polygon intersection, are optional and not installed by default. To enable all of them, run:

```console
$ pip install pypoman[all]
```

## Examples

### Vertex enumeration

We can compute the list of vertices of a polytope described in halfspace representation by $A x \leq b$:

```python
import numpy as np
from pypoman import compute_polytope_vertices

A = np.array([
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
b = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 2, 2, 1, 2, 3])

vertices = compute_polytope_vertices(A, b)
```

### Halfspace enumeration

The other way round, assume we know the vertices of a polytope, and want to get its halfspace representation $A x \leq b$.

```python
import numpy as np
from pypoman import compute_polytope_halfspaces

vertices = map(
    np.array,
    [[1, 0, 0], [0, 1, 0], [1, 1, 0], [0, 0, 1], [0, 1, 1]],
)

A, b = compute_polytope_halfspaces(vertices)
```

### Polytope projection

Let us project an $n$-dimensional polytope $A x \leq b$ over $x = [x_1\ \ldots\ x_n]$ onto its first two coordinates $proj(x) = [x_1 x_2]$:

```python
from numpy import array, eye, ones, vstack, zeros
from pypoman import plot_polygon, project_polytope

n = 10  # dimension of the original polytope
p = 2   # dimension of the projected polytope

# Original polytope:
# - inequality constraints: \forall i, |x_i| <= 1
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

vertices = project_polytope(proj, ineq, eq, method='bretl')
```

We can then plot the projected polytope:

```python
import pylab

pylab.ion()
pylab.figure()
plot_polygon(vertices)
```

## See also

- A short introduction to [Polyhedra and polytopes](https://scaron.info/blog/polyhedra-and-polytopes.html)
- Komei Fukuda's [Frequently Asked Questions in Polyhedral Computation](https://www.inf.ethz.ch/personal/fukudak/polyfaq/polyfaq.html)
- The [Polyhedron](http://doc.sagemath.org/html/en/reference/discrete_geometry/sage/geometry/polyhedron/constructor.html) class in [Sage](http://www.sagemath.org/)
- [StabiliPy](https://github.com/haudren/stabilipy): a Python package implementing a more general recursive method for polytope projection
