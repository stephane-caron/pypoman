# Polyhedron manipulation in Python

[![Build](https://img.shields.io/github/actions/workflow/status/stephane-caron/pypoman/ci.yml?branch=main)](https://github.com/stephane-caron/pypoman/actions)
[![Documentation](https://img.shields.io/github/actions/workflow/status/stephane-caron/pypoman/docs.yml?branch=main&label=docs)](https://stephane-caron.github.io/pypoman/)
[![Coverage](https://coveralls.io/repos/github/stephane-caron/pypoman/badge.svg?branch=main)](https://coveralls.io/github/stephane-caron/pypoman?branch=main)
[![PyPI version](https://img.shields.io/pypi/v/pypoman?color=blue)](https://pypi.org/project/pypoman/)
[![PyPI downloads](https://img.shields.io/pypi/dm/pypoman?color=blue)](https://pypistats.org/packages/pypoman)

This library interfaces common operations over [convex polyhedra](https://en.wikipedia.org/wiki/Convex_polyhedron) such as [polytope projection](https://stephane-caron.github.io/pypoman/index.html#module-pypoman.projection) and [vertex enumeration](https://stephane-caron.github.io/pypoman/index.html#module-pypoman.duality). Check out the [documentation](https://stephane-caron.github.io/pypoman/) for details.

## Installation

```console
$ pip install pypoman
```

This installs the core library, covering for instance polytope projection with the Bretl method, Chebyshev centers and polygon operations. The following features rely on optional dependencies:

| Extra | Dependency | Features |
|-------|------------|----------|
| `cdd` | [pycddlib](https://pypi.org/project/pycddlib/) | Vertex and halfspace enumeration, polytope projection with `method="cdd"` |
| `all` | All of the above, plus [pyclipper](https://pypi.org/project/pyclipper/) and [qpsolvers](https://pypi.org/project/qpsolvers/) | Everything, including polygon intersection and point-to-polytope projection |

### Installing cdd

The `cdd` and `all` extras depend on pycddlib, which is built from source against [cddlib](https://github.com/cddlib/cddlib) and [GMP](https://gmplib.org/). The simplest way is to get these libraries from conda-forge:

```console
$ conda install cddlib
$ pip install pypoman[cdd]
```

Alternatively, install them as system packages, for instance on Debian-based Linux distributions:

```console
$ sudo apt-get install libcdd-dev libgmp3-dev
$ pip install pypoman[cdd]
```

Replace `pypoman[cdd]` with `pypoman[all]` to enable all optional features.

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
