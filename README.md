# Polyhedron Manipulation in Python

This library implements common operations over [convex polyhedra](https://en.wikipedia.org/wiki/Convex_polyhedron) such as [polytope projection](https://scaron.info/doc/pypoman/index.html#module-pypoman.projection), [double description](https://scaron.info/doc/pypoman/index.html#module-pypoman.duality) (conversions between halfspace and vertex representations), [Chebyshev center](https://scaron.info/doc/pypoman/index.html#chebyshev-center) computation, ... See the complete [API documentation](https://scaron.info/doc/pypoman/) for details.

## Installation

Make sure you have all system-wide dependencies by:
```
sudo apt-get install cython libglpk-dev python python-dev python-pip python-scipy
```
Then, install the module itself:
```
pip install pypoman
```

## Examples

### Vertex enumeration

We can compute the list of vertices of a polytope described in halfspace
representation by ``A * x <= b``:

```python
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
```

### Polytope projection

Let us project an n-dimensional polytope over ``x = [x_1 ... x_n]`` onto its first two coordinates ``proj(x) = [x_1 x_2]``:

```python
import pypoman
from numpy import array, eye, ones, vstack, zeros

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

vertices = pypoman.project_polytope(proj, ineq, eq, method='bretl')

if __name__ == "__main__":   # plot projected polytope
    import pylab
    pylab.ion()
    pylab.figure()
    pypoman.plot_polygon(vertices)
```

## See also

- A short introduction to [Polyhedra and
  polytopes](https://scaron.info/teaching/polyhedra-and-polytopes.html)
- Komei Fukuda's [Frequently Asked Questions in Polyhedral Computation](http://www.cs.mcgill.ca/~fukuda/soft/polyfaq/)
- The
  [Polyhedron](http://doc.sagemath.org/html/en/reference/discrete_geometry/sage/geometry/polyhedron/constructor.html) class in [Sage](http://www.sagemath.org/)
- [StabiliPy](https://github.com/haudren/stabilipy): a Python package
  implementing a more general recursive method for polytope projection
