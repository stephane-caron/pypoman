# Polyhedron Manipulation in Python

This library implements common operations over [convex
polyhedra](https://en.wikipedia.org/wiki/Convex_polyhedron):

- Conversion between halfspace (H-rep) and vertex (V-rep) representations
- Polytope projection

## Installation

Install dependencies by:
```
sudo apt-get install cython libglpk-dev python python-dev python-pip python-scipy
CVXOPT_BUILD_GLPK=1 pip install cvxopt --user
pip install pycddlib --user
```
You can remove all ``--user`` arguments to install these Python modules system-wide.

Finally, clone this repository and run its setup script:
```
git clone https://github.com/stephane-caron/pypoman.git && cd pypoman
python setup.py build
python setup.py install --user
```

## Examples

### Polytope projection

Let us project an n-dimensional polytope over ``x = [x_1 ... x_n]`` onto its first two coordinates ``y = [x_1 x_2]``:

```python
import pylab
import pypoman
from numpy import array, eye, ones, vstack, zeros

n = 10

# Output is y = [x_0 x_1]
E = zeros((2, n))
E[0, 0] = 1.
E[1, 1] = 1.
f = zeros(2)
proj = (E, f)  # y = E * x + f

# Inequality constraints: \forall i, |x_i| <= 1
A = vstack([+eye(n), -eye(n)])
b = ones(2 * n)
ineq = (A, b)  # A * x <= b

# Equality constraint: sum_i x_i = 0
C = ones(n).reshape((1, n))
d = array([0])
eq = (C, d)  # C * x == d

vertices = pypoman.project_polytope(proj, ineq, eq, method='bretl')
pylab.ion()
pylab.figure()
pypoman.plot_polygon(vertices)
```

## See also

- Komei Fukuda's [Frequently Asked Questions in Polyhedral Computation](http://www.cs.mcgill.ca/~fukuda/soft/polyfaq/polyfaq.html)
- The
  [Polyhedron](http://doc.sagemath.org/html/en/reference/discrete_geometry/sage/geometry/polyhedron/constructor.html) class in [Sage](http://www.sagemath.org/)
- The [StabiliPy](https://github.com/haudren/stabilipy) package provides a more
  general recursive projection method
