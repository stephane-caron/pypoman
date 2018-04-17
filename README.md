# Convex Polyhedron Manipulation in Python

- Conversion between halfspace (H-rep) and vertex (V-rep) representations
- Polytope projection

## Installation

Install dependencies by:
```
sudo apt-get install cython libglpk-dev python python-dev python-pip python-scipy
CVXOPT_BUILD_GLPK=1 pip install cvxopt --user
pip install pycddlib --user
```
Remove ``--user`` arguments to install Python modules system-wide.

Finally, clone this repository and run the setup script:
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
