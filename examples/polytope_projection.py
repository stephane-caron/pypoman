# SPDX-License-Identifier: GPL-3.0-or-later

import pylab
from numpy import array, eye, ones, vstack, zeros

import pypoman

n = 10  # dimension of the original polytope
p = 2  # dimension of the projected polytope

# Original polytope:
# - inequality constraints: \forall i, |x_i| <= 1
# - equality constraint: sum_i x_i = 0
A = vstack([+eye(n), -eye(n)])
b = ones(2 * n)
C = ones(n).reshape((1, n))
d = array([0])
ineq = (A, b)  # A * x <= b
eq = (C, d)  # C * x == d

# Projection is proj(x) = [x_0 x_1]
E = zeros((p, n))
E[0, 0] = 1.0
E[1, 1] = 1.0
f = zeros(p)
proj = (E, f)  # proj(x) = E * x + f

vertices = pypoman.project_polytope(proj, ineq, eq, method="bretl")

if __name__ == "__main__":  # plot projected polytope
    pylab.ion()
    pylab.figure()
    pypoman.plot_polygon(vertices, resize=True)
    pylab.show(block=True)
