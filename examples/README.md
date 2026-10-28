# Examples

You can run examples from the root of the repository using [uv](https://docs.astral.sh/uv/):

```console
$ uv run examples/polytope_projection.py
```

Some examples rely on optional dependencies, which you can enable as extras. See the readme for installation instructions.

## Polytope projection

<img src="https://github.com/user-attachments/assets/e809b4a9-0c60-47c6-8df2-aa011e674892" height=300 align="right">

In the `polytope_projection.py` example, we project the polytope $\mathcal{P} = \\{ x \ | \ A x \leq b, C x = d \\}$ defined by:

```math
\begin{align*}
A & = \begin{bmatrix} I_10 \\ -I_10 \end{bmatrix} & b = 1_{20 \times 1} \\
C & = 1_{1 \times 10} & d = [0]
\end{align*}
```

Our projection $\pi(\mathcal{P})$, defined by \$\pi(x) = E x + f$, is the projection on the first two coordinates:

```math
\begin{align*}
E & = \begin{bmatrix} I_2 \\ 0_{2 \times 8} \end{bmatrix} & f = 0_{2 \times 1}
\end{align*}
```

The output of the example is a list of vertices $V = \\{ v_i \\}$ such that $\pi(\mathcal{P}) = \mathrm{conv}(V)$.
