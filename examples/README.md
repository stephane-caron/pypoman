# Examples

## Polytope projection

<img src="https://github.com/user-attachments/assets/e809b4a9-0c60-47c6-8df2-aa011e674892" height=200 align="right">

In this example, we project the polytope $\mathcal{P} = \\{ x \ | \ A x \leq b, C x = d \\}$ defined by:

```math
\begin{align*}
A & = \begin{bmatrix} I_10 \\ -I_10 \end{bmatrix} & b = 1_{20 \times 1} \\
C & = 1_{1 \times 10} & d = [0]
\end{align*}
```

We apply the projection $y = \pi(x) = E x + f$ along the first two coordinates:

```math
\begin{align*}
E & = \begin{bmatrix} I_2 \\ 0_{2 \times 8} \end{bmatrix} & f = 0_{2 \times 1}
\end{align*}
```
