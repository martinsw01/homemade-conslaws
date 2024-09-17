```@raw html
# Finite volume schemes

We will now design efficient schemes for the conservation law


$$
\begin{equation}
    U_t + f(U)_x = 0 \label{eq:conservation_law}
\end{equation}
$$


We have already seen that finite difference schemes works poorly, even for linear transport:


- In this case, we had to "upwind" wind scheme, taking the derivative in the direction of  information propagation. This is not possible a-priori for non-linear equations like $\eqref{eq:conservation_law}$.
- It further requires smoothness and that $\eqref{eq:conservation_law}$ is satisfied pointwise. However, we are looking for entropy solutions, which may be discontinuous.

## Finite volume schemes

### The grid


For simplicity, we use uniform discretization of $[x_L, x_R]$:


- points:

    $$
    \begin{aligned}
        x_j &= x_L + \qty(j+\frac{1}{2})\Delta x \\
        \Delta x &= \frac{x_R - x_L}{N+1}
    \end{aligned}
    $$

    for $j=0,\dots, N$

- midpoints:

    $$x_{j-\frac{1}{2}} = x_L - \frac{1}{2}\Delta x = x_L + j\Delta x$$

    for $j=0,\dots, N+1$

- control volumes:

    $$\Cell_j = [x_{j-1/2}, x_{j+1/2}]$$
    
    for $j=1,\dots, N$

- uniform discretization of time:

    $$t^n = n\Delta t$$

### Cell averages

The finite difference schemes approximates the point values, whereas the finite volume schemes aims to approximate the cell averages

$$U_j^n \approx \frac{1}{\Delta x}\int_{\Cell_j} U(x, t^n) \dd x$$

starting with

$$U_j^0 = \frac{1}{\Delta x}\int_{\Cell_j} U_0(x) \dd x$$

### Integral form

Assuming $U_j^n$ is known for som time $t^n$, we can integrate $\eqref{eq:conservation_law}$ over the swuare $\Cell_j \times [t^n, t^{n+1})$ to get the next averages:

$$\int_{t^n}^{t^{n+1}} \int_{\Cell_j} U_t + f(U)_x \dd x \dd t = 0.$$

Defining the flux $\overline F_{j+1/2}^n$ across the boundary $x_{j+1/2}$ over the time interval $[t^n, t^{n+1})$:

$$\overline F_{j+1/2}^n = \frac{1}{\Delta t}\int_{t^n}^{t^{n+1}} f(U(x_{j+1/2}, t)) \dd t$$

we can rewrite the above equation as

$$U_j^{n+1} = U_j^n + \frac{\Delta t}{\Delta x}\qty(\overline F_{j+1/2}^n - \overline F_{j-1/2}^n).$$

In other words, the change in cell averages is given by difference in fluxes across the cell boundaries. It now remains to approximate the fluxes.
```