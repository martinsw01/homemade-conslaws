# Linear transport equations

The simplest case of the one-dimensional version of the transport equation,

```math
\begin{equation}
    U_t + a(x, t) U_x = 0 \quad \text{ in } \R \times \R_+, \label{eq:transport}
\end{equation}
```

is when the velocity field ``a`` is constant. The Cauchy problem is solved using the method of characteristics, where the solution is constant along the characteristic curves ``x = x_0 + a t``. The solution is then

```math
\begin{equation}
    U(x, t) = U_0(x - a t) \label{eq:transport_const}
\end{equation}
```

for any ``(x,t) \in \R \times \R_+``. We see that the initial data is transported with the velocity ``a``. For more general cases, the characteristics may not be possible to solve explicitly. However, we can obtain some information of the sulutions with the following *a priori* energy estimate:

!!! lemma
    Assume ``U`` is a smooth solution of the transport equation decaying to zero at infinity for all ``t \in R_+`` and ``a \in C^1(\R, \R_+)``. Then, ``U`` satisfies the energy bound

    ```math
        \int_{\R} U^2(x, t) \dd x \leq e^{\norm{a}_{C^1}t} \int_\R U_0^2(x) \dd x.
    ```

    !!! proof
        Follows by multiplying the transport equation by ``U`` and integrating over space.
        Then, use that ``U`` decays to zero at infinity and apply [Grönwall's inequality](https://en.wikipedia.org/wiki/Gr%C3%B6nwall%27s_inequality).


The lemma shows that the energy is bounded. Using another functional, the assumptions on ``U`` can be relaxed:

!!! lemma
    Assume ``U`` is a smooth bounded solution of the transport equation. Then, we have

    ```math
        \sup_{x \in \R} \abs{U(x, t)} \leq \norm{U_0}_{L^\infty(\R)}
    ```

    for any ``t > 0``. 
    
    !!! proof

        For any ``(x, t) \in \R \times \R_+``, there exists ``\xi \in \R`` such that ``U(x, t) = U_0(\xi)``


## Finite difference schemes for the transport equation

For some velocity fields $a(x,t)$, it may not be possible to derive an explicit formula for the characteristic equation. We thus use numerical methods to approximate the solutions of the transport equation.

### Discretization

For simplicity, assume that the velocity field is positive. AS $\R$ is unbounded, we truncate the domain into $\Omega = [x_L, X_R]$. Thus, we must impose boundary conditions, which will be discussed below. For simplicity, we use a uniform mesh of mesh size $\Delta x$ with $N+1$ points $x_j$:

```math
    x_L = x_0 < x_1 < \cdots < x_N = x_R
```

We further choose some terminal time $T$ and divide into $M+1$ points $t^n = n \Delta t$. We set the initial approximaition $U_j^0 := U_0(x_j)$ and update the next approximation $U_j^{n+1}$ using a finite difference scheme.

### Centered finite difference scheme

One such scheme is the forward difference in time and cetral difference in space approximating $(\ref{eq:transport_const})$:

```math
    \frac{U_j^{n+1} - U_j^n}{\Delta t} + \frac{a(U_{j+1}^n - U_{j-1}^n)}{2 \Delta x} = 0,
    \quad 0 < j < N.
```

### Example

We consider the domain $[0, 1]$ with initial data

```math
    U_0(x) = \sin(2\pi x)
```

and $a = 1$. Since the data is periodic, we impose periodic boundary conditions. Numerically, we implement this by setting


```math
    U_{0}^n = U_N^n, \quad U_{N+1}^n = U_1^n
```

Thus, on the boundary, $j = 0, N$, we have

```math
\begin{align*}
    \frac{U_1^{n+1} - U_1^n}{\Delta t} + \frac{a(U_2^n - U_N^n)}{2 \Delta x} = 0, \\
    \frac{U_N^{n+1} - U_N^n}{\Delta t} + \frac{a(U_1^n - U_{N-1}^n)}{2 \Delta x} = 0.
\end{align*}
```

In addition, reshape $U_j^n$ into a column vector $U_j^n = U_{j+nN}$. Then, then above equations can be formulated as $A \bm U = \bm f$.