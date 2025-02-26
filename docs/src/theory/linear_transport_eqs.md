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

```@raw html
!!! lemma
    Assume $U$ is a smooth solution of the transport equation decaying to zero at infinity for all $t \in R_+$ and $a \in C^1(\R, \R_+)$. Then, $U$ satisfies the energy bound

    $$
        \int_{\R} U^2(x, t) \dd x \leq e^{\norm{a}_{C^1}t} \int_\R U_0^2(x) \dd x.
    $$

    ??? proof
        Follows by multiplying the transport equation by $U$ and integrating over space.
        Then, use that $U$ decays to zero at infinity and apply [Grönwall's inequality](https://en.wikipedia.org/wiki/Gr%C3%B6nwall%27s_inequality).
```


The lemma shows that the energy is bounded. Using another functional, the assumptions on ``U`` can be relaxed:

```@raw html
!!! lemma
    Assume $U$ is a smooth bounded solution of the transport equation. Then, we have

    $$
        \sup_{x \in \R} \abs{U(x, t)} \leq \norm{U_0}_{L^\infty(\R)}
    $$

    for any $t > 0$. 
    
    ??? proof

        For any $(x, t) \in \R \times \R_+$, there exists $\xi \in \R$ such that $U(x, t) = U_0(\xi)$
```


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

```@example 1
u0(x) = sin(2*pi*x)
nothing # hide
```

```@setup 1
ENV["GKSwstype"]="nul" #https://discourse.julialang.org/t/deactivate-plot-display-to-avoid-need-for-x-server/19359/2
using homemade_conslaws, homemade_conslaws.Viz

function upwind_scheme(x, t, a::Function, u0::Function)
    u = zeros(length(t) + 1, length(x))
    u[1,:] = u0.(x)

    dt = t[2] - t[1]
    dx = x[2] - x[1]
    nu = dt/dx

    for n in eachindex(t)
        a_eval = a.(x, t[n])
        u[n+1,:] = nu * a_eval .* circshift(u[n,:], 1) + (1 .- nu * a_eval) .* u[n,:]
    end

    return @view u[2:end, :]
end

function central_difference_scheme(x, t, u0::Function)
    u = zeros(length(t) + 1, length(x))
    u[1,:] = u0.(x)

    dt = t[2] - t[1]
    dx = x[2] - x[1]
    nu = dt/(2*dx)

    for n in 1:length(t)
        u[n+1,:] = u[n,:] - nu * (circshift(u[n,:], -1) - circshift(u[n,:], 1))
    end

    return @view u[2:end,:]
end
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

For the first time step ``n = 1``, we have

```math
\begin{equation}
    U_j^1 = U_0(x_j) - \frac{\Delta t}{2 \Delta x} (U_0(x_{j+1}) - U_0(x_{j-1})) \label{eq:transport_first}
\end{equation}
```

Using a grid of ``50`` mesh points, simulating to time ``T = 3``, we get the following result:

```@example 1
x_L, x_R = 0, 1 # hide
T = 3 # hide
N = 100 # hide
dx = (x_R - x_L)/N # hide
dt = 1 * dx # hide

x = x_L+dx:dx:x_R # hide
t = dt:dt:T # hide

U = central_difference_scheme(x, t, u0) # hide

Viz.animate_solution(U, # hide
                     (x, t) -> u0(x-t), # hide
                     x, t, # hide
                     6) # hide
```

After some time, the solution diverges. Intuitively, the information should propagate from
left to right. However, the scheme uses information from both sides. This can be explained
rigorously using the the discrete energy

```math
    E^n = \frac{1}{2} \Delta x \sum_j (U_j^n)^2.
```

We know that the exact solution have bounded energy. We say that a scheme is *energy stable*
if ``E^n \leq E^0`` for all ``n``.

```@raw html
!!! lemma
    Let $U_j^n$ be the solutions computed with the central difference scheme. Then, we have the 
    following dicsrete energy estimate:

    $$
        E^{n+1} = E^n + \frac{\Delta x}{2} \sum_j (U_j^{n+1} - U_j^n)^2
    $$

    Thus, the energy grows at each time step for any choice of $\Delta x$ and $\Delta t$.

    ??? proof
        Similar to the proof of the continuous energy estimate and using the identity

        $$
            d_2(d_1 - d_2) = \frac{1}{2}\qty(d_1^2 - d_2^2)
        $$
```


### Upwind scheme

To respect the flow of information, we can use forward and backward differences in space
depending on the direction of propagation of information, i.e.

```math
    \frac{U_j^{n+1} - U_j^n}{\Delta t}
    + \frac{a^+ (U_j^n - U_{j-1}^n)}{\Delta x} + \frac{a^- (U_{j+1}^n - U_j^n)}{\Delta x} = 0.
```

Information is "carried by the wind", hence the name upwind. The above equations can be written as

```math
    \frac{U_j^{n+1} - U_j^n}{\Delta t} + \frac{a (U_{j+1}^n - U_{j-1}^n)}{2 \Delta x}
    = \frac{\abs{a}}{2\Delta x} (U_{j+1}^n - 2 U_j^n + U_{j-1}^n).
```

We see we have the central difference scheme with a diffusion term; the right hand side approximates
``\frac{\Delta x \abs{a}}{2} U_xx``. Thus, it adds *numerical viscosity* to the unstable central
difference scheme, which will play a crucial role later.

Now, we do the same numerical experiment as previously with ``a = 1``:

```@example 1
x_L, x_R = 0, 1 # hide
T = 1 # hide
N = 50 # hide
dx = (x_R - x_L)/N # hide
dt = 1.3 * dx # hide
x = x_L + dx:dx:x_R # hide
t = dt:dt:T # hide

a(x, t) = 1 # hide

U = upwind_scheme(x, t, a, u0) # hide


Viz.animate_solution(U, # hide
                    (x, t) -> u0(x-t), # hide
                    x, t, # hide
                    2) # hide
```

```@example 1
x_L, x_R = 0, 1 # hide
T = 1 # hide
N = 100 # hide
dx = (x_R - x_L)/N # hide
dt = 0.7 * dx # hide
x = x_L+dx:dx:x_R # hide
t = dt:dt:T # hide

a(x, t) = 1 # hide

U = upwind_scheme(x, t, a, u0) # hide

Viz.animate_solution(U, # hide
                     (x, t) -> u0(x-t), # hide
                     x, t, # hide
                     2) # hide
```

We see that stability depends on the relation ``\frac{\Delta t}{\Delta x}``. 

```@raw html
!!! lemma
    If the mesh parameters satisfy the condition

    $$\frac{\abs{a} \Delta t}{\Delta x} \leq 1,$$

    then the upwind solution satisfies the estimate

    $$E^{n+1} \leq E^n$$

    so the scheme is conditionally stable.

    ??? proof
        Start similarly to the proof of the unconditional unstability of the central difference scheme and use the mentioned identity several times.
```

The above condition is called A *CFL condition*. We also have ``L^1`` and ``L^\infty`` stability:

```@raw html
!!! lemma
    Assume the above CFL condition holds. Then, the solutions of the upwind scheme satisfy

    $$
        \norm{U^{n+1}}_{L^1} \le \norm{U^n}_{L^1}, \quad
        \norm{U^{n+1}}_{L^\infty} \le \norm{U^n}_{L^\infty}.
    $$

    ??? proof
        The first inequality follows from $U_j^{n+1}$being a convex combination of $U_{j-1}^n$,
        $U_j^n$ and $U_{j+1}^n$.
```


### Discontinuous initial data

We now consider the transport equation with ``a = 1`` in the domain ``[0, 1]`` with initial data

```math
    U_0(x) = \begin{cases}
        2, & x < 0.5, \\
        1, & x \geq 0.5.
    \end{cases}
```

```@example 1
u0(x) = x .< 0.5 ? 2 : 1 # hide
a(x, t) = 1 # hide

nothing # hide
```

This yields the discontinuous solution ``U(x, t) = U_0(x - t)``. We use the upwind scheme with ``N = 50`` and ``N = 200``
mesh points.

```@example 1
x_L, x_R = 0, 1 # hide
T = 0.25 # hide
N = 50 # hide

dx = (x_R - x_L)/N # hide
dt = 0.9 * dx # hide
x = x_L+dx:dx:x_R # hide
t = dt:dt:T # hide

U = upwind_scheme(x, t, a, u0) # hide

Viz.animate_solution(U, # hide
                     (x, t) -> u0(x-t), # hide
                     x, t, # hide
                     1) # hide
```

```@example 1
x_L, x_R = 0, 1 # hide
T = 0.25 # hide
N = 200 # hide
dx = (x_R - x_L)/N # hide
dt = 0.9 * dx # hide

x = x_L+dx:dx:x_R # hide
t = dt:dt:T # hide

U = upwind_scheme(x, t, a, u0) # hide

Viz.animate_solution(U, # hide
                     (x, t) -> u0(x-t), # hide
                     x, t, # hide
                     1) # hide
```
