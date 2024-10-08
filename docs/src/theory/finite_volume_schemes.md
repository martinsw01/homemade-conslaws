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

    $$x_{j-\frac{1}{2}} = x_j - \frac{1}{2}\Delta x = x_L + j\Delta x$$

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

$$
\begin{equation}    
    \overline F_{j+1/2}^n = \frac{1}{\Delta t}\int_{t^n}^{t^{n+1}} f(U(x_{j+1/2}, t)) \dd t
    \label{eq:flux}
\end{equation}
$$

we can rewrite the above equation as

$$
\begin{equation}
    U_j^{n+1} = U_j^n - \frac{\Delta t}{\Delta x}\qty(\overline F_{j+1/2}^n - \overline F_{j-1/2}^n).
    \label{eq:cell_average_update}
\end{equation}
$$


In other words, the change in cell averages is given by difference in fluxes across the cell boundaries. It now remains to approximate the fluxes.

### Godunov method

Godunov noticed that the cell averages are constant in each cell $\Cell_j$ at each time. We thus get a Riemann problem at each cell interface $x_{j+1/2}$:

$$
\begin{equation}
    \left\{\begin{aligned}
        U_t + f(U)_x & = 0 \\
        U(x, t^n) & = \begin{cases}
            U_j^n, & x < x_{j+1/2} \\
            U_{j+1}^n, & x > x_{j+1/2}
        \end{cases}
    \end{aligned}\right. \label{eq:riemann_problem} \tag{RP}
\end{equation}
$$

So at each time, we get a superposition of Riemann problems of the form $\eqref{eq:riemann_problem}$, at each interface. We have seen that the entropy solution is a combination of shocks, rarefactions and compound waves.

As previously discussed, the solution $\overline U_j(x,t)$ of each Riemann problem $\eqref{eq:riemann_problem}$ is self-similar:

$$
\begin{equation}
    \overline U_j(x,t) = \overline U_j\qty(\frac{x-x_{j+1/2}}{t-t^n}) = \overline U_j(\xi). \label{eq:self_similar}
\end{equation}    
$$

However, the waves intersect after some time. The maximal speed of the waves is given by $\displaystyle\max_j \abs{f'(U_j^n)}$. Thus, to ensure that two neighboring waves do not intersect, they cannot travel in total more that $\Delta x$ in time $\Delta t$. This yields the CFL-condition

$$
\begin{equation}
    \max_j \abs{f'(U_j^n)} \frac{\Delta t}{\Delta x} \leq \frac{1}{2} \label{eq:CFL_condition} \tag{CFL}
\end{equation}
$$

Assume that this condition is satisfied, and consider the curve given by $\xi = 0$, i.e. the cell interface $x_{j+1/2}$.

??? claim "Claim: $f(\overline U_j(0))$ is well defined"
    By $\eqref{eq:self_similar}$, we have that $\overline U_j$ is constant whenever $\xi$ is, and particularly at $\xi = 0$:

    $$f(U(x_{j+1/2}, t)) = f(\overline U_j(0))$$

    Along this curve, $\overline U_j$ is either continuous or not.

    === "Continuous"
        In this case, we obviously have that

        $$f(\overline U_j(0^-)) = f(\overline U_j(0^+)),$$

    === "Discontinuous"
        For $\overline U_j$ to be an entropy solution, the shock must satisfy the Rankine-Hugoniot condition

        $$
        \begin{equation}
            \jump{f(U)} = s\jump{U}. \label{eq:rankine_hugoniot} \tag{RH}
        \end{equation}
        $$

        with $s=0$. We immidiately get that

        $$f(\overline U_j(0^-)) = f(\overline U_j(0^+))$$

    so $f(\overline U_j(0)) = f(\overline U_j(0^-)) = f(\overline U_j(0^+))$ is well defined.


Now, we can define the (Riemann) edge-centered fluxes

$$
\begin{equation}
    F_{j+1/2}^n := f(\overline U_j(0)).
    \label{eq:edge_centered_flux}
\end{equation}

Since this is constant in time, we can insert into $\eqref{eq:flux}$ and explicitly calculate

$$\overline F_{j+1/2}^n = F_{j+1/2}^n.$$

Then, substituting into $\eqref{eq:cell_average_update}$, we get

$$
\begin{equation}
    U_j^{n+1} = U_j^n - \frac{\Delta t}{\Delta x}\qty(F_{j+1/2}^n - F_{j-1/2}^n).
    \label{eq:godunov}
\end{equation}
$$

### Godunov flux

It can be shown that the Riemann fluxes at the interfaces $x_{j+1/2}$ is

$$
\begin{equation}
    F_{j+1/2}^n = F(U_j^n, U_{j+1}^n)
    = \begin{cases}
        \displaystyle
        \min_{U_j^n \leq \theta \leq U_{j+1}^n} f(\theta), & U_j^n \leq U_{j+1}^n \\
        \displaystyle
        \max_{U_{j+1}^n \leq \theta \leq U_j^n} f(\theta), & U_j^n > U_{j+1}^n.
    \end{cases}
    \label{eq:godunov_flux}
\end{equation}
$$

This also holds for non-convex flux functions. Further, for flux functions $f$ with a unique minimizer $\omega$, we can simplify further to

$$
\begin{equation}
    F_{j+1/2}^n = F(U_j^n, U_{j+1}^n) = \max\qty{f\qty(\max\{U_j^n, \omega\}), f\qty(\min\{U_{j+1}^n, \omega\})}.
    \label{eq:godunov_flux_simplified}
\end{equation}
$$

### Numerical example

Now, using $\eqref{eq:godunov_flux_simplified}$, one can easily implement the Godunov method. Testing against Burgers' equation with initial conditions

$$
\begin{equation}
    U_0(x) = \begin{cases}
        1, & x < 0 \\
        0, & x \geq 0
    \end{cases} \label{eq:initial_condition_shock}
\end{equation}
$$

we know the solution is a shock wave at $x = 1/2t$.
```

```@setup 1
ENV["GKSwstype"]="nul" #https://discourse.julialang.org/t/deactivate-plot-display-to-avoid-need-for-x-server/19359/2
using homemade_conslaws.Godunov, homemade_conslaws.Viz
```

```@example 1
x_L, x_R = -1, 1
f(U) = U^2/2+1
df(U) = U
ω = 0
N = 50
dx = (x_R - x_L) / (N + 1)
x_mid = x_L + 0:dx:x_R
U0 = x_mid .< 0
BC(t, U) = [1;; 0]
dt = 0.1 # max time step
T = 1
U, t = Godunov.godunovmethod(f, df, ω, U0, BC, dx, dt, T)

Viz.animate_solution((U,),
                    "Approximate solution",
                    x_mid, t)
```

The numerical solution approxumates this well, having a sharp shock at ``x = 1/2t``. Further, it seems stable without any oscillations. Next, for the initial condition

```math
\begin{equation}
    U_0(x) = \begin{cases}
        -1, & x < 0 \\
        1, & x \geq 0
    \end{cases} \label{eq:initial_condition_rarefaction}
\end{equation}
```

we expect a rarefaction wave between the curves ``x = -t`` and ``x = t``:

```@example 1
U0 = (x_mid .> 0) * 2 .-1
BC(t, U) = [-1;; 1]
U, t = Godunov.godunovmethod(f, df, ω, U0, BC, dx, dt, T)

Viz.animate_solution((U,),
                    "Approximate solution",
                    x_mid, t)
```

For a more complicated example, ``U_0 = \sin(4\pi x)``, we expect a compound wave:

```@example 1
U0 = sin.(4π*x_mid)
BC(t, U) = [0;; 0]
U, t = Godunov.godunovmethod(f, df, ω, U0, BC, dx, dt, T)

Viz.animate_solution((U,),
                    "Approximate solution",
                    x_mid, t)
```

```@raw html
### Limitations

Despite having many disirable properties, the Godunov method has some limitations:

- It requires having an explicit formula for the solutions of the Riemann problems $\eqref{eq:riemann_problem}$. This is not always the case for systems of conservation laws.
- In the numerical flux $\eqref{eq:edge_centered_flux}$, the only required information is the flux at the interface. It may thus be unnecessary to solve an entire Riemann problem for this.
- For more complicated flux functions with multiple extremas, one cannot use $\eqref{eq:godunov_flux_simplified}$. Instead, one must solve the optimization problem $\eqref{eq:godunov_flux}$, which may be computationally expensive.


## Approximate Riemann solvers

Instead of solving $\eqref{eq:riemann_problem}$ exactly, one can approximate them. These solutions can be used to define the numerical flux $F$. Such schems are called *approximate Riemann solvers*.

### Linearized (Roe) solvers

One simple method is to linearize the conservation law $\eqref{eq:conservation_law}$:

$$f(U)_x = f'(U)U_x \approx \hat A_{j+1/2} U_x.$$

For example, one may use the *Roe average*

$$\hat A_{j+1/2} = \begin{cases}
    \frac{f(U_{j+1}^n) - f(U_j^n)}{U_{j+1}^n - U_j^n}, & U_{j+1}^n \neq U_j^n \\
    f'(U_j^n), & U_{j+1}^n = U_j^n.
\end{cases}$$

Now, we can obtain the numerical flux $F$ by replacing $\eqref{eq:riemann_problem}$ with the linearized problem

$$
\begin{equation}
    \left\{\begin{aligned}
        U_t + \hat A_{j+1/2} U_x & = 0 \\
        U(x, t^n) & = \begin{cases}
            U_j^n, & x < x_{j+1/2} \\
            U_{j+1}^n, & x > x_{j+1/2}
        \end{cases}
    \end{aligned}\right. \label{eq:linearized_riemann_problem} \tag{LRP}
\end{equation}
$$

Notice that this is just the linear transport equation with the explicit solution

$$
F_{j+1/2}^n = F^\text{Roe}(U_j^n, U_{j+1}^n) = \begin{cases}
    f(U_j^n), & \hat A_{j+1/2} \ge 0 \\
    f(U_{j+1}^n), & \hat A_{j+1/2} < 0.
\end{cases}
$$

The numerical scheme $\eqref{eq:cell_average_update}$ with this flux is called the *Roe* or *Murman-Roe* scheme. This preserves shocks well, but fails at rarefaction waves, for example our usual test case with initial data $\eqref{eq:initial_condition_rarefaction}$.

### Central schemes

Due to the linearization, we only get a single wave travelling in a single direction depending on the sign of the Roe average $\hat A$. This is also the case for shocks for the exact solution, but not for rarefaction waves. Therefore, the Roe scheme fails for rarefaction waves. Hence, instead of linearizing the conservation law, one can approximate the solution of the Riemann problem with two waves travelling in opposite directions from the interface with speeds $s_{j+1/2}^l$ and $s_{j+1/2}^r$. Different such speeds yields different schemes. The approximate solution is then

$$
U(x, t) = \begin{cases}
    U_j^n, & x < s_{j+1/2}^l t \\
    U_{j+1}^*n, & s_{j+1/2}^l t < x < s_{j+1/2}^r t \\
    U_{j+1}^n, & x > s_{j+1/2}^r t.
\end{cases}
$$

The middle state $U_{j+1/2}^*$ can be determined by using the Rankine Hugoniot condition:

$$
\begin{aligned}
    f(U_{j+1}^n) - f_{j+1/2}^* &= s_{j+1/2}^r(U_{j+1}^n - U_{j+1/2}^*) \\
    f_{j+1/2}^* - f(U_j^n) &= s_{j+1/2}^l(U_{j+1/2}^* - U_j^n),
\end{aligned}
$$

where $f_{j+1/2}^*$ is the middle flux. Now, solving for $f_{j+1/2}^*$, we get

$$f_{j+1/2}^* = \frac{s_{j+1/2}^r f(U_j^n) - s_{j+1/2}^l f(U_{j+1}^n) + s_{j+1/2}^r s_{j+1/2}^l(U_{j+1}^n - U_j^n)}{s_{j+1/2}^r - s_{j+1/2}^l}.$$

This yields the numerical flux

$$F_{j+1/2}^n = F(U_j^n, U_{j+1}^n) = f_{j+1/2}^*.$$

### Lax-Friedrichs scheme

By using the maximum speed the waves can travel without intersecting at the interface,

$$s_{j+1/2}^l = -\frac{\Delta x}{\Delta t}, \quad s_{j+1/2}^r = \frac{\Delta x}{\Delta t},$$

we get the Lax-Friedrichs flux

$$F_{j+1/2}^n = F^\text{LxF}(U_j^n, U_{j+1}^n) = \frac{f(U_j^n) + f(U_{j+1}^n)}{2} - \frac{\Delta x}{2\Delta t}(U_{j+1}^n - U_j^n).$$

The solutions are stable and non-oscillatory. Unlike the Roe scheme, it also approximates the entropy solution. However, shocks are not well preserved. See for example the usual Burgers' test case with initial data $\eqref{eq:initial_condition_shock}$:
```

```@example 1
using homemade_conslaws.LaxFriedrichs   # hide
U0 = x_mid .< 0     # hide
BC(t, U) = [1;; 0]  # hide
dt = 0.1 # max time step    # hide
T = 1   # hide
U, t = LaxFriedrichs.lax_friedrichs_method(f, df, U0, BC, dx, dt, T)    # hide

Viz.animate_solution(U, (x, t) -> x < 0.5t, x_mid, t) #hide
```

```@raw html

### Rusanov scheme

If we instead of using the maximum speeds use the speeds of propagation, we get the Rusanov scheme:

$$s_{j+1/2}^l = -s_{j+1/2}, \quad s_{j+1/2}^r = s_{j+1/2}$$

with

$$s_{j+1/2} = \max\qty{\abs{f'(U_j^n)}, \abs{f'(U_{j+1}^n)}}.$$

Then, we get the Rusanov (or Local Lax-Friedrichs) flux

$$
\begin{aligned}
    F_{j+1/2}^n
    &= F^\text{Rus}(U_j^n, U_{j+1}^n) \\
    &= \frac{f(U_j^n) + f(U_{j+1}^n)}{2}
    - \frac{\max\qty{\abs{f'(U_j^n)}, \abs{f'(U_{j+1}^n)}}}{2}(U_{j+1}^n - U_j^n).
\end{aligned}
$$

For the same problem as above with get
```

```@example 1
using homemade_conslaws.Rusanov   # hide
U, t = Rusanov.rusanov_method(f, df, U0, BC, dx, dt, T)    # hide

Viz.animate_solution(U, (x, t) -> x < 0.5t, x_mid, t) #hide
```

and with the initial data $\eqref{eq:initial_condition_rarefaction}$, we get

```@example 1
U0 = (x_mid .> 0) * 2 .-1
BC(t, U) = [-1;; 1]
U, t = Rusanov.rusanov_method(f, df, U0, BC, dx, dt, T)

Viz.animate_solution((U,), "Approximation", x_mid, t)
```