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
using homemade_conslaws.FiniteVolumes: godunov_scheme, lax_friedrichs_scheme, rusanov_scheme
using homemade_conslaws
using homemade_conslaws.Viz: animate_solution, animate_solutions
```

```@example 1
x_L, x_R = -1, 1
N = 50
u0(x) = x < 0 ? 1. : 0.

grid = UniformGrid1D(N, NeumannBC(), u0, (x_L, x_R))
system = ConservedSystem(BurgersEQ(), NoReconstruction(grid), GodunovFlux(0.), ForwardEuler(grid))
simulator = Simulator(system, grid, 0.)

x_mid = cell_centers(grid)

dt = 0.1 # max time step
T = 1
U, t = simulate_and_aggregate!(simulator, T, dt)

animate_solution(U', (x, t) -> u0(x - 0.5t),
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
u0(x) = (x > 0) * 2 - 1
grid = UniformGrid1D(N, NeumannBC(), u0, (x_L, x_R))
simulator = Simulator(system, grid, 0.)
U, t = simulate_and_aggregate!(simulator, T, dt)

animate_solution(U',
                 "Approximate solution",
                 x_mid, t)
```

For a more complicated example, ``U_0 = \sin(4\pi x)``, we expect a compound wave:

```@example 1
u0(x) = sin(4π*x)
grid = UniformGrid1D(N, PeriodicBC(), u0, (x_L, x_R))
simulator = Simulator(system, grid, 0.)
U, t = simulate_and_aggregate!(simulator, T, dt)
# U, t = godunov_scheme(f, df, ω, U0, BC, dx, dt, T)

animate_solution(U',
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
\begin{equation}
    F_{j+1/2}^n = F^\text{Roe}(U_j^n, U_{j+1}^n) = \begin{cases}
        f(U_j^n), & \hat A_{j+1/2} \ge 0 \\
        f(U_{j+1}^n), & \hat A_{j+1/2} < 0.
    \end{cases}
    \label{eq:roe_flux}
\end{equation}
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
u0(x) = x < 0     # hide
grid = UniformGrid1D(N, NeumannBC(), u0, (x_L, x_R)) # hide
system = ConservedSystem(BurgersEQ(), NoReconstruction(grid), LaxFriedrichsFlux(), ForwardEuler(grid)) # hide
simulator = Simulator(system, grid, 0.) # hide
dt = 0.1 # max time step    # hide
T = 1   # hide
U, t = simulate_and_aggregate!(simulator, T, dt)    # hide

animate_solution(U', (x, t) -> x < 0.5t, x_mid, t) #hide
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
U0 = u0.(x_mid) # hide
ω = 0. # hide
f(u) = 0.5u^2 # hide
df(u) = u # hide
dx = grid.dx # hide
BC(t, U) = [U[1]; U[end]] # hide

U, t = rusanov_scheme(f, df, U0, BC, dx, dt, T)    # hide

animate_solution(U, (x, t) -> x < 0.5t, x_mid, t) #hide
```

and with the initial data $\eqref{eq:initial_condition_rarefaction}$, we get

```@example 1
U0 = (x_mid .> 0) * 2 .-1 # hide
BC(t, U) = [-1;; 1] # hide
U, t = rusanov_scheme(f, df, U0, BC, dx, dt, T) # hide

animate_solution(U, "Approximation", x_mid, t) # hide
```

### Comparison

```@example 1
U0 = x_mid .< 0     # hide
BC(t, U) = [1;; 0] # hide

U_Rus, t = rusanov_scheme(f, df, U0, BC, dx, dt, T) # hide
U_LxF, _ = lax_friedrichs_scheme(f, df, U0, BC, dx, dt, T) # hide
U_God, _ = godunov_scheme(f, df, ω, U0, BC, dx, dt, T) # hide

animate_solutions((U_Rus, U_LxF, U_God, x_mid' .< 0.5t), # hide
                  ["Rusanov" "Lax-Friedrichs" "Godunov" "Exact"], # hide
                  x_mid, t) # hide
```


```@raw html
## Consistent, conservative and monotone schemes

### Conservative schemes

A numerical scheme approximating $\eqref{eq:conservation_law}$ can be formulated as

$$
\begin{equation}
    U_j^{n+1} = H(U_{j-p}^n, \dots, U_{j+p}^n) \label{eq:general_scheme}
\end{equation}
$$

for some update function $H$ depending on the $2p+1$ points stencil $\{U_{j-p}^n, \dots, U_{j+p}^n\}$ for the scheme. Until now, we have worked with $3$-stencil schemes, i.e. $p=1$. 

!!! definition "Definition (Conservative shceme)"
    The generic scheme $\eqref{eq:general_scheme}$ approximating $\eqref{eq:conservation_law}$ is called *conservative* if

    $$\sum_j U_j^{n+1} = \sum_j U_j^n,$$

    ignoring the boundary conditions.


!!! theorem
    Assume $H(0, \dots, 0) = 0$. 

    $$
    \begin{gather*}
        \eqref{eq:general_scheme} \text{ is conservative } \\
        \Updownarrow \\
        \text{There is a function } F_{j+1/2}^n = F(U_{j-p}^n, \dots, U_{j+p}^n) \text{ such that \eqref{eq:general_scheme} can be written as } \eqref{eq:cell_average_update}.
    \end{gather*}
    $$

    ??? proof
        === "$\implies$"
            Define

            $$G(U_{-p}, \dots, U_p) := \frac{\Dx}{\Dt}(U_0 - H(U_{-p}, \dots, U_p)).$$

            By conservation, we have that

            $$\sum_j G(U_{j-p}^n, \dots, U_{j+p}^n) = \frac{\Dx}{\Dt}\sum_j (U_j^n - H(U_{j-p}^n, \dots, U_{j+p}^n)) = 0$$

            for any sequence $\{U_j^n\}_j$.
            
            === "$F_{i+1/2}$"
                Now, define the sequence

                $$V_j^n :=  \begin{cases}
                    0, & j -i \le -p \lor j-i < p \\
                    U_j^n, & \text{otherwise}.
                \end{cases}$$

                Then, we have

                $$
                \begin{aligned}
                    0 &= \sum_j G(V_{j-p}^n, \dots, V_{j+p}^n) \\
                    &= \underbrace{G(0, \dots, 0, U_{1-p-i}^n) + \dots + G(0, U_{1-p-i}^n, \dots, U_{p-i}^n)}_{F_{i+1/2}} \\
                    &\quad + \underbrace{G(U_{1-p-i}^n, \dots, U_{p-i}^n, 0) + \dots + G(U_{-p-i}^n, 0, \dots, 0)}_B
                \end{aligned}
                $$

            === "$F_{i-1/2}$"
                Now, define the sequence

                $$W_j^n :=  \begin{cases}
                    0, & j -i < -p \lor j-i < p \\
                    U_j^n, & \text{otherwise}.
                \end{cases}$$

                Then, we have

                $$
                \begin{aligned}
                    0 &= \sum_j G(W_{j-p}^n, \dots, W_{j+p}^n) \\
                    &= \underbrace{G(0, \dots, 0, U_{-p-i}^n) + \dots + G(0, U_{-p-i}^n, \dots, U_{p-i}^n)}_{F_{i-1/2}} \\
                    &\quad + G(U_{-p-i}^n, \dots, U_{p-i}^n) \\
                    &\quad + \underbrace{G(U_{1-p-i}^n, \dots, U_{p-i}^n, 0) + \dots + G(U_{-p-i}^n, 0, \dots, 0)}_B
                \end{aligned}
                $$

            ---

            This results in the form $\eqref{eq:godunov}$:

            $$
            \begin{aligned}
                F_{i+1/2} + B &= F_{i-1/2} + G(U_{-p-i}^n, \dots, U_{p-i}^n) + B \\
                F_{i+1/2} - F_{i-1/2} &= \frac{\Dx}{\Dt}(U_i^n - H(U_{i-p}^n, \dots, U_{i+p}^n)) \\
                U_i^{n+1} &= U_i^n - \frac{\Dt}{\Dx}(F_{i+1/2} - F_{i-1/2}).
            \end{aligned}
            $$

        === "$\impliedby$"
            Then, we have

            $$U_j^{n+1} = H(U_{j-p}^n, \dots, U_{j+p}^n) = U_j^n - \frac{\Delta t}{\Delta x}\qty(F_{j+1/2}^n - F_{j-1/2}^n).$$

            This yields a telescoping sum

            $$
            \begin{aligned}
                \sum_j U_j^{n+1} &= \sum_j \qty[U_j^n - \frac{\Delta t}{\Delta x}\qty(F_{j+1/2}^n - F_{j-1/2}^n)] \\
                &= \sum_j U_j^n
            \end{aligned}
            $$



### Consistent schemes

We now consider numerical fluxes $F_{j+1/2}^n$ defined on $2p+1$ stencils

$$
\begin{equation}
    F_{j+1/2}^n = F(U_{j-p}^n, \dots, U_{j+p}^n) \label{eq:numerical_flux}
\end{equation}
$$

!!! definition "Definition (Consistency)"
    The numerical flux $\eqref{eq:numerical_flux}$ is *consistent* if $F(U, \dots, U) = f(U)$ for all $U \in \R$.

!!! example
    All the numerical fluxes disussed above are consistent. Take for example Lax-Friedrichs:

    $$F^\text{LxF}(U, U) = \frac{f(U) + f(U)}{2} - \frac{\Delta x}{2\Delta t}(U - U) = f(U).$$

However, conservation and consistency does not imply stability or convergence. Take for example the Roe flux $\eqref{eq:roe_flux}$, which is trivially consistent but does not converge for certain problems.


### Monotone schemes

Recall that entropy solutions to the conservation law $\eqref{eq:conservation_law}$ are monotinicity preserving. It would therefore be desirable to have numerical schemes that have a discrete version of this property.

!!! definition "Definition (Monotone scheme)"
    The numericals scheme $\eqref{eq:general_scheme}$ is *monotone* if it is non-decreasing in each argument, i.e. $(\nabla H)_i \ge 0$ for all $i$.


!!! lemma
    Let $F$ be a locally Lipschitz continuous two-point numerical flux. Then we have

    $$\eqref{eq:godunov} \text{ is monotone} \iff \begin{aligned}
        a & \mapsto F(a, b) \text{ is non-decreasing for fixed } b \\
        b & \mapsto F(a, b) \text{ is non-increasing for fixed } a.
    \end{aligned}$$

    and we get the following CFL type condition

    $$
    \begin{equation}
        \abs{\pdv{F}{a}(v,w)} + \abs{\pdv{F}{b}(u,v)} \le \frac{\Delta x}{\Delta t} \quad \forall u, v, w. \label{eq:monotone_CFL}
    \end{equation}
    $$


One cannowuse monotoicity to differentiate between robust and possibly non-robust schemes. For example, the Roe scheme is not monotone.

## Stability properties of monotone schemes

Recall that the entropy solutions of $\eqref{eq:conservation_law}$ have

1. $L^\infty$ bound
2. $L^p$ bounds
3. $TV$ bound

and satisfy time constinuity. We will show that monotone schemes have these properties.

### L<sup>∞</sup> bound

!!! lemma
    Let $U_j^n$ be the approximate solution using a consistent, monotone scheme of the form $\eqref{eq:general_scheme}$. Then we have

    $$\min\{U_{j-p}^n, \dots, U_{j+p}^n\} \le U_j^{n+1} \le \max\{U_{j-p}^n, \dots, U_{j+p}^n\},$$

    and in particular

    $$\min_i U_j^0 \le U_j^n \le \max_i U_j^0 \quad \forall n,j$$

    ??? proof
        Let $\overline U_j^n = \max\{U_{j-p}^n, \dots, U_{j+p}^n\}$ be the maximum of the stencil. Now, as $\eqref{eq:general_scheme}$ is monotone, we have

        $$
        \begin{aligned}
            U_j^{n+1} &= H(U_{j-p}^n, \dots, U_{j+p}^n) \\
            &\le H(\overline U_{j-p}^n, U_{j-p+1}^n, \dots, U_{j+p}^n) \\
            &\vdots \\
            &\le H(\overline U_{j-p}^n, \dots, \overline U_{j+p}^n) \\
            &= \overline U_j^n.
        \end{aligned}
        $$

        Similarly, we get the minimum principle.

        Propagating this through time, we get the last claim too.


### Entropy-inequalities and Lᵖ bounds

In the continuous case, we used the entropy inequality to obtain $L^p$ bounds. We will use discrete counterpart to the Kruzkov entropy inequality. The *Crandall-Majda numerical entropy flux* is given by

$$Q_{j+1/2}^n = Q(U_j^n, U_{j+1}^n) = F(U_j^n \lor k, U_{j+1}^n \lor k) - F(U_j^n \land k, U_{j+1}^n \land k),$$

where $a \lor b = \max\{a, b\}$ and $a \land b = \min\{a, b\}$. For consistent numerical fluxes, we see that this entropy flux is consistent with the Kruskov entropy flux:

$$
\begin{aligned}
    Q(U, U) &= f(U \lor k) - f(U \land k) \\
    &= \sign(U - k)\qty[f(U) - f(k)] \\
    &= q(U; k).
\end{aligned}
$$

!!! lemma "Lemma (Crandall-Majda [CM80])"
    Let $U_j^n$ be the approximate solution using a consistent, conservative, monotone scheme of the form $\eqref{eq:general_scheme}$. Then we have the discrete entropy inequality

    $$
    \begin{equation}
        \abs{U_j^{n+1} - k}  - \abs{U_j^n - k} + \frac{\Delta t}{\Delta x}\qty(Q_{j+1/2}^n - Q_{j-1/2}^n) \le 0 \quad \forall n,j. \label{eq:discrete_entropy_inequality}
    \end{equation}
    $$

    Moreover, if $U_0 \in L^1(\R)$, then

    $$\sum_j \abs{U_j^n} \Delta x \le \norm{U_0}_{L^1(\R)}.$$

    ??? proof
        The scheme is conservative, so we have

        $$
        \begin{aligned}
            H(U_{j-p}^n \lor k, &\dots, U_{j+p}^n \lor k) \\
            & = U_j^n + \frac{\Delta t}{\Delta x}\qty(F(U_{j-p}^n \lor k, \dots, U_{j+p}^n \lor k) - F(U_{j-p}^n \land k, \dots, U_{j+p}^n \land k)). \\
        \end{aligned}
        $$

        and similarly for the minimum. Taking their difference, we get

        $$
        \begin{aligned}
            H(U_{j-p}^n \lor k, &\dots, U_{j+p}^n \lor k)
            - H(U_{j-p}^n \land k, \dots, U_{j+p}^n \land k) \\
            & = U_j^n \lor k - U_j^n \land k
            - \frac{\Delta t}{\Delta x}\qty(Q_{j+1/2}^n - Q_{j-1/2}^n) \\
            &= \abs{U_j^{n+1} - k}  - \frac{\Delta t}{\Delta x}\qty(Q_{j+1/2}^n - Q_{j-1/2}^n)
        \end{aligned}
        $$

        Further, by monotonicity, we have

        $$H(U_{j-p}^n \lor k, \dots, U_{j+p}^n \lor k) \ge \begin{matrix}
            U_j^{n+1} \\ k
        \end{matrix}
        \ge H(U_{j-p}^n \land k, \dots, U_{j+p}^n \land k),$$

        so inserting into the above equation, we get the first claim.

        Now, let $k=0$ and sum over $j$ to get the second claim.


### TV bound

If we interpret the cell averages as step functions, then the total variation reduces to the sum of the jumps:

$$\norm{U^n}_{TV(\R)} = \sum_j \abs{U_j^n - U_{j-1}^n}.$$

We can rewrite the finite volume scheme $\eqref{eq:godunov}$ in the incremental form

$$
\begin{equation}
    U_j^{n+1} = U_j^n + C_{j+1/2} \qty(U_{j+1}^n - U_j^n) - D_{j-1/2} \qty(U_j^n - U_{j-1}^n)
    \label{eq:incremental_form}
\end{equation}
$$

where

$$C_{j+1/2}^n = \frac{\Delta t}{\Delta x} \frac{f(U_j^n) - F_{j+1/2}^n}{U_{j+1}^n - U_j^n},
\quad D_{j+1/2}^n = \frac{\Delta t}{\Delta x} \frac{f(U_{j+1}^n) - F_{j+1/2}}{U_{j+1}^n - U_j^n}.$$

Then, wecanuse the following result:

!!! lemma "Hartem's Lemma [Har83]"
    Let $U_j^n$ be computed with $\eqref{eq:incremental_form}$.

    1. If $C_{j+1/2}^n, D_{j+1/2}^n \ge 0$ and $C_{j+1/2}^n + D_{j+1/2}^n \le 1$ for all $n, j$, then

        $$\norm{U^{n+1}}_{TV(\R)} \le \norm{U^n}_{TV(\R)}.$$

    2. If $C_{j+1/2}^n, D_{j-1/2}^n \ge 0$ and $C_{j+1/2}^n + D_{j-1/2}^n \le 1$ for all $n, j$, then

        $$\norm{U^{n+1}}_{L^\infty} \le \norm{U^n}_{L^\infty}.$$

    ??? proof
        === "1."
            Using $\eqref{eq:incremental_form}$, we get the difference

            $$
            \begin{aligned}
                U_{j+1}^{n+1} - U_j^{n+1}
                \le& \qty(1-C_{j+1/2}^n - D_{j+1/2}^n)\qty(U_{j+1}^n - U_j^n) \\
                &+ C_{j+3/2}^n \qty(U_{j+2}^n - U_{j+1}^n) \\
                &+ D_{j-1/2}^n \qty(U_j^n - U_{j-1}^n).
            \end{aligned}
            $$

            Taking the absolute value, we get the inequality

            $$
            \begin{aligned}
                \abs{U_{j+1}^{n+1} - U_j^{n+1}}
                &\le \qty(1-C_{j+1/2}^n - D_{j+1/2}^n)\abs{U_{j+1}^n - U_j^n} \\
                &\quad + C_{j+3/2}^n \abs{U_{j+2}^n - U_{j+1}^n} \\
                &\quad + D_{j-1/2}^n \abs{U_j^n - U_{j-1}^n}.
            \end{aligned}
            $$

            Summing over $j$, we get a telescoping sum, resulting in

            $$\sum_j \abs{U_{j+1}^{n+1} - U_j^{n+1}} \le \sum_j \abs{U_{j+1}^n - U_j^n}.$$

        === "2."
            Rewriting $\eqref{eq:incremental_form}$ as

            $$U_j^{n+1} = C_{j+1/2}^n U_{j+1}^n
            + \qty(1 - C_{j+1/2}^n - D_{j-1/2}^n) U_j^n
            + D_{j-1/2}^n U_{j-1}^n,$$

            we see that $U_j^{n+1}$ is a convex combination of $U_{j-1}^n, U_j^n, U_{j+1}^n$. Thus, $U_j^{n+1}$ is bounded by the maximum and minimum of these values, so we get the claim.


It turns out that consistent, conservative, monotone schemes satisfies these criteria, so we can use this lemma to show that the schemes are TVD.

!!! lemma
    Monotone, consistent, conservative 3-point schemes are TVD under the CFL condition $\eqref{eq:monotone_CFL}$.

    ??? proof
        The coefficients are positive:
        
        $$
        \begin{align*}
            C_{j+1/2}^n &= \frac{\Delta t}{\Delta x} \frac{f(U_j^n) - F_{j+1/2}^n}{U_{j+1}^n - U_j^n} \\
            &= \frac{\Delta t}{\Delta x} \frac{F(U_j^n, U_j^n) - F_{j+1/2}^n}{U_{j+1}^n - U_j^n}
            \tag*{(Consistency)} \\
            &\ge 0
            \tag*{(Monotonicity)}
        \end{align*}
        $$

        Similarly for $D_{j-1/2}^n$. Further, we have

        $$
        \begin{align*}
            C_{j+1/2}^n &= \frac{\Delta t}{\Delta x} \frac{F(U_j^n, U_j^n) - F_{j+1/2}^n}{U_{j+1}^n - U_j^n} \\
            & \le \frac{\Delta t}{\Delta x} \norm{\pdv{F}{a}}_{L^\infty}
            \tag*{(Lipschitz flux)}
        \end{align*}
        $$

        and similarly for $D_{j-1/2}^n$. Then, under the CFL condition $\eqref{eq:monotone_CFL}$, we get 

        $$1 - C_{j+1/2}^n - D_{j-1/2}^n \ge 1 - \frac{\Delta x}{\Delta t} \qty(\abs{\pdv{F}{a}} + \abs{\pdv{F}{b}}) \ge 0.$$


### Time continuity

In the continuous case, we have the time continuity property

$$\norm{U(\cdot, t) - U(\cdot, s)}_{L^1(\R)} \le |t-s| \norm{f}_{\text{Lip}} \norm{U_0}_{TV(\R)}.$$

We can show that the discrete solution has a similar property.

!!! lemma
    Let $U_j^n$ computed using  consistent, conservative and monotone scheme with Lipschitz continuous numerical fluxe $F=F(a,b)$ whose Lipschitz constant is $C_F$. Then we get

    $$\Delta x \sum_j \abs{U_j^n - U_j^m} \le 2C_F \abs{t^n - t^m} \norm{U_0}_{TV(\R)}.$$

    ??? proof
        Taking the absolute value of $\eqref{eq:godunov}$ and summing over $j$, we get

        $$
        \begin{aligned}
            \Delta x \sum_j \abs{U_j^{n+1} - U_j^n}
            &= \Delta t \sum_j \abs{F_{j+1/2}^n - F_{j-1/2}^n} \\
            &\le \Delta t \sum_j \abs{F(U_j^n, U_{j+1}^n) - F(U_j^n, U_j^n)} \\
            & \quad + \Delta t \sum_j \abs{F(U_j^n, U_j^n) - F(U_{j-1}^n, U_j^n)} \\
            & \le \Delta t \sum_j C_F \abs{U_{j+1}^n - U_j^n} + \Delta t \sum_j C_F \abs{U_j^n - U_{j-1}^n} \\
            &= 2C_F \Delta t \sum_j \abs{U_{j+1}^n - U_j^n}.
        \end{aligned}
        $$

        Now, taking the difference from $k=n$ to $k=m$ and using the TVD property, we get

        $$
        \begin{aligned}
            \Delta x \sum_j \abs{U_j^n - U_j^m}
            & \le \Delta x \sum_{k=n}^{m-1} \qty(\sum_j \abs{U_j^{k+1} - U_j^k}) \\
            & \le 2C_F \Delta t \sum_{k=n}^{m-1} \norm{U^k}_{TV(\R)} \\
            & \le 2C_F \abs{t^n - t^m} \norm{U_0}_{TV(\R)}.
        \end{aligned}
        $$


## Convergence of monotone schemes

We have now shown that conserative, consistent and monotone schemes are stable in $L^\infty$ and $TV$, are TVD and satisfy the discrete entropy inequality. We will use this to prove convergence. Further, Lax-Wendroff theorem shows that the limit is the entropy solution. To simplify notation, we will set

$$
\newcommand{\UDx}{U^{\Delta x}}
\newcommand{\pDx}{\phi^{\Delta x}}
\UDx(x, t) = U_j^n \quad \text{for } (x, t) \in \Cell_j \times [t^n, t^{n+1}).
$$

!!! theorem "Lax-Wendroff theorem"
    Let $U_j^n$ an approxiamate solution of $\eqref{eq:conservation_law}$ using a consistent and conservative finite volume scheme. Assume that the numerical flux is differentiable and that $U_j^0$ are the cell averages of the initial data.

    Assume further that

    - $U_0 \in L^\infty(\R)$,
    - $\UDx$ is uniformly bounded in $\Delta x$,
    - there exists a function $U \in L_{\text{loc}}^1(\R \times R_+)$ such that $\UDx \to U$ in this space as $\Delta x, \Delta t \to 0$.

    Then, $U$ is a weak solution of $\eqref{eq:conservation_law}$ with initial data $U_0$.

    ??? proof
        Let $\phi \in C_c^\infty(\R \times \R_+)$, be a test function. Next, sum $\eqref{eq:godunov}$ over $j, n$:

        $$
        \begin{aligned}
            0 &= \Dx \Dt \sum_j \sum_{n=0}^\infty
            \qty(\frac{U_j^{n+1} - U_j^n}{\Dt} + \frac{F_{j+1/2}^n - F_{j-1/2}^n}{\Dx}) \phi_j^n \\
            &= - \Dx \sum_j U_j^0 \phi_j^0 - \Dx \sum_j \sum_{n=0}^\infty
            \qty(U_j^{n+1} \frac{\phi_j^{n+1} - \phi_j^n}{\Dt} + F_{j+1/2}^n \frac{\phi_{j+1}^n - \phi_j^n}{\Dx})
        \end{aligned}
        $$

        As $\phi$ has compact support, the above sum is finite. Define $\phi_j^n = \phi(x_j, t^n)$ and denote $\pDx(x, t)$ similarly to $\UDx$. Then, we can express the above sums as integrals:

        $$
        \begin{aligned}
            0 &= \underbrace{\int_\R \UDx(x, 0) \pDx(x, 0) \dd x
            + \int_\R \int_0^\infty \UDx(x, t + \Dt) \frac{\pDx(x, t + \Dt) - \pDx(x, t)}{\Dt} \dd t \dd x}_{I_1^{\Delta x}} \\
            &\quad + \underbrace{\int_\R \int_0^\infty F\qty(\UDx(x, t), \UDx(x + \Dx, t)) \frac{\pDx(x + \Dx, t) - \pDx(x, t)}{\Dx} \dd t \dd x}_{I_2^{\Delta x}}.
        \end{aligned}
        $$

        We have that $\abs{\UDx} \le \abs{U}$ and $\abs{\pDx} \le \abs{\phi}$, so we can use the dominated convergence theorem for the first two integrals $I_1^{\Delta x}$ yielding

        $$I_1^{\Delta x} \xrightarrow{\Delta x, \Delta t \to 0} \int_\R U_0(x) \phi(x, 0) \dd x + \int_\R \int_0^\infty U \phi_t \dd t \dd x.$$

        Next, rewrite $I_2^{\Delta x}$ as

        $$
        \begin{aligned}
            I_2^{\Delta x} &= \underbrace{\int_\R \int_0^\infty f(\UDx(x, t)) \frac{\pDx(x + \Dx, t) - \pDx(x, t)}{\Dx} \dd t \dd x}_{I_{2,1}^{\Delta x}} \\
            & \quad\quad + \underbrace{\int_\R \int_0^\infty \qty[F\qty(\UDx(x, t), \UDx(x + \Dx, t)) - f(\UDx(x, t))]
            \frac{\pDx(x + \Dx, t) - \pDx(x, t)}{\Dx} \dd t \dd x}_{I_{2,2}^{\Delta x}}.
        \end{aligned}
        $$

        ??? claim "Claim: $I_{2,2}^{\Delta x} \xrightarrow{\Delta x, \Delta t \to 0} 0$"
            Recall that $F$ is differentiable and $U \in L^\infty(\R \times \R_+)$, so we have

            $$
            \begin{aligned}
                & \abs{F\qty(\UDx(x, t), \UDx(x + \Dx, t)) - f(\UDx(x, t))} \\
                & \quad\quad = \abs{F(\UDx(x, t), \UDx(x + \Dx, t)) - F(\UDx(x, t), \UDx(x, t))} \\
                & \quad\quad \le \tilde C \abs{\UDx(x + \Dx, t) - \UDx(x, t)}.
            \end{aligned}
            $$

            we can use the above inequality for $I_{2,2}^{\Delta x}$ to get

            $$\abs{I_{2,2}^{\Delta x}} \le \tilde C \norm{\phi_x}_{L^\infty}
            \iint_{\supp \phi_x} \abs{\UDx(x + \Dx, t) - \UDx(x, t)} \dd x \dd t$$

            Using approximation by $C^\infty$ functions, we get thedesired result.

            Let $\psi_l \in C^\infty(\R)$ be a sequence of functions converging to $U$ in $L_{\text{loc}}^1(\R\times \R_+)$. The integrand can be bounded by

            $$
            \begin{aligned}
                & {\abs{\UDx(x + \Dx, t) - \psi_l(x,t)}}
                + \abs{\UDx(x, t) - \psi_l(x,t)} \\
                & \quad \quad \le \underbrace{\abs{\UDx(x + \Dx, t) - U(x+\Dx, t)}}_{J_1}
                + \underbrace{\abs{U(x+\Dx, t) - \psi_l(x,t)}}_{J_2} \\
                & \quad \quad \quad + \underbrace{\abs{\UDx(x, t) - U(x, t)}}_{J_3}
                + \underbrace{\abs{U(x, t) - \psi_l(x,t)}}_{J_4}.
            \end{aligned}
            $$

            By dominated convergence, $\iint_{\supp \phi_x} J_i \dd x \dd t \to 0$ as $\Delta x, \Delta t \to 0$ for $i=1,3$. By definition of $\psi_l$, the same holds for $i=4$. Convince yourself that it stays true for $i=2$.


        Again, by diminated convergence, we get that

        $$I_{2,1}^{\Delta x} \xrightarrow{\Delta x, \Delta t \to 0}
        \int_\R \int_0^\infty f(U(x, t)) \phi_t \dd t \dd x.$$

        Combining the results, we get

        $$0 = \int_\R U_0(x) \phi(x, 0) \dd x + \int_\R \int_0^\infty U \phi_t + f(U) \phi_t \dd t \dd x.$$


!!! lemma
    Assume the same conditions as in the Lax-Wendroff theorem, including that the discrete entropy inequality $\eqref{eq:discrete_entropy_inequality}$ holds. Then,

    $$U = \lim_{\Delta x, \Delta t \to 0} \UDx$$

    is the entropy solution of $\eqref{eq:conservation_law}$.

    ??? proof
        The proof is almost identical to that of the Lax-Wendroff theorem. We start with the discrete entropy inequality $\eqref{eq:discrete_entropy_inequality}$ instead of $\eqref{eq:godunov}$ and use a non-negative test function $\phi$. Then, following the same steps as in the proof of the Lax-Wendroff theorem, we get

        $$0 \ge I_1^{\Delta x} + I_2^{\Delta x}$$

        where

        $$
        \begin{aligned}
            I_1^{\Delta x} &= \int_\R \abs{\UDx(x, 0) - k} \pDx(x, 0) \dd x \\
            & \quad\quad+ \int_\R \int_0^\infty \abs{\UDx(x, t + \Dt) - k} \frac{\pDx(x, t + \Dt) - \pDx(x, t)}{\Dt} \dd t \dd x, \\
            I_2^{\Delta x} &= \int_\R \int_0^\infty Q\qty(\UDx(x, t), \UDx(x + \Dx, t)) \frac{\pDx(x + \Dx, t) - \pDx(x, t)}{\Dx} \dd t \dd x.
        \end{aligned}
        $$

        As with in the proof of the Lax-Wendroff theorem, we have that

        $$I_1^{\Delta x} \xrightarrow{\Delta x, \Delta t \to 0} \int_\R \abs{U_0(x)-k} \phi(x, 0) \dd x + \int_\R \int_0^\infty \abs{U-k} \phi_t \dd t \dd x.$$

        We also rewrite $I_2^{\Delta x} = I_{2,1}^{\Delta x} + I_{2,2}^{\Delta x}$ as in Lax-Wendroff. Again,

        $$I_{2,1}^{\Delta x} \xrightarrow{\Delta x, \Delta t \to 0}
        \int_\R \int_0^\infty q(U(x, t); k) \phi_t \dd t \dd x.$$

        Doing as in the proof of the Lax-Wendroff theorem, we see that $I_{2,2}^{\Delta x}$ vanishes. Combining the results, we get

        $$0 \ge \int_\R \abs{U_0(x)-k} \phi(x, 0) \dd x + \int_\R \int_0^\infty \abs{U-k} \phi_t + q(U; k) \phi_t \dd t \dd x.$$


Now, we can use monotinicity, conservation, consistency and the following two known results to show convergence.

??? theorem "Ascoli's Theorem"
    Let $(x, d_X)$ be a metric space and $K \subset X$ be relatively compact. Further, let

    $$\qty{u_k : [0, T] \to K}_k$$

    be a sequence of uniformly Lipschitz continuous functions, i.e.

    $$d_X(u_k(t), u_k(s)) \le C \abs{t-s} \quad \forall k\in \N, t, s \in [0, T].$$

    Then, there exists a subsequence $u_{k_\ell}$ converging to a Lipschitz continuous function $u$ uniformly on $[0, T]$.


??? theorem "Helly's Theorem [Giu84, T.1.19]"
    Let $a, b \in \R$ and $K \subset L^1([a, b])$. If there exists a number $M>0$ such that

    $$\sup_{U \in K} \|U\|_{BV([a, b])} \le M,$$

    then $K$ is relatively compact in $L^1([a, b])$.

    ---

    If additionally, $K \subset L^1(\R)$ and

    $$\lim_{R\to\infty} \sup_{U\in K} \int_{\R \setminus [-R, R]} \abs{U(x)} \dd x = 0,$$

    then $K$ is relatively compact in $L^1(\R)$.


For the convergence, we define the piecewise affine (in time) function

$$\UDx(x, t) = \frac{t^{n+1} - t}{\Dt} U_j^n + \frac{t - t^n}{\Dt} U_j^{n+1}, \quad x \in \Cell_j, t \in [t^n, t^{n+1}).$$

!!! theorem
    Consider the conservation law $\eqref{eq:conservation_law}$ with $f \in C^1(\R)$ and $U_0 \in BV(\R)$, and a conmsistent, conservative, monotone finite volumes scheme under the CFL condition $\eqref{eq:monotone_CFL}$. Assume further that there is some $c>0$ such that $\frac{\Dt}{\Dx} \ge c$.

    Then, $\UDx \to U$ in $L^1(\R \times [0,T])$ and $U$ is the entropy solution.

    ??? proof
        Let $K = \qty{U = \UDx(\cdot, t) \mid t \in [0, T], \Dx > 0}$ be functions $U : \R \to \R$ attained at some time $t$ by the scheme.

        ??? claim "Claim: The conditions of Helly's Theorem are satisfied"
            By Crandall-Majda, we have that $K$ is bounded in $L^1(\R)$ by $\norm{U_0}_{L^1(\R)}$. We have also shown that the scheme is TVD, so $K \subset BV(\R)$.

            For the second part of the theorem, we have that $U_0 \in L^1(\R)$, so

            $$\forall \eps > 0 \exists\, r > 0 : \left\{\begin{aligned}
                & \abs{U_0(x)} \le \eps \quad \forall \abs{x} > r, \\
                & \int_{\R \setminus [-r, r]} \abs{U_0(x)} \dd x \le \eps.
            \end{aligned}\right.$$

            Additionally, by monotinicity, we have that

            $$\min\qty{U_{j-1}^n, U_j^n, U_{j+1}^n} \le U_j^{n+1} \le \max\qty{U_{j-1}^n, U_j^n, U_{j+1}^n}, \quad \forall n, j.$$

            So, propagating at the speed $1/c$, we get that

            $$\UDx(x, t) \le \eps \quad \forall \abs{x} > R, t\in[0, T],$$

            where $R = r + T/c$. Next, let $J \in \N$ be such that $R \in \Cell_J$. Then, we sum the discrete entropy inequality $\eqref{eq:discrete_entropy_inequality}$ over $\abs{j} > J$ and over $n=0,\dots,N$ for $k=0$:

            $$
            \begin{aligned}
                \sum_{\abs{j} > J} \sum_{n=0}^N \qty(\abs{U_j^{n+1}} - \abs{U_j^n})
                + \frac{\Dt}{\Dx} \sum_{\abs{j} > J} \sum_{n=0}^N \qty(Q_{j+1/2}^n - Q_{j-1/2}^n)
                & \le 0 \\
                \sum_{\abs{j} > J} \qty(\abs{U_j^{N+1}} - \abs{U_j^0})
                + \frac{\Dt}{\Dx} \sum_{n=0}^N \qty(Q_{-J-1/2}^n - Q_{J+1/2}^n)
                &\le 0
            \end{aligned}
            $$

            Now, $U_j^0$ is the cell average of $U_0$, so we further get the bound

            $$\Dx\sum_{\abs{j} > J} \abs{U_j^{N+1}}
            \le \int_{\R \setminus [-R, R]} \abs{U_0} \dd x
            + \Dt \sum_{n=0}^N \qty(\abs{Q_{-J-1/2}^n} + \abs{Q_{J+1/2}^n}).$$

            Finally, for $t \in [t^m, t^{m+1})$, we have

            $$
            \begin{aligned}
                \int_{\R \setminus [-R, R]} \abs{\UDx} \dd x
                &\le \Dx \sum_{\abs{j} > J} \qty(\abs{U_j^m} + \abs{U_j^{m+1}}) \\
                &\le 2\int_{\R \setminus [-R, R]} \abs{U_0} \dd x
                + 2\Dt \sum_{n=0}^m \qty(\abs{Q_{-J-1/2}^n} + \abs{Q_{J+1/2}^n}) \\
                &\le 2\eps + 2\Dt \sum_{n=0}^N \qty(\abs{Q_{-J-1/2}^n} + \abs{Q_{J+1/2}^n}).
            \end{aligned}
            $$

            ??? claim "Claim: $\abs{Q_{\pm J \pm 1/2}^n} \le 2C_F\eps$"
                Recall that $F$ is Lipschitz continuous, so
                === "$\sign{U_J^n} = \sign{U_{J+1}^n}$"
                    $$
                    \begin{aligned}
                        \abs{Q_{J+1/2}^n} &= \abs{F(U_J^n, U_{J+1}^n) - F(0, 0)} \\
                        & \le \abs{F(U_J^n, U_{J+1}^n) - F(U_J^n, 0)} + \abs{F(U_J^n, 0) - F(0, 0)} \\
                        & \le C_F \qty(\abs{U_{J+1}^n} + \abs{U_J^n})
                    \end{aligned}
                    $$

                === "$\sign{U_J^n} \neq \sign{U_{J+1}^n}$"
                    $$
                    \begin{aligned}
                        \abs{Q_{J+1/2}^n} &= \abs{F(U_J^n, 0) - F(0, U_{J+1}^n)} \\
                        & \le \abs{F(U_J^n, 0) - F(0, 0)} + \abs{F(0, 0) - F(0, U_{J+1}^n)} \\
                        & \le C_F \qty(\abs{U_J^n} + \abs{U_{J+1}^n}).
                    \end{aligned}
                    $$

                Similarly for $Q_{-J-1/2}^n$. Recall further that $J$ is chosen such that $U_J^n, U_{J+1}^n \le \eps$ for all $n$. Hence, we get the desired result.

            Then, we get the bound

            $$\int_{\R \setminus [-R, R]} \abs{\UDx} \dd x \le 2\eps \qty(1 + 4C_F T)$$

            independent of $t$ and $\Dx$. This proves the claim.

        So $K$ is relatively compact in $L^1(\R)$. Then, by Ascoli's theorem, there exists a subsequence $\Dx'$ and $U \in L^1(\R \times [0, T])$ such that $\UDx{}^{'} \to U$. By Lax-Wendroff, $U$ is the entropy solution.

        

```