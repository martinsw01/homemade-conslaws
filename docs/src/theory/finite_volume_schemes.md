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
using homemade_conslaws.Viz: animate_solution
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
U, t = godunov_scheme(f, df, ω, U0, BC, dx, dt, T)

animate_solution((U,),
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
U, t = godunov_scheme(f, df, ω, U0, BC, dx, dt, T)

animate_solution((U,),
                 "Approximate solution",
                 x_mid, t)
```

For a more complicated example, ``U_0 = \sin(4\pi x)``, we expect a compound wave:

```@example 1
U0 = sin.(4π*x_mid)
BC(t, U) = [0;; 0]
U, t = godunov_scheme(f, df, ω, U0, BC, dx, dt, T)

animate_solution((U,),
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
U0 = x_mid .< 0     # hide
BC(t, U) = [1;; 0]  # hide
dt = 0.1 # max time step    # hide
T = 1   # hide
U, t = lax_friedrichs_scheme(f, df, U0, BC, dx, dt, T)    # hide

animate_solution(U, (x, t) -> x < 0.5t, x_mid, t) #hide
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
U, t = rusanov_scheme(f, df, U0, BC, dx, dt, T)    # hide

animate_solution(U, (x, t) -> x < 0.5t, x_mid, t) #hide
```

and with the initial data $\eqref{eq:initial_condition_rarefaction}$, we get

```@example 1
U0 = (x_mid .> 0) * 2 .-1 # hide
BC(t, U) = [-1;; 1] # hide
U, t = rusanov_scheme(f, df, U0, BC, dx, dt, T) # hide

animate_solution((U,), "Approximation", x_mid, t) # hide
```

### Comparison

```@example 1
U0 = x_mid .< 0     # hide
BC(t, U) = [1;; 0] # hide

U_Rus, t = rusanov_scheme(f, df, U0, BC, dx, dt, T) # hide
U_LxF, _ = lax_friedrichs_scheme(f, df, U0, BC, dx, dt, T) # hide
U_God, _ = godunov_scheme(f, df, ω, U0, BC, dx, dt, T) # hide

animate_solution((U_Rus, U_LxF, U_God, x_mid' .< 0.5t), # hide
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

???+ definition "Definition (Conservative shceme)"
    The generic scheme $\eqref{eq:general_scheme}$ approximating $\eqref{eq:conservation_law}$ is called *conservative* if

    $$\sum_j U_j^{n+1} = \sum_j U_j^n,$$

    ignoring the boundary conditions.


???+ theorem
    Assume $H(0, \dots, 0) = 0$. 

    $$
    \begin{gather*}
        \eqref{eq:general_scheme} \text{ is conservative } \\
        \Updownarrow \\
        \text{There is a function } F_{j+1/2}^n = F(U_{j-p}^n, \dots, U_{j+p}^n) \text{ such that \eqref{eq:general_scheme} can be written as } \eqref{eq:cell_average_update}.
    \end{gather*}
    $$

    ??? proof
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

We now cinsider numerical fluxes $F_{j+1/2}^n$ defined on $2p+1$ stencils

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

    $$\abs{\pdv{F}{a}(v,w)} + \abs{\pdv{F}{b}(u,v)} \le \frac{\Delta x}{\Delta t} \quad \forall u, v, w.$$


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
    Let $U_j^n$ be the approximate solution using a consistent, conservative, monotone scheme of the form $\eqref{eq:general_scheme}$. Then we have

    $$\abs{U_j^{n+1} - k}  - \abs{U_j^n - k} + \frac{\Delta t}{\Delta x}\qty(Q_{j+1/2}^n - Q_{j-1/2}^n) \le 0 \quad \forall n,j.$$

    Moreover, if $U_0 \in L^1(\R)$, then

    $$\sum_j \abs{U_j^n} \Delta x \le \norm{U_0}_{L^1{\R}}.$$

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

If we We define the discrete version of the total variation as

$$\norm{U^n}_{TV(\R)} = \sum_j \abs{U_j^n - U_{j-1}^n}.$$
```