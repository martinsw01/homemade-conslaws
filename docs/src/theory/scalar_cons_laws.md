# Scalar conservation laws

In the previous chapter, we considered the scalar conservation law

```math
U_t + a(x,t) U_x = 0
```

This was linear, as ``a`` was a given function. However, the velocity field can depend on the solution,
and the equation becomes nonlinear. For example, we can have

```math
a(x, t) = U(x, t)
```

and we get the *inviscid Burgers' equation*, which can be written on the conservative form

```math
\begin{equation}
    U_t + \qty(\frac{1}{2} U^2)_x = 0 \label{eq:burgers}
\end{equation}
```

whenever ``U`` is smooth. More generally, we may have

```math
\begin{equation}
    U_t + f(U)_x = 0 \label{eq:scalar_cons_law}
\end{equation}
```

where ``f`` is a given *flux* function. This is a *nonlinear scalar conservation law*.

## Charachteristics for Burgers' equation

The characteristics for the Burgers' equation are given by

```math
\begin{aligned}
    x'(t) &= U(x(t), t) \\
    x(0) &= x_0
\end{aligned}
```

We consider first the initial data

```math
U_0(x) = \begin{cases}
    U_L, & x < 0 \\
    U_R, & x > 0
\end{cases}
```

The soluiton is constant along the charachteristics, so

```math
x(t) = U_0(x_0) t + x_0.
```

Consider ``U_L = 1`` and ``U_R = 0``. Then, we clearly see that the characteristics will intersect. Thus,
we get a discontinuity in the solution. This can also happen for smooth initial data.

## Weak solutions

We therefore need to use the weak formulation

```math
\begin{equation}
    \int_{\R_+} \int_\R U \phi_t + f(U) \phi_x \dd x \dd t + \int_\R U_0(x) \phi(x, 0) \dd x = 0 \label{eq:weak_scalar_cons_law}
\end{equation}
```

for any test function ``\phi \in C_c^1(\R \times \R_+)``.

### The Rankine-Hugoniot condition

In nature, the discontinuities appear as *shock waves*. These must saisfy some conditions. Assume that we have a weak solution ``U`` that is smooth on two regions ``\Omega^-, \Omega^+`` separated by a shock wave. Let the shock wave be defined by the curve ``x = \gamma(t)``. Lastly, assume that ``U \in C^1(\Omega^-) \cap C^1(\Omega^+)``. Then, we have

```math
\begin{aligned}
    \int_\Omega U \phi_t + f(U) \phi_x \dd \Omega
    =& -\int_{\Omega^+} (U_t + f(U)_x) \phi \dd \Omega
    + \int_{\partial \Omega^+} \bm F(U^+) \cdot \bm \nu \phi \dd \Omega \\
    &- \int_{\Omega^-} (U_t + f(U)_x) \phi \dd \Omega
    + \int_{\partial \Omega^-} \bm F(U^-) \cdot \bm \nu \phi \dd \Omega \\
    =& 0
\end{aligned}
```

where ``\bm F(U) = (U, f(U))``. Up to normalizaton, we have that ``\bm \nu = (1, -s(t))``. Where ``s(t) = \gamma'(t)`` is the speed of the shock curve. Then we have

```math
\int_{\Omega^+\cup\Omega^-} \underbrace{\qty(U_t + f(U)_x)}_0 \phi \dd \Omega
+ \int_{\gamma} \qty(s(t)\qty(U^+ - U^-) - \qty(f(U^+) - f(U^-))) \phi \dd \Omega = 0
```

This needs to hold for all test functions and ``U^+, U^-`` are smooth on ``\gamma``, so we get the *Rankine-Hugoniot condition*

```math
\begin{equation}
    s(t) \qty(U^+ - U^-) = f(U^+) - f(U^-) \label{eq:rankine_hugoniot}
\end{equation}
```

where ``U^\pm`` are evaluated on ``(\gamma(t),t)``.

!!! theorem
    Let ``\gamma \in C^1(\R_+)`` and ``U \in L^\infty(\R \times \R_+)`` be of the form

    ```math
    U(x, t) = \begin{cases}
        U^+(x,t), & x < \gamma(t) \\
        U^-(x,t), & x > \gamma(t)
    \end{cases}
    ```

    where ``U^\pm \in C^1(\R \times \R_+)``. Then, ``U`` solves ``(\ref{eq:scalar_cons_law})`` weakly if and only if
    
    - ``U^\pm`` solves it in the classical sense, and
    - the shock speed ``s(t)`` satisfies the Rankine-Hugoniot condition ``(\ref{eq:rankine_hugoniot})`` at ``x = \gamma(t)``.


### Lax entropy condition

We consider again the Burgers' equation ``(\ref{eq:burgers})`` on the Riemann problem with ``U_L = 1`` and ``U_R = 0``. The yields the shock speed

```math
    s(t) = \frac{\frac{U_R^2}{2} - \frac{U_L^2}{2}}{U_R - U_L} \equiv \frac{1}{2}
```

The shock starts at ``x = 0``, so we get the weak solution

```math
\begin{equation}
    U(x, t) = \begin{cases}
        1, & x < \frac{1}{2} t \\
        0, & x > \frac{1}{2} t.
    \end{cases} \label{eq:riemann_solution_1}
\end{equation}
```

This satisfies ``(\ref{eq:weak_scalar_cons_law})`` and we see that the characteristsics flow into the shock. 

Next, consider the Riemann problem with ``U_L = 0`` and ``U_R = 1``. Then, the characteristics do not fill the plane. However, we see that 

```math
\begin{equation}
    U(x, t) = \begin{cases}
        0, & x < \frac{1}{2} t \\
        1, & x > \frac{1}{2} t
    \end{cases} \label{eq:riemann_solution_2}
\end{equation}
```

is a solution of the weak formulation ``(\ref{eq:weak_scalar_cons_law})``. However, one can easily construct another solution

```math
\begin{equation}
    U(x, t) = \begin{cases}
        0, & x < \frac{1}{3} t \\
        \frac{2}{3}, & \frac{1}{3} t < x < \frac{5}{6} t \\
        0, & x > \frac{5}{6} t.
    \end{cases} \label{eq:riemann_solution_3}
\end{equation}
```

This nonuniqueness is implicit in the weak formulation, so we need to impose extra conditions. Observe that the charachteristics of ``(\ref{eq:riemann_solution_1})`` flow into the shock, while the characteristics of ``(\ref{eq:riemann_solution_2})`` and ``(\ref{eq:riemann_solution_3})`` flow out from the shock. Intuitively, information should always flow from the initial data. This is the case for ``(\ref{eq:riemann_solution_1})``, but not for the other two. For the Burgers' equation, this requirement can be expressed as

```math
U^- > s > U^+.
```

For convex ``f``, this can be generalized to the *Lax entropy condition*

```math
\begin{equation}
    f'(U^-) \leq s \leq f'(U^+).
\end{equation}
```

!!! exercise
    We know that

    ```math
    U(x, t) = \begin{cases}
        U_L, & x < s t \\
        U_R, & x > s t
    \end{cases}
    ```

    is a weak solution of the Riemann problem.

    Assume that ``f`` is strictly convex and that ``U_L > U_R``. Show that ``U`` satisfies the Lax entropy condition.

    Similarly, show that if ``U_L < U_R``, then ``U`` does not satisfy the Lax entropy condition.

It turns out that in the last case, a continuous solution exists.


### Rarefaction waves

We will assume that ``f`` is strictly convex. Note that scaling ``(x, t) \mapsto (\lambda x, \lambda t)`` does not change the solution and that the initial data is invariant under this scaling. Thus, it is natural to assume that the solution only depends on the ratio ``\xi:=\frac{x}{t}``:

```math
U(x, t) = V\qty(\frac{x}{t}).
```

Substituting this into the conservation law ``\ref{eq:scalar_cons_law}``, we get

```math
(f'(V(\xi)) - \xi) V'(\xi) = 0.
```

As ``f'`` is strictly increasing, the inverse function ``f'^{-1}`` exists. Thus, we get the *rarefaction wave*

```math
V'(x/t) = \qty(f')^{-1}(x/t).
```

Then we can construct a weak solution satisfying the Lax entropy condition by

```math
U(x, t) = \begin{cases}
    U_L, & x \le f'(U_L) t \\
    (f')^{-1}(x/t), & f'(U_L) t < x \le f'(U_R) t \\
    U_R, & x > f'(U_R) t.
\end{cases}
```

## Entropy solutions

The entropy condition is a local condition at shocks, so it may be difficoult to apply it in proofs of global stability estimates. We will here derive an equivalent global condition.

### The entropy condition

We add a viscous term to the conservation law:
    
```math
\begin{equation}
    U_t^\eps + f(U^\eps)_x = \eps U_{xx}^\eps \label{eq:convection-diffusion}
\end{equation}
```

and later let ``\eps \to 0``. We now have a parabolic *convection-diffusion equation*. These equations have ``C^\infty`` solutions. Now, let ``\eta : \R \to \R`` be a strictly convex function and define

```math
q(U) = \int_0^U f'(v) \nu'(v) \dd v.
```

We call ``(\eta, q)`` an *entropy pair* and respectively the *entropy function* and the *entropy flux*. Note that ``q' = f' \nu'``. Then, multiplying ``(\ref{eq:convection-diffusion})`` with ``\eta'(U^\eps)``, we get

```math
\eta(U^\eps)_t + q(U^\eps)_x = \eps \eta(U^\eps)_{xx} - \eps \eta''(U^\eps) (U^\eps_x)^2.
```

By the conevxity of ``\eta``, the last term is nonnegative, so we have

```math
\eta(U^\eps)_t + q(U^\eps)_x \leq \eps \eta(U^\eps)_{xx}.
```

Then, for the vanishing viscosity solution ``\displaystyle U = \lim_{\eps \to 0} U^\eps``, we get the *entropy condition*

```math
\begin{equation}
    \eta(U)_t + q(U)_x \leq 0. \label{eq:entropy_condition} \tag{EI}
\end{equation}
```

Again, this must be interpreted in the sense of distributions:

```math
\int_{\R_+} \int_\R \eta(U) \phi_t + q(U) \phi_x \dd x \dd t + \int_\R \eta(U_0(x)) \phi(x, 0) \dd x \leq 0
```

!!! definition "Definition (Entropy solution)"
    A function ``U \in L^\infty(\R \times \R_+)`` is an *entropy solution* of ``(\ref{eq:scalar_cons_law})`` if both

    - ``U`` is a weak solution of ``(\ref{eq:scalar_cons_law})``, and
    - ``U`` satisfies the entropy condition ``(\ref{eq:entropy_condition})`` for all entropy pairs ``(\eta, q)``.


### Kruzkov Entropy Condition

Of special importance is the family of entropy pairs

```math
\begin{aligned}
    \eta(U; c) &= \abs{U - c} \\
    q(U; c) &= \sign(U - c) [f(U) - f(c)].
\end{aligned}
```

for ``c \in \R``. ``U``satisfies the *Kruzkov entropy condition* if

```math
\begin{equation}
    \abs{U - c}_t + \sign(U - c) [f(U) - f(c)]_x \leq 0 \label{eq:kruzkov_entropy_condition} \tag{KEC}
\end{equation}
```

in the weak sence for all ``c \in \R``. This is important because of the following lemma:

!!! lemma
    Let ``U \in L^\infty(\R \times \R_+)``.

    ```math
    U \text{ is an entropy solution of } (\ref{eq:scalar_cons_law}) \iff U \text{ satisfies } (\ref{eq:kruzkov_entropy_condition})
    ```

    !!! proof
        === "``\implies``"
            Follows from the definition of entropy solutions.
        === "``\impliedby``"
            === "Weak solution"
                Let ``a, b \in \R``. such that ``a \le U \le b`` a.e. in ``\R \times \R_+``. Then, selecting ``c=a``, we get by ``(\ref{eq:kruzkov_entropy_condition})``

                ```math
                U_t + f(U)_x \leq 0 \quad \text{a.e. in } \R \times \R_+.
                ```

                Similarly, selecting ``c=b``, we get

                ```math
                U_t + f(U)_x \geq 0 \quad \text{a.e. in } \R \times \R_+.
                ```

                This yields

                ```math
                U_t + f(U)_x = 0 \quad \text{a.e. in } \R \times \R_+.
                ```

                so ``U`` is a weak solution. 
            === "Entropy condition"
                Any convex combination can be approximated by linear combinations of ``\eta(U; a)``. Thus, the entropy condition holds for all entropy pairs.


!!! theorem
    Let ``\gamma \in C^1(\R_+)`` be a curve and ``s = \gamma'`` be its speed. Assume ``U \in L^\infty(\R \times \R_+)`` is a weak solution of ``(\ref{eq:scalar_cons_law})`` of the form

    ```math
    U(x, t) = \begin{cases}
        U^-(x, t), & x < \gamma(t) \\
        U^+(x, t), & x > \gamma(t)
    \end{cases}
    ```

    where ``U^\pm \in C^1``. Then, the following are equivalent:

    1. ``U`` is an entropy solution.
    2. At ``x = \gamma(t)``, ``U`` satisfies

        
        ``\jump{q(U)} - s\jump{U} \le 0``

        for all entropy pairs ``(\eta, q)``.
    3. For all ``v \in (U^-, U^+)``, we have

        ``\frac{f(v) - f(U^-)}{v - U^-} \le s \le \frac{f(v) - f(U^+)}{v - U^+}``

        along ``s``.
