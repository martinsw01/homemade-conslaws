```@raw html
# Scalar conservation laws

In the previous chapter, we considered the scalar conservation law

$$
U_t + a(x,t) U_x = 0
$$

This was linear, as $a$ was a given function. However, the velocity field can depend on the solution,
and the equation becomes nonlinear. For example, we can have

$$
a(x, t) = U(x, t)
$$

and we get the *inviscid Burgers' equation*, which can be written on the conservative form

$$
\begin{equation}
    U_t + \qty(\frac{1}{2} U^2)_x = 0 \label{eq:burgers}
\end{equation}
$$

whenever $U$ is smooth. More generally, we may have

$$
\begin{equation}
    U_t + f(U)_x = 0 \label{eq:scalar_cons_law}
\end{equation}
$$

where $f$ is a given *flux* function. This is a *nonlinear scalar conservation law*.

## Charachteristics for Burgers' equation

The characteristics for the Burgers' equation are given by

$$
\begin{aligned}
    x'(t) &= U(x(t), t) \\
    x(0) &= x_0
\end{aligned}
$$

We consider first the initial data

$$
U_0(x) = \begin{cases}
    U_L, & x < 0 \\
    U_R, & x > 0
\end{cases}
$$

The soluiton is constant along the charachteristics, so

$$
x(t) = U_0(x_0) t + x_0.
$$

Consider $U_L = 1$ and $U_R = 0$. Then, we clearly see that the characteristics will intersect. Thus,
we get a discontinuity in the solution. This can also happen for smooth initial data.

## Weak solutions

We therefore need to use the weak formulation

$$
\begin{equation}
    \int_{\R_+} \int_\R U \phi_t + f(U) \phi_x \dd x \dd t + \int_\R U_0(x) \phi(x, 0) \dd x = 0 \label{eq:weak_scalar_cons_law}
\end{equation}
$$

for any test function $\phi \in C_c^1(\R \times \R_+)$.

### The Rankine-Hugoniot condition

In nature, the discontinuities appear as *shock waves*. These must saisfy some conditions. Assume that we have a weak solution $U$ that is smooth on two regions $\Omega^-, \Omega^+$ separated by a shock wave. Let the shock wave be defined by the curve $x = \gamma(t)$. Lastly, assume that $U \in C^1(\Omega^-) \cap C^1(\Omega^+)$. Then, we have

$$
\begin{aligned}
    \int_\Omega U \phi_t + f(U) \phi_x \dd \Omega
    =& -\int_{\Omega^+} (U_t + f(U)_x) \phi \dd \Omega
    + \int_{\partial \Omega^+} \bm F(U^+) \cdot \bm \nu \phi \dd \Omega \\
    &- \int_{\Omega^-} (U_t + f(U)_x) \phi \dd \Omega
    + \int_{\partial \Omega^-} \bm F(U^-) \cdot \bm \nu \phi \dd \Omega \\
    =& 0
\end{aligned}
$$

where $\bm F(U) = (U, f(U))$. Up to normalizaton, we have that $\bm \nu = (1, -s(t))$. Where $s(t) = \gamma'(t)$ is the speed of the shock curve. Then we have

$$
\int_{\Omega^+\cup\Omega^-} \underbrace{\qty(U_t + f(U)_x)}_0 \phi \dd \Omega
+ \int_{\gamma} \qty(s(t)\qty(U^+ - U^-) - \qty(f(U^+) - f(U^-))) \phi \dd \Omega = 0
$$

This needs to hold for all test functions and $U^+, U^-$ are smooth on $\gamma$, so we get the *Rankine-Hugoniot condition*

$$
\begin{equation}
    s(t) \qty(U^+ - U^-) = f(U^+) - f(U^-) \label{eq:rankine_hugoniot} \tag{RH}
\end{equation}
$$

where $U^\pm$ are evaluated on $(\gamma(t),t)$.

!!! theorem
    Let $\gamma \in C^1(\R_+)$ and $U \in L^\infty(\R \times \R_+)$ be of the form

    $$
    U(x, t) = \begin{cases}
        U^+(x,t), & x < \gamma(t) \\
        U^-(x,t), & x > \gamma(t)
    \end{cases}
    $$

    where $U^\pm \in C^1(\R \times \R_+)$. Then, $U$ solves $\eqref{eq:scalar_cons_law}$ weakly if and only if
    
    - $U^\pm$ solves it in the classical sense, and
    - the shock speed $s(t)$ satisfies the Rankine-Hugoniot condition $\eqref{eq:rankine_hugoniot}$ at $x = \gamma(t)$.


### Lax entropy condition

We consider again the Burgers' equation $\eqref{eq:burgers}$ on the Riemann problem with $U_L = 1$ and $U_R = 0$. The yields the shock speed

$$
    s(t) = \frac{\frac{U_R^2}{2} - \frac{U_L^2}{2}}{U_R - U_L} \equiv \frac{1}{2}
$$

The shock starts at $x = 0$, so we get the weak solution

$$
\begin{equation}
    U(x, t) = \begin{cases}
        1, & x < \frac{1}{2} t \\
        0, & x > \frac{1}{2} t.
    \end{cases} \label{eq:riemann_solution_1}
\end{equation}
$$

This satisfies $\eqref{eq:weak_scalar_cons_law}$ and we see that the characteristsics flow into the shock. 

Next, consider the Riemann problem with $U_L = 0$ and $U_R = 1$. Then, the characteristics do not fill the plane. However, we see that 

$$
\begin{equation}
    U(x, t) = \begin{cases}
        0, & x < \frac{1}{2} t \\
        1, & x > \frac{1}{2} t
    \end{cases} \label{eq:riemann_solution_2}
\end{equation}
$$

is a solution of the weak formulation $\eqref{eq:weak_scalar_cons_law}$. However, one can easily construct another solution

$$
\begin{equation}
    U(x, t) = \begin{cases}
        0, & x < \frac{1}{3} t \\
        \frac{2}{3}, & \frac{1}{3} t < x < \frac{5}{6} t \\
        0, & x > \frac{5}{6} t.
    \end{cases} \label{eq:riemann_solution_3}
\end{equation}
$$

This nonuniqueness is implicit in the weak formulation, so we need to impose extra conditions. Observe that the charachteristics of $\eqref{eq:riemann_solution_1}$ flow into the shock, while the characteristics of $\eqref{eq:riemann_solution_2}$ and $\eqref{eq:riemann_solution_3}$ flow out from the shock. Intuitively, information should always flow from the initial data. This is the case for $\eqref{eq:riemann_solution_1}$, but not for the other two. For the Burgers' equation, this requirement can be expressed as

$$
U^- > s > U^+.
$$

For convex $f$, this can be generalized to the *Lax entropy condition*

$$
\begin{equation}
    f'(U^-) \leq s \leq f'(U^+).
\end{equation}
$$

!!! exercise
    We know that

    $$
    U(x, t) = \begin{cases}
        U_L, & x < s t \\
        U_R, & x > s t
    \end{cases}
    $$

    is a weak solution of the Riemann problem.

    Assume that $f$ is strictly convex and that $U_L > U_R$. Show that $U$ satisfies the Lax entropy condition.

    Similarly, show that if $U_L < U_R$, then $U$ does not satisfy the Lax entropy condition.

It turns out that in the last case, a continuous solution exists.


### Rarefaction waves

We will assume that $f$ is strictly convex. Note that scaling $(x, t) \mapsto (\lambda x, \lambda t)$ does not change the solution and that the initial data is invariant under this scaling. Thus, it is natural to assume that the solution only depends on the ratio $\xi:=\frac{x}{t}$:

$$
U(x, t) = V\qty(\frac{x}{t}).
$$

Substituting this into the conservation law $\eqref{eq:scalar_cons_law}$, we get

$$
(f'(V(\xi)) - \xi) V'(\xi) = 0.
$$

As $f'$ is strictly increasing, the inverse function $f'^{-1}$ exists. Thus, we get the *rarefaction wave*

$$
V'(x/t) = \qty(f')^{-1}(x/t).
$$

Then we can construct a weak solution satisfying the Lax entropy condition by

$$
\begin{equation}
    U(x, t) = \begin{cases}
        U_L, & x \le f'(U_L) t \\
        (f')^{-1}(x/t), & f'(U_L) t < x \le f'(U_R) t \\
        U_R, & x > f'(U_R) t.
    \end{cases} \label{eq:rarefaction_solution}
\end{equation}
$$

## Entropy solutions

The entropy condition is a local condition at shocks, so it may be difficoult to apply it in proofs of global stability estimates. We will here derive an equivalent global condition.

### The entropy condition

We add a viscous term to the conservation law:
    
$$
\begin{equation}
    U_t^\eps + f(U^\eps)_x = \eps U_{xx}^\eps \label{eq:convection-diffusion}
\end{equation}
$$

and later let $\eps \to 0$. We now have a parabolic *convection-diffusion equation*. These equations have $C^\infty$ solutions. Now, let $\eta : \R \to \R$ be a strictly convex function and define

$$
q(U) = \int_0^U f'(v) \nu'(v) \dd v.
$$

We call $(\eta, q)$ an *entropy pair* and respectively the *entropy function* and the *entropy flux*. Note that $q' = f' \nu'$. Then, multiplying $\eqref{eq:convection-diffusion}$ with $\eta'(U^\eps)$, we get

$$
\eta(U^\eps)_t + q(U^\eps)_x = \eps \eta(U^\eps)_{xx} - \eps \eta''(U^\eps) (U^\eps_x)^2.
$$

By the conevxity of $\eta$, the last term is nonnegative, so we have

$$
\eta(U^\eps)_t + q(U^\eps)_x \leq \eps \eta(U^\eps)_{xx}.
$$

Then, for the vanishing viscosity solution $\displaystyle U = \lim_{\eps \to 0} U^\eps$, we get the *entropy condition*

$$
\begin{equation}
    \eta(U)_t + q(U)_x \leq 0. \label{eq:entropy_condition} \tag{EI}
\end{equation}
$$

Again, this must be interpreted in the sense of distributions:

$$
\int_{\R_+} \int_\R \eta(U) \phi_t + q(U) \phi_x \dd x \dd t + \int_\R \eta(U_0(x)) \phi(x, 0) \dd x \leq 0
$$

!!! definition "Definition (Entropy solution)"
    A function $U \in L^\infty(\R \times \R_+)$ is an *entropy solution* of $\eqref{eq:scalar_cons_law}$ if both

    - $U$ is a weak solution of $\eqref{eq:scalar_cons_law}$, and
    - $U$ satisfies the entropy condition $\eqref{eq:entropy_condition}$ for all entropy pairs $(\eta, q)$.


### Kruzkov Entropy Condition

Of special importance is the family of entropy pairs

$$
\begin{aligned}
    \eta(U; c) &= \abs{U - c} \\
    q(U; c) &= \sign(U - c) [f(U) - f(c)].
\end{aligned}
$$

for $c \in \R$. $U$satisfies the *Kruzkov entropy condition* if

$$
\begin{equation}
    \abs{U - c}_t + \sign(U - c) [f(U) - f(c)]_x \leq 0 \label{eq:kruzkov_entropy_condition} \tag{KEC}
\end{equation}
$$

in the weak sence for all $c \in \R$. This is important because of the following lemma:

!!! lemma
    Let $U \in L^\infty(\R \times \R_+)$.

    $$
    U \text{ is an entropy solution of } \eqref{eq:scalar_cons_law} \iff U \text{ satisfies } \eqref{eq:kruzkov_entropy_condition}
    $$

    ??? proof
        === "$\implies$"
            Follows from the definition of entropy solutions.
        === "$\impliedby$"
            === "Weak solution"
                Let $a, b \in \R$. such that $a \le U \le b$ a.e. in $\R \times \R_+$. Then, selecting $c=a$, we get by $\eqref{eq:kruzkov_entropy_condition}$

                $$
                U_t + f(U)_x \leq 0 \quad \text{a.e. in } \R \times \R_+.
                $$

                Similarly, selecting $c=b$, we get

                $$
                U_t + f(U)_x \geq 0 \quad \text{a.e. in } \R \times \R_+.
                $$

                This yields

                $$
                U_t + f(U)_x = 0 \quad \text{a.e. in } \R \times \R_+.
                $$

                so $U$ is a weak solution. 
            === "Entropy condition"
                Any convex combination can be approximated by linear combinations of $\eta(U; a)$. Thus, the entropy condition holds for all entropy pairs.


!!! theorem
    Let $\gamma \in C^1(\R_+)$ be a curve and $s = \gamma'$ be its speed. Assume $U \in L^\infty(\R \times \R_+)$ is a weak solution of $\eqref{eq:scalar_cons_law}$ of the form

    $$
    U(x, t) = \begin{cases}
        U^-(x, t), & x < \gamma(t) \\
        U^+(x, t), & x > \gamma(t)
    \end{cases}
    $$

    where $U^\pm \in C^1$. Then, the following are equivalent:

    1. $U$ is an entropy solution.
    2. At $x = \gamma(t)$, $U$ satisfies

        
        $$\jump{q(U)} - s\jump{U} \le 0$$

        for all entropy pairs $(\eta, q)$.
    3. For all $v \in (U^-, U^+)$, we have

        $$
        \begin{equation}
            \frac{f(v) - f(U^-)}{v - U^-} \ge s \ge \frac{f(v) - f(U^+)}{v - U^+} \label{eq:oleinik-E} \tag{OCD}
        \end{equation}
        $$

        along $x = \gamma(t)$. Also called the *Oleinik codition* $E$.
    4. If $f$ is convex/concave, then

        $$f'(U^-) \ge s \ge f'(U^+)$$

        along $x = \gamma(t)$.

    ??? proof
        The other direction is left as an exercise for the reader.
        === "1. $\implies$ 2."
            $U$ is an entropy solution, so we have that

            $$\abs{U - c}_t + \sign(U - c) [f(U) - f(c)]_x \le 0 \quad \forall c \in \R.$$

            Now, define

            $$\Omega_{>c} := \qty{(x, t) \in \Omega \mid U(x, t) > c}$$

            and similarly for $\Omega_{<c}$. Now, we can integrate the above inequality, split it into the two regions and use the divergence theorem for both regions. All the terms will cancel except for those at the shock curve:

            $$0 \ge \int_\gamma \qty(U^+ \nu^t + f(U^+) \nu^x) \phi \dd t + \int_\gamma \qty(U^- \nu^t + f(U^-) \nu^x) \phi \dd t$$

            Now, this holds for all test functions $\phi$, so we get

            $$\jump{q(U)} - s\jump{U} \le 0.$$

        === "2. $\implies$ 3."
            We have

            $$
            \begin{aligned}
                \jump{\eta(U)} &= \int_{U^-}^{U^+} \eta'(v) \dd v \\
                &= \eval[\nu'(v) (v-U^-)|_{U^-}^{U^+} - \int_{U^-}^{U^+} \nu''(v) (v-U^-) \dd v \\
                &= \eta' \jump{U} - \int_{U^-}^{U^+} \nu''(v) (v-U^-) \dd v.
            \end{aligned}
            $$

            Similiarly, and with a few more steps, we get

            $$\jump{q(U)} = \eta'(U^+) \jump{f(U)} - \int_{U^-}^{U^+} \nu''(v) (f(v) - f(U^-)) \dd v.$$

            Now, we can use 2. and $\eqref{eq:rankine_hugoniot}$ to get

            $$
            \begin{aligned}
                0 & \ge \jump{q(U)} - s\jump{U} \\
                & = \eta'(U^+) \Big(\underbrace{\jump{f(U)} - s\jump{U}}_0\Big) - \int_{U^-}^{U^+} \nu''(v) \Big(s(v - U^-) - (f(v) - f(U^-))\Big) \dd v.
            \end{aligned}
            $$


            Now, this must hold for all entrpy pairs $(\eta, q)$, we end up with

            $$s \le \frac{f(v) - f(U^-)}{v - U^-} \quad \forall\ v$$

            and similarly for the other inequality with $U^+$.

        === "3. $\implies$ 4."
            Follows by taking the limit $v \to U^\pm$.

### Stability estimates

We can now integrate the entropy condition $\eqref{eq:entropy_condition}$ over space to get stability estimates on solutions. We get that the entropy decreases in time:

$$0 \ge \int_\R \eta(U)_t + q(U)_x \dd x = \dv{t} \int_\R \eta(U) \dd x$$

Thus, we get the bound

$$\int_\R \eta(U) \dd x \le \int_\R \eta(U_0) \dd x.$$

Notably, for $\eta(U) = \abs{U}^p$, we get $L^p$ bounds

$$\norm{U(\cdot, t)}_{L^p(\R)} \le \norm{U_0}_{L^p(\R)}.$$

This holds for all $p \ge 1$, so taking the limit $p \to \infty$, we get the maximum principle

$$
\begin{equation}
    \norm{U(\cdot, t)}_{L^\infty(\R)} \le \norm{U_0}_{L^\infty(\R)} \label{eq:max_principle} \tag{MP}
\end{equation}
$$

We now have a bound on the amplitude, but we can also derive a bound on the derivate, i.e. the oscillations in $U$.

!!! definition "Definition (Total variation)"
    Let $g$ be defined on $[a, b]$. The *total variation* of $g$ is defined as

    $$\norm{g}_{TV([a,b])} = \sup_{\mathcal P} \sum_{j=1}^{N-1} \abs{g(x_{j+1}) - g(x_j)}$$

    For differentiable $g$, this reduces to the $L^1$ norm of the derivative:

    $$\norm{g'}_{TV([a,b])} = \int_a^b \abs{\dv{g}{x}} \dd x.$$

    This is inly a seminorm (all constant functions have zero total variation). The *bounded variation*, however is a norm:

    $$\norm{g}_{BV([a,b])} = \norm{g}_{L^1([a,b])} + \norm{g'}_{TV([a,b])}.$$

    The space of functions of bounded variation on $\R$ is

    $$BV(\R) = \qty{g \in L^1(\R) \mid \norm{g}_{BV(\R)} < \infty}.$$


!!! theorem
    Assume $f \in C^1(\R)$ and $U_0 \in L^1(\R) \cap L^\infty(\R)$.

    Then, there exists a unique entropy soluiton of $\eqref{eq:scalar_cons_law}$ such that

    1. $L^1$-bound:

        $$\norm{U(\cdot, t)}_{L^1(\R)} \le \norm{U_0}_{L^1(\R)}$$

    2. $L^\infty$-bound:

        $$\norm{U(\cdot, t)}_{L^\infty(\R)} \le \norm{U_0}_{L^\infty(\R)}$$

    3. Total variation bound:

        $$U_0 \in BV(\R) \implies \norm{U(\cdot, t)}_{TV(\R)} \le \norm{U_0}_{TV(\R)}$$

    4. Time continuity:

        $$U_0 \in BV(\R) \implies \norm{U(\cdot, t) - U(\cdot, s)}_{L^1(\R)} \le |t-s| M \norm{U_0}_{TV(\R)}$$

        - $M = M(U_0) := \max_{\underline{u} \le u \le \overline{u}} \abs{f'(u)}$
        - $\displaystyle \underline{u} = \essinf_\R U_0$
        - $\displaystyle \overline{u} = \esssup_\R U_0$.

    5. $L^1$-stability:

        $$\norm{U(\cdot, t) - V(\cdot, t)}_{L^1(\R)} \le \norm{U_0 - V_0}_{L^1(\R)} \quad \forall t \ge 0$$

    6. Local $L^1$-stability:

        $$\int_a^b \abs{U(x, t) - V(x, t)} \dd x \le \int_{a-Mt}^{b+Mt} \abs{U_0(x) - V_0(x)} \dd x \quad \forall t \ge 0, a < b$$

        - $M = \max\{M(U_0), M(V_0)\}$.

    7. Monotonicity: if we have $U_0 \le V_0$ almost everywhere in $\R$, then

    $$U(\cdot, t) \le V(\cdot, t) \quad \text{a.e. in } \R$$

    ??? proof
        === "Existence"
            Let $U^\eps$ solve $\eqref{eq:convection-diffusion})$ and let $U = \lim_{\eps \to 0} U^\eps$ be a limit. Then, we have seen that $U$ satisfies the entropy condition $\eqref{eq:entropy_condition}$. Additionally, from the maximum principle $\eqref{eq:max_principle}$, we get that $U$ is an entropy solution.
        === "Uniqueness"
            Follows from 5. Set $V_0 = U_0$. Then, we get that $U(\cdot, t) = V(\cdot, t)$ almost everywhere in $\R$ for all $t \ge 0$.
        === "1."
            Follows from 5. by setting $V_0 = 0$. Then, $V$ must be zero.
        === "2."
            This is the maximum principle $\eqref{eq:max_principle}$.
        === "3."
            Follows from 5. by setting $V_0(x) = U_0(x+h)$. Then, we have that

            $$\int_\R \abs{U(x, t) - U(x+h, t)} \dd x \le \int_\R \abs{U_0(x) - U_0(x+h)} \dd x$$

            for all $h > 0$. Then, dividing by $h$, and taking the limit $h \to 0$, we get the desired result.
        === "4."
            A discrete version of time continuity will be shown later for monotone finite volume schemes. The convergence of these to the entropy solution will then yield the result.
        === "5."
            Follows from 5. by taking the limit $(a, b) \to (-\infty, \infty)$.
        === "6."
            Let $U, V$ be two entropy solutions of $\eqref{eq:scalar_cons_law}$ with initial data $U_0, V_0$. Then, it follows that

            $$
            \begin{aligned}
                \partial_t \eta(U; c) + \partial_x q(U; c) &\le 0 && \forall c \in \R \\
                \partial_t \eta(d; V) + \partial_x q(d; V) &\le 0 && \forall d \in \R.
            \end{aligned}
            $$

            Thus, using the chain rule, we get

            $$
            \begin{aligned}
                \partial_t \eta(U; V) & = \partial_t \eval{\eta(U; c)}_{c=V} + \partial_t \eval{\eta(d; V)}_{d=U} \\
                & \le - \partial_x \eval{q(U; c)}_{c=V} - \partial_x \eval{q(d; V)}_{d=U} \\
                & = - \partial_x q(U; V)
            \end{aligned}
            $$

            Equivalently,

            $$\partial_t \eta(U; V) + \partial_x q(U; V) \le 0.$$

            Now, integrating over the trapezoid

            $$\{(x, t) : 0 \le t \le T, a - M(T-t) \le x \le b + M(T-t)\}$$

            we get

            $$
            \begin{aligned}
                0 \le & \int_0^T \int_{x_L(t)}^{x_R(t)} \partial_t \eta(U; V) + \partial_x q(U; V) \dd x \dd t \\
                =& -\int_{a-MT}^{b+MT} |U_0 - V_0| \dd x + \int_{a}^{b} |U(x, T) - V(x, T)| \dd x \\
                & - \int_0^T [q(U; V) - M\eta(U; V)](x_L(t), t) \dd t
                - \int_0^T [q(U; V) - M\eta(U; V)](x_R(t), t) \dd t
            \end{aligned}
            $$

            Notice that $q(U; V) - M\eta(U; V) = |U-V| \qty(\frac{q(U; V)}{|U-V|} - M) \le 0$, so we have

            $$\int_{a}^{b} |U(x, T) - V(x, T)| \dd x \le \int_{a-MT}^{b+MT} |U_0 - V_0| \dd x.$$

        === "7."
            By $\eqref{eq:scalar_cons_law}$, we have the conservation

            $$\int_\R U(x, t) \dd x = \int_\R U_0(x) \dd x \quad \forall\ t > 0$$

            Therefore, we have that

            $$\int_\R (u - V)^+ \dd x \le \int_\R (U_0 - V_0)^+ \dd x.$$

            Then, we can conclude with

            $$
            \begin{aligned}
                U_0(x) &\le V_0(x) & \implies (U_0(x) - V_0(x))^+ = 0 \\
                & \implies 0 \le \int_\R (U - V)^+ \dd x \le 0 \\
                & \implies U(\cdot, t) \le V(\cdot, t) \quad \text{a.e. in } \R.
            \end{aligned}
            $$


## Solutions for the Riemann problem for general f

We can now use $\eqref{eq:oleinik-E}$ to solve the Riemann problem for general $f \in C^1(\R)$ not necessarily convex. This condition is equivalent to 

> The chord joining $(U_L, f(U_L))$ and $(U_R, f (U_R))$ must lie below the graph of the function $f$ between these points when $U_L < U_R$. Similarly, it must lie above the graph when $U_L > U_R$.

Assume without loss of generality that $U_L < U_R$. Let $f_c$ be the lower convex envelope of $f$:

$$f_c(x):= \sup_g \qty{g(x) \mid g \text{ convex}, g \le f}$$

Now, the interval $[U_L, U_R]$ can be divided into two parts. One part is where $f_c = f$, and the other is where $f_c < f$. In the second part, $f_c$ is affine. We now use shocks in the affine region and rarefaction waves in the complement. Then, the solution is $\eqref{eq:rarefaction_solution}$ with $f$ replaced by $f_c$. For $U_L > U_R$, we instead use the upper convex envelope. Now, we often get rarefaction waves followed by shocks and vice versa. This is called a *compound wave*.

## Summary [^1]

Solutions of the conservation law $\eqref{eq:scalar_cons_law}$ may develop discontinuities or shock waves, even for smooth initial data. Consequently, weak solutions are sought. Shock speeds are computed with the Rankine–Hugoniot condition $\eqref{eq:rankine_hugoniot}$.

- Weak solutions are not necessarily unique. Entropy conditions like Oleinik’s condition E $\eqref{eq:oleinik-E}$
 have to be imposed. Self-similar continuous solutions or rarefaction waves have to be considered.
- Explicit solutions for the Riemann problem (even for non-convex fluxes) can be constructed in terms of shocks, rarefaction waves and compound shocks.
- Entropy solutions exist, are unique and are stable in $L^1$ with respect to the initial data. Furthermore, the entropy solutions satisfy an $L^\infty$ estimate, $L^p$ estimates and are Total Variation Diminishing (TVD)—that is, the total variation decreases in time.

[^1]: https://www.uio.no/studier/emner/matnat/math/MAT-IN9240/h17/pensumliste/numcl_notes.pdf
```