# Second order schemes

```@raw html

$$
\begin{equation}
    U_t + f(U)_x = 0 \label{eq:conservation_law}
\end{equation}
$$


The schemes we have seen this far has slow convergence. We can improve the convergence by using higher order schemes.

## Order of convergence

For a $(2p+1)$-point scheme with update function $H$ and $\lambda \equiv \frac{\Dt}{\Dx}$, we define the truncation error as usual as

$$\tau_j^n := U(x_j, t_{n+1}) - H(U_{j-p}^n, \dots, U_{j+p}^n).$$

We call it $q$-th order accurate if $q \in \N$ is the largest integer such that

$$\tau_j^n = \O(\Dt^{q+1}) \quad \forall j, n.$$


!!! lemma
    If $U \in C^2$ is the exact solution and $H$ is a consistent, conservative $3$-point FVM, then the scheme is at least first order accurate.

    ??? proof
        Follows as usual from Taylor expansion.

!!! example
    - Both Lax-Friedrichs and Engquist-Osher use $C^2$ numerical fluxes, and are thus first order accurate.
    - Godunov is not, so we cannot use the lemma on it. However, one can still show that it is of first order accuracy.


Kuznetsov [GR91] shows that monotone schemes has an order of convergence of at least $1/2$:

$$\norm{U^{\Dx} - U}_{L^1(\R)} \le C \Dx^{1/2}.$$


### Lax-Wendroff scheme

For weak solutions, one cannot use Taylor expansions to find the error. However, we can use it to design schemes that formally are better than first order. Let $U$ be a smooth solution. Then, we have

$$U_tt = \qty(f'(U)U_x)_x.$$

Inserting into the taylor expansion yields

$$U(x_j, t_{n+1}) = U(x_j, t_n) + \Dt U_t(x_j, t_n) + \frac{\Dt^2}{2} \qty(f'(U)U_x)_x + \O(\Dt^3).$$

The last two term are approximated by

$$
\begin{aligned}
    f(U)_x &= \frac{f(U_{j+1}) - f(U_{j-1})}{2\Dx} + \O(\Dx^3) \\
    \qty(f'(U)U_x)_x &= \frac{1}{\Dx} \qty(a_{j+1/2}^n \frac{f(U_{j+1}^n) - f(U_j^n)}{\Dx}
    - a_{j-1/2}^n \frac{f(U_j^n) - f(U_{j-1}^n)}{\Dx}) + \O(\Dx^3),
\end{aligned}
$$

yielding the Lax-Wendroff scheme

$$U_j^{n+1} = U_j^n - \frac{\Dt}{2\Dx} \Big(f(U_{j+1}^n) - f(U_{j-1}^n)\Big)
+ \frac{\Dt^2}{2\Dx^2} \qty(a_{j+1/2}^n \frac{f(U_{j+1}^n) - f(U_j^n)}{\Dx}
- a_{j-1/2}^n \frac{f(U_j^n) - f(U_{j-1}^n)}{\Dx}),$$

where $a_{j+1/2}^n = f'\qty(\frac{U_j^n + U_{j+1}^n}{2})$ is an approximation of $f'(U_{j+1/2}^n)$. Alternatively, the numerical flux can be written as

$$F_{j+1/2}^n = F^{\text{LxW}}(U_j^n, U_{j+1}^n) = \frac{f(U_j^n) + f(U_{j+1}^n)}{2} - a_{j+1/2}^n \frac{\Dt}{2\Dx} \qty(f(U_{j+1}^n) - f(U_j^n)).$$

Numerical experiments show that the Lax-Wendroff scheme is second order accurate for smooth soluitons. Additionally, the discontinuities are sharper than for first order schemes. However, there are heavy oscillations near the shocks.


## REA Algorithm

The derivation of the schemes presented in the last chapter went as follows:

- **R**econstruction: Recreate $U(x,t^n)$ from the cell averages $U_j^n$. Prevoiusly, we used piecewise constant functions.
- **E**volution: Solve the superposition of Riemann problems. This can be done approximately (LxF) or exactly (Godunov).
- **A**veraging: Average the solutions for each cell.

The difference in many of the schemes presented is the evolution step. However, one can also increase the accuracy of the reconstruction step. One can for example use higher order interpolation.

### Second order reconstruction

To determine the piecewise affine function, we need some requirements. One can show that the last two steps are conservative. It would be natural to require that the reconstruction too. For the interpolation function $p_j$ in the cell $\Cell_j$, we thus require that

$$p_j(x) = U_j^n + \sigma_j^n(x - x_j).$$

or instance, one can apply the standard slope approximation using central, backward, or forward difference methods. However, in the case of the linear advection equation with constant velocity, it can be demonstrated that using the forward difference approach leads to the Lax-Wendroff scheme. Oscillations also arise using the two other methods.

### Source of oscillations

The solutions of $\eqref{eq:conservation_law}$ are TVD. Thus, they are not oscillatory whenever the initial data is not. Finite volume methods should therefore not introduce oscillations. The averaging operator and suitable Riemann solvers respects TVD. Thus, the source must be the reconstruction step. Consider for example the cell averages

$$U_j^n = \begin{cases}
    1 & j \le 0 \\
    0 & j > 0.
\end{cases}$$

It has total variation $\norm{U^n}_{TV} = 1$. Using the downwind reconstruction, we get

$$\sigma_j^n = \begin{cases}
    0 & j \neq 0 \\
    -1/2 & j = 0.
\end{cases}$$

Thus, the reconstruction overshoots and reaches $3/2$. The total variation is now $\norm{p^n}_{TV} = 3/2$. Similar arguments can be made for the other reconstruction methods. Therefore, we should require that the reconstruction operator is TVD: $\norm{p}_{BV} \le \norm{U^{\Dx}}_{BV}$.


## Minmod limiter

One such TVD reconstruction operator is the minmod limiter. The problem above at discontinuities occurs whenever the downwind and upwind slopes have opposite signs. One can therefore default to constant interpolation, that is, a slope of zero. This results in the minmod limiter:

$$\sigma_j^n = \minmod\qty(\frac{U_{j+1}^n - U_j^n}{\Dx}, \frac{U_j^n - U_{j-1}^n}{\Dx}),$$

where

$$\minmod\{a_k\}_k = \begin{cases}
    \sign(a_1) \min_k\abs{a_k}, & \sign(a_1) = \dots = \sign(a_n) \\
    0, & \text{otherwise}.
\end{cases}$$

```