```@raw html
# Nonlinear hyperbolic systems in one dimension

Most equations of the form

$$
\newcommand{\U}{\mathbf{U}}
\newcommand{\f}{\mathbf{f}}
\newcommand{\r}{\mathbf{r}}
\begin{equation}
    \begin{aligned}
        \partial_t \U + \partial_x \f(\U) &= 0 \\
        \U(x, 0) &= \U_0(x)
    \end{aligned} \label{eq:nonlinear_hyperbolic_system}
\end{equation}
$$

from physics are nonlinear hyperbolic systems. We will now look at 

- structural properties and
- well-posedness,

in particular, for the Riemann problem

$$
\begin{equation}
    \begin{aligned}
        \partial_t \U + \partial_x \f(\U) &= 0 \\
        \U(x, 0) &= \begin{cases}
            \U_L, & x < 0 \\
            \U_R, & x > 0,
        \end{cases}
    \end{aligned} \label{eq:riemann_problem} \tag{RP}
\end{equation}
$$

an determine the general form of the entropy solution. We will assume $\U \in L_\text{loc}^1(\R \times \R_+, \mathcal U)$ for some $\mathcal U \subset \R^m$, for which $\eqref{eq:nonlinear_hyperbolic_system}$ makes sense.

As with the scalar case, $\U$ can develop discontinuities in finite time, so we need to interprate $\eqref{eq:nonlinear_hyperbolic_system}$ in the weak sense:

$$
\begin{equation}
    \int_{\R_+} \int_\R \U \partial_t \phi + \f(\U) \partial_x \phi \, \dd x \, \dd t
    + \int_\R \U_0(x) \phi(x, 0) \, \dd x = 0
    \label{eq:weak_nonlinear_hyperbolic_system}
\end{equation}
$$ 

for all test functions $\phi \in C_c^\infty(\R \times \R_+)$. Similarly, one can derive the rankine-hugoniot condition: If $\U$ is piecewise $C^1$ with only jump discontinuities, then $\U$ is a weak solution of $\eqref{eq:nonlinear_hyperbolic_system}$ if and only if it is a classical solution wherever it is $C^1$ and satisfies

$$
\begin{equation}
    s = \frac{[\![\f(\U)]\!]}{[\![\U]\!]}
    \label{eq:rankine_hugoniot}
    \tag{RH}
\end{equation}
$$

where $s = \gamma'(t)$ is the speed of the discontinuity. We now have $2m + 1$ unknowns, $\U_L$, $\U_R$, and $s$, but only $m$ equations. The entropy condition will deal with the last $m+1$ unknowns.

## Structural properties

To develop existence and uniqueness results, we will need to impose some assumptions on $\f$. 

!!! definition "Definition (Hyperbolic system)"
    The system of equations $\eqref{eq:nonlinear_hyperbolic_system}$ is called hyperbolic if the Jacobian matrix $\f'(\U)$ is real diagonizable for all $\U \in \mathcal U$:

    $$\f'(\U) = R(\U) \Lambda(\U) R(\U)^{-1}$$

    It is called strictly hyperbolic if

    $$\lambda_1(\U) < \lambda_2(\U) < \ldots < \lambda_m(\U).$$


!!! definition
    The $j$-th wave-family is **genuinely nonlinear** if for all $\U \in \mathcal U$,

    $$\nabla \lambda_j(\U) \cdot \r_j(\U) \neq 0.$$

    ---

    It is called **linearly degenerate** if for all $\U \in \mathcal U$,

    $$\nabla \lambda_j(\U) \cdot \r_j(\U) = 0.$$


??? example "Example (Scalar conservation law)"
    For $m=1$, the system reduces to a scalar conservation law. Then, we only have a single eigenvalue $\lambda(U) = f'(U) \in \R$ and can choose $r(U) = 1$. Thus, $\eqref{eq:nonlinear_hyperbolic_system}$ is always hyperbolic.

    Additionally,

    $$\nabla \lambda(U) \cdot r(U) = f''(U).$$

    We therefore have that it is

    - genuinely nonlinear if $f$ is strictly convex or concave,
    - linearly degenerate if $f$ is affine.

??? example "Example (Shallow water equations)"
    A model for water waves in a shallow body of water is given by

    $$\begin{aligned}
        \partial_t h + \partial_x (hv) &= 0, \\
        \partial_t (hv) + \partial_x \qty(\frac{1}{2} g h^2 + hv^2) &= 0
    \end{aligned}$$

    where the height $h$ and the momentum $m:=hv$ are the conserved quantities. We can write this as $\eqref{eq:nonlinear_hyperbolic_system}$ with

    $$\U = \begin{pmatrix} h \\ m \end{pmatrix}, \quad \f(\U) = \begin{pmatrix} m \\ \frac{1}{2} g h^2 + \frac{m^2}{h} \end{pmatrix}.$$

    This gives the Jacobian

    $$\f'(\U) = \begin{pmatrix} 0 & 1 \\ gh - \frac{m^2}{h^2} & \frac{2m}{h} \end{pmatrix},$$

    and the eigenvalues and eiegenvectors

    $$\lambda_\pm = v \pm c,\quad \r_\pm = \begin{pmatrix} 1 \\ v \pm c \end{pmatrix}.$$

    The values if $\U$ for which the equation makes sense are

    $$\mathcal U = \qty{ (h, m) \in \R^2 \mid h > 0 }.$$

    The families of waves are genuinely nonlinear, as

    $$\nabla \lambda_\pm \cdot \r_\pm = \mp \frac{3}{2} \sqrt{\frac{g}{h}} \neq 0.$$


## Simple solutions

From now on, we will assume that all the families of waves are either genuinely nonlinear or linearly degenerate. Additionally, if the $j$-th family is genuinely nonlinear, we will normalize $\r_j$ such that

$$\nabla \lambda_j(\U) \cdot \r_j(\U) = 1 \quad \forall \U \in \mathcal U$$

We will now fix $\U_L \in \mathcal U$ and find all $U_R$ such that we get a rarefaction wave or a shock wave in the Riemann problem $\eqref{eq:riemann_problem}$.

### Rarefaction waves

A rarefaction wave is a smooth solution of $\eqref{eq:nonlinear_hyperbolic_system}$ connecting $\U_L$ and $\U_R$. Recall that the solution is invariant under the transformation $(x, t) \mapsto (\alpha x, \alpha t)$ so it takes the form $\newcommand{\u}{\mathbf{u}}\U(x, t) = \u(x/t)$, where $\u \in C^1(\R, \R^m)$. Inserting into $\eqref{eq:nonlinear_hyperbolic_system}$ gives

$$\u'(\xi) f'(u(\xi)) = \xi \u'(\xi).$$

Then, we either have $\u'(\xi) = 0$ or $\xi$ is an eigenvalue of the Jacobian $f'(\u(\xi))$. In the latter case, we have that

$$\u'(\xi) = \r_j(\u(\xi)), \quad \xi = \lambda_j(\u(\xi))$$

for some $j$. Now, if there are some $\xi_L, \xi_R \in \R$ such that $\u(\xi_L) = \U_L$ and $\u(\xi_R) = \U_R$, they must be $\xi_L = \lambda_j(\U_L)$ and $\xi_R = \lambda_j(\U_R)$. The rarefaction solution of $\eqref{eq:riemann_problem}$ is then given by

$$\U(x, t) = \left\{\begin{aligned}
    \U_L, &&& \frac{x}{t} < \lambda_j(\U_L) \\
    \u\qty(\frac{x}{t}), && \lambda_j(\U_L) <  &\frac{x}{t} < \lambda_j(\U_R) \\
    \U_R, &&  \lambda_j(\U_R) < & \frac{x}{t}
\end{aligned}\right.$$

```

