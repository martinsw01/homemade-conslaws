# Homemade Conslaws

## Conservation laws

Let ``\boldsymbol U`` be a quantity defined on a domain ``\Omega \subset \R^n``. For any subdomain ``\omega \subset \Omega``, the temporal rate of change of ``\bm U`` is equal to the amount of ``\bm U`` created or destroyed and the flux going through the boundary ``\partial \omega``. It can be described mathematically by

```math
    \dv{t} \int_{\omega} \bm U \dd \x
    = - \int_{\partial \omega} \bm F \cdot \bm \nu \dd \s
    + \int_{\omega} \bm S \dd \x,
```

where ``\bm F`` is the flux and ``\bm S`` is the source term. By the Gauss divergence theorem, we can rewrite this as

```math
    \dv{t} \int_{\omega} \bm U \dd \x
    + \int_{\omega} \div {\bm F} \dd \x
    = \int_{\omega} \bm S \dd \x.
```

Since this equation holds for all subdomains ``\omega \subseteq \Omega``, we can write

```math
    \bm U_t + \div {\bm F} = \bm S \quad \text{ in } \Omega \times \R_+
```

We call this equation a conservation law.