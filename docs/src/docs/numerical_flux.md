# Numerical flux

```@docs
NumericalFlux
```

## Dispatch functions

An implementation of `NumericalFlux` should provide the following functions:

```@autodocs
Modules = [homemade_conslaws]
Pages = ["numerical_fluxes/numerical_flux.jl"]
Order = [:function]
```

## Implementations

### Lax-Friedrichs flux

```@docs
LaxFriedrichsFlux
```
---
### Rusanov flux

```@docs
RusanovFlux
```
---
### Godunov flux

```@docs
GodunovFlux
```
---
### Central upwind flux

```@docs
CentralUpwind
```
