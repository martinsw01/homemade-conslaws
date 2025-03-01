# Boundary condition

```@docs
BoundaryCondition
```

## Dispatch functions

It is done by implementing [`homemade_conslaws.for_each_boundary_cell`](@ref) for subtypes of `BoundaryCondition`.

## Implementations

### Neumann boundary condition

```@docs
NeumannBC
```
---
### Periodic boundary condition

```@docs
PeriodicBC
```
---
### Wall boundary condition

```@docs
WallBC
```
---
### Multiple walls boundary condition

```@docs
WallsBC
```


