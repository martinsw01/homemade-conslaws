# Docs

In order to make a solve a conservation law, need to construct instances of the following types:

| Abstract type | Implementations |
| :------------ | :-------------- |
| [`Equation`](@ref)      | [`BurgersEQ`](@ref), [`ShallowWater1D`](@ref) |
| [`BoundaryCondition`](@ref) | [`PeriodicBC`](@ref), [`NeumannBC`](@ref), [`WallBC`](@ref), [`WallsBC`](@ref) |
| [`Grid`](@ref)          | [`UniformGrid1D`](@ref), [`UniformGrid1DWalls`](@ref) |
| [`Reconstruction`](@ref) | [`NoReconstruction`](@ref), [`LinearReconstruction`](@ref) |
| [`NumericalFlux`](@ref)  | [`LaxFriedrichsFlux`](@ref), [`RusanovFlux`](@ref), [`GodunovFlux`](@ref), [`CentralUpwind`](@ref) |
| [`TimeStepper`](@ref)    | [`ForwardEuler`](@ref), [`RK2`](@ref) |

gather them in a [`ConservedSystem`](@ref) and [`Simulator`](@ref) object:

```julia
system = ConservedSystem(equation, reconstruction, numerical_flux, timestepper)
simulator = Simulator(system, grid, t0)
```

and solve in place using [`simulate!`](@ref) or [`simulate_and_aggregate!`](@ref):

```julia
simulate!(simulator, T, max_dt, callbacks)
```

## Example

For example, we can solve Burgers' equation with Neumann boundary conditions and Riemann initial condition on a uniform
grid with $N=20$ cells spanning $(x_L, x_R) = (0,1)$. The simplest way to do this is is with no reconstruction, a forward
Euler time-stepper, and for example the Lax-Friedrichs numerical flux.

```@setup 1
using homemade_conslaws, homemade_conslaws.Viz
```


```@example 1
eq = BurgersEQ()
bc = NeumannBC()

N = 20
x_L, x_R = 0.0, 1.0
x = cell_centers(N, x_L, x_R)
u0(x) = x .< 0.5 ? 1.0 : 0.0
U0 = u0.(x)
grid = UniformGrid1D(N, bc, U0, (x_L, x_R))

F = LaxFriedrichsFlux()
reconstruction = NoReconstruction()
timestepper = ForwardEuler(grid)
nothing # hide
```

we gather these into a `System` object, and solve using a `Simulator` object.

```@example 1
system = ConservedSystem(eq, reconstruction, F, timestepper)
simulator = Simulator(system, grid, 0.)

U, t = simulate_and_aggregate!(simulator, 1., 0.1)
nothing # hide
```

Next, we can animate the solution

```@example 1
u(x, t) = u0(x-0.5t)
Viz.animate_solution(U', u, grid, t, 2)
```

## Extending the functionallity

The package `homemade_conslaws` is designed to be easily extensible, making heavy use of multiple dispatch. This is done by implementing the abstract types in the [table](#docs) above. 
