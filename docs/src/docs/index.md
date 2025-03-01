# Docs

In order to make a solve a conservation law, we need to specify the following:

- What equation we are solving,
- on what kind of grid,
- what initial conditions to use,
- what boundary conditions to use,
- what numerical flux to use,
- how to reconstruct the solution,
- how to evolve the solution in time.


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
reconstruction = NoReconstruction(grid)
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

The `homemade_conslaws` package is designed to be easily extensible. This is done by extending the abstract types


- [`homemade_conslaws.Equation`](@ref)
- [`homemade_conslaws.Grid`](@ref)
- [`homemade_conslaws.Reconstruction`](@ref)
- [`homemade_conslaws.NumericalFlux`](@ref)
- [`homemade_conslaws.TimeStepper`](@ref)
- [`homemade_conslaws.BoundaryCondition`](@ref)
