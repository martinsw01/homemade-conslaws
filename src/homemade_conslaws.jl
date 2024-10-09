module homemade_conslaws
# include("finite_volumes.jl")
# include("central_difference.jl")
# include("upwind.jl")
include("viz.jl")

include("types.jl")
include("timesteppers/forwardeuler.jl")
include("equations/burgers.jl")
include("bc/periodic.jl")
include("bc/neumann.jl")
include("grids/uniform1d.jl")
include("numerical_fluxes/laxfriedrichs.jl")
include("system.jl")
include("simulator.jl")


# using ElasticArrays

# function main()

#     u0(x) = x < 0

#     N = 50
#     eq = BurgersEQ()
#     grid = UniformGrid1D(N, NeumannBC(), u0, (-1., 1.))

#     F = LaxFriedrichsFlux()
#     system = ConservedSystem(eq, grid, F)

#     simulator = Simulator(system, grid, 0.)

#     U = ElasticMatrix(reshape(grid.cells .+ 0., :, 1))
#     t = ElasticVector([0.])

#     function add_state(simulator)
#         append!(U, simulator.grid.cells)
#         append!(t, simulator.t[])
#     end


#     function print_time(simulator)
#         # println("t = ", simulator.t[])
#     end

#     function print_state(simulator)
#         # println("u = ", simulator.grid.cells)
#     end

#     simulate!(simulator, 1., 0.5, [print_state, print_time, add_state])

#     U, t
#     # Viz.animate_solution(U, (x, t) -> u0(x - 0.5t),
#     #                      range(-1, 1, length=N+1), t)
# end

# # main()

end