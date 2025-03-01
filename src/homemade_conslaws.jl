module homemade_conslaws

include("equations/equation.jl")
include("bc/boundary_condition.jl")
include("grids/grid.jl")
include("reconstruction/reconstruction.jl")
include("timesteppers/timestepper.jl")

include("numerical_fluxes/numerical_flux.jl")
include("system.jl")

include("simulator.jl")

include("autodiff/dual_numbers.jl")
include("autodiff/forward_diff.jl")

include("viz/viz.jl")

end