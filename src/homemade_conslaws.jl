module homemade_conslaws

include("types.jl")
include("bc/periodic.jl")
include("bc/neumann.jl")
include("bc/wall.jl")
include("bc/walls.jl")
include("equations/burgers.jl")
include("equations/linear_advection.jl")
include("equations/shallow_water.jl")
include("grids/uniform1d.jl")
include("grids/uniform1dwalls.jl")
include("system.jl")
include("timesteppers/forwardeuler.jl")
include("timesteppers/rk2.jl")
include("numerical_fluxes/laxfriedrichs.jl")
include("numerical_fluxes/central_upwind.jl")
include("numerical_fluxes/godunov.jl")
include("numerical_fluxes/rusanov.jl")
include("reconstruction/noreconstruction.jl")
include("reconstruction/linear_reconstruction.jl")
include("simulator.jl")

include("autodiff/dual_numbers.jl")
include("autodiff/forward_diff.jl")

include("viz/viz.jl")

end