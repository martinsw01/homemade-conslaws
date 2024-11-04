module homemade_conslaws
include("finite_volumes.jl")
include("central_difference.jl")
include("upwind.jl")
include("viz.jl")

include("types.jl")
include("bc/periodic.jl")
include("bc/neumann.jl")
include("bc/wall.jl")
include("equations/burgers.jl")
include("equations/linear_advection.jl")
include("equations/shallow_water.jl")
include("grids/uniform1d.jl")
include("system.jl")
include("timesteppers/forwardeuler.jl")
include("timesteppers/rk2.jl")
include("numerical_fluxes/laxfriedrichs.jl")
include("numerical_fluxes/central_upwind.jl")
include("reconstruction/noreconstruction.jl")
include("reconstruction/linear_reconstruction.jl")
include("simulator.jl")



end