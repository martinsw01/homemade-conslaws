module homemade_conslaws
include("finite_volumes.jl")
include("central_difference.jl")
include("upwind.jl")
include("viz.jl")

include("types.jl")
include("bc/periodic.jl")
include("bc/neumann.jl")
include("grids/uniform1d.jl")
include("timesteppers/forwardeuler.jl")
include("timesteppers/rk2.jl")
include("equations/burgers.jl")
include("numerical_fluxes/laxfriedrichs.jl")
include("system.jl")
include("simulator.jl")



end