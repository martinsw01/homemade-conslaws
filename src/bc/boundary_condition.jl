export BoundaryCondition

"""
    BoundaryCondition

Describes the boundary conditions for a conservation law. It is passed to the grid, determining how to handle the boundaries.
"""
abstract type BoundaryCondition end


include("periodic.jl")
include("neumann.jl")
include("wall.jl")
include("walls.jl")