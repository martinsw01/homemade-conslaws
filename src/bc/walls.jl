using StaticArrays

export WallsBC

"""
    WallsBC{Float <: AbstractFloat}

Simulates walls starting and stopping at cell interfaces. Equivalent to [`homemade_conslaws.WallBC`](@ref) when no walls are given.
"""
struct WallsBC{Float <: AbstractFloat} <: BoundaryCondition
    walls::Vector{SVector{2, Float}}
    function WallsBC(walls::AbstractVector{Wall}) where Wall
        new{eltype(Wall)}([SVector(wall...) for wall in walls])
    end
    function WallsBC()
        WallsBC(
            Vector{SVector{2, Float64}}(undef, 0), 
        )
    end
end