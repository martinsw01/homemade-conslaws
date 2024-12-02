using StaticArrays

export WallsBC

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