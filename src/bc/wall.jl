export WallBC

"""
    WallBC

Implements a wall boundary condition, in the sense that the water height is reflected and
momentum is reflected and negating, simulating an opposite and equal wave.
"""
struct WallBC <: BoundaryCondition end