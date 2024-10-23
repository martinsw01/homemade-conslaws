export ConservedSystem

struct ConservedSystem{EquationType <: Equation, GridType <: Grid} <: System
    eq::EquationType
    grid::GridType
    numerical_flux::NumericalFlux
end