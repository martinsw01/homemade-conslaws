export ConservedSystem

struct ConservedSystem{EquationType <: Equation, ReconstructionType <: Reconstruction, NumericalFluxType <: NumericalFlux, TimeStepperType <: TimeStepper} <: System
    eq::EquationType
    reconstruction::ReconstructionType
    numerical_flux::NumericalFluxType
    timestepper::TimeStepperType

    function ConservedSystem(equation, reconstruction, numerical_flux, timestepper)
        new{typeof(equation), typeof(reconstruction), typeof(numerical_flux), typeof(timestepper)}(
            equation, reconstruction, numerical_flux, timestepper)
    end
end