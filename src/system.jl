export ConservedSystem

"""
    ConservedSystem{EquationType <: Equation, ReconstructionType <: Reconstruction, NumericalFluxType <: NumericalFlux, TimeStepperType <: TimeStepper}

A wrapper around the objects (equation, reconstruction, numerical flux, and time stepper)
needed to solve a conservation law.
"""
struct ConservedSystem{EquationType <: Equation, ReconstructionType <: Reconstruction, NumericalFluxType <: NumericalFlux, TimeStepperType <: TimeStepper}
    eq::EquationType
    reconstruction::ReconstructionType
    numerical_flux::NumericalFluxType
    timestepper::TimeStepperType

    function ConservedSystem(equation, reconstruction, numerical_flux, timestepper)
        new{typeof(equation), typeof(reconstruction), typeof(numerical_flux), typeof(timestepper)}(
            equation, reconstruction, numerical_flux, timestepper)
    end
end