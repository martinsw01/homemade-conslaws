using homemade_conslaws
using Test

@testset "Conslaws tests" begin
    @testset "UniformGrid1D tests" begin
        include("uniformgrid1d_tests.jl")
    end
    @testset "Simulator tests" begin
        include("simulator_test.jl")
    end
    @testset "Runge kutta 2 tests" begin
        include("rk2_test.jl")
    end
    @testset "Linear reconstruction tests" begin
        include("linear_reconstruction_test.jl")
    end
    @testset "Linear system tests" begin
        include("linear_system_test.jl")
    end
    @testset "Test multiple walls" begin
        include("multiple_walls_test.jl")
    end
end

@testset "Autodiff tests" begin
    @testset "Test autodiff" begin
        include("autodiff_test.jl")
    end
end

nothing