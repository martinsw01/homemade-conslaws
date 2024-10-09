using homemade_conslaws
using Test

@testset "Conslaws tests" begin
    @testset "UniformGrid1D tests" begin
        include("uniformgrid1d_tests.jl")
    end
    @testset "Simulator tests" begin
        include("simulator_test.jl")
    end
end

nothing