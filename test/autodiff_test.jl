using Random

@testset "Test promotion" begin
    a = DualNumber(1, [1])
    b = DualNumber(2., [2])
    c = 1

    abc = [a, b, c]
    @test eltype(abc) == DualNumber{1, Float64}
end


generate_random_pair(d) = ((rand(), rand(d)), (rand(), rand(d)))
numbers = [generate_random_pair(4) for _ in 1:10]

@testset "Test arithmetic" begin
    for op in (+, -)
        for (a, b) in numbers
            a_real, a_dual = a
            b_real, b_dual = b
            ab = op(DualNumber(a...), DualNumber(b...))
            ab2 = DualNumber(op(a_real, b_real), op(a_dual, b_dual))
            @test ab == ab2

            @test op(DualNumber(a...), b_real) == DualNumber(op(a_real, b_real), a_dual) == op(b_real, DualNumber(a...))
        end
    end

    @testset "Test multiplication" begin
        @test 3 * DualNumber(1, [1]) == DualNumber(3, [3])
        for (a, b) in numbers
            a_real, a_dual = a
            b_real, b_dual = b
            ab = DualNumber(a...) * DualNumber(b...)
            ab2 = DualNumber(a_real * b_real, a_real * b_dual + a_dual * b_real)
            @test ab == ab2
        end
    end

    @testset "Test division" begin
        @test 3 / DualNumber(1, [1]) == DualNumber(3, [-3])
        for (a, b) in numbers
            a_real, a_dual = a
            b_real, b_dual = b
            ab = DualNumber(a...) / DualNumber(b...)
            ab2 = DualNumber(a_real / b_real, (a_dual * b_real - a_real * b_dual) / b_real^2)
            @test ab == ab2
        end
    end

    @testset "Test powers" begin
        for _ in 1:4
            a, b = rand(1:5), rand(1:5)

            @test a^b ≈
                DualNumber(a, [0])^b ≈
                a^DualNumber(b, [0]) ≈
                DualNumber(a^b, [0]) ≈
                DualNumber(a, [0])^Float64(b)
        end
    end
end

@testset "Test functions" begin
    f(x) = sin(x + cos(x))
    df(x) = (1-sin(x)) * cos(x + cos(x))
    for (a, b) in numbers
        a_real, a_dual = a
        b_real, b_dual = b
        ab = f(DualNumber(a...))
        ab2 = DualNumber(f(a_real), df(a_real) * a_dual)
        @test ab ≈ ab2
    end
end

@testset "Test forward_diff" begin
    f(x) = sum(sin, x)
    df(x) = cos.(x)

    x = rand(4)
    fx, grad_fx = gradient(f, x)
    @test fx ≈ f(x)
    @test grad_fx ≈ df(x)
end

@testset "Test on simulation" begin
    N = 10
    gradient(ones(N)) do u0
        bc = NeumannBC()
        x_L, x_R = -1, 1
        grid = UniformGrid1D(N, bc, u0, (x_L, x_R))
        dt = grid.dx
        T = 1
        eq = BurgersEQ()
        F = LaxFriedrichsFlux()
        reconstruction = NoReconstruction()
        timestepper = ForwardEuler(grid)
        system = ConservedSystem(eq, reconstruction, F, timestepper)
        simulator = Simulator(system, grid, 0.)

        simulate!(simulator, T, dt)

        sum(cells(grid))
    end
end

@testset "Test on water simulation" begin
    N = 20
    gradient([0.]) do (u0, )
        bc = WallsBC([[0.,0.]])
        q0 = [[exp(-(x - 0.5)^2 / 0.001) .+ 1.5, u0] for x in range(-1, 1, length=N)] 
        grid = UniformGrid1DWalls(N, bc, q0, (-1, 1))
        dt = grid.dx
        T = 1

        F = CentralUpwind()

        eq = ShallowWater1D(1.)
        reconstruction = NoReconstruction()
        timestepper = ForwardEuler(grid)
        system = ConservedSystem(eq, reconstruction, F, timestepper)
        simulator = Simulator(system, grid, 0.)

        simulate!(simulator, T, dt)

        maximum(uh for (h, uh) in cells(grid))
    end
end