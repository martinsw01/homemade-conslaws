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
            a, b = rand(1:10), rand()
            @test DualNumber(a, [0])^b ≈ DualNumber(a^b, [0])
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