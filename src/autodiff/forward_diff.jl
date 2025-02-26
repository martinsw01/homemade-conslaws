export gradient, gradient!

gradient(f, x) = gradient!(f, similar(x), x)


function gradient!(f, grad_out, x)
    d = length(grad_out)
    dual_numbers = similar(x, DualNumber{d, eltype(x)})
    for i in eachindex(x)
        dual_numbers[i] = DualNumber(x[i], 1:d .== i)
    end
    result = f(dual_numbers)
    for i in eachindex(x)
        grad_out[i] = result.dual[i]
    end
    result.real, grad_out
end