using LinearAlgebra

function w_to_a(w, T, n)
    a = zeros(n)
    if T == 0.0
        a .= sqrt(1.0 ./ (2.0 .* w))
    else
        a .= sqrt.((1.0 ./ tanh.(0.5 .* w .* 315774.65221921849 ./ T)) ./ (2.0 .* w))
    end
    return a
end

function w_to_da(w, T, n)
    da = zeros(n)
    if T == 0.0
        da .= -sqrt.(1.0 ./ (8.0 .* w.^3.0))
    else
        beta = 315774.65221921849 / T
        a = w .* beta .+ sinh.(w .* beta)
        b = sqrt.(1.0 ./ (32.0 .* w.^3.0 .* (sinh.(0.5 .* w .* beta).^3.0) .* cosh.(0.5 .* w .* beta)))
        da .= -a .* b
    end
    return da
end

function dW_f0_u0(w, T)
    if T == 0.0
        return 0.25
    else
        return 0.25 * (2.0 * nb(w, T) + 1.0 + 2.0 * (315774.65221921849 / T) * w * exp(w * 315774.65221921849 / T) * nb(w, T)^2.0)
    end
end

function nb(w, T)
    if T == 0.0
        return 0.0
    else
        return 1.0 ./ (exp.(w * 315774.65221921849 / T ) .- 1.0)
    end
end
