function fibonacci_julia(n)
    #= Genera la secuencia de Fibonacci hasta el n-ésimo término =#
    secuencia = Int64[] # Definimos un array de enteros vacío
    a, b = 0, 1

    for i in 1:n
        push!(secuencia, a)
        # Asignación múltiple (similar a Python)
        a, b = b, a + b
    end

    return secuencia
end

# Ejecución
n_terminos = 10
resultado = fibonacci_julia(n_terminos)
println("Julia (n=$n_terminos): $resultado")
