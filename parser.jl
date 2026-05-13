#!/usr/bin/env julia
using ArgParse

function main()
    s = ArgParseSettings(description = "Este código es una prueba para introducir argumentos desde la terminal y dar ayuda.")

    @add_arg_table! s begin
        "integers"
            help = "un entero en el rango 0..9"
            arg_type = Int
            range_tester = (x -> 0 <= x <= 9) # Equivalente a choices=range(10)
            nargs = '+'
            required = true
        "--sum"
            help = "suma los enteros (por defecto: busca el máximo)"
            dest_name = "accumulate"
            action = :store_const
            constant = sum
            default = max
    end

    args = parse_args(s)

    # En Julia, los resultados se guardan en un Diccionario (Dict)
    # Llamamos a la función (sum o max) pasando el array de enteros
    println(args["accumulate"](args["integers"]))
end

main()
