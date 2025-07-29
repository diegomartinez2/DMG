using Terminals # Necesitarás este módulo para has_buffered_input

"""
    non_blocking_getc()

Intenta leer un carácter del teclado sin bloquear el hilo principal.
Devuelve el carácter si se ha pulsado uno, o `nothing` si no hay entrada.
"""
function non_blocking_getc()
    # Asegurarse de que la terminal esté en modo raw para lectura inmediata
    ret_enter_raw = ccall(:jl_tty_set_mode, Int32, (Ptr{Cvoid}, Int32), stdin.handle, true)
    ret_enter_raw == 0 || @warn "No se pudo cambiar a modo raw. La lectura podría seguir bloqueando."

    c = nothing
    if Terminals.has_buffered_input(stdin)
        c = read(stdin, Char)
    end

    # Devolver la terminal a su modo normal
    ret_exit_raw = ccall(:jl_tty_set_mode, Int32, (Ptr{Cvoid}, Int32), stdin.handle, false)
    ret_exit_raw == 0 || @warn "No se pudo restaurar el modo normal de la terminal."

    return c
end

# --- Ejemplo de uso ---
function main_loop_non_blocking()
    println("Presiona una tecla (o 'q' para salir). El programa seguirá funcionando.")

    # Asegurarse de que stdin esté disponible para la lectura en el ciclo.
    # Esto es crucial para que Terminals.has_buffered_input funcione correctamente.
    Base.start_reading(stdin)

    while true
        tecla = non_blocking_getc()

        if tecla !== nothing
            println("\n¡Se pulsó la tecla: '$tecla'!")
            if tecla == 'q'
                println("Saliendo del bucle.")
                break
            end
        end

        # Aquí es donde tu programa haría otras cosas importantes
        print(".")
        sleep(0.1) # Pausa breve para no consumir el 100% de la CPU
    end

    # Detener la lectura de stdin cuando hayas terminado.
    Base.stop_reading(stdin)
    println("Programa terminado.")
end

# Llama a la función principal para ver el ejemplo
main_loop_non_blocking()
