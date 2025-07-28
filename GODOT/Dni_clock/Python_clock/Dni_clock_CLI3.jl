import Pkg; # Ya no necesitas Pkg.add("Tk") si ya no lo usas.
using Dates # Para fechas y horas
using Printf # Para formato de cadenas
using Terminals # Para control avanzado de terminal (modo raw, etc.)

# -------
# Constantes D'ni
# -------
const DNISH_VAILEE_DAYS = 29
const DNISH_HAHR_VAILEE = 10 # Number of Vailee in a Hahr
const DNISH_EPOCH_HAHR = 9647

const DNISH_YR_TO_GREGORIAN_DAY_RATIO = 1.25945582758621
const GREGORIAN_TO_DNISH_YR_RATIO = 0.793993705929756
const JULIAN_DAY_OFFSET_DNISH_EPOCH = 727249.704166666

# D'ni time units per Yahr
const DNISH_YAHR_UNITS = 78125.0
const GARHTAHVO_UNIT_VAL = 15625.0
const TAHVO_UNIT_VAL = 625.0
const GORAHN_UNIT_VAL = 25.0
const PRORAHN_UNIT_VAL = 1.0

const VAILEE_NAMES = Dict(
    1 => "Leefo", 2 => "Leebro", 3 => "Leesahn", 4 => "Leetar", 5 => "Leevot",
    6 => "Leevofo", 7 => "Leevobro", 8 => "Leevosahn", 9 => "Leevotar", 10 => "Leenovoo"
)

# ----------
# Funciones de Utilidad
# ----------

"""
    to_digits_base_n(number::Int, base::Int)

Convert a positive integer `number` to its digit representation in `base`.
"""
function to_digits_base_n(number::Int, base::Int)::Vector{Int}
    if number == 0
        return [0]
    end
    if number < 0
        throw(ArgumentError("Number must be non-negative."))
    end
    if base < 2
        throw(ArgumentError("Base must be an integer greater than or equal to 2."))
    end

    digits = Int[]
    temp_n = number
    while temp_n > 0
        pushfirst!(digits, temp_n % base) # pushfirst! adds to the beginning
        temp_n = div(temp_n, base)
    end
    return digits
end

"""
    from_digits_base_n(digits::Vector{Int}, base::Int)

Compute the number given by `digits` in `base`.
"""
function from_digits_base_n(digits::Vector{Int}, base::Int)::Int
    if any(d < 0 || d >= base for d in digits)
        throw(ArgumentError("Digits must be integers within the range [0, base-1]."))
    end
    if base < 2
        throw(ArgumentError("Base must be an integer greater than or equal to 2."))
    end

    number = 0
    for d in digits
        number = base * number + d
    end
    return number
end

# ----------
# Conversiones de Fecha y Hora
# ----------

"""
    gregorian_to_julian(day::Int, month::Int, year::Int, hour::Int=0, minute::Int=0, second::Int=0)

Converts a Gregorian date and time to a Julian Day number.
"""
function gregorian_to_julian(day::Int, month::Int, year::Int, hour::Int=0, minute::Int=0, second::Int=0)::Float64
    adj_year = year
    adj_month = month

    if adj_year < 0
        adj_year = -(adj_year - 1)
    end

    if adj_month < 3
        adj_month += 12
        adj_year -= 1
    end

    # Whole days
    wd = day + floor(Int, ((153 * adj_month) - 457) / 5) +
             floor(Int, 365.25 * adj_year) - floor(Int, 0.01 * adj_year) + floor(Int, 0.0025 * adj_year)

    # Fractional day
    fd = ((hour * 3600) + (minute * 60) + second) / 86400.0

    return Float64(wd) + fd
end

"""
    julian_to_gregorian(jd::Float64)

Converts a Julian Day number to a Gregorian date and time.
Returns (day, month, year, hour, minute, second).
"""
function julian_to_gregorian(jd::Float64)::Tuple{Int, Int, Int, Int, Int, Int}
    z = floor(Int, jd)
    f = jd - z # Fractional part for time calculation

    g = z - 0.25
    a = floor(Int, g / 36524.25)
    b = a - floor(Int, (0.25 * a))
    year = floor(Int, (g + b) / 365.25)
    c = z + b - floor(Int, 365.25 * year)

    month = floor(Int, ((5 * c) + 456) / 153)
    day = c - floor(Int, ((153 * month) - 457) / 5)

    if month > 12
        year += 1
        month -= 12
    end

    if year < 1
        year = 1 - year
    end

    # Time calculation
    total_seconds = f * 86400.0
    hour = floor(Int, total_seconds / 3600.0)
    remainder_seconds = total_seconds - (hour * 3600.0)
    minute = floor(Int, remainder_seconds / 60.0)
    second = round(Int, remainder_seconds - (minute * 60.0))

    return Int(day), Int(month), Int(year), Int(hour), Int(minute), Int(second)
end

"""
    cavernian_to_atrian_yahr(yahr::Int, vailee::Int, hahr::Int, gahrtahvo::Int=0, tahvo::Int=0, gorahn::Int=0, prorahn::Int=0)

Converts a Cavernian (D'ni) date and time to Atrian Yahr (decimal).
"""
function cavernian_to_atrian_yahr(yahr::Int, vailee::Int, hahr::Int, gahrtahvo::Int=0, tahvo::Int=0, gorahn::Int=0, prorahn::Int=0)::Float64
    # Whole yahrtee (days)
    whole_yahrtee = yahr + ((vailee - 1) * DNISH_VAILEE_DAYS) +
                    ((hahr - DNISH_EPOCH_HAHR) * (DNISH_VAILEE_DAYS * DNISH_HAHR_VAILEE))

    # Fractional yahr (time of yahr)
    fractional_yahr = ((gahrtahvo * GARHTAHVO_UNIT_VAL) +
                       (tahvo * TAHVO_UNIT_VAL) +
                       (gorahn * GORAHN_UNIT_VAL) +
                       (prorahn * PRORAHN_UNIT_VAL)) / DNISH_YAHR_UNITS

    return Float64(whole_yahrtee) + fractional_yahr
end

"""
    atrian_yahr_to_cavernian(ay::Float64)

Converts Atrian Yahr (decimal) to a Cavernian (D'ni) date and time.
Returns (hahr, vailee, yahr, gahrtahvo, tahvo, gorahn, prorahn).
"""
function atrian_yahr_to_cavernian(ay::Float64)::Tuple{Int, Int, Int, Int, Int, Int, Int}
    z = floor(Int, ay)

    # Date calculation
    g = z - 0.25 # Similar offset as in Julian conversion
    a = floor(Int, g / (DNISH_VAILEE_DAYS * DNISH_HAHR_VAILEE))
    hahr = DNISH_EPOCH_HAHR + a
    c = z - (a * (DNISH_VAILEE_DAYS * DNISH_HAHR_VAILEE))

    vailee = floor(Int, (c - 0.25) / DNISH_VAILEE_DAYS) + 1
    yahr = c - ((vailee - 1) * DNISH_VAILEE_DAYS)

    # Time calculation
    total_time_units = round(Int, (ay - floor(ay)) * DNISH_YAHR_UNITS)
    gahrtahvo = floor(Int, total_time_units / GARHTAHVO_UNIT_VAL)
    remainder = total_time_units - (gahrtahvo * GARHTAHVO_UNIT_VAL) # <-- CORREGIDO: GARHTAHVO_UNIT_VAL
    tahvo = floor(Int, remainder / TAHVO_UNIT_VAL)
    remainder = remainder - (tahvo * TAHVO_UNIT_VAL)
    gorahn = floor(Int, remainder / GORAHN_UNIT_VAL)
    prorahn = remainder - (gorahn * GORAHN_UNIT_VAL)

    return Int(hahr), Int(vailee), Int(yahr), Int(gahrtahvo), Int(tahvo), Int(gorahn), Int(prorahn)
end

"""
    cavernian_to_gregorian(ay::Float64)

Converts Atrian Yahr to Gregorian date and time.
Returns (day, month, year, hour, minute, second).
"""
function cavernian_to_gregorian(ay::Float64)::Tuple{Int, Int, Int, Int, Int, Int}
    ayd = ay - 1.0 # Atrian Yahr Day since epoch start

    jdd = ayd * DNISH_YR_TO_GREGORIAN_DAY_RATIO
    jd = jdd + JULIAN_DAY_OFFSET_DNISH_EPOCH

    return julian_to_gregorian(jd)
end

"""
    gregorian_to_cavernian(day::Int, month::Int, year::Int, hour::Int=0, minute::Int=0, second::Int=0)

Converts Gregorian date and time to Cavernian (D'ni) date and time.
Returns (hahr, vailee, yahr, gahrtahvo, tahvo, gorahn, prorahn).
"""
function gregorian_to_cavernian(day::Int, month::Int, year::Int, hour::Int=0, minute::Int=0, second::Int=0)::Tuple{Int, Int, Int, Int, Int, Int, Int}
    jd = gregorian_to_julian(day, month, year, hour, minute, second)
    jdd = jd - JULIAN_DAY_OFFSET_DNISH_EPOCH
    ayd = jdd * GREGORIAN_TO_DNISH_YR_RATIO
    ay = ayd + 1.0 # Adjust to make 1 AY the epoch start

    return atrian_yahr_to_cavernian(ay)
end

# ----------
# Interfaz de Terminal
# ----------

function main_julia()
    base_display = 25 # Base de visualización para números D'ni

    # Códigos de escape ANSI:
    # \033[2J   - Borra toda la pantalla
    # \033[H    - Mueve el cursor a la posición Home (0,0)
    # \033[<N>A - Mueve el cursor N líneas hacia arriba
    # \033[<N>B - Mueve el cursor N líneas hacia abajo
    # \033[<N>;<M>H - Mueve el cursor a la fila N, columna M
    # \033[K    - Borra desde el cursor hasta el final de la línea

    # Configuración del terminal para lectura de un solo carácter
    term = Terminals.TTYTerminal("xterm", stdin, stdout, stderr)
    Terminals.raw!(term, true) # Entra en modo raw para lectura de caracteres directos
    Base.exit_on_sigint(false) # Deshabilita Ctrl+C para que 'q' sea la única forma de salir

    try
        # Borra la pantalla y mueve el cursor al inicio para un arranque limpio
        print("\033[2J\033[H")

        # Imprime encabezados estáticos
        println("------------------------------------")
        println("  Convertidor de Reloj D'ni (Terminal)   ")
        println("------------------------------------")
        println("\nFecha y Hora Actual (Gregoriana):")
        println("\nFecha y Hora D'ni:")
        println("\nPresiona 'q' para salir.") # Mensaje de salida actualizado

        # Definir las filas para la salida dinámica
        const GREG_DATE_ROW = 6
        const GREG_TIME_ROW = 7
        const DNI_LINE1_ROW = 10
        const DNI_LINE2_ROW = 11
        const EXIT_MESSAGE_ROW = 13 # Fila del mensaje "Presiona 'q' para salir."

        # Columna inicial para la mayoría del texto dinámico
        const START_COL = 1

        should_exit = false
        while !should_exit
            now = Dates.now()
            day, month, year = Dates.day(now), Dates.month(now), Dates.year(now)
            hour, minute, second = Dates.hour(now), Dates.minute(now), Dates.second(now)

            hahr, vailee, yahr, gahrtahvo, tahvo, gorahn, prorahn =
                gregorian_to_cavernian(day, month, year, hour, minute, second)

            # Formatear la salida Gregoriana
            greg_date_str = @sprintf("Fecha: %d-%02d-%02d", year, month, day)
            greg_time_str = @sprintf("Hora: %02d:%02d:%02d", hour, minute, second)

            # Formatear la salida D'ni
            hahr_dni_digits_str = join(to_digits_base_n(hahr, base_display))
            yahr_dni_digits_str = join(to_digits_base_n(yahr, base_display))
            gahrtahvo_dni_digits_str = join(to_digits_base_n(gahrtahvo, base_display))
            tahvo_dni_digits_str = join(to_digits_base_n(tahvo, base_display))
            gorahn_dni_digits_str = join(to_digits_base_n(gorahn, base_display))
            prorahn_dni_digits_str = join(to_digits_base_n(prorahn, base_display))

            vailee_name = get(VAILEE_NAMES, vailee, "Desconocido")

            dni_output_line1 = @sprintf("Hahr: %s / Yahr: %s / Vailee: %d (%s)",
                                        hahr_dni_digits_str, yahr_dni_digits_str, vailee, vailee_name)
            dni_output_line2 = @sprintf("Gahrtahvo: %s / Tahvo: %s / Gorahn: %s / Prorahn: %s",
                                        gahrtahvo_dni_digits_str, tahvo_dni_digits_str, gorahn_dni_digits_str, prorahn_dni_digits_str)

            # --- Actualizar la pantalla ---
            # Fecha Gregoriana
            print("\033[$(GREG_DATE_ROW);$(START_COL)H") # Mueve el cursor a la fila 6, columna 1
            print("\033[K") # Borra la línea
            print(greg_date_str)

            # Hora Gregoriana
            print("\033[$(GREG_TIME_ROW);$(START_COL)H") # Mueve el cursor a la fila 7, columna 1
            print("\033[K") # Borra la línea
            print(greg_time_str)

            # Línea 1 D'ni
            print("\033[$(DNI_LINE1_ROW);$(START_COL)H") # Mueve el cursor a la fila 10, columna 1
            print("\033[K") # Borra la línea
            print(dni_output_line1)

            # Línea 2 D'ni
            print("\033[$(DNI_LINE2_ROW);$(START_COL)H") # Mueve el cursor a la fila 11, columna 1
            print("\033[K") # Borra la línea
            print(dni_output_line2)

            flush(stdout) # Fuerza la actualización de la pantalla

            # Comprobar si se ha presionado 'q'
            if bytesavailable(stdin) > 0
                char_input = read(stdin, Char)
                if char_input == 'q'
                    should_exit = true
                end
            end

            sleep(1) # Espera 1 segundo
        end

    finally
        # Restaurar el terminal a su estado normal al salir
        Terminals.raw!(term, false)
        Base.exit_on_sigint(true) # Reactivar Ctrl+C

        # Limpiar la pantalla y dejar un mensaje final
        print("\033[2J\033[H") # Borra toda la pantalla y mueve el cursor al inicio
        println("------------------------------------")
        println("  Convertidor de Reloj D'ni (Terminal)   ")
        println("------------------------------------")
        println("\nPrograma terminado limpiamente.                   ")
        flush(stdout)
    end
end

# Ejecutar la aplicación
main_julia()
