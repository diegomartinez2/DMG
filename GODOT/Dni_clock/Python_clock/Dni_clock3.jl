import Pkg; Pkg.add("Tk")

# Dni_clock.jl
#
# Copyright 2025 Diego Martinez Gutierrez <diego.martinez@ehu.eus>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY and FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
# MA 02110-1301, USA.
#
# D'ni clock converter from Gregorian to D'ni
# ---------------------------

using Dates # Para fechas y horas
using Printf # Para formato de cadenas
using Tk # Para la interfaz gráfica

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
    remainder = total_time_units - (gahrtahvo * GARHTAHVO_UNIT_VAL)
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
# Interfaz Gráfica (Tk.jl)
# ----------

mutable struct DniClockApp
    root::Tk.Widget # Changed from Tk.Toplevel to Tk.Widget
    gregorian_date_label::Tk.Widget
    gregorian_time_label::Tk.Widget
    dni_date_label::Tk.Widget
    base_display::Int

    function DniClockApp(master::Tk.Toplevel) # This remains Tk.Toplevel for the argument type
        t = new()
        t.root = master
        Tk.title(master, "D'ni Clock Converter")

        t.base_display = 25 # Default display base for D'ni numbers

        # Labels for Gregorian time
        Tk.pack(Tk.label(master, text="Fecha y Hora Actual (Gregoriana)", font=("Arial", 14, "bold")), pady=5)
        t.gregorian_date_label = Tk.label(master, text="", font=("Arial", 12))
        Tk.pack(t.gregorian_date_label)
        t.gregorian_time_label = Tk.label(master, text="", font=("Arial", 12))
        Tk.pack(t.gregorian_time_label)

        # Labels for D'ni time
        Tk.pack(Tk.label(master, text="\nFecha y Hora D'ni", font=("Arial", 14, "bold")), pady=5)
        t.dni_date_label = Tk.label(master, text="", font=("Arial", 12))
        Tk.pack(t.dni_date_label)

        update_time(t) # Call the update function for the first time
        return t
    end
end

function update_time(app::DniClockApp)
    now = Dates.now()
    day, month, year = Dates.day(now), Dates.month(now), Dates.year(now)
    hour, minute, second = Dates.hour(now), Dates.minute(now), Dates.second(now)

    hahr, vailee, yahr, gahrtahvo, tahvo, gorahn, prorahn =
        gregorian_to_cavernian(day, month, year, hour, minute, second)

    # Update Gregorian labels
    Tk.configure(app.gregorian_date_label, text=@sprintf("Fecha: %d-%02d-%02d", year, month, day))
    Tk.configure(app.gregorian_time_label, text=@sprintf("Hora: %02d:%02d:%02d", hour, minute, second))

    # Format D'ni numbers in base 25 for display
    # Join digits as string
    hahr_dni_digits_str = join(to_digits_base_n(hahr, app.base_display))
    yahr_dni_digits_str = join(to_digits_base_n(yahr, app.base_display))
    gahrtahvo_dni_digits_str = join(to_digits_base_n(gahrtahvo, app.base_display))
    tahvo_dni_digits_str = join(to_digits_base_n(tahvo, app.base_display))
    gorahn_dni_digits_str = join(to_digits_base_n(gorahn, app.base_display))
    prorahn_dni_digits_str = join(to_digits_base_n(prorahn, app.base_display))

    vailee_name = get(VAILEE_NAMES, vailee, "Desconocido") # Handle unknown Vailee

    Tk.configure(app.dni_date_label,
        text=@sprintf("Hahr: %s / Yahr: %s / Vailee: %d (%s)\n" * "Gahrtahvo: %s / Tahvo: %s / Gorahn: %s / Prorahn: %s",
                     hahr_dni_digits_str, yahr_dni_digits_str, vailee, vailee_name,
                     gahrtahvo_dni_digits_str, tahvo_dni_digits_str, gorahn_dni_digits_str, prorahn_dni_digits_str)
    )

    Tk.after(app.root, 1000, () -> update_time(app)) # Schedule next update
end

function main_julia()
    root = Tk.toplevel() # This is the function call
    app = DniClockApp(root)
    Tk.wait_until_closed(root)
end

# Run the application
main_julia()
