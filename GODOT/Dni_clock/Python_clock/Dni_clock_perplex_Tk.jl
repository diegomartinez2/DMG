import Pkg; Pkg.add("Tk")

using Dates
using Printf
using Tk

# -------------------------
# (Aquí va TODO el código de definición de constantes y funciones que ya tenías,
# incluyendo: to_digits_base_n, from_digits_base_n,
# gregorian_to_julian, julian_to_gregorian,
# cavernian_to_atrian_yahr, atrian_yahr_to_cavernian,
# cavernian_to_gregorian, gregorian_to_cavernian
# etc. (como en tu script original)
# -------------------------

# Por brevedad no repito las funciones aquí, pero en tu script completo debes incluirlas.

# ---- Interfaz gráfica corregida ----

mutable struct DniClockApp
    root::Tk.TkWindow
    gregorian_date_label::Any
    gregorian_time_label::Any
    dni_date_label::Any
    base_display::Int

    function DniClockApp(master::Tk.TkWindow)
        t = new()
        t.root = master
        Tk.title(master, "D'ni Clock Converter")

        t.base_display = 25 # Base para mostrar números D'ni

        # Labels para Gregorian time
        Tk.pack(Tk.Label(master, text="Fecha y Hora Actual (Gregoriana)", font=("Arial", 14, "bold")), pady=5)
        t.gregorian_date_label = Tk.Label(master, text="", font=("Arial", 12))
        Tk.pack(t.gregorian_date_label)
        t.gregorian_time_label = Tk.Label(master, text="", font=("Arial", 12))
        Tk.pack(t.gregorian_time_label)

        # Labels para D'ni time
        Tk.pack(Tk.Label(master, text="\nFecha y Hora D'ni", font=("Arial", 14, "bold")), pady=5)
        t.dni_date_label = Tk.Label(master, text="", font=("Arial", 12))
        Tk.pack(t.dni_date_label)

        update_time(t)

        return t
    end
end

function update_time(app::DniClockApp)
    now = Dates.now()
    day, month, year = Dates.day(now), Dates.month(now), Dates.year(now)
    hour, minute, second = Dates.hour(now), Dates.minute(now), Dates.second(now)

    hahr, vailee, yahr, gahrtahvo, tahvo, gorahn, prorahn =
        gregorian_to_cavernian(day, month, year, hour, minute, second)

    # Actualizar etiquetas Gregorianas
    Tk.configure(app.gregorian_date_label, text=@sprintf("Fecha: %d-%02d-%02d", year, month, day))
    Tk.configure(app.gregorian_time_label, text=@sprintf("Hora: %02d:%02d:%02d", hour, minute, second))

    # Formatear valores D'ni en base 25 para mostrar
    hahr_dni_digits_str = join(to_digits_base_n(hahr, app.base_display))
    yahr_dni_digits_str = join(to_digits_base_n(yahr, app.base_display))
    gahrtahvo_dni_digits_str = join(to_digits_base_n(gahrtahvo, app.base_display))
    tahvo_dni_digits_str = join(to_digits_base_n(tahvo, app.base_display))
    gorahn_dni_digits_str = join(to_digits_base_n(gorahn, app.base_display))
    prorahn_dni_digits_str = join(to_digits_base_n(prorahn, app.base_display))

    vailee_name = get(VAILEE_NAMES, vailee, "Desconocido") # Nombre Vailee

    Tk.configure(app.dni_date_label,
        text=@sprintf("Hahr: %s / Yahr: %s / Vailee: %d (%s)\n" *
            "Gahrtahvo: %s / Tahvo: %s / Gorahn: %s / Prorahn: %s",
            hahr_dni_digits_str, yahr_dni_digits_str, vailee, vailee_name,
            gahrtahvo_dni_digits_str, tahvo_dni_digits_str, gorahn_dni_digits_str, prorahn_dni_digits_str)
    )

    # Llama a update_time cada 1000 ms para refrescar la hora
    Tk.after(app.root, 1000, () -> update_time(app))
end

function main_julia()
    root = Tk.Tk() # Ventana principal, no Toplevel
    app = DniClockApp(root)
    Tk.mainloop()
end

# Ejecutar aplicación
main_julia()
