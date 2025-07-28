import Pkg; Pkg.add("Gtk")

using Gtk
using Dates
using Printf

# Aquí debes definir las funciones de conversión D'ni que ya tienes,
# pero para el ejemplo corto sólo uso la fecha actual para mostrar.

function dniclock_window()
    win = GtkWindow("D'ni Clock Converter", 400, 200)

    vbox = GtkBox(:v)
    push!(win, vbox)

    label_greg_date = GtkLabel("")
    label_greg_time = GtkLabel("")
    label_dni = GtkLabel("")

    push!(vbox, GtkLabel("Fecha y Hora Actual (Gregoriana)"))
    push!(vbox, label_greg_date)
    push!(vbox, label_greg_time)
    push!(vbox, GtkSeparator())
    push!(vbox, GtkLabel("Fecha y Hora D'ni"))
    push!(vbox, label_dni)

    function update_labels()
        now = Dates.now()
        day, month, year = Dates.day(now), Dates.month(now), Dates.year(now)
        hour, minute, second = Dates.hour(now), Dates.minute(now), Dates.second(now)

        set_gtk_property!(label_greg_date, :label,
            @sprintf("Fecha: %d-%02d-%02d", year, month, day))
        set_gtk_property!(label_greg_time, :label,
            @sprintf("Hora: %02d:%02d:%02d", hour, minute, second))

        ## Aquí llamas a tus funciones de conversión gregorian_to_cavernian y formateas igual

        # Por ejemplo, llamo a tu función (debes copiarla aquí)
        hahr, vailee, yahr, gahrtahvo, tahvo, gorahn, prorahn =
            gregorian_to_cavernian(day, month, year, hour, minute, second)

        texto_dni = @sprintf("Hahr: %d / Yahr: %d / Vailee: %d\nGahrtahvo: %d / Tahvo: %d / Gorahn: %d / Prorahn: %d",
            hahr, yahr, vailee, gahrtahvo, tahvo, gorahn, prorahn)

        set_gtk_property!(label_dni, :label, texto_dni)
    end

    # Refrescar cada segundo
    GLib.timeout_add_seconds(GLib.PRIORITY_DEFAULT_IDLE, 1) do
        update_labels()
        true
    end

    update_labels()

    showall(win)
    return win
end

function main()
    dniclock_window()
    Gtk.GtkMain()
end

main()
