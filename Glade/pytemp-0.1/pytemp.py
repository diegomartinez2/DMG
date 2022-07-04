#! /usr/bin/env python
# -*- coding: UTF-8 -*-

# pytemp - A temperature converter written in python using PyGTK
# Copyright (C) 2009 Daniel Fuentes Barría <dbfuentes@gmail.com>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

# Importamos los módulos necesarios
import sys
try:
    import pygtk
    pygtk.require('2.0') # Intenta usar la versión2
except:
    # Algunas distribuciones vienen con GTK2, pero no con pyGTK (o pyGTKv2)
    pass

try:
    import gtk
    import gtk.glade
except:
    print ("You need to install pyGTK or GTKv2 or set your PYTHONPATH correctly")
    sys.exit(1)

# Definimos varias funciones que convierten la temperatura.

def fahrenheit2celsius(temp):
    "Convert Fahrenheit to celsius"
    celsius = (temp - 32) * 5.0 / 9.0
    return celsius

def celsius2fahrenheit(temp):
    "Convert Celsius to Fahrenheit"
    fahrenheit = (9.0 / 5.0) * temp + 32
    return fahrenheit

def kelvin2celsius(temp):
    "Convert kelvin to Celsius"
    celsius = temp - 273.15
    return celsius

def celsius2kelvin(temp):
    "Convert Celsius to kelvin"
    kelvin = temp + 273.15
    return kelvin

# Función que obtiene el texto de la opción seleccionada en un ComboBox
# Se usa para obtener que unidad de temperatura seleccionada de la lista
def valor_combobox(combobox):
    model = combobox.get_model()
    activo = combobox.get_active()
    if activo <0:
        return None
    return model[activo][0]

# Creamos una clase que almacena la información del programa (después se usara)
class Info:
    "Store the program information"
    name = "pyTemp"
    version = "0.1"
    copyright = "Copyright © 2009 Daniel Fuentes B."
    authors = ["Daniel Fuentes Barría <dbfuentes@gmail.com>"]
    website = "http://pythonmania.wordpress.com/tutoriales/"
    description = "A temperature converter written in python using PyGTK"
    license = "This program is free software; you can redistribute it and/or \
modify it under the terms of the GNU General Public License as published by \
the Free Software Foundation; either version 2 of the License, or (at your \
option) any later version. \n\nThis program is distributed in the hope that \
it will be useful, but WITHOUT ANY WARRANTY; without even the implied \
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. \
See the GNU General Public License for more details. \n\nYou should have \
received a copy of the GNU General Public License along with this program; \
if not, write to the Free Software Foundation, Inc., 51 Franklin Street, \
Fifth Floor, Boston, MA 02110-1301, USA."


# Interfaz gráfica (gtk-glade), Clase para el Loop principal (de la GUI)
class MainGui:
    "GTK/Glade User interface. This is a pyGTK window"
    def __init__(self):
        # Le indicamos al programa que archivo XML de glade usar
        self.widgets = gtk.glade.XML("pytemp.glade")

        # se definen las signals
        signals = { "on_button1_clicked" : self.on_button1_clicked,
                    "on_about1_activate" : self.on_about1_activate,
                    "gtk_main_quit" : gtk.main_quit }

        # y se autoconectan las signals.
        self.widgets.signal_autoconnect(signals)

        # Del archivo glade obtenemos los widgets a usar
        self.entry1 = self.widgets.get_widget("entry1")
        self.textview1 = self.widgets.get_widget("textview1")
        self.combobox1 = self.widgets.get_widget("combobox1")
        self.combobox2 = self.widgets.get_widget("combobox2")

        # Para el ComboBox1 se fija por defecto la primera opción de la lista
        self.combobox1.set_active(0)
        # y en el ComboBox2 se fija por defecto la segunda opción de la lista
        self.combobox2.set_active(1)

    # A continuación se definen/crean las ventanas especiales (about, dialogo,
    # etc) y las acciones a realizar en ellas.

    # Ventana genérica de error (se le pasan los mensajes de error y los
    # muestra en una ventana de dialogo)
    def error(self, message):
        "Display the error dialog "
        dialog_error = gtk.MessageDialog(parent=None, flags=0, buttons=gtk.BUTTONS_OK)
        dialog_error.set_title("Error")
        label = gtk.Label(message)
        dialog_error.vbox.pack_start(label, True, True, 0)
        # Con show_all() mostramos el contenido del cuadro de dialogo (en este
        # caso solo tiene la etiqueta) si no se hace el dialogo aparece vacío
        dialog_error.show_all()
        # El run y destroy hace que la ventana se cierre al apretar el botón
        dialog_error.run()
        dialog_error.destroy()

    # Ventana About (conocida como Acerca de).
    def about_info(self, data=None):
        "Display the About dialog "
        about = gtk.AboutDialog()
        about.set_name(Info.name)
        about.set_version(Info.version)
        about.set_comments(Info.description)
        about.set_copyright(Info.copyright)

        def openHomePage(widget,url,url2): # Para abrir el sitio
            import webbrowser
            webbrowser.open_new(url)

        gtk.about_dialog_set_url_hook(openHomePage,Info.website)
        about.set_website(Info.website)
        about.set_authors(Info.authors)
        about.set_license(Info.license)
        about.set_wrap_license(True) # Adapta el texto a la ventana
        about.run()
        about.destroy()


    # Ahora declaramos las acciones a realizar (por menús, botones, etc.):

    # Definimos la ventana about (help > About)
    def on_about1_activate(self, widget):
        "Open the About windows"
        self.about_info()

    # Definimos las acciones a realizar al apretar el botón de convertir
    def on_button1_clicked(self, widget):
        "Convert button"
        # Se crea un buffer en donde se guardaran los resultados
        text_buffer = gtk.TextBuffer()
        # Se obtiene el valor para convertir desde la entrada
        valor = self.entry1.get_text()
        # Obtiene la opción escogida en los 2 ComboBoxs
        selec1 = valor_combobox(self.combobox1)
        selec2 = valor_combobox(self.combobox2)
        try:
            # Intenta transformar el valor ingresado en un numero. En caso
            # de fallar (por ejemplo falla si lo ingresado son letras) se
            # lanza la excepción, si es exitoso se continua con la conversión
            temp_ini = float(valor)

            # Inicia la conversión adecuada dependiendo de la opción
            # escogida en los 2 ComboBoxs (selec1 y selec2)
            if selec1 == "Celsius" and selec2 == "Fahrenheit":
                text_buffer.set_text(str(celsius2fahrenheit(temp_ini)))
            elif selec1 == "Celsius" and selec2 == "Kelvin":
                text_buffer.set_text(str(celsius2kelvin(temp_ini)))
            elif selec1 == "Fahrenheit" and selec2 == "Celsius":
                text_buffer.set_text(str(fahrenheit2celsius(temp_ini)))
            elif selec1 == "Fahrenheit" and selec2 == "Kelvin":
                # Pasamos primero de F a Celsius, luego de Celsius a Kelvin
                conversion1 = fahrenheit2celsius(temp_ini)
                text_buffer.set_text(str(celsius2kelvin(conversion1)))
            elif selec1 == "Kelvin" and selec2 == "Celsius":
                text_buffer.set_text(str(kelvin2celsius(temp_ini)))
            elif selec1 == "Kelvin" and selec2 == "Fahrenheit":
                 # Pasamos primero de Kelvin a Celsius, luego de Celsius a F
                 conversion1 = kelvin2celsius(temp_ini)
                 text_buffer.set_text(str(celsius2fahrenheit(conversion1)))
            else:
                # Se produce cuando las dos selecciones son iguales
                self.error("The initial and target units are the same")

            # Luego se fija (muestra) el buffer (que contiene la temperatura
            # convertida) en textview1 (el cuadro la lado del botón)
            self.textview1.set_buffer(text_buffer)

        except:
            if (len(valor) == 0):
                # Se produce si no se ingresa nada en la entry1
                self.error("Please enter a value")
            else:
                # Se produce si no se ingresa un numero ( por ejemplo se
                # ingresan letras o símbolos)
                self.error("The value entered is not valid.\nEnter a \
valid number")


if __name__== "__main__":
    MainGui()
    gtk.main()
