'''
---------------------------------------------------------------
Estos son trozos de código para hacer un widget GTK
---------------------------------------------------------------
'''
#!/usr/bin/env python

import sys
try:
    import pygtk
      pygtk.require("2.0")
except:
      pass
try:
    import gtk
      import gtk.glade
except:
    sys.exit(1)
'''
-----------------------------------------------------------------
(2) Se carga el archivo 'loquesea'.glade con:
-----------------------------------------------------------------
'''
class HellowWorldGTK:
    """This is an Hello World GTK application"""

    def __init__(self):

        #Set the Glade file
        self.gladefile = "loquesea.glade"
            self.wTree = gtk.glade.XML(self.gladefile)

        #Get the Main Window, and connect the "destroy" event
        self.window = self.wTree.get_widget("MainWindow")
        if (self.window):
            self.window.connect("destroy", gtk.main_quit)
'''
-----------------------------------------------------
(3) Se hace que 'rule':
-----------------------------------------------------
'''
if __name__ == "__main__":
    hwg = HellowWorldGTK()
    gtk.main()
'''
-----------------------------------------------------
(4) lo que queda es añadir algo de funcionalidad (que haga algo), primero con el diccionario de acciones: (en este caso el boton "btnHelloWorld" y el destruir la ventana:
-----------------------------------------------------
'''
#Create our dictionay and connect it
dic = { "on_btnHelloWorld_clicked" : self.btnHelloWorld_clicked,
    "on_MainWindow_destroy" : gtk.main_quit }
self.wTree.signal_autoconnect(dic)
'''
------------------------------------------------------
(5) y que el botón haga algo (definimos lo que hace el boton, en este caso dice "¡hola mundo!" cuando se clikea):
------------------------------------------------------
'''
    def btnHelloWorld_clicked(self, widget):
        print "Hello World!"
'''
------------------------------------------------------
'''
