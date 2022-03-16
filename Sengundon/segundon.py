#!/usr/bin/python
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
	print("GTK Not Availible")
	sys.exit(1)

class adder:

	result = 0
	
	def __init__( self, number1, number2,number3,number4,number5,number6 ):
		self.result = int( number1 )*32140800 + int( number2 )*2678400+int( number3 )*86400+int( number4 )*3600+int( number5 )*60+int( number6 )
		
	def giveResult( self ):
		return str(self.result)
		
class es_genial:

	wTree = None

	def __init__( self ):
		self.wTree = gtk.glade.XML( "main.glade" )
		
		dic = { 
			"on_buttonQuit_clicked" : self.quit,
			"on_calcular_clicked" : self.add,
			"on_windowMain_destroy" : self.quit,
		}
		
		self.wTree.signal_autoconnect( dic )
		
		gtk.main()

	def add(self, widget):
		try:
			thistime = adder( self.wTree.get_widget("in_anos").get_text(), self.wTree.get_widget("in_meses").get_text(), self.wTree.get_widget("in_dias").get_text(), self.wTree.get_widget("in_horas").get_text(), self.wTree.get_widget("in_minutos").get_text(), self.wTree.get_widget("in_segundos").get_text() )
		except ValueError:
			self.wTree.get_widget("hboxWarning").show()
			self.wTree.get_widget("out_segundos").set_text("ERROR")
			return 0
		self.wTree.get_widget("hboxWarning").hide()
		self.wTree.get_widget("out_segundos").set_text(thistime.giveResult())
	
	def quit(self, widget):
		sys.exit(0)
		
letsdothis = es_genial()
