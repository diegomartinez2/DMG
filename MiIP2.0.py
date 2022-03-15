#!/usr/bin/env python
# -*- coding: cp1252 -*-
#Desarrollado por RadicalEd
import re, urllib2
#obtenemos el fuente de la pagina
print "\nObteniendo IP."

s = urllib2.urlopen('http://www.showmyip.com').read()
#Con expresiones regulares se obtiene todo lo que venga despues del 'displaycopy'
m = re.search('(?<=displaycopy).*', s) #Esta línea es un pedazo de JS que imprime en la página la IP
ip = m.group(0)
#Reemplazamos por nada los datos que no sirven
ip = ip.replace('("', '')
ip = ip.replace('");', '')
print ip
print "\nOk ^-^."

print "\nEscribiendo."

text_file = open("IPs.txt","w")

text_file.write(ip)
text_file = open("IPs.txt","r")

old_ip=text_file.read()
print "\nOk ^-^."

if old_ip != ip:
	print "\nEnviando."
	#envia correo
	#from email.MIMEText import MIMEText
	#msg = MIMEText(ip)
	#msg['Subject'] = 'Esto es una prueba para el envio del IP'
	#msg['From'] = "La_Python-isa"
	#msg['Reply-to'] = "Nadie_aqui"
	#msg['To'] = "root"
	
	#import getpass
	#import smtplib
	#sender = smtplib.SMTP('smtp.gmail.com')
	#sender.ehlo()
	#sender.starttls()
	#usuario_gmail = 'usuario_gmail'
	#sender.login(usuario_gmail, getpass.getpass())
	#origen = usuario_gmail + '@gmail.com'
	#destino = 'ssmm_sofi@yahoo.es'
	#sender.sendmail(origen, destino, msg.as_string())
	#sender.close()
	
	#otra forma
	# Import smtplib for the actual sending function

	import smtplib

	# Import the email modules we'll need
	from email.mime.text import MIMEText

	# Open a plain text file for reading.  For this example, assume that
	# the text file contains only ASCII characters.
	#fp = open(textfile, 'rb')
	# Create a text/plain message
	#msg = MIMEText(fp.read())
	msg = MIMEText(ip)
	#fp.close()

	# me == the sender's email address
	me = "usuario_gmail@gmail.com"
	# you == the recipient's email address
	you = "usuario_yahoo@yahoo.es"
	msg['Subject'] = 'La IP: %s' % ip
	msg['From'] = me
	msg['To'] = you
	
	# Send the message via our own SMTP server, but don't include the
	# envelope header.
	s = smtplib.SMTP()
	#s.connect("www.wellho.net") #cambiarlo por el servidor de correo adecuado
	s.connect("smtp.gmail.com")
	s.sendmail(me, [you], msg.as_string())
	s.quit()
	
	
	text_file.write(old_ip)
	print "\nOk ^-^."
else:
	print "La IP no ha cambiado."

text_file.close()
