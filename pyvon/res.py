from Tkinter import *
import sys,os

root = Tk()

larg = root.winfo_screenwidth()
haut = root.winfo_screenheight()
largeur=int(larg/6.2)
hauteur=int(haut/60)
pos_vert=(haut/4)*3


os.system("xterm -cr blue -bg black -fg white -C -title streams -geometry %dx%d+0+%d -e python -OO /usr/local/pyvon/start.py" % (largeur,hauteur,pos_vert) )
#os.system("Eterm -font1 -n streams -b 'rgb:FF/FF/BA' -f black -g %dx%d+0+%d -e python c51.py" % (largeur,hauteur,pos_vert) )
# Eterm -b "rgb:f/f/b" -f black -g 165x13 -n streams # pour 1024x768

sys.exit
