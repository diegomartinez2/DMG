#!/usr/bin/env python

import sys, os, string,re, signal, popen2 # à verifier
from os import spawnv, P_NOWAIT
from string import split
from Tkinter import *
from tkFont import families
import stat





thisPyFile = os.path.join(os.getcwd(), sys.argv[0])
thisDir = os.path.normpath(os.path.dirname(thisPyFile))
langDir = os.path.join(thisDir, 'Lang_kw')
PmwDir = os.path.join(thisDir, 'Pmw')
IconsDir = os.path.join(thisDir, 'Icons')
InsertDir = os.path.join(thisDir, 'Insert')
sys.path.append(langDir)
sys.path.append(PmwDir)



import Pmw
from kws35 import *



try :
	homedir = os.path.expanduser('~')
	fd = open(homedir+"/.pyvonrc", 'r')
	for line in fd.readlines():
		if line.startswith("Language :"):
			line=line[11:]
			if line.endswith('\n'):
				line=line[:-1]
			Lang='menu_'+line
			break
except:
	Lang="menu_english"


evalstr = 'from ' + Lang + ' import *'
exec evalstr


###############################################################
#--                 Loads pyvonrc in home                   --#

#--          Checks if the path for the executable          --#
#--        found in pyvonrc exist on the hard drive         --#
###############################################################



def loadini():
	try:
		homedir = os.path.expanduser('~')
		fd=open(homedir+"/.pyvonrc")
		for line in fd.readlines():
			for val in ini_list:
				if line.find(val[0])==0:
					val[1]=line[(val[2]):]
					if val[1].endswith("\n"):
						val[1]=val[1][:-1]
					if val[1].endswith("\r"):
						val[1]=val[1][:-1]
		fd.close()
	finally:
		if ini_list[3][1]=="":ini_list[3][1]=" "
		if ini_list[5][1]=="":ini_list[5][1]="arial"
		if ini_list[6][1]=="":ini_list[6][1]="10"
		if ini_list[7][1]=="":ini_list[7][1]="Yes"
		if ini_list[8][1]=="":ini_list[8][1]="Yes"
		if ini_list[9][1]=="":ini_list[9][1]="english"
		if ini_list[10][1]=="":ini_list[10][1]="#840084"
		if ini_list[11][1]=="":ini_list[11][1]="#5D00A4"
		if ini_list[12][1]=="":ini_list[12][1]="#008284"
		if ini_list[13][1]=="":ini_list[13][1]="#008200"
		if ini_list[14][1]=="":ini_list[14][1]="#FF0000"
		if ini_list[15][1]=="":ini_list[15][1]="#FFFF11"
		if ini_list[16][1]=="":ini_list[16][1]="#895511"
		if ini_list[17][1]=="":ini_list[17][1]="#753B93"
		if ini_list[18][1]=="":ini_list[18][1]="#0000FF"
		if ini_list[19][1]=="":ini_list[19][1]="#FFFFFF"
		if ini_list[20][1]=="":ini_list[20][1]="#000000"
		if ini_list[21][1]=="":ini_list[21][1]="40"
		if ini_list[22][1]=="":ini_list[22][1]="Yes"
		if ini_list[23][1]=="":ini_list[23][1]="Yes"




plat=sys.platform
if plat=="win32":police = ('helvetica',10)
else:police = '-Adobe-Helvetica-Medium-R-Normal-*-120-*-*-*-*-*-*'


#root=Pmw.initialise(root = None, size = None, fontScheme = 'pmw2', useTkOptionDb = 0, noBltBusy = 1, disableKeyboardWhileBusy = 1)
Pmw.clearbusycursor()
root = Tk()
nouveau = PhotoImage(file=(os.path.join(IconsDir,"new.gif")))
ouvrir = PhotoImage(file=(os.path.join(IconsDir,"open.gif")))
closee = PhotoImage(file=(os.path.join(IconsDir,"close.gif")))
sauver = PhotoImage(file=(os.path.join(IconsDir,"save.gif")))
lancer = PhotoImage(file=(os.path.join(IconsDir,"run.gif")))
stop = PhotoImage(file=(os.path.join(IconsDir,"stop.gif")))
pausee = PhotoImage(file=(os.path.join(IconsDir,"pause.gif")))
resumee = PhotoImage(file=(os.path.join(IconsDir,"resume.gif")))




#build the keywords pattern string
keypat=[r"(?P<%s>" % name + string.join(keywords[name], "|")+")"
        for name in keywords.keys()]
keypat=r"\b(?:"+string.join(keypat, "|")+r")\b"
rkeypat=re.compile(keypat)


#build the directives pattern string
dirpat=[r"#\s*(?P<%s>" % name + string.join(directives[name], "|")+")"
         for name in directives.keys()]
dirpat=r"(?:"+string.join(dirpat, "|")+r")\b"

#patterns not coverd by the above
mlcompat=r'(?P<mlcom>/\*)'      #match start of multiline comment
slcompat=r'(?P<slcom>//[^\n]*)' #single line comment
strpat=r'(?P<str>(?<!\\)")'     #matches quot (") but not escaped quot (\")
squigpat=r'(?P<squig>}|{)'      # { }
numpat=r'(?P<num>[.\d]+\d*([eE][+-]\d+)?)'          #numbers
mathpat=r'(?P<math>\^|\*|/|\+|-|=|!|>|<|&|\||\?)'   #math-symbols
identpat=r"(?P<ident>[A-Za-z_][A-Za-z_0-9]*)"       #identifier
end_of_linepat=r'(?P<end_of_line>[\n])'


#join it all to one string and compile to regular expression
patlist=[mlcompat,slcompat,strpat,squigpat,dirpat,keypat,identpat,numpat,mathpat,end_of_linepat]

pat=string.join(patlist,"|")
global pattern
pattern=re.compile(pat)         #match all keywords
#some patterns for separate matching steps
global mlcom
global mlcompattern
mlcompattern=re.compile(r'(?P<mlcom>\*/|/\*)')       #multi-line comment
strpattern=re.compile(r'(?P<str>(?<!\\)")')          #string



global filename
filename='untitled.pov'
global MLC
MLC = "\*/|/\*"
global scene_contents
global ident_flag
global tabu
tabu='[ \t]'
global fin_de_ligne
fin_de_ligne='[ \n]'
global sel
sel=""
global General_Buffer
General_Buffer=[]
global Buffer
Buffer=[0,0,0,0]
scene_contents=""
scene_color=[]
global scene_idx
scene_idx=0
global Files_list
Files_list=['untitled.pov']
global includedir
includedir=[]
global fin
global sts_pause
sts_pause=0
global sts_rendu
sts_rendu=0
global re_nice
re_nice=0
global preset_resol
preset_resol=""
global SC,EC,SR,ER
global start_column,end_column,start_row,end_row
start_column=0.0
end_column=1.0
start_row=0.0
end_row=1.0
global forbidden_keys
forbidden_keys=[65364,65361,65363,65362,65288,65507,
		65508,65379,65433,65430,65432,65431,
		65549,65407,65434,65435,65366,65365,
		65505,65506,65407,65513,65514,65307,
		65535,65387,65299,65377]




global Oldfiles			# list to store previously opened files
Oldfiles=[]


color_tags


def save_b_render(s):
        if s==1:
            ini_list[7][1]="Yes"
        if s==0:
            ini_list[7][1]="No"


def show_hlp_b(s):
	if s==1:
		ini_list[22][1]="Yes"
	if s==0:
		ini_list[22][1]="No"


def auto_exkw(s):
	if s==1:
		ini_list[23][1]="Yes"
	if s==0:
		ini_list[23][1]="No"




def color_onload(i):
	if i==1:
		ini_list[8][1]="Yes"
	if i==0:
		ini_list[8][1]="No"



########################################################	
## ---------------------- Expand word ------------------
########################################################

def change_color():
	global rv,gv,bv,Display_color_box,etiquette_size,font_editor_listbox,Change_color_TL
	global SRED,SGREEN,SBLUE,drop_down_color_class,TEXT_FONT
	rv=0
	gv=0
	bv=0
	Change_color_TL=Toplevel(width=80)
	Change_color_TL.resizable(width=0,height=0)
	Change_color_TL.geometry(newGeometry="+250+150")
	Change_color_TL.title( Trans_Font_Options )
	F_LEFT=Frame(Change_color_TL)
	F_RIGHT=Frame(Change_color_TL)


##	- SCALES -       ##
	SRED=Scale(F_LEFT, orient='horizontal', borderwidth=1, from_=0, to=1, resolution=0.01, sliderlength=20,width=10,length=150,command=RED_COMM)
	SRED.grid(column=1,row=0,sticky=NW)
	SGREEN=Scale(F_LEFT, orient='horizontal', borderwidth=1, from_=0, to=1, resolution=0.01, sliderlength=20,width=10,length=150,command=GREEN_COMM)
	SGREEN.grid(column=1,row=1,sticky=NW)
	SBLUE=Scale(F_LEFT, orient='horizontal', borderwidth=1, from_=0, to=1, resolution=0.01, sliderlength=20,width=10,length=150,command=BLUE_COMM)
	SBLUE.grid(column=1,row=2,sticky=NW)

##	- LABELS -       ##
	Label(F_LEFT,text= Trans_Red ).grid(column=0,row=0,sticky=NW)
	Label(F_LEFT,text= Trans_Green ).grid(column=0,row=1,sticky=NW)
	Label(F_LEFT,text= Trans_Blue ).grid(column=0,row=2,sticky=NW)


##	- COLORED BOX -       ##
	Display_color_box=Canvas(F_LEFT,height=50,width=130)
	Display_color_box.grid(row=3,column=0,columnspan=2)
	F_Tab=Frame(F_LEFT)
	Label(F_Tab,text= Trans_Indentation ).grid(row=0,column=0,sticky=W)
	global Tab_Entry
	Tab_Entry=Entry(F_Tab,width=5)
	Tab_Entry.insert(INSERT,ini_list[21][1])
	Tab_Entry.grid(row=0,column=1)
	Label(F_Tab,text=" pixels").grid(row=0,column=2,sticky=W)
	F_Tab.grid(row=4,column=0,columnspan=2,sticky=W,pady=5)
	F_LEFT.grid(row=0,column=0,sticky=N)

	Label(F_RIGHT,text= Trans_Select_Word_Class).grid(row=0,column=0,sticky=N)
	KWDS=["keywords","bold keywords","comments","numbers & math","strings","patterns",
			"directives","braces","identifiers","editor background","highlighting"]
	drop_down_color_class=Pmw.ComboBox(F_RIGHT,scrolledlist_items=KWDS,
			selectioncommand=drop_down_choose_type)
	drop_down_color_class.selectitem(KWDS[0],1)
	drop_down_choose_type("keywords")
	drop_down_color_class.grid(row=1,column=0,sticky=N)

	Button(Change_color_TL,text= Trans_Apply ,command=Apply_color).grid(row=6,column=0,sticky=E+W)
	Button(Change_color_TL,text= Trans_Exit ,command=Exit_change_colors).grid(row=6,column=1,sticky=E+W)





###############################################################
#--                         Font menu                       --#
###############################################################


	global fonts
	fonts=list(families(root))
	fonts.sort()
	sizes=[4,5,6,7,8,9,10,11,12,13,14,15,16,18,20,22,24,26,30]


## etiquette Font size :

	Label(F_RIGHT,text='').grid(row=2,column=0)
	Label(F_RIGHT,text= Trans_Size ).grid(row=3,column=0)
	etiquette_size=Pmw.ComboBox(F_RIGHT,scrolledlist_items=sizes,listheight=100,selectioncommand=show_size)
	etiquette_size.selectitem(ini_list[6][1],1)
	etiquette_size.grid(row=4,column=0,padx=2,pady=2)

## bouton pour valider

	Label(F_RIGHT,text='').grid(row=5,column=0)
	Label(F_RIGHT,text= Trans_Font ,anchor=W).grid(row=6,column=0)
	font_editor_listbox =  Pmw.ComboBox(F_RIGHT, scrolledlist_items=fonts,selectioncommand=show_font)

	try:
		font_editor_listbox.selectitem(ini_list[5][1],1)
	except:
		font_editor_listbox.selectitem(fonts[0],1)
	font_editor_listbox.grid(row=7,column=0,padx=2,pady=2)
	Label(F_RIGHT,text='').grid(row=8,column=0)
	TEXT_FONT=Label(F_RIGHT,text='{ [ Test Font " >',background='white')
	TEXT_FONT.grid(row=9,column=0,padx=2,pady=2)
	F_RIGHT.grid(row=0,column=1,sticky=N)


def show_font(ft):
	ini_list[6][1] = etiquette_size.get()
	TEXT_FONT.configure(font=(ft, ini_list[6][1]))
	font_editor_listbox.focus_force()

def show_size(sz):
	ini_list[5][1] = font_editor_listbox.get()
	TEXT_FONT.configure(font=(ini_list[5][1],sz))
	etiquette_size.focus_force()
	

def click_event(event=None):
	text.tag_config('keyword', font=(ini_list[5][1], ini_list[6][1], 'normal' ),foreground=ini_list[10][1])
	text.tag_config('keywordbold', font=(ini_list[5][1], ini_list[6][1], 'bold' ),foreground=ini_list[11][1])
	text.tag_config('num', font=(ini_list[5][1], ini_list[6][1], 'normal' ),foreground=ini_list[12][1])
	text.tag_config('slcomment', font=(ini_list[5][1], ini_list[6][1], 'italic' ),foreground=ini_list[13][1])
	text.tag_config('mlcomment', font=(ini_list[5][1], ini_list[6][1], 'italic' ),foreground=ini_list[13][1])
	text.tag_config('math', font=(ini_list[5][1], ini_list[6][1], 'normal' ),foreground=ini_list[12][1])
	text.tag_config('squig', font=(ini_list[5][1], ini_list[6][1], 'normal' ),foreground=ini_list[18][1])
	text.tag_config('dirc', font=(ini_list[5][1], ini_list[6][1], 'normal' ),foreground=ini_list[17][1])
	text.tag_config('str', font=(ini_list[5][1], ini_list[6][1], 'normal' ),foreground=ini_list[14][1])
	text.tag_config('pattern', font=(ini_list[5][1], ini_list[6][1], 'normal' ),foreground=ini_list[16][1])
	text.tag_config('highl',background=ini_list[15][1])
	text.configure(font=(ini_list[5][1], ini_list[6][1]),background=ini_list[19][1],
			foreground=ini_list[20][1],insertbackground = ini_list[20][1])
	ini_list[21][1]=Tab_Entry.get()
	text.configure(tabs=ini_list[21][1])
	Change_color_TL.focus_force()


def Exit_change_colors():
	ini_list[21][1]=Tab_Entry.get()
	text.configure(tabs=ini_list[21][1])
	Change_color_TL.destroy()


def Apply_color():
	rv=int(SRED.get()*255)
	gv=int(SGREEN.get()*255)
	bv=int(SBLUE.get()*255)
	ln=drop_down_color_class.get()
	count=0
	for val in ini_list:
		if val[0][:-3]==ln:
			ini_list[count][1]= "#%02x%02x%02x" % (rv, gv, bv)
		count=count+1

	ini_list[5][1] = font_editor_listbox.get()
	ini_list[6][1]=etiquette_size.get()

	click_event()


def drop_down_choose_type(type_chosen):
	count=0
	for val in ini_list:
		if val[0][:-3]==type_chosen:
			rouge = ini_list[count][1][1:3]
			vert  = ini_list[count][1][3:5]
			bleu  = ini_list[count][1][5:7]
		count=count+1
	val_rouge = float(int("0x"+rouge, 16))/255
	val_vert  = float(int("0x"+vert, 16))/255
	val_bleu  = float(int("0x"+bleu, 16))/255
	SRED.set(val_rouge)
	SGREEN.set(val_vert)
	SBLUE.set(val_bleu)


def RED_COMM(rv):
	rv=int(SRED.get()*255)
	gv=int(SGREEN.get()*255)
	bv=int(SBLUE.get()*255)
	canvas_color="#%02x%02x%02x" % (rv, gv, bv)
	Display_color_box.create_rectangle(2,2,150,50,fill=canvas_color)

def GREEN_COMM(gv):
	rv=int(SRED.get()*255)
	gv=int(SGREEN.get()*255)
	bv=int(SBLUE.get()*255)
	canvas_color="#%02x%02x%02x" % (rv, gv, bv)
	Display_color_box.create_rectangle(2,2,150,50,fill=canvas_color)

def BLUE_COMM(bv):
	rv=int(SRED.get()*255)
	gv=int(SGREEN.get()*255)
	bv=int(SBLUE.get()*255)
	canvas_color="#%02x%02x%02x" % (rv, gv, bv)
	Display_color_box.create_rectangle(2,2,150,50,fill=canvas_color)



########################################################
## -------------------- Auto indent --------------------
########################################################


def auto_indent(event=None):
	try:
		a0=str(string.join((Ligne,str(0)),'.'))
		a1=str(string.join((Ligne,END),'.'))
		line_cont=text.get(a0,a1)
		b=line_cont.count('\t')
		if b != 0 :
			c=b*'\t'
			text.insert(INSERT,c)
	except:
		pass


########################################################
## ---------------------Delete line --------------------
########################################################

def delete_line(event=None):
	Ligne, column = text.index(INSERT).split('.', 1)
	a0=str(string.join((Ligne,str(0)),'.'))
	a1=str(string.join((Ligne,'end'),'.'))
	if text.get(a1)=='\n':
		a1=str(string.join((Ligne,'end+1chars'),'.'))
	text.delete(a0,a1)
	suppr_mlcomment()


########################################################
## ---------------------Delete word --------------------
########################################################

def delete_word(event=None):
	Ligne, column = text.index(INSERT).split('.', 1)
	col_0=str(string.join((str(Ligne),str(0)),'.'))
	espace_avant=text.search(tabu,INSERT,backwards='yes',regexp='yes',stopindex=col_0)
	if espace_avant=="":
		text.delete(col_0,INSERT)
	else:
		text.delete(espace_avant,INSERT)
	suppr_mlcomment()


########################################################
## ---------------------- Expand word ------------------
########################################################


def auto_expand_word():
	global kw_loc_start
	try:
		pos_actu=str(string.join((Ligne, column),'.')) # 1.7
		espace_avant=text.search(tabu,pos_actu,backwards='yes',regexp='yes',stopindex=str(string.join((Ligne,str(0)),'.')))
		if text.get(pos_actu) == '{' or text.get(pos_actu) == '}':match_brace()
		else:text.tag_remove('accolade',0.1,END)


		if espace_avant != "" : # si il y a un espace avant dÃ©but de ligne
			kw_loc_start=str(string.join((Ligne, str(int(string.split(espace_avant,'.')[1])+1)),'.'))

			to_find=text.get(kw_loc_start,INSERT)


		elif espace_avant == "" : # si il n'y a pas d'espace avant dÃ©but de ligne -> col=0
			kw_loc_start=str(string.join((Ligne,str(0)),'.')) # ('4', '0')
			to_find=text.get(kw_loc_start,INSERT) # 1.7->1.11

		m=0
		lst=[]
		for i in kws:    # teste toutes les occurences possibles
			if i[0:len(to_find)]==to_find:
				m=i
				lst.append(m)

		global longueur_liste
		longueur_liste=len(lst)
		if longueur_liste==1:    #  si il n'y a qu'une occurence
			text.delete(kw_loc_start,INSERT)
			text.insert(INSERT,m)     #  -> insere le mot trouvÃ©
			if m in keywords['objbold']:
				text.insert(INSERT," {")
			kw_loc_start_x=int(string.split(kw_loc_start,'.')[1])
			kw_loc_end_x=kw_loc_start_x+len(m)
			kw_loc_end=string.join((str(int(string.split(kw_loc_start,'.')[0])),INSERT),'.')
			to_find=text.get(kw_loc_start,INSERT)
			p=pattern.search(to_find)
			color(p,kw_loc_start,INSERT)
	except:
		pass


def expand_word(event=None):
	global kw_loc_start
#	try:
	pos_actu=str(string.join((Ligne, column),'.')) # 1.7
	espace_avant=text.search(tabu,pos_actu,backwards='yes',regexp='yes',stopindex=str(string.join((Ligne,str(0)),'.')))
	if text.get(pos_actu) == '{' or text.get(pos_actu) == '}':match_brace()
	else:text.tag_remove('accolade',0.1,END)


	if espace_avant != "" : # si il y a un espace avant début de ligne
		kw_loc_start=str(string.join((Ligne, str(int(string.split(espace_avant,'.')[1])+1)),'.'))

		to_find=text.get(kw_loc_start,INSERT)


	elif espace_avant == "" : # si il n'y a pas d'espace avant début de ligne -> col=0
		kw_loc_start=str(string.join((Ligne,str(0)),'.')) # ('4', '0')
		to_find=text.get(kw_loc_start,INSERT) # 1.7->1.11


	m=0
	lst=[]
	for i in kws:    # teste toutes les occurences possibles
		if i[0:len(to_find)]==to_find:
			m=i
			lst.append(m)

	global longueur_liste
	longueur_liste=len(lst)

	if longueur_liste==1:
#		print kw_loc_start
#		print len(to_find)
		text.delete(kw_loc_start,'insert')
		text.insert('insert',lst[0])

	if longueur_liste>1:    # si il y a plus d'une occurence
		global Top_level_kw_sel
		Top_level_kw_sel=Toplevel()
		global kw_listbox
		Label(Top_level_kw_sel,text= Trans_Phrase, bg="white",justify=LEFT,width=30).pack()
		kw_listbox = Pmw.ScrolledListBox(Top_level_kw_sel, items=lst)
		Top_level_kw_sel.geometry(newGeometry="+250+150")
		Top_level_kw_sel.title( Trans_Word_Expansion )
		kw_listbox.pack(expand=1, fill="both")
		kw_listbox.bind("<Return>",sel_lst)
		kw_listbox.bind("<Escape>",escape)
		kw_listbox.bind("<KeyPress-Down>",sel_down)
		kw_listbox.bind("<KeyPress-Up>",sel_up)
		global selected_it
		selected_it=0
		kw_listbox.selection_set(0)
		kw_listbox.focus_set()

	if longueur_liste==0:  # si il n'y a aucune occurence
		text.configure(bg=ini_list[15][1])
		root.after(100,r_white)
#	except:
#		pass



def escape(cur):
	Top_level_kw_sel.destroy()



def r_white():
	text.configure(bg=ini_list[19][1])



def sel_down(cur):
	global selected_it
	kw_listbox.selection_clear()
	if longueur_liste>selected_it+1:
		selected_it=selected_it+1
	kw_listbox.selection_set(selected_it)
	kw_listbox.see(selected_it)



def sel_up(cur):
	global selected_it
	kw_listbox.selection_clear()
	if selected_it>0:
		selected_it=selected_it-1
	kw_listbox.selection_set(selected_it)
	kw_listbox.see(selected_it)



def sel_lst(cur):
	kw_sel=kw_listbox.getcurselection()
	text.delete(kw_loc_start,'insert')
	text.insert('insert',kw_sel)
	if kw_sel[0] in keywords['objbold']:
		text.insert(INSERT," {")
			
	Top_level_kw_sel.destroy()
	kw_loc_start_x=int(string.split(kw_loc_start,'.')[1])
	kw_loc_end_x=kw_loc_start_x+len(kw_sel[0])
	kw_loc_end=string.join((str(int(string.split(kw_loc_start,'.')[0])),str(kw_loc_end_x)),'.')
	to_find=text.get(kw_loc_start,kw_loc_end)
	p=pattern.search(to_find)
	color(p,kw_loc_start,kw_loc_end)



########################################################
## ---------------------- Fontify ----------------------
########################################################


def couleur(event=None):
	for i in color_tags:
		text.tag_remove (i,'%d.0' % int(Ligne),'%d.end' % int(Ligne))   # remove all the tags in the text
		

	for tag in taglist:							# set the tags in the text
		if tag[2]=="obj":
			text.tag_add ('keyword','%d.%d' % (int(Ligne),int(tag[0])),'%d.%d' % (int(Ligne),int(tag[1])))
		elif tag[2]=="objbold":
			text.tag_add ('keywordbold','%d.%d' % (int(Ligne),int(tag[0])),'%d.%d' % (int(Ligne),int(tag[1])))
		elif tag[2]=="pattern":
			text.tag_add ('pattern','%d.%d' % (int(Ligne),int(tag[0])),'%d.%d' % (int(Ligne),int(tag[1])))
		elif tag[2]=="num":
			text.tag_add ('num','%d.%d' % (int(Ligne),int(tag[0])),'%d.%d' % (int(Ligne),int(tag[1])))
		elif tag[2]=="str":
			text.tag_add ('str','%d.%d' % (int(Ligne),int(tag[0])),'%d.%d' % (int(Ligne),int(tag[1])))
		elif tag[2]=="squig":
			text.tag_add ('squig','%d.%d' % (int(Ligne),int(tag[0])),'%d.%d' % (int(Ligne),int(tag[1])))
		elif tag[2]=="math":
			text.tag_add ('math','%d.%d' % (int(Ligne),int(tag[0])),'%d.%d' % (int(Ligne),int(tag[1])))
		elif tag[2]=="dirc":
			text.tag_add ('dirc','%d.%d' % (int(Ligne),int(tag[0])),'%d.%d' % (int(Ligne),int(tag[1])))
		elif tag[2]=="slcom":
			text.tag_add ('slcomment','%d.%d' % (int(Ligne),int(tag[0])),'%d.%d' % (int(Ligne),int(tag[1])))



def fontify():
	if text.get('insert') == '{' or text.get('insert') == '}':match_brace()   # searches for braces
	else:text.tag_remove('parenthese',0.1,END)

	if not 'mlcomment' in text.tag_names('insert'):
		line_contents=text.get('%d.0' % int(Ligne),'%d.end' % int(Ligne))
		tag_it(line_contents)
		couleur()


	# searches for multiline comments
	st=text.search( '/*' ,'%d.0' % int(Ligne),stopindex='%d.end' % int(Ligne))
	if st:ml_comment_start(st)
	en=text.search( '*/' ,'%d.0' % int(Ligne),stopindex='%d.end' % int(Ligne))
	if en:ml_comment_end(en)



def ml_comment_start(deb):
	end_ml_com_1 = text.search( "*/" , 'insert' , backwards = 'yes' ,stopindex= '1.0')
	end_ml_com_2 = text.search( "*/" , 'insert' , stopindex= 'end')
	if end_ml_com_1:
		end_ml_com_11=string.split(end_ml_com_1,'.')[0]+'.'+str(int(string.split(end_ml_com_1,'.')[1])+2)
		text.tag_remove('mlcomment', end_ml_com_11 , deb )
	else:
		text.tag_remove('mlcomment', '1.0' , deb )

	if end_ml_com_2:
		end_ml_com_22=string.split(end_ml_com_2,'.')[0]+'.'+str(int(string.split(end_ml_com_2,'.')[1])+2)
		text.tag_add('mlcomment', deb , end_ml_com_22 )
		text.tag_raise('mlcomment' , aboveThis = None)
	else:
		text.tag_add('mlcomment', deb ,  'end' )
		text.tag_raise('mlcomment' , aboveThis = None)		


def ml_comment_end(fin):
	start_ml_comment_1 = text.search( "/*" , 'insert' , backwards = 'yes' ,stopindex= '1.0')
	start_ml_comment_2 = text.search( "/*" , 'insert' , stopindex='end')	
	end_ml_comment=string.split(fin,'.')[0]+'.'+str(int(string.split(fin,'.')[1])+2)
	if start_ml_comment_1: # si /* trouvé -> comment de start_ml_comment_1 à fin
		text.tag_add('mlcomment',start_ml_comment_1,end_ml_comment)
		text.tag_raise('mlcomment' , aboveThis = None)
	else:
		text.tag_add('mlcomment','1.0',end_ml_comment)
		text.tag_raise('mlcomment' , aboveThis = None)

	if start_ml_comment_2:
		text.tag_remove('mlcomment',end_ml_comment,start_ml_comment_2)
	else:
		text.tag_remove('mlcomment',end_ml_comment, 'end')



def suppr_mlcomment(event=None):
	deb=text.search( MLC ,'insert',backwards = 'yes' , stopindex = '1.0',regexp='yes')
	fin=text.search( MLC ,'insert', stopindex = 'end', regexp='yes')
	if deb and fin:
		fin_2=string.split(fin,'.')[0]+'.'+str(int(string.split(fin,'.')[1])+2)
		if text.get(deb)=="/":
			text.tag_add('mlcomment',deb,fin_2)
			text.tag_raise('mlcomment' , aboveThis = None)
		else:
			text.tag_remove('mlcomment',deb,fin_2)
			text.tag_lower('mlcomment' , belowThis = None)

	elif deb and not fin:
		if text.get(deb)=="/":
			text.tag_add( 'mlcomment', deb , 'end' )
			text.tag_raise('mlcomment' , aboveThis = None)
		else:
			text.tag_remove( 'mlcomment', deb , 'end' )
			text.tag_lower('mlcomment' , belowThis = None)

	elif not deb and fin:
		fin_2=string.split(fin,'.')[0]+'.'+str(int(string.split(fin,'.')[1])+2)
		if text.get(fin)=="/":
			text.tag_remove( 'mlcomment', '1.0' , fin_2 )
			text.tag_lower('mlcomment' , belowThis = None)
		else:
			text.tag_add( 'mlcomment', '1.0' , fin_2 )
			text.tag_raise('mlcomment' , aboveThis = None)

	elif not deb and not fin:
		text.tag_remove( 'mlcomment', '1.0' , 'end' )
		text.tag_lower('mlcomment' , belowThis = None)


def color(p,kw_loc_start,kw_loc_end):
	try:
		if p.lastgroup=='obj':
			text.tag_add ('keyword',kw_loc_start,kw_loc_end)
		elif p.lastgroup=='objbold':
			text.tag_add ('keywordbold',kw_loc_start,kw_loc_end)
		elif p.lastgroup=="num":
			text.tag_add ('num',kw_loc_start,kw_loc_end)
		elif p.lastgroup=="slcom":
			text.tag_add ('slcomment',kw_loc_start,kw_loc_end)
		elif p.lastgroup=="mlcom":
			text.tag_add ('mlcomment',kw_loc_start,kw_loc_end)
		elif p.lastgroup=="math":
			text.tag_add ('math',kw_loc_start,kw_loc_end)
		elif p.lastgroup=="squig":
			text.tag_add ('squig',kw_loc_start,kw_loc_end)
		elif p.lastgroup=="dirc":
			text.tag_add ('dirc',kw_loc_start,kw_loc_end)
		elif p.lastgroup=="str":
			text.tag_add ('str',kw_loc_start,kw_loc_end)
		elif p.lastgroup=="pattern":
			text.tag_add ('pattern',kw_loc_start,kw_loc_end)
	except AttributeError:
		pass




def match_brace(arg=None):
	old = text.index(INSERT)
	n = PositionParentheseAssociee(text)
	text.tag_remove('parenthese',0.1,END)
	if n != None :
		text.tag_add('parenthese',old)
		text.tag_config("parenthese", background = 'red',foreground='blue',relief = GROOVE )
		text.tag_add('parenthese',text.index(n))
		text.tag_config("parenthese", background = 'red',foreground='blue',relief = GROOVE )


def PositionParentheseAssociee(event=None):
	counter=0
	size = len(text.get(1.0,END))
	posini = text.index(INSERT)
	if text.get('insert') == '{':
		index = 0
		for i in range(size):
			counter=counter+1
			if counter>1500:
				break
			if text.get('insert + %d chars'%i) == '}': index = index - 1
			if text.get('insert + %d chars'%i) == '{': index = index + 1
			if index == 0 :
				return 'insert + %d chars'%i
	if  text.get('insert') == '}':
		index = 0
		ok = 0
		for i in range(size):
			counter=counter+1
			if counter>1500:
				break
			if  text.get('insert - %d chars'%i) == '}': index = index - 1;ok=1
			if  text.get('insert - %d chars'%i) == '{': index = index + 1;ok=1
			if index == 0 and ok == 1:
				return 'insert - %d chars'%i
	return None



##########################################################
## ---------------------- Open file ----------------------
##########################################################


def get_script():
	global scene_idx
	scene_contents=text.get('1.0','end')				# checks if text is empty
	if len(scene_contents)==1 : 					# if text is empty
		pass
	else:            						# if text is not empty
		scene_contents=text.get("1.0","end")
		filename=General_Buffer[scene_idx][0]
		store_scene(scene_contents,scene_idx,filename)



def add_path(workingdir):						# used to add path in '-l/path/to/the/current/file
	if len(includedir)==0:
		includedir.append(workingdir)
	else:
		if workingdir in includedir:
			pass
		else:
			includedir.append(workingdir)
			ini_list[4][1]=workingdir


def warning_already_open():
	W_already_open=Pmw.MessageDialog(title= Trans_Warning ,message_text = Trans_Already_Open)
	W_already_open.resizable(width=0,height=0)
	W_already_open.geometry(newGeometry="+250+150")


def info_nb_of_changes(counter,string_to_find,replace_string):
	Pmw.MessageDialog(title='Info....',message_text = Trans_Replaced_DLG % (string_to_find,replace_string,counter))



def warning_large_file(filename):
	W_large_file_TL=Pmw.MessageDialog(title= Trans_Warning,
			message_text =Trans_Too_Large %int(os.stat(filename)[6]/1000))

	W_large_file_TL.resizable(width=0,height=0)
	W_large_file_TL.geometry(newGeometry="+250+150")
	i=0
	color_onload(i)


def open_homepage(event=None):
	com=ini_list[2][1]+" http://pyvon.sourceforge.net"
	popen2.Popen3(com)




def destroy_About():
	About_Tl.destroy()


def Keys_help():
	global About_Tl
	About_Tl=Toplevel()
	About_Tl.resizable(width=0,height=0)
	About_Tl.geometry(newGeometry="+50+150")
	About_Tl.title( Trans_KB_SC )
	F_text1=Frame(About_Tl,relief=RAISED)
	W_help=Text(F_text1,width=70,font=('Lucida Console',11),background=ini_list[19][1])
#text.tag_config('keyword', font=(ini_list[5][1], ini_list[6][1], 'normal' ),foreground=ini_list[10][1])
	W_help.tag_configure('bold_red',font=('Lucida Console', 11, 'bold' ),foreground='red')
	W_help.tag_configure('blue',font=('Lucida Console', 11, 'normal' ),foreground='blue') #,background='yellow')
	chem=os.path.join(langDir,'help.txt')
	fd=open(chem,'r')
	for line in fd.readlines():
		if line.endswith('\r\n'):
			line = line[:-2] + "\n"

		W_help.insert('end',line)


	W_help.configure(relief='raised')
	fd.close()
	yscroll_i = Scrollbar(F_text1,orient=VERTICAL, command=W_help.yview)
	yscroll_i.pack(side=RIGHT, fill=Y)
	xscroll_i.pack(side=BOTTOM, fill=X)
	W_help['yscrollcommand']=yscroll_i.set

	W_help.tag_add('bold_red','4.0','4.10')
	W_help.tag_add('bold_red','6.0','6.10')
	W_help.tag_add('bold_red','8.0','8.10')
	W_help.tag_add('bold_red','10.0','10.10')
	W_help.tag_add('bold_red','12.0','12.16')
	W_help.tag_add('bold_red','14.0','15.10')
	W_help.tag_add('bold_red','17.0','17.20')
	W_help.tag_add('bold_red','20.0','20.10')
	W_help.tag_add('bold_red','22.0','22.10')
	W_help.tag_add('bold_red','24.0','24.10')
	W_help.tag_add('bold_red','26.0','26.10')
	W_help.tag_add('bold_red','28.0','28.10')
	W_help.tag_add('bold_red','30.0','30.10')



	W_help.tag_add('blue','4.26','4.end')
	W_help.tag_add('blue','6.26','6.end')
	W_help.tag_add('blue','8.26','8.end')
	W_help.tag_add('blue','10.26','10.end')
	W_help.tag_add('blue','12.26','12.end')
	W_help.tag_add('blue','13.26','13.end')
	W_help.tag_add('blue','15.26','15.end')
	W_help.tag_add('blue','17.26','17.end')
	W_help.tag_add('blue','18.26','18.end')
	W_help.tag_add('blue','20.26','20.end')
	W_help.tag_add('blue','22.26','22.end')
	W_help.tag_add('blue','24.26','24.end')
	W_help.tag_add('blue','26.26','26.end')
	W_help.tag_add('blue','28.26','28.end')
	W_help.tag_add('blue','30.26','30.end')

	W_help.configure(state='disabled')

	W_help.pack()
	F_text1.pack()


def about():
	global About_Tl
	About_Tl=Toplevel()
	About_Tl.resizable(width=0,height=0)
	About_Tl.geometry(newGeometry="+250+150")
	About_Tl.title('About ....')
	F_text1=Frame(About_Tl,relief=RAISED)
	W_about=Canvas(F_text1,relief=RAISED,cursor="circle",height=270)
	W_about.create_text(150,30,text='P Y V O N  Version : 1.3',font=("helvetica",14,'bold'),fill='red')
	W_about.create_text(150,60,text='Linux-Unix GUI scene editor for POV-RAY',fill='blue')
	W_about.create_text(150,80,text='send bugs report, suggestions to :')
	W_about.create_text(150,100,text='povray.unix forum or fabien.henon@caramail.com')
	W_about.create_text(150,120,text='Credits to Ingo for his syntax highlighting code')
	W_about.create_text(150,140,text='Credits to GP Guignot for his brace-matching code')
	W_about.create_text(150,160,text='Thanks to the translators : Bonsai, Ingo,')
	W_about.create_text(150,180,text='John Coppens, Alfonso,')
	W_about.create_text(150,200,text='Alessandro Falappa and Slawomir Szczyrba')


	W_about.create_text(150,240,text='Special Thanks to Céline for her patience....')


#	W_about.create_text(150,260,text='Homepage :')
	W_about.grid(column=0)
#	W_about.tag_config('link',foreground='yellow')
	W_about_link=Canvas(F_text1,relief=RAISED,cursor="hand2",height=30)
	W_about_link.create_text(150,15,text='http://pyvon.sourceforge.net',font=("arial",15,'normal'),tags="highl",fill='blue')
#	W_a.bind('highl',"<ButtonPress-1>",quitter)
	W_about_link.tag_bind('highl',"<Button-1>",open_homepage)
	W_about_link.grid()
	Button(F_text1,text= Trans_Close , command=destroy_About).grid(padx=3,pady=3)
	F_text1.grid()




def Partial_Render():
	global SC_scale,EC_scale,SR_scale,ER_scale
	global drawing_box,Partial_Render_TL
	global Activate_partial_render
	Partial_Render_TL=Toplevel(height=460,width=460)

	Partial_Render_TL.resizable(width=0,height=0)
	Partial_Render_TL.geometry(newGeometry="+150+150")
	Partial_Render_TL.title(Trans_Partial_Render)
	
	drawing_box=Canvas(Partial_Render_TL,width=320,height=240,relief=FLAT)
	drawing_box.grid(column=1,row=1)

	Text_SC=Label(Partial_Render_TL,text="sc:").grid(column=0,row=2,sticky=S)
	Text_EC=Label(Partial_Render_TL,text="ec:").grid(column=0,row=3,sticky=S)
	Text_SC=Label(Partial_Render_TL,text="sr:").grid(column=2,row=0,sticky=SE)
	Text_EC=Label(Partial_Render_TL,text="er:").grid(column=3,row=0,sticky=SE)


	SC_scale=Scale(Partial_Render_TL,orient='horizontal',length=350,from_=0,to=1,resolution=0.001,command=SC_line)
	SC_scale.set(start_column)
	SC_scale.grid(row=2,column=1,sticky=NW)
	EC_scale=Scale(Partial_Render_TL,orient='horizontal',length=350,from_=0,to=1,resolution=0.001,command=EC_line)
	EC_scale.set(end_column)
	EC_scale.grid(row=3,column=1,sticky=NW)
	SR_scale=Scale(Partial_Render_TL,orient='vertical',length=270,from_=0,to=1,resolution=0.001,command=SR_line)
	SR_scale.set(start_row)
	SR_scale.grid(column=2,row=1,rowspan=2,sticky=NE)
	ER_scale=Scale(Partial_Render_TL,orient='vertical',length=270,from_=0,to=1,resolution=0.001,command=ER_line)
	ER_scale.set(end_row)
	ER_scale.grid(column=3,row=1,rowspan=2,sticky=NE)


	Button_partial_render_OK=Button(Partial_Render_TL,text='OK',command=get_partial)
	Button_partial_render_OK.grid(columnspan=2,row=3,column=2)



def SC_line(ki):
	try:
		drawing_box.delete('Front')
		drawing_box.delete('Back')
	except:
		pass
	SC=int(SC_scale.get()*320+1)
	EC=int(EC_scale.get()*320)
	SR=int(SR_scale.get()*240+1)	
	ER=int(ER_scale.get()*240)
	Back=drawing_box.create_rectangle(0,0,320,240,fill='gray')
	Front=drawing_box.create_rectangle(SC,SR,EC,ER,outline='blue')



def EC_line(ki):
	try:
		drawing_box.delete('Front')
		drawing_box.delete('Back')
	except:
		pass
	SC=int(SC_scale.get()*320+1)
	EC=int(EC_scale.get()*320)
	SR=int(SR_scale.get()*240+1)	
	ER=int(ER_scale.get()*240)
	Back=drawing_box.create_rectangle(0,0,320,240,fill='gray')
	Front=drawing_box.create_rectangle(SC,SR,EC,ER,outline='blue')




def SR_line(ki):
	try:
		drawing_box.delete('Front')
		drawing_box.delete('Back')
	except:
		pass
	SC=int(SC_scale.get()*320+1)
	EC=int(EC_scale.get()*320)
	SR=int(SR_scale.get()*240+1)	
	ER=int(ER_scale.get()*240)
	Back=drawing_box.create_rectangle(0,0,320,240,fill='gray')
	Front=drawing_box.create_rectangle(SC,SR,EC,ER,outline='blue')




def ER_line(ki):
	try:
		drawing_box.delete('Front')
		drawing_box.delete('Back')
	except:
		pass
	SC=int(SC_scale.get()*320+1)
	EC=int(EC_scale.get()*320)
	SR=int(SR_scale.get()*240+1)	
	ER=int(ER_scale.get()*240)
	Back=drawing_box.create_rectangle(0,0,320,240,fill='gray')
	Front=drawing_box.create_rectangle(SC,SR,EC,ER,outline='blue')





def get_partial():

	global start_column
	global end_column
	global start_row
	global end_row		
	start_column=SC_scale.get()
	end_column=EC_scale.get()
	start_row=SR_scale.get()
	end_row=ER_scale.get()
	if start_column > end_column:
		import tkMessageBox
		MSG=tkMessageBox.showwarning(Trans_Warning,"The start column value is larger\nthan the end column value")
		Partial_Render_TL.update()
		Partial_Render_TL.focus_force()	        	
	elif start_row > end_row:
		import tkMessageBox
		tkMessageBox.showwarning(Trans_Warning,"The start row value is larger\nthan the end row value")
		Partial_Render_TL.update()
		Partial_Render_TL.focus_force()	        
	else:		
		Partial_Render_TL.destroy()






def check_file_size(filename):
	if os.stat(filename)[6]>500000:
		warning_large_file(filename)

		
def open_file(event=None):
	global includedir
	get_script()

	scene_contents=""
	import tkFileDialog
	global filename
	filename_bk=filename

	filetypes = ( ( "Pov files","*.pov"),( "Pov files","*.POV"), ( "Include files","*.inc"), ( "Include files","*.INC"),( "All","*"))
	filename = tkFileDialog.askopenfilename(filetypes=filetypes,initialdir=ini_list[4][1])
	if filename:
		check_file_size(filename)
		if len(General_Buffer)>0:
			for n in range(len(General_Buffer)):
				if General_Buffer[n][0]==filename:
					warning_already_open()
					filename=filename_bk
					return
	    
		text.delete(1.0, END)
		fd = open(filename)
		scene_contents=string.replace(fd.read(),'\r\n','\n')
		fd.close()


		if len(Oldfiles)>6:
			del Oldfiles[0]
		if filename not in Oldfiles:
			Oldfiles.append(filename)

            
		ini_list[4][1]=os.path.split(filename)[0] # add path to workingdir
		workingdir=ini_list[4][1]
		add_path(workingdir)
            
		Files_list.append(os.path.split(filename)[1])
		drop_down_file.setlist(Files_list) # updates the combo-box
		num=len(Files_list)
		drop_down_file.selectitem(num-1,1) # show which file is selected in the combo-box

		if ini_list[8][1]=="Yes":
			tag_it(scene_contents)
			paste=0
			resyntax=0
			colorize(scene_contents,taglist,paste,resyntax)
			store_scene_open(scene_contents,filename)
		if ini_list[8][1]=="No":
			text.insert('1.0',scene_contents)
			store_scene_open(scene_contents,filename)
	else:
		try: 
			filename=filename_bk

		except NameError:
			pass


def open_old_file(oldfile):
	global includedir
	get_script()
	scene_contents=""
	import tkFileDialog

	check_file_size(oldfile)

	if len(General_Buffer)>0:
		for n in range(len(General_Buffer)):
			if General_Buffer[n][0]==oldfile:
				warning_already_open()
				return

	filename_bk=oldfile
	text.delete(1.0, END)
	fd = open(oldfile)
	scene_contents=string.replace(fd.read(),'\r\n','\n')
	fd.close()


	ini_list[4][1]=os.path.split(oldfile)[0] 	# add path to workingdir
	workingdir=ini_list[4][1]
	add_path(workingdir)

	Files_list.append(os.path.split(oldfile)[1])
	drop_down_file.setlist(Files_list) 		# updates the combo-box
	num=len(Files_list)
	drop_down_file.selectitem(num-1,1) 		# show which file is selected in the combo-box

	if ini_list[8][1]=="Yes":
		tag_it(scene_contents)
		paste=0
		resyntax=0
		colorize(scene_contents,taglist,paste,resyntax)
		store_scene_open(scene_contents,oldfile)
	if ini_list[8][1]=="No":
		text.insert('1.0',scene_contents)
		store_scene_open(scene_contents,oldfile)
	else:
		try: 
			oldfile=filename_bk
		except NameError:
			pass




def tag_it(scene_contents):
	ident_flag = 0
	global taglist
	taglist=[]
	start=0
	while 1:
		p=pattern.search(scene_contents, start)
		try:
			p.group()
		except:
			break                           #no more patterns that match.
		start=p.end()
		if p.lastgroup == 'mlcom':          #find end of comment and
			lastgr=p.lastgroup              #deal with nested comments or
			s=p.start()                     #missing end of comment.
			end=find_end_comment(scene_contents, start)
			taglist.append((s,end,lastgr))
			start=end
		elif p.lastgroup == 'str':          #find end of string or deal
			s=p.start()                     #with missing end of string.
			lastgr=p.lastgroup
			p=strpattern.search(scene_contents, start)
			if not p:
#                     '\n#warning: Looks like there is a " missing'
				end=len(scene_contents)
				taglist.append((s,end,lastgr))
				break
			taglist.append((s,p.end(),p.lastgroup))
			start=p.end()
		else:
			taglist.append((p.start(),p.end(),p.lastgroup))
			if p.lastgroup == 'dirc':
				if p.group() == '#declare' or p.group() == '#macro':
					ident_flag=1
	return taglist




def colorize(scene_contents,taglist,paste,resyntax):
	if resyntax == 1 :
		try:
			text.delete(SEL_FIRST,SEL_LAST)
			text.insert('insert',scene_contents)
		except:
			text.delete('1.0','end')
			text.insert('insert',scene_contents)

	elif resyntax==0 :
		text.insert('insert',scene_contents)


	if paste==0:
		line=1
		col=0
		end_line=0
	elif paste==1:
		line=int(l_sel)
		col=int(c_sel)
		end_line=0
	resyntax=0
	paste=0
	for tag in taglist:
		if tag[2]=="end_of_line":
			line=line+1
			col=0
			end_line=tag[1]
		elif tag[2]=="obj":
			tag_deb=string.join((str(line),str(int(tag[0]-end_line+col))),'.')
			tag_fin=string.join((str(line),str(int(tag[1]-end_line+col))),'.')
			text.tag_add ('keyword',tag_deb,tag_fin)
		elif tag[2]=="objbold":
			tag_deb=string.join((str(line),str(int(tag[0]-end_line+col))),'.')
			tag_fin=string.join((str(line),str(int(tag[1]-end_line+col))),'.')
			text.tag_add ('keywordbold',tag_deb,tag_fin)
		elif tag[2]=="pattern":
			tag_deb=string.join((str(line),str(int(tag[0]-end_line+col))),'.')
			tag_fin=string.join((str(line),str(int(tag[1]-end_line+col))),'.')
			text.tag_add ('pattern',tag_deb,tag_fin)
		elif tag[2]=="num":
			tag_deb=string.join((str(line),str(int(tag[0]-end_line+col))),'.')
			tag_fin=string.join((str(line),str(int(tag[1]-end_line+col))),'.')
			text.tag_add ('num',tag_deb,tag_fin)
		elif tag[2]=="slcom":
			tag_deb=string.join((str(line),str(int(tag[0]-end_line+col))),'.')
			tag_fin=string.join((str(line),str(int(tag[1]-end_line+col))),'.')
			text.tag_add ('slcomment',tag_deb,tag_fin)
		elif tag[2]=="mlcom":
			tag_deb=string.join((str(line),str(int(tag[0]-end_line+col))),'.')
			contents=text.get(tag_deb,'end')
			com_fin=text.search("*/",tag_deb)
			col_fi=int(string.split(com_fin,'.')[1])+2
			li_fi=int(string.split(com_fin,'.')[0])
			tag_fin=str(string.join((str(li_fi),str(col_fi)),'.'))
			text.tag_add ('mlcomment',tag_deb,tag_fin)
			line=int(string.split(tag_fin,'.')[0])
			col=0
		elif tag[2]=="str":
			tag_deb=string.join((str(line),str(int(tag[0]-end_line+col))),'.')
			tag_fin=string.join((str(line),str(int(tag[1]-end_line+col))),'.')
			text.tag_add ('str',tag_deb,tag_fin)
		elif tag[2]=="squig":
			tag_deb=string.join((str(line),str(int(tag[0]-end_line+col))),'.')
			tag_fin=string.join((str(line),str(int(tag[1]-end_line+col))),'.')
			text.tag_add ('squig',tag_deb,tag_fin)
		elif tag[2]=="math":
			tag_deb=string.join((str(line),str(int(tag[0]-end_line+col))),'.')
			tag_fin=string.join((str(line),str(int(tag[1]-end_line+col))),'.')
			text.tag_add ('math',tag_deb,tag_fin)
		elif tag[2]=="dirc":
			tag_deb=string.join((str(line),str(int(tag[0]-end_line+col))),'.')
			tag_fin=string.join((str(line),str(int(tag[1]-end_line+col))),'.')
			text.tag_add ('dirc',tag_deb,tag_fin)



def store_scene_untitled():
	global General_Buffer
	global Buffer
	global scene_idx
	Buffer=['untitled.pov', '', [('keyword', ()), ('keywordbold', ()), ('pattern', ()),
	 ('num', ()), ('slcomment', ()), ('mlcomment', ()),('str', ()), ('squig', ()), ('math', ()), ('dirc', ())], 'status']
	General_Buffer.append(Buffer)
	Buffer=[0,0,0,0]



def store_scene_open(scene_contents,filename):
	global scene_idx
	global General_Buffer
	global Buffer
	global scene_idx
	scene_color=[]
	for tagname in color_tags:
		scene_color.append((tagname,text.tag_ranges(tagname)))
	Buffer[0]=filename
	Buffer[1]=scene_contents
	Buffer[2]=scene_color
	Buffer[3]="status"
	General_Buffer.append(Buffer)
	scene_idx=len(Files_list)-1
	Buffer=[0,0,0,0]
	ini_list[8][1]='Yes'





############################################################
#--       Store all the information about a scene        --#
############################################################

def store_scene(scene_contents,scene_idx,filename):
	global General_Buffer
	global Buffer
	scene_color=[]
	for tagname in color_tags:
		scene_color.append((tagname,text.tag_ranges(tagname)))
	Buffer[0]=filename
	Buffer[1]=scene_contents
	Buffer[2]=scene_color
	Buffer[3]="status"
	General_Buffer[scene_idx]=[Buffer[0],Buffer[1],Buffer[2],Buffer[3]]
	Buffer=[0,0,0,0]


############################################################
#--   Select the resolution chosen in the drop down box  --#
############################################################

def click_Resolution(Resol):
	global preset_resol
	if   Resol=="1280 x 1024":preset_resol=" +w1280 +h1024"
	elif Resol=="1024 x 768" :preset_resol=" +w1024 +h768"
	elif Resol== "800 x 600" :preset_resol=" +w800 +h600"
	elif Resol== "640 x 480" :preset_resol=" +w640 +h480"
	elif Resol== "400 x 300" :preset_resol=" +w400 +h300"
	elif Resol== "320 x 240" :preset_resol=" +w320 +h240"
	elif Resol== "128 x 96"  :preset_resol=" +w128 +h96"
	else :preset_resol=""
	return preset_resol




############################################################
#--   Select which file is chosen in the file selector   --#
############################################################

def drpdwn_choose_file_error(Lin,i):
	global scene_idx
	scene_contents=text.get("1.0","end")
	if scene_contents.endswith("\n"):
		scene_contents=scene_contents[:-1]
	filename=General_Buffer[scene_idx][0]
	store_scene(scene_contents,scene_idx,filename)
	scene_idx=i #drop_down_file.index('active')
	restore_scene(scene_idx)
	drop_down_file.selectitem(scene_idx,1)


def error_in_file(file_error,line):
	global tl
	tl=Toplevel()
	tl.title("Error !")
	tl.resizable(width=0,height=0)
	tl.geometry(newGeometry="+250+150")
	warninglb=Label(tl, font=('arial', 12,'bold'),text= Trans_File_Error % (file_error,line),justify='left')
	ok=Button(tl, text= Trans_Close , fg="red", command=destroywarning)
	ok.pack(side=BOTTOM,padx=4,pady=10)
	warninglb.pack()




def drpdwn_choose_file(Lin):
	global scene_idx
	scene_contents=text.get("1.0","end")
	if scene_contents.endswith("\n"):
		scene_contents=scene_contents[:-1]
	filename=General_Buffer[scene_idx][0]
	store_scene(scene_contents,scene_idx,filename)
	scene_idx=drop_down_file.index('active')
	restore_scene(scene_idx)

def switch_file(Lin):
	global scene_idx
	global new_idx

	scene_contents=text.get("1.0","end")
	if scene_contents.endswith("\n"):
		scene_contents=scene_contents[:-1]
	filename=General_Buffer[scene_idx][0]
	store_scene(scene_contents,scene_idx,filename)

	restore_scene(new_idx)
	scene_idx=new_idx
	drop_down_file.selectitem(scene_idx,1)


####################################################
#--   Paste scene & color into the editor        --#
####################################################

def restore_scene(scene_idx):
	global filename
	tag_deb=''
	tag_fin=''
	for i in color_tags:
		text.tag_remove(i,'1.0','end')
	text.delete('1.0','end')
	text.insert('insert',General_Buffer[scene_idx][1])
	filename=General_Buffer[scene_idx][0]
	for tagname in General_Buffer[scene_idx][2]:
		for word_type in color_tags:
			if tagname[0]==word_type:
				for tags in xrange( 0, len(tagname[1]), 2 ):
					tag_deb=tagname[1][0+tags]
					tag_fin=tagname[1][1+tags]
					text.tag_add (word_type,tag_deb,tag_fin)




####################################################
#--                   New File                   --#
####################################################--

def clr(event=None):
	global scene_idx
	if len(General_Buffer)==1 and General_Buffer[0][0]=='untitled.pov':
		pass
	else:
		filename=General_Buffer[scene_idx][0]
		store_scene_untitled()
		scene_contents=text.get('1.0','end')
		store_scene(scene_contents,scene_idx,General_Buffer[scene_idx][0])				
		filename="untitled.pov"
		store_scene_untitled()
	        text.delete(1.0, END)
	        Files_list.append("untitled.pov")
		scene_idx=len(Files_list)-1
	        drop_down_file.setlist(Files_list)
	        drop_down_file.selectitem(len(Files_list)-1,1)	



####################################################
#--                   Close File                 --#
####################################################

def close_checking():
	global Close_check_TL
	Close_check_TL=Toplevel()
	Close_check_TL.resizable(width=0,height=0)
        Close_check_TL.geometry(newGeometry="+250+150")
	Close_check_TL.title(Trans_Warning)
	Label(Close_check_TL,text= Trans_Really_close_File ).grid(row=0,columnspan=2)
	Button(Close_check_TL,text= Trans_Yes ,width=10,command=Close_File_Yes).grid(row=1,column=0)
	Button(Close_check_TL,text= Trans_No ,width=10,command=Close_File_No).grid(row=1,column=1)

def Close_File_Yes():
	close()





def Close_File_No():
	Close_check_TL.destroy()
	return


def close():
        global scene_idx
	global General_Buffer
	global Files_list
	Close_check_TL.destroy()
	if len(General_Buffer)==1:
	        for i in color_tags:
	            text.tag_remove(i,'1.0','end')
	        text.delete('1.0','end')
	        General_Buffer=[]
	        Files_list=[]
		scene_idx=0
	        store_scene_untitled()
		Files_list.append('untitled.pov')
		drop_down_file.setlist(Files_list)
		drop_down_file.selectitem(0)

	elif len(General_Buffer)>1:
	        for i in color_tags:
	            text.tag_remove(i,'1.0','end')
	        text.delete('1.0','end')
	        del General_Buffer[scene_idx]
	        del Files_list[scene_idx]
		scene_idx=0
	        restore_scene(scene_idx)
		drop_down_file.setlist(Files_list)
		drop_down_file.selectitem(scene_idx)


####################################################
#--                       Save                   --#
####################################################

def save(event=None):
	global filename
	global scene_idx


	scene_contents = text.get(1.0, END)
        if scene_contents:
		try:
			fd = open(General_Buffer[scene_idx][0], 'w')
	                for line in string.split(scene_contents, '\n'):
			    fd.write(line)
			    fd.write('\n')
	                fd.close
		except:
			write_permissions()



def write_permissions():
	global Check_w
	Check_w=Toplevel()
	Check_w.resizable(width=0,height=0)
        Check_w.geometry(newGeometry="+250+150")
	Check_w.title('Warning !')

	Label(Check_w,text='File could not be written.\nCheck your writing permissions').grid(row=0,column=0)
	Button(Check_w,text='OK', command=destroy_Check_w).grid(row=1,column=0)

def destroy_Check_w():
	Check_w.destroy()


####################################################
#--                  Save As                     --#
####################################################

def saveas(event=None):
	global filename

	scene_contents = text.get(1.0, END)
        import tkFileDialog
        filetypes = ( ( "Pov files","*.pov"), ( "Include files","*.inc"), ( "All","*"))

        filename = tkFileDialog.asksaveasfilename(filetypes=filetypes)
	if filename:
		General_Buffer[scene_idx][0]=filename
		Files_list[scene_idx]=os.path.split(filename)[1]
		drop_down_file.setlist(Files_list)
		drop_down_file.selectitem(scene_idx,1)
		scene_contents=text.get("1.0","end")

		fd = open(filename, 'w')
	        for line in string.split(scene_contents, '\n'):
		    fd.write(line)
		    fd.write('\n')
	        fd.close
	else:
		pass


##############################################################
#--       Checks if the renderer is running or not         --#
##############################################################


def render(event=None):
        global sts_rendu
	if os.path.isfile(ini_list[0][1])==0 or ini_list[0][1]=="":
		loadinierror()

	else:
		if ini_list[7][1]=="Yes":
		    save()
		if sts_rendu==0:
		     runpov()
		     button_render.configure(image=stop)
		     sts_rendu=1
		     rendermenu.delete(0)
		     rendermenu.insert_command(0,label= Trans_Stop,underline=12, command=render)
		elif sts_rendu==1:
		     button_render.configure(image=lancer)
		     sts_rendu=0
		     rendermenu.delete(0)
		     rendermenu.insert_command(0,label= Trans_Render ,underline=12, command=render)	     
		     kill()
	     

##############################################################
#--                     Render section                     --#
##############################################################


def runpov():
        global includedir
	global preset_resol

	global start_column
        totalpath=''
	if sys.platform!="win32":
		os.system('rm -f ${HOME}/err.txt')
	else:
		pass
        for i in includedir:
            totalpath=totalpath+' -l'+i
        ini_list[3][1] = fill_params.get()
	QS=Quality_scale.get()
	Resol=drop_down_resolution.get()
	click_Resolution(Resol)
	if APR.get()==1:Partial_render= " +sc"+str(start_column)+" +ec"+str(end_column)+" +sr"+str(start_row)+" +er"+str(end_row)
	else : Partial_render=""
	if FP.get()==1:Fast_preview=" +sp64 +ep2"
	else:Fast_preview=""
  	c = ['Input_File_Name="',str(General_Buffer[scene_idx][0]),'" ',ini_list[3][1],totalpath," +GF${HOME}/err.txt +q",str(QS),Fast_preview,preset_resol,Partial_render]
        cmd=string.replace((string.join(c,sep=',')),',','')
        total=ini_list[0][1]+" "+cmd
        global run
        run=popen2.Popen3(total)
        global PID
        PID= run.pid
        tempo()
        text.tag_remove('highl','1.0','end')




##############################################################
#--               Rendering Error Section                  --#
##############################################################

def tempo():
     if run.poll()!=-1:
        finished()
	try:
		homedir = os.path.expanduser('~')
	        fd=open(homedir+"/err.txt")
	        for line in fd.readlines():
			if line.startswith("File"):  # POV 3.5
				pos_lin=string.rfind(line,"Line: ")
				pos_file=string.find(line,"File:")
				file_error=line[pos_file+6:pos_lin-2]
				line_err=line[pos_lin+6:]


		try:
			if General_Buffer[scene_idx][0]==file_error: # si c'est la scene courante
				try:
					text.tag_add ('highl','%d.0'% int(line_err),'%d.end'% int(line_err) )
					text.see('%d.0'% (int(line_err)+5))
				except ValueError:
				        pass

			else:
				l=[]
				for i in range(len(General_Buffer)):
					l.append(General_Buffer[i][0])
				if file_error in l:
					f=os.path.split(General_Buffer[i][0])[1]
					ind=l.index(file_error)
					drpdwn_choose_file_error(f,ind)
					try:
						text.tag_add ('highl','%d.0'% int(line_err),'%d.end'% int(line_err) )
						text.see('%d.0'% (int(line_err)+5))
					except ValueError:
						pass
				else:
					error_in_file(file_error,line_err)
		except:
			pass


	except IOError:
		pass

     else:
		root.after(200,tempo)



###############################################################
#--                    Kill Render                          --#
###############################################################

def kill(event=None):
	     os.system("killall -9 %s"%ini_list[0][1])
	     PID=0




###############################################################
#--                      Pause Render                       --#
###############################################################

def pause(event=None):
   global sts_rendu
   if sts_rendu==0:
      pass
   else:      
      global sts_pause
      if sts_pause==0:
	   rendermenu.delete(1)
           rendermenu.insert_command(1,label= Trans_Resume , command=pause)
           os.system("kill -19 %d" % PID)
	   sts_pause=1
           button_pause.configure(image=resumee,relief=SUNKEN)

      elif sts_pause==1:
	   rendermenu.delete(1)
           rendermenu.insert_command(1,label= Trans_Pause , command=pause)
           os.system("kill -18 %d" % PID)
	   sts_pause=0
           button_pause.configure(image=pausee,relief=FLAT)



###############################################################
#--       When the raytrace is finished or stopped-          -#
###############################################################

def finished():
        global sts_rendu
	global re_nice	
        sts_rendu=0
	re_nice=0
        button_render.configure(image=lancer)




###############################################################
#--                    Replace/Find box                     --#
###############################################################

def box_replace_word(event=None):
	global on_off
	on_off=IntVar()
	off_off=0
	global find_word_TL
	find_word_TL=Toplevel()
	find_word_TL.resizable(width=0,height=0)
	find_word_TL.geometry(newGeometry="+250+150")			     
	find_word_TL.title( Trans_Find_Replace )
	Label(find_word_TL,text='').grid(row=0,column=0)
	Label(find_word_TL,text= Trans_Find ).grid(row=1,column=0,sticky=W)
	Label(find_word_TL,text= Trans_Replace_BY ).grid(row=2,column=0,sticky=W)	
	global To_find_word
	To_find_word=Entry(find_word_TL,width=30)
	To_find_word.focus_set()
	To_find_word.grid(row=1,column=1)
	global To_replace_word
	To_replace_word=Entry(find_word_TL,width=30)
	To_replace_word.grid(row=2,column=1)
	Label(find_word_TL,text= Trans_Match_Case ).grid(row=4,column=0)
	Match_CASE=Checkbutton(find_word_TL,variable=on_off)
	Match_CASE.grid(row=4,column=1,sticky=W)
	Button(find_word_TL,text= Trans_Next,command=find_next_word).grid(row=1,column=3,sticky=E+W,padx=2,pady=2)	
	Button(find_word_TL,text= Trans_Replace ,command=replace_word).grid(row=2,column=3,sticky=E+W,padx=2,pady=2)	
	Button(find_word_TL,text= Trans_Replace_All ,command=replace_all_words).grid(row=3,column=3,sticky=E+W,padx=2,pady=2)		
	Label(find_word_TL,text='').grid(row=3,column=0)
	Button(find_word_TL,text= Trans_Cancel ,command=destroy_find_box).grid(row=4,column=3,padx=2,pady=2,sticky=E+W)
	find_word_TL.bind("<Escape>", destroy_find_box)
	find_word_TL.bind("<Return>", find_next_word)



def find_next_word(event=None):
	global To_find_word
	global pos_ini_to_find
	global start_pos
	global end_pos
	global on_off
	if on_off.get()==1:MAJ=0
	else:MAJ=1
	text.tag_remove(SEL,'1.0','end')
	string_to_find=To_find_word.get()
	if string_to_find=='':return
	try:
		start_pos=text.search(string_to_find,pos_ini_to_find,nocase=MAJ)
	except:
		pos_ini_to_find='1.0'		
		start_pos=text.search(string_to_find,pos_ini_to_find,nocase=MAJ)
	try:
		end_col=int(string.split(start_pos,'.')[1])+len(string_to_find)
		end_pos=string.join(((string.split(start_pos,'.')[0]),str(end_col)),'.')
		text.tag_add(SEL,start_pos,end_pos)
		text.see(start_pos)
		pos_ini_to_find=end_pos	
	except:
		root.bell()
		return



def replace_word():
	global To_find_word
	global pos_ini_to_find
	global To_replace_word
	global start_pos
	global end_pos
	find_next_word()
	try:
		text.delete(start_pos,end_pos)
		text.insert(start_pos,To_replace_word.get())
	except:
		root.bell()

def replace_all_words():
	global To_find_word
	global To_replace_word
	global on_off
	if on_off.get()==1:MAJ=0
	else:MAJ=1
	string_to_find=To_find_word.get()
	pos_ini_to_find='1.0'
	counter=0
	replace_string=To_replace_word.get()
	while text.search(string_to_find,pos_ini_to_find,nocase=MAJ)!='':
		start_pos=text.search(string_to_find,pos_ini_to_find,nocase=MAJ)
		end_col=int(string.split(start_pos,'.')[1])+len(string_to_find)
		end_pos=string.join(((string.split(start_pos,'.')[0]),str(end_col)),'.')
		text.delete(start_pos,end_pos)
		text.insert(start_pos,replace_string)
		text.tag_add(SEL,start_pos,end_pos)
		text.see(start_pos)
		pos_ini_to_find=end_pos	
		counter=counter+1
	info_nb_of_changes(counter,string_to_find,replace_string)







###############################################################
#--           Called with the cancel                        --#
###############################################################

def destroy_find_box(event=None):
	global find_word_TL
	find_word_TL.destroy()	





def undo(event=None):
	global B_undo
	if B_undo!=[]:
		text.delete('1.0','end')
		text.insert('1.0',B_undo[len(B_undo)-2])




def copy(event=None):
	global B_undo
	global taglist_paste
        try:
            global scene_contents
            scene_contents = text.selection_get()
            tag_it(scene_contents)
	    taglist_paste=taglist
        except:
            pass



def cut(event=None):
	global B_undo
	global taglist_paste
        try:
            global sel
            global scene_contents
            scene_contents = text.selection_get()
            tag_it(scene_contents)
	    taglist_paste=taglist
            sel= ""
            sel = text.selection_get()
            text.delete(SEL_FIRST, SEL_LAST)
        except:
            pass




def paste(event=None):
	global B_undo
	global taglist_paste
	pos_ini=text.index(INSERT)
        try:
            global c_sel,l_sel
	    global scene_contents
            l_sel,c_sel = text.index(INSERT).split('.', 1)#3 4
            paste=1
            resyntax=0
            colorize(scene_contents,taglist_paste,paste,resyntax)
	    text.see(pos_ini)        
	except:
            pass
	text.see(pos_ini)



###############################################################
#--              Selects all the text                       --#
###############################################################

def selectall(event=None):
        sel= ""
        text.focus_force()
	sel = text.tag_add("sel",1.0,"end")




###############################################################
#--                   Find /* in mutiline comments          --#
###############################################################

def find_end_comment(scene_contents, start):
        end=start
        end_comment=0
        while end_comment ==0:
            p=mlcompattern.search(scene_contents, end)
            if not p:
                end=len(scene_contents)
                end_comment=1
            elif p.group() == '*/':
                end=p.end()
                end_comment=1
            else:
                end=p.end()
                end=find_end_comment(scene_contents, end)
        return end




###############################################################
#--                Povray executable options                --#
###############################################################



def show_location():
        global tl_exe_loc
	global re_nice
	nice_value=IntVar()
	nice_value.set(re_nice)		
        tl_exe_loc=Toplevel()
	tl_exe_loc.resizable(width=0,height=0)
	Label(tl_exe_loc,text= Trans_POV_Location ).grid(row=0,column=0,sticky=W,padx=2,pady=2)
	Label(tl_exe_loc,text='').grid(row=2,column=0)	
	Label(tl_exe_loc,text= Trans_Render_Niceness ).grid(row=3,column=0,padx=0,pady=0)
        tl_exe_loc.title( Trans_Exe_Opt )
        tl_exe_loc.geometry(newGeometry="+50+150")	
        pov_location=Entry(tl_exe_loc,width=60)
        if ini_list[0][1].endswith("\n"):
            ini_list[0][1] = ini_list[0][1][:-1]
        pov_location.insert(INSERT,ini_list[0][1])
	global scale_exe
	scale_exe=Scale(tl_exe_loc,orient=HORIZONTAL,from_=-20,to=20,
			length=280,
			resolution=1,
			tickinterval=5,
			variable=15
			)
        Search=Button(tl_exe_loc, text= Trans_Search_Location ,command=ask_location)
        Ok=Button(tl_exe_loc, text= Trans_Close,command=quit_ask_location,fg='red')
        tl_exe_loc.focus_force()
	pov_location.grid(row=1,column=0,sticky=W,padx=2,pady=2)
        Ok.grid(row=4,column=3,padx=2,pady=2)
        Search.grid(row=1,column=3,sticky=W,padx=2,pady=2)
        scale_exe.grid(row=4,column=0)



def renice_exe():
        global scale_exe
	global re_nice
	try:
		re_nice=scale_exe.get()
		os.system("renice %d %d" % (int(re_nice),PID)) 
	except NameError:
	        pass



def ask_location():

        import tkFileDialog
        tl_exe_loc.destroy()
        ini_list[3][1]=fill_params.get()
        ini_list[0][1] = tkFileDialog.askopenfilename()
        if ini_list[0][1].endswith("\n"):
            ini_list[0][1] = ini_list[0][1][:-1]
        if ini_list[3][1].endswith("\n"):
            ini_list[3][1] = ini_list[3][1][:-2]
        show_location()




def quit_ask_location():
        renice_exe() 
        tl_exe_loc.destroy()




###############################################################
#--                    Save before exit                     --#
###############################################################

def save_on_exit():
	homedir = os.path.expanduser('~')
        ini_list[3][1]=fill_params.get()

	try:
		fd=open(homedir+'/.pyoldfiles','w')
		for file in Oldfiles:
			fd.write(file+'\n')
		fd.close	
	except:
		pass


        fd = open(homedir+"/.pyvonrc", 'w')
        for (a,b,c) in ini_list:
            fd.write(a)
            fd.write(b+"\n")            
        fd.close




###############################################################
#--                 Exit the application                    --#
###############################################################

def quitter(event=None):
	save_on_exit()
	sys.exit()



###############################################################
#--                 When a key is pressed                   --#
#--                Shows the line & column                  --#
###############################################################


def keyPress(event):
        global Ligne, column
        Ligne, column = text.index(INSERT).split('.', 1)
        status_line.config(text="%s" % (Ligne))
        status_column.config(text="%s" % (column))
	if event.keysym_num not in forbidden_keys and ini_list[23][1]=="Yes":
			auto_expand_word()
        fontify()




###############################################################
#--              When a mouse-key is pressed                --#
#--                Shows the line & column                  --#
###############################################################

def buttonPress(event):
        Ligne, column = text.index(INSERT).split('.', 1)
        status_line.config(text="%s" % (Ligne))
        status_column.config(text="%s" % (column))








###############################################################	
#--                       Message                           --#
###############################################################	


def loadinierror():
	global tl
        tl=Toplevel()
        tl.title("Warning !")
	tl.resizable(width=0,height=0)
        warninglb=Label(tl, font=('arial', 12, 'bold')
                             ,text= Trans_Where_POV ,fg="blue")
        ok=Button(tl, text="Close", fg="red", command=destroywarning)
        ok.pack(side=BOTTOM,padx=4,pady=10)
        root.iconify()
        warninglb.pack()



def destroywarning():
	global tl
        tl.destroy()
        root.deiconify()





###############################################################	
#--                    Insertions menu                      --#
###############################################################	

def importedir(dirname,menu):
    listeobjets = os.listdir(dirname)
    for x in listeobjets :
        longname = os.path.join(dirname,x) # 
        menu1 = Menu(menu,tearoff=0)
        if os.path.isdir(longname) :
            menu.add_cascade(label=x,menu=menu1)
            importedir(longname,menu1)
        if os.path.isfile(longname) :
            if x[-3:] == 'txt' or  x[-3:] == 'pov' :
                menu.add_command(label=x[:-4],command= lambda nom=longname : insert_example(nom))
	    else:
                menu.add_command(label=x,command= lambda nom=longname : insert_example(nom))


def insert_example(nom):
        scene_contents=''
	text.insert('insert','\n')
        fd=open(nom)
        for line in fd.readlines():
            if line.endswith('\r\n'):
                line = line[:-2] + "\n"
            scene_contents=scene_contents+line
        tag_it(scene_contents)
        global c_sel,l_sel
        l_sel,c_sel = text.index(INSERT).split('.', 1)#3 4
        paste=1
        resyntax=0
        colorize(scene_contents,taglist,paste,resyntax)
        fd.close()
	text.insert('insert','\n')

###############################################################	
#--                     No help viewer                      --#
###############################################################	


def no_help_options():
	No_help_msg=Pmw.MessageDialog(title= Trans_Warning ,message_text = Trans_Where_Help )




###############################################################
#--                 Choose help viewer                      --#
###############################################################


def choose_viewer():
	global help_location_entry
	global choose_viewer
	global hlp_view
	global LANG
	hlp_view=StringVar()
	hlp_view.set(ini_list[2][1])
	choose_viewer=Toplevel(relief=SUNKEN)
	choose_viewer.resizable(width=0,height=0)
	choose_viewer.geometry(newGeometry="+100+150")
	choose_viewer.focus_force()
	choose_viewer.title(Trans_Choose_HLP_Viewer)
	Label(choose_viewer,text='').grid(row=0,column=0)
	Radiobutton(choose_viewer,text="Galeon",variable=hlp_view,  # galeon
      				value="galeon",
				command=lambda viewer="galeon":select_viewer(viewer)).grid(row=1,column=0,sticky=W)
	Radiobutton(choose_viewer,text="Mozilla",variable=hlp_view,  # mozilla
      				value="mozilla",
				command=lambda viewer="mozilla":select_viewer(viewer)).grid(row=2,column=0,sticky=W)
	Radiobutton(choose_viewer,text="Konqueror",variable=hlp_view, # konqueror
      				value="konqueror",
				command=lambda viewer="konqueror":select_viewer(viewer)).grid(row=3,column=0,sticky=W)
	Radiobutton(choose_viewer,text="Netscape",variable=hlp_view,  # netscape
      				value="netscape",
				command=lambda viewer="netscape":select_viewer(viewer)).grid(row=4,column=0,sticky=W)
	Radiobutton(choose_viewer,text="Ghostview",variable=hlp_view,  # gv
				value="gv",
				command=lambda viewer="gv":select_viewer(viewer)).grid(row=5,column=0,sticky=W)
	Radiobutton(choose_viewer,text="Xpdf",variable=hlp_view,  # xpdf
      				value="xpdf",
				command=lambda viewer="xpdf":select_viewer(viewer)).grid(row=6,column=0,sticky=W)
	help_location_label=Label(choose_viewer,text= Trans_HLP_PATH ).grid(row=7,column=0,sticky=W)
	help_location_entry=Entry(choose_viewer,width=60)
	help_location_entry.insert(INSERT,ini_list[1][1])
	help_location_entry.grid(row=8,column=0,sticky=W,padx=2,pady=2)
	Button(choose_viewer,text= Trans_Close ,command=exit_help_viewer).grid(row=10,column=1,padx=4,pady=4)
	Button(choose_viewer,text= Trans_Search_HLP_LOC ,command=search_help_location).grid(row=8,column=1,padx=2,pady=2)
	Label(choose_viewer,text='').grid(row=9,column=0)
	Button(choose_viewer,foreground="red",text= Trans_View_POV_HLP ,command=show_help).grid(row=10,column=0,padx=4,pady=4)

	LANG=StringVar()
	LANG.set(ini_list[9][1])
	F_LANG=Frame(choose_viewer,relief=RIDGE,borderwidth=2)
	Label(F_LANG,text= Trans_Choose_Lang ).grid(row=1,column=0)
	global Choose_Lang

	Choose_Lang=Pmw.ComboBox(F_LANG,
		scrolledlist_items=[Trans_FRENCH,Trans_ENGLISH,Trans_GERMAN,Trans_ITALIAN,Trans_DUTCH,Trans_SPANISH,Trans_POLISH],
		entry_width=14,
		listheight=150,
		selectioncommand=select_LANG)
      
	if ini_list[9][1]=="francais":Choose_Lang.selectitem(0,1)
	elif ini_list[9][1]=="english":Choose_Lang.selectitem(1,1)
	elif ini_list[9][1]=="deutsch":Choose_Lang.selectitem(2,1)      
	elif ini_list[9][1]=="italiano":Choose_Lang.selectitem(3,1)
	elif ini_list[9][1]=="nederlands":Choose_Lang.selectitem(4,1)
	elif ini_list[9][1]=="espanol":Choose_Lang.selectitem(5,1)
	elif ini_list[9][1]=="polski":Choose_Lang.selectitem(6,1)      
	else :Choose_Lang.selectitem(1,1)
#      if L.selectitem(0,1)
	Choose_Lang.grid(row=2,column=0,sticky=W)
	F_LANG.grid(rowspan=2,row=1,column=1)





def select_viewer(viewer):
	ini_list[2][1]=viewer


def select_LANG(Language):
	global Choose_Lang
	t=Choose_Lang.curselection()[0]
	if t=="0":ini_list[9][1]="francais"
	elif t=="1":ini_list[9][1]="english"
	elif t=="2":ini_list[9][1]="deutsch"
	elif t=="3":ini_list[9][1]="italiano"
	elif t=="4":ini_list[9][1]="nederlands"
	elif t=="5":ini_list[9][1]="espanol"
	elif t=="6":ini_list[9][1]="polski"
	else : ini_list[9][1]="english"

def exit_help_viewer():
	ini_list[1][1]=help_location_entry.get()
	choose_viewer.destroy()


def search_help_location():
	import tkFileDialog
	global help_location_entry
	global choose_viewer
	filetypes = (( "index","index.html"), ( "Intro","intro.html"), ("Pdf files","*.pdf") , ("Poscript files","*.ps") , ( "All","*"))
	ini_list[1][1]=tkFileDialog.askopenfilename(filetypes=filetypes)
	help_location_entry.delete(0,END)
	help_location_entry.insert(INSERT,ini_list[1][1])
	choose_viewer.focus_force()


def show_help(event=None):
	global choose_viewer
	global help_location_entry      
	if ini_list[2][1]=="" or ini_list[1][1]=="":no_help_options()	
	try:
		ini_list[1][1]=help_location_entry.get()
		choose_viewer.destroy()
	except:
		pass	
	cmd_hlp=ini_list[2][1]+" "+ini_list[1][1]+" &"
	h_in,h_out=os.popen4(cmd_hlp)




###################################################
#--        right-click on text                  --#
###################################################


def right_click_text(event=None):
	menuZ.tk_popup(event.x_root, event.y_root)



menuZ = Menu(root,font=police)
menuZ.option_add('*font',police)
menuZ.add_command(label=Trans_Cut, command=cut)
menuZ.add_command(label=Trans_Copy, command=copy)
menuZ.add_command(label=Trans_Paste, command=paste)
menuZ.add_command(label= Trans_Delete_line , command=delete_line)
menuZ.add_command(label= Trans_Delete_word , command=delete_word)
menuZ['tearoff'] = 0
menu_insert_popup = Menu(menuZ,font=police)
menu_insert_popup['tearoff'] = 0
menuZ.add_separator()
menuZ.add_cascade(label = Trans_Insert,menu=menu_insert_popup,font=police)


listeobjets = os.listdir(InsertDir)
for x in listeobjets :
	longname = os.path.join(InsertDir,x)
	menu1 = Menu(menu_insert_popup,tearoff=0,font=police)
	if os.path.isdir(longname) :
		menu_insert_popup.add_cascade(label=x,menu=menu1)
		importedir(longname,menu1)
	if os.path.isfile(longname) :
		if x[-3:] == 'txt' or  x[-3:] == 'pov' :
			menu1.add_command(label=x[:-4],command= lambda nom=longname : insert_example(nom),font=police)
		else:
			menu1.add_command(label=x,command= lambda nom=longname : insert_example(nom),font=police)



###################################################
#--           Pyvonrc loading                   --#
###################################################


try:
	loadini()
except IOError:
	loadinierror
except IndexError:
	loadinierror
store_scene_untitled()




###################################################
#--            GUI starts here                  --#
###################################################

mframe = Frame(root)
mframe.option_add('*font',police)
mframe.pack(expand=1, fill=BOTH)
mframe.focus_force()
frame_input=Frame(mframe)
try:
	frame_input.option_add('*font',police)
except:
	pass




###################################################
#--                    Menu                     --#
###################################################

menubar=Menu(frame_input)
filemenu=Menu(menubar,tearoff=0)
try:
	filemenu.add_command(label= Trans_New ,command=clr)
except:
	from menu_english import *
	ini_list[9][1]="english"
	LANG="english"
	filemenu.add_command(label= Trans_New ,command=clr)
filemenu.add_command(label= Trans_Open,command=open_file)
filemenu.add_command(label= Trans_Close,command=close_checking)
filemenu.add_command(label= Trans_Save,command=save)
filemenu.add_command(label= Trans_SaveAs,command=saveas)
filemenu.add_separator()
filemenu.add_command(label= Trans_Quit,command=quitter)
filemenu.add_separator()
try:
	homedir = os.path.expanduser('~')
	fd=open(homedir+'/.pyoldfiles','r')
	for line in fd.readlines():
		if line.endswith('\n'):
			line=line[:-1]
		filemenu.add_command(label=os.path.split(line)[1],command=lambda oldfile=line:open_old_file(oldfile))

		Oldfiles.append(line)
	fd.close()

except:
	pass

menubar.add_cascade(label=Trans_File,menu=filemenu)

editmenu=Menu(menubar,tearoff=0)

editmenu.add_command(label=Trans_Undo,command=undo,state="disabled")
editmenu.add_separator()
editmenu.add_command(label=Trans_Cut,command=cut)
editmenu.add_command(label=Trans_Copy, command=copy)
editmenu.add_command(label=Trans_Paste,command=paste)
editmenu.add_command(label=Trans_SelAll,command=selectall)
editmenu.add_separator()
editmenu.add_command(label= Trans_Delete_line , command=delete_line)
editmenu.add_command(label= Trans_Delete_word , command=delete_word)
editmenu.add_separator()
editmenu.add_command(label= Trans_Replace_Menu,command=box_replace_word)



menubar.add_cascade(label=Trans_Edit,menu=editmenu)
insertmenu=Menu(menubar,tearoff=0)        

importedir(InsertDir,insertmenu)

menubar.add_cascade(label=Trans_Insert,menu=insertmenu)


auto_save=StringVar()
auto_save.set(ini_list[7][1])

show_help_bubbles=StringVar()
show_help_bubbles.set(ini_list[22][1])



if ini_list[23][1]=="":ini_list[23][1]="yes"

auto_expand_keyword=StringVar()
auto_expand_keyword.set(ini_list[23][1])


#	ini_list[22][1]="Yes"	
#	show_help_bubbles.set(ini_list[22][1])

optionsmenu=Menu(menubar,tearoff=0)
optionsmenu.add_command(label=Trans_Exe_Opt, command=show_location)
options_auto_save=Menu(optionsmenu,tearoff=0)
options_help_bubbles=Menu(optionsmenu,tearoff=0)
options_auto_expand=Menu(optionsmenu,tearoff=0)
optionsmenu.add_cascade(label=Trans_Auto_Save, menu=options_auto_save)


ASY=options_auto_save.add_radiobutton(label=Trans_Yes,variable=auto_save,
                                             value= "Yes",
                                             command=lambda s=1:save_b_render(s))
ASN=options_auto_save.add_radiobutton(label=Trans_No,variable=auto_save,
                                             value= "No" ,
                                             command=lambda s=0:save_b_render(s))

optionsmenu.add_cascade(label= Trans_Display_HLP_B , menu=options_help_bubbles)
ASY=options_help_bubbles.add_radiobutton(label=Trans_Yes,variable=show_help_bubbles,
                                             value= "Yes",
                                             command=lambda s=1:show_hlp_b(s))
ASN=options_help_bubbles.add_radiobutton(label=Trans_No,variable=show_help_bubbles,
                                             value= "No" ,
                                             command=lambda s=0:show_hlp_b(s))




optionsmenu.add_command(label= Trans_Expand_Keyword ,command=expand_word)

optionsmenu.add_cascade(label= Trans_Expand_Kd , menu=options_auto_expand,command=render)

AEKW_Y=options_auto_expand.add_radiobutton(label=Trans_Yes,variable=auto_expand_keyword,
                                             value= "Yes",
                                             command=lambda s=1:auto_exkw(s))


AEKW_N=options_auto_expand.add_radiobutton(label=Trans_No,variable=auto_expand_keyword,
                                             value= "No" ,
                                             command=lambda s=0:auto_exkw(s))

optionsmenu.add_command(label= Trans_Font ,command=change_color)
optionsmenu.add_command(label= Trans_HLP_Lang_Menu ,command=choose_viewer)

menubar.add_cascade(label=Trans_Options,menu=optionsmenu)

rendermenu=Menu(menubar,tearoff=0)
rendermenu.add_command(label=Trans_Render,command=render)
rendermenu.add_command(label=Trans_Pause,command=pause)
rendermenu.add_separator()
rendermenu.add_command(label= 'Queue rendering',command=pause,state="disabled")
menubar.add_cascade(label=Trans_Render_Menu,menu=rendermenu)

helpmenu=Menu(menubar,tearoff=0)
helpmenu.add_command(label= Trans_About ,command=about )
helpmenu.add_command(label= Trans_Show_Help ,command=show_help )
helpmenu.add_command(label= Trans_KB_SC , command=Keys_help )
menubar.add_cascade(label= Trans_Help ,menu=helpmenu )
root.config(menu=menubar)



###################################################
#--                  Buttons                    --#
###################################################
balloon = Pmw.Balloon(root)


global frame_buttons
frame_buttons=Frame(mframe)
button_new = Button(frame_buttons, image=nouveau, width=37, relief=FLAT,command=clr)
if ini_list[22][1]=="Yes":
	balloon.bind(button_new, Trans_New)
button_new.pack(side=LEFT)
button_open = Button(frame_buttons, image=ouvrir, width=37, relief=FLAT,command=open_file)
if ini_list[22][1]=="Yes":
	balloon.bind(button_open, Trans_Open)
button_open.pack(side=LEFT)
button_save = Button(frame_buttons, image=sauver, width=37, relief=FLAT,command=save)
if ini_list[22][1]=="Yes":
	balloon.bind(button_save, Trans_Save)
button_save.pack(side=LEFT)
button_close = Button(frame_buttons, image=closee, width=37, relief=FLAT,command=close_checking)
if ini_list[22][1]=="Yes":
	balloon.bind(button_close, Trans_Close)
button_close.pack(side=LEFT)
button_render = Button(frame_buttons, image=lancer, width=37, relief=FLAT,command=render)
if ini_list[22][1]=="Yes":
	balloon.bind(button_render, Trans_Render)
button_render.pack(side=LEFT)
button_pause = Button(frame_buttons, image=pausee,relief=FLAT, width=37, command=pause)
if ini_list[22][1]=="Yes":
	balloon.bind(button_pause, Trans_Pause)
button_pause.pack(side=LEFT)



###################################################
#--               File selector                 --#
###################################################



F_file_selector=Frame(frame_buttons,borderwidth=2,relief=FLAT)

Label(F_file_selector,text= Trans_File_Selector ).pack()
drop_down_file=Pmw.ComboBox(F_file_selector,
	scrolledlist_items=Files_list,
	entry_width=14,
	selectioncommand=drpdwn_choose_file)

if ini_list[22][1]=="Yes":
	balloon.bind(F_file_selector, Trans_File_Select_B )

drop_down_file.setlist(Files_list)
drop_down_file.selectitem(0,1)
drop_down_file.pack(side=LEFT,padx=2)
F_file_selector.pack(side=LEFT,padx=2)



###################################################
#--         Resolution settings                 --#
###################################################

Resolution_list=[ Trans_Default ,"1280 x 1024","1024 x 768","800 x 600","640 x 480","400 x 300","320 x 240","128 x 96"]
F_Resolutions=Frame(frame_buttons,borderwidth=2,relief=FLAT)
Label(F_Resolutions,text= Trans_Fast_Resol ).pack()
drop_down_resolution=Pmw.ComboBox(F_Resolutions,scrolledlist_items=Resolution_list,
	entry_width=11,
	listheight=165)
if ini_list[22][1]=="Yes":
	balloon.bind(F_Resolutions,  Trans_Resol_Select_B )

drop_down_resolution.setlist(Resolution_list)
drop_down_resolution.selectitem(0,1)
drop_down_resolution.pack()
F_Resolutions.pack(side=LEFT,padx=2)



###################################################
#--                Quality scale                --#
###################################################

global Q_setting
Q_setting=IntVar()
Q_setting.set(9)
F_Quality=Frame(frame_buttons,borderwidth=2,relief=GROOVE)
T_Quality=Label(F_Quality,text= "+Q",anchor=CENTER,height=1)
Quality_scale=Scale(F_Quality,orient='vertical',from_=0,to=9,width=15,variable=Q_setting,length=55,sliderlength=20)

if ini_list[22][1]=="Yes":
	balloon.bind(F_Quality, Trans_Quality_B )

T_Quality.pack(side=LEFT)
Quality_scale.pack(side='left')
F_Quality.pack(side='left',padx=4,pady=0)



###################################################
#--                Fast preview                 --#
###################################################


FP=IntVar()
F_Fast_Preview=Frame(frame_buttons,borderwidth=2,relief=GROOVE)
Preview_Button=Checkbutton(F_Fast_Preview,variable=FP).pack(expand=0,anchor=CENTER,side=BOTTOM)
Label(F_Fast_Preview,text= Trans_Fast_Prev ,width=8).pack(side=BOTTOM)

if ini_list[22][1]=="Yes":
	balloon.bind(F_Fast_Preview, Trans_Preview_B )

F_Fast_Preview.pack(side='left',padx=4,pady=0,fill='y')
frame_buttons.pack(expand=0, fill=BOTH)


###################################################
#--                Partial Render               --#
###################################################


APR=IntVar()
F_Partial_render=Frame(frame_buttons,borderwidth=2,relief=GROOVE,height=4)
Sub_F1=Frame(F_Partial_render)
Sub_F2=Frame(F_Partial_render)
Button_Partial_render=Button(Sub_F1,command=Partial_Render,text= "+sc +ep")
if ini_list[22][1]=="Yes":
	balloon.bind(Button_Partial_render, Trans_Part_Butt_B )
Button_Partial_render.pack(side='left',expand=0)
L_ACTIVATE=Label(Sub_F2,text= Trans_activate ).pack(side='left')
Activate_partial_render=Checkbutton(Sub_F2,variable=APR,anchor=E)
if ini_list[22][1]=="Yes":
	balloon.bind(Activate_partial_render,  Trans_Part_Check_B)
Activate_partial_render.pack(side='right')
Sub_F1.pack(side='top')
Sub_F2.pack(side='bottom')
F_Partial_render.pack(side='left')


###################################################
#--                Status Bar                   --#
###################################################

frame_status_bar=Frame(mframe,height=1,relief=GROOVE,borderwidth=3)
if ini_list[22][1]=="Yes":
	balloon.bind(frame_status_bar, Trans_Parameter_Bar )
parameters = Label(frame_status_bar,text=Trans_Parameters, anchor="w",bg='gray70')
parameters.pack(side='left')
fill_params=Entry(frame_status_bar,bg='gray100')
fill_params.insert(INSERT,ini_list[3][1])
fill_params.pack(side='left', expand=1, fill=BOTH)


status_column = Label(frame_status_bar,padx=1,
                         relief=FLAT,bg="gray90",
                         width=8,
                         anchor="w",
                          text="0"
                         )
status_column.pack(side='right')
status_coltext = Label(frame_status_bar,padx=1,bg='gray70',
                         width=4,
                         anchor="w",
                         text="Cl : ")
status_coltext.pack(side='right')


status_line = Label(frame_status_bar,padx=1,
                         relief=FLAT,bg="gray90",
                         width=8,
                         anchor="w",
                          text="1"
                         )
status_line.pack(side='right')


status_linetext = Label(frame_status_bar,padx=1,bg='gray70',
                         width=4,
                         anchor="w",
                         text="Ln : ")
status_linetext.pack(side='right')


frame_status_bar.pack(side='top',expand=0, fill='x')


###################################################
#--              The editor itself              --#
###################################################


if ini_list[21][1]=="":ini_list[21][1]="40"
text = Text(frame_input,tabs=ini_list[21][1],font=(ini_list[5][1], ini_list[6][1]),
	    wrap=NONE,bg=ini_list[19][1],foreground=ini_list[20][1],insertbackground = ini_list[20][1],selectborderwidth=2)
yscroll_i = Scrollbar(frame_input,orient=VERTICAL, command=text.yview)
xscroll_i = Scrollbar(frame_input,orient=HORIZONTAL, command=text.xview)
yscroll_i.pack(side=RIGHT, fill=Y)
xscroll_i.pack(side=BOTTOM, fill=X)
text['yscrollcommand']=yscroll_i.set
text['xscrollcommand']=xscroll_i.set
text.pack(side=LEFT,expand=1, fill=BOTH, )

frame_input.pack(expand=1, fill=BOTH)





###################################################
#--         The color tags definition           --#
###################################################


text.tag_config('keyword', font=(ini_list[5][1], ini_list[6][1], 'normal' ),foreground=ini_list[10][1])
text.tag_config('keywordbold', font=(ini_list[5][1], ini_list[6][1], 'bold' ),foreground=ini_list[11][1])
text.tag_config('num', font=(ini_list[5][1], ini_list[6][1], 'normal' ),foreground=ini_list[12][1])
text.tag_config('slcomment', font=(ini_list[5][1], ini_list[6][1], 'italic' ),foreground=ini_list[13][1])
text.tag_config('mlcomment', font=(ini_list[5][1], ini_list[6][1], 'italic' ),foreground=ini_list[13][1])
text.tag_config('math', font=(ini_list[5][1], ini_list[6][1], 'normal' ),foreground=ini_list[12][1]) #0000FF
text.tag_config('squig', font=(ini_list[5][1], ini_list[6][1], 'normal' ),foreground=ini_list[18][1])
text.tag_config('dirc', font=(ini_list[5][1], ini_list[6][1], 'normal' ),foreground=ini_list[17][1])
text.tag_config('str', font=(ini_list[5][1], ini_list[6][1], 'normal' ),foreground=ini_list[14][1])
text.tag_config('pattern', font=(ini_list[5][1], ini_list[6][1], 'normal' ),foreground=ini_list[16][1])
text.tag_config('ident', font=(ini_list[5][1], ini_list[6][1], 'normal' ),foreground="#000000")  ###########
text.tag_config('brace', font=(ini_list[5][1], ini_list[6][1], 'normal' ),foreground=ini_list[18][1],background="#FF0000",relief=GROOVE)
text.tag_config('highl',background=ini_list[15][1])




###################################################
#--            Key shortcuts                    --#
###################################################


text.bind("<Button-3>",right_click_text)
root.bind("<KeyRelease>",keyPress)
root.bind("<Return>",auto_indent)
root.bind("<Control-space>",expand_word)
root.bind("<Control-c>",copy)
root.bind("<Control-x>",cut)
root.bind("<Control-v>",paste)
root.bind("<Control-y>",delete_line)
root.bind("<Control-BackSpace>",delete_word)
root.bind("<Control-Home>",text.see('1.0'))
root.bind("<Control-End>",text.see('end'))
root.bind("<Delete>",suppr_mlcomment)
root.bind("<BackSpace>",suppr_mlcomment)
root.bind("<ButtonPress>",buttonPress)
root.bind("<Control-a>",selectall)
root.bind("<Control-f>",box_replace_word)
root.bind("<Control-g>",render)
root.bind("<Control-n>",clr)
root.bind("<Control-o>",open_file)
root.bind("<Control-p>",pause)
root.bind("<Control-q>",quitter)
root.bind("<Control-s>",save)


###################################################
#--                    ROOT                     --#
###################################################


larg = root.winfo_screenwidth()
haut = root.winfo_screenheight()

haut_vert=(haut/4)*2.45
pos_vert=(haut/4)





root.geometry("%dx%d+0+0" % (larg,haut_vert))
root.title("PYVON - POVRAY editor")
root.protocol("WM_DELETE_WINDOW",quitter)
root.mainloop()
