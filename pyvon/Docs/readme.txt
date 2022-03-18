24 May 2002

Pyvon is a Graphic User Interface (GUI) for the POV-RAY raytracer.
It uses the graphical interface of Python (www.python.org) and I tried to clone the windows version of POV.



*******************************************
* * * * * * HOW TO START IT * * * * * * * *
*******************************************

To start this program just type ./pyvon in a shell




***********************************************
* * * * * * KEYBOARDS SHORTCUTS * * * * * * * *
***********************************************

CRTL + c		-> copy
CRTL + v		-> paste
CRTL + a		-> select all
CRTL + f		-> find/replace
CRTL + <space>		-> Show a list of keywords that match what you have started typing
CRTL + y		-> delete the current line
CRTL + <backspace>	-> delete the current line
CRTL + n		-> create new file
CRTL + o		-> open a file
CRTL + s		-> save the current file
CRTL + q		-> exit pyvon
CRTL + g		-> starts/stops the render
CRTL + p		-> pause/resume the render

___________________________________________________

The first time you use Pyvon, use the menu /options/Pov-ray location to tell where the pov
executable is. The location of the executable, the font, the font size, the parameters are stored in a hidden file '.pyvonrc' located in your home directory.

Python is a script language ( like Perl, javascript,...). On files larger than 500 Kb, the syntax highlighting is automatically disabled.





CHANGELOG :
d44 - 10/04/2003 - Pyvon 1.3.2
Add a button to enable/disable auto expansion keyword
Fixed bug with moving text when copying text


d42 - 28/01/2003 - Pyvon 1.3
'Beautify' the code : Replaced spaces with tabs inside the code
Much remains to be done
Use the Python optimization starting Pyvon -> It should be somewhat faster



d38 - 20/01/2003 - Pyvon 1.3
Added small tips when moving the mouse over the interface
Added the translations for the tips.
Added the names of the translators in the language files and in the help menu
Show a warning if the +sc or +sr value in the partial render are greater than the +ec or +er values



d37 - 14/01/2003 - Pyvon 1.3
Changed the structure of the Insert directory
When selecting a font, pyvon searches all the available fonts in the system and displays them.
A small label shows what the font looks like if you want to apply it in the editor.
Added an entry in the menu for a pyvon help window for the keyboard shortcuts



d36 - 08/01/2003 - Pyvon 1.2a
Added the possibility to change the color of the identifiers
Changed the way the Insert directory is sorted


d35 - 03/01/2003 - Pyvon 1.2
Automatic keyword : 
	For example, When you type 'hei', Pyvon automatically types the rest of the keyword : height_field.
The CTRL+Space key combination still works :
	If you type 'ab' the CTRL+Space, Pyvon will display a box the words that begin with 'ab' -> abs, absorption



d31 - 19/12/2002 - Pyvon 1.2
Reworked the rendering error section :
	Pyvon switch to the erronous file if it is in the file selector or displays an error message.


d29 - 17/12/2002 - Pyvon 1.1
Improved the multicomments highlight syntax -> it works like WinPov now


d23 - 01/12/2002 - Pyvon 1.0
Added a shortkey combination from emacs : CTRL+Y -> delete all the current line
Added a shortkey combination from emacs : CTRL+BACKSAPCE -> delete all the word behind the cursor
Added a shortkey combination from windows : CTRL+Home
Added a shortkey combination from windows : CTRL+End
Added PMW 1.1


d21 - 28/11/2002 - Pyvon 1.0
Rework all the indentation :
  - As in WinPov, Pyvon remembers the tab of the previous written 	    line
  - Indentation can be changed and is stored in .pyvonrc


d19 - 19/11/2002 - Pyvon 1.0
Added a small warning when the user writes a file where (s)he has not writing permissions


d17 - 19/11/2002 - Pyvon 1.0
Added right-click feature : Cut, copy, paste, insertions


d16 - 18/11/2002 - Pyvon 1.0
Enlightenment of the code. It is a bit faster now.



d13 - 16/11/2002 - Pyvon 1.0
Insert menu : The insert folder can now be changed to whatever you want. Pyvon reads
and includes all the folders and sub-folders found.
Mega pack can be pasted instead for those you wish.



d12 - 15/11/2002 - Pyvon 1.0
Fixed a few bugs ( editor background color)

d11 - 14/11/2002 - Pyvon 1.0
Syntax highlighting colors and font can now be changed interactively in the interface.
The settings are saved in .pyvonrc


d06 - 10/11/2002 - Pyvon 1.0
Pyvon does not have to be started from a directory with writing permissions any more.
The error file is written to the home of the user


d04 - 03/11/2002 - Pyvon 1.0
Fixed the cut/paste bug
Changed the 'Disabled' keyword in resolutions by 'Default'


d03 - 03/11/2002 - Pyvon 1.0
Added a warning message box when closing a file. Some people close instead of saving.


d02 - 20/10/2002 - Pyvon 1.0
Added a partial render button + an activation checkbutton


c99 - 03/10/2002 - Pyvon 1.0
Moved the status line below the buttons

c98 - 01/10/2002 - Pyvon 1.0
Removed the -f command line when no previous file was found

c97 - 16/09/2002 - Pyvon 0.99
Created install script for the tar.gz archive

    
c96 - 12/09/2002 - Pyvon 0.99
Fixed the problem with space in script names
Compiled rpm package


c95 - 10/09/2002 - Pyvon 0.99
Added the Polish language. Pb with fonts


c94 - 09/09/2002 - Pyvon 0.99
Added the GPL2 licence in the archive 


c93 - 09/09/2002 - Pyvon 0.99
Pyvon now renders files with space in their name (thanks to Roz), even if
the normal command line doesn't.
The optimized version is included in the archive.


c92 - 06/09/2002 - Pyvon 0.98
Added Dutch language (thanks to Ingo)
Added Spanish language (thanks to John Coppens) 


c91 - 04/09/2002 - Pyvon 0.97
If an error occurs because of an include that is already open, Pyvon
automatically swtich to that file and highlight the error


c90 - 04/09/2002 - Pyvon 0.97
If an error occurs because of an include file that is not open, Pyvon
automatically opens that file and highlight the error 


c89 - 03/09/2002 - Pyvon 0.96
Pyvon Will never crash again because .pyvonrc or pyvon.ini are not found
I included default parameters in case those file miss


c88 - 03/09/2002 - Pyvon 0.96
Fixed minor bug with font
Added Italian language (Thanks to Alfonso Martone)
Added German language (Thanks to Bonsai)


c87 - 30/09/2002 - Pyvon 0.95
Fixed small bug with word expansion
Added the possibility to choose between French-English (more to be added)


c85 - 29/08/2002 - Pyvon 0.92
Optimized the code. It is now about 10 % faster


c82 - 29/08/2002 - Pyvon 0.92
Completely changed the About menu (now with working link)


c81 - 21/08/2002 - Pyvon 0.92
Rebuild the help menu


c80 - 28/08/2002 - Pyvon 0.92
Fix a bug about the insertions. Now a line is added at the beginning and at the end


c79 - 28/08/2002 - Pyvon 0.92
Added the possiblity for multilanguage support (menu_eng.py, menu_fra.py)


c78 - 27/08/2002 - Pyvon 0.92
Added comments for easy reading of source


c77 - 27/08/2002 - Pyvon 0.92
Rebuild the povray executable options 

c76 - 27/08/2002 - Pyvon 0.92
Rebuild the font menu 

c75 - 26/08/2002 - Pyvon 0.92
Added a few shortcut keys (Ctrl+o : Open, Ctrl+q : Quit,....)


c72 - 26/08/2002 - Pyvon 0.92
Changed the interface font


c70 - 25/08/2002 - Pyvon 0.92
Added file/replace/replace all menu 


C67 - 22/08/2002 - Pyvon 0.92
Added file history in the file menu ( up to 7 files )


C62->c65 - 22/08/2002 - Pyvon 0.92
Fixed minor bugs


C61 - 21/08/2002 - Pyvon 0.92
Syntax coloring is by default turned on, but when opening files larger than 500 Kb, 
it is automatically turned off.
Line syntax coloring (where the cursor is located) is still effective.s



C60 - 20/08/2002 - Pyvon 0.92
pyvon.ini is now replaced by pyvonrc which is located in your home directory.


C58 - 20/08/2002 - Pyvon 0.92
Fixed a the slowdown when matching brace inside a very large file


C57 - 19/08/2002 - Pyvon 0.92
Improved the opening file speed by ~ 100. Especially useful when opening large file.
It now takes 25 seconds to open, tag and colorize a 1 Mb large file on a PIII 1 Mghz.
(compared to minutes with a 200 Ko file)


C56 - 18/08/2002 - Pyvon 0.92
Added a scale for the quality setting


C55 - 18/08/2002 - Pyvon 0.91a
Fixed stupid bug in the res.py 

C54 - 09/08/2002 - Pyvon 0.91
Checks for writing permissions in the current working directory
If not, it exit back to the OS.
Going on holidays for a week

C53 - 09/08/2002
The render and pause are changed dynamically in the menu bar according  
(as in WINPOV 3.5) to the state of the renderer


C52 - 08/08/2002
Pyvon checks whether the file found in Pyvon.ini exists on the hard drive
If not, it displays an warning message

C51 - 07/08/2002
 - Pyvon 0.9
Added syntax highlighting support 'while typing' for multiline comments

C46 - 06/08/2002
Added an 'About' box
Fixed the bug in the help viewer selector

C42 - 04/08/2002
Added syntax highlighting support 'while typing' for the single line comments

C40 - 03/08/2002
Corrected the size of the xterm console when the display is 1024x768
or 800x600.
Not tested for the other resolutions

C38 - 02/08/2002
Added message dialog if the same file is opened more than once

C37 - 01/08/2002
Disabled the undo feature. For the time being, it does more harm
than good

C36 - 27/07/2002
Fixed the evil bug that caused loss of data when switching from one file to another 
with the file selector

C32 - 24/07/2002
Added highlight feature in the editor when an error occurs in the script
(it works as in WINPOV 3.5)

c30 - 23/07/2002
Added the CPU priority (from CPU HOG to IDLE) in the POV executable options


c23 - 21/07/2002
Added the pause-resume icons
The render-stop and pause-resume buttons now switches from one the another

C19 - 19/07/2002
Added the new POV 3.5 windows icons

C18 - 18/07/2002
I removed the output window at the bottom. It was too much hassle for too
little result.

c14 - I am going abroad on holidays for a week -> there won't be much update.


c13 - 03/07/2002
Pyvon now 'communicates' with POV : The line being raytraced is shown at the
bottom left of Pyvon and updated every second.
BUG : For the 'slow-to-render-scenes', the line number is shown in the output.
I'll fiw that later.
Added the POVRAY 3.5 keywords highlighting and insertions.


c10 - 01/07/2002
The internal is 60% rewritten and lightened (it should be faster now)


b39 - 20/06/2002
Added the expand keyword feature:
Ctrl+space or right-click auto-expands a keyword if only one matches.
If more than one matches, a list is displayed with all the possible matches. Then use the up and down 
arrows to select (with the Return key) the wanted keyword.
If no word matches, the screen blinks for 1/10 second.


b37 - 18/06/2002
Cleaned up the code, and splitting it into several files (to improve speed)


b35 - 11/06/2002
Added help menu :
- The user have the choice between several browsers
- The path to the index.html, intro.html, pdf or ps files has to be defined first
- The choices are stored in the pyvon.ini file


b29 - 07/06/2002
Added combo-box for the font choice and font size


b28 - 05/06/2002
The autosave and syntax coloring on loading are now stored in Pyvon.ini


b27 - 04/06/2002
Added close feature


b23 - 01/06/2002
Added multifile opening feature


b22 - 31/05/2002
Added the Pwm module
Removed the file name at the statusbar and added a drop down-box to choose which scene to render
Added syntax highlighting when switching from one scene to another 


b20 - 30/05/2002
Added a fixed font to the stream out window


b18 - 29/05/2002
Pyvon.ini is updated when exiting with the quit menu or the
"kill" button at the top right of pyvon


b16 - 29/05/2002
Now keeps track of the working directory when exiting
(updates pyvon.ini)


b15 - 28/05/2002 -
Added brace matching function
Added the scene menu ( still disabled)


b13 - 26/05/2002
Added the insert menu with syntax coloring
Syntax coloring :
    Added the force coloring menu
    Improved and cleaned the code


b08 - 25/05/2002
Added the current directory to the path


b05
Changed the font for the stream out



Please send me comments, flames, whatever you want.


Regards,

Fabien HENON
