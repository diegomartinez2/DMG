#!/bin/bash
# Date : (2019-12-19)
# Last revision : see the changelog below
# Wine version used : see the changelog below
# Distribution used to test : XUbuntu 19.04 x64
# Author : Dadu042
# Licence : Retail
# Only For : http://www.playonlinux.com
#
# TESTED Editions: local 7Z archive.
#
# Middlewares used by this software : OpenBOR.
#
# CHANGELOG
# [Dadu042] (2019-12-19)
#   First script.
# [Dadu042] (2019-12-28)
#   POL_RequiredVersion 4.2.12 -> 4.0.0
#
# KNOWN ISSUES:
#
# Ideas to improve this script: select archive, then decide if extension is RAR or ZIP or 7Z...

[ "$PLAYONLINUX" = "" ] && exit 0
source "$PLAYONLINUX/lib/sources"

TITLE="Golden Axe Myth"
PREFIX="Golden_Axe_Myth"
EDITOR="Collaborative"
GAME_URL="https://archive.org/details/Golden_Axe_Myth"
AUTHOR="Dadu042"
STEAM_ID=""
WORKING_WINE_VERSION="3.0.3"
GAME_VMS="256"
SHORTCUT_FILENAME="OpenBOR.exe"
SOFTWARE_CATEGORIES="Game;Arcade;"
# http://wiki.playonlinux.com/index.php/Scripting_-_Chapter_9:_Standardization#Advanced_Standardization

# Starting the script
POL_SetupWindow_Init

# Starting debugging API
POL_Debug_Init

# Open dialogue box
POL_SetupWindow_presentation "$TITLE" "$EDITOR" "$GAME_URL" "$AUTHOR" "$PREFIX"

# POL_SetupWindow_message "$(eval_gettext 'WARNING: this software does exist in Linux native version.\n\nThis script only allow to run the Windows version on Linux, please prefer the Linux edition for better 3D speed.')" "$TITLE"

# POL_SetupWindow_message "$(eval_gettext 'This game requires a fast 3D GPU (FYI: a iGPU Intel HD Graphics 530 does only become comfortable with low details).')" "$TITLE"

POL_RequiredVersion "4.0.0" || POL_Debug_Fatal "$APPLICATION_TITLE $VERSION is required to install $TITLE"

# Setting prefix path
POL_Wine_SelectPrefix "$PREFIX"

# Determine Architecture
POL_System_SetArch "amd64"
# POL_System_SetArch "x86"

# Downloading wine if necessary and creating prefix
POL_Wine_PrefixCreate "$WORKING_WINE_VERSION"

Set_OS "winxp"

# Installing mandatory dependencies
# POL_Call POL_Install_corefonts
# POL_Call POL_Instatl_directmusic
# POL_Call POL_Install_d3dx9_43
# POL_Call POL_Install_msxml4
# POL_Call POL_Install_riched30
# POL_Call POL_Install_phzysx
# POL_Call POL_Install_corefonts
# POL_Call POL_Install_d3dx11
# POL_Call POL_Install_mono210

## Sound problem fix - pulseaudio related
# [ "$POL_OS" = "Linux" ] && Set_SoundDriver "alsa"
# [ "$POL_OS" = "Linux" ] && Set_SoundEmulDriver "Y"
## End Fix

# Choose between Steam and other Digital Download versions
# POL_SetupWindow_InstallMethod "STEAM,DVD,LOCAL,DOWNLOAD"
POL_SetupWindow_InstallMethod "LOCAL,DOWNLOAD"

# POL_SetupWindow_message "Note: at the end of the installation, please do not run the game, and do not install DirectX 9." "$TITLE"

# Begin game installation
if [ "$INSTALL_METHOD" == "STEAM" ]; then
        POL_Call POL_Install_steam
        # Mandatory pre-install fix for steam
        POL_Call POL_Install_steam_flags "$STEAM_ID"
        # Shortcut done before install for steam version
        POL_Shortcut "steam.exe" "$TITLE" "" "steam://rungameid/$STEAM_ID"
        # Steam install
        POL_SetupWindow_message "$(eval_gettext 'When $TITLE download by Steam is finished,\nDo NOT click on Play.\n\nClose COMPLETELY the Steam interface, \nso that the installation script can continue')" "$TITLE"
        cd "$WINEPREFIX/drive_c/$PROGRAMFILES/Steam"
        POL_Wine start /unix "steam.exe" steam://install/$STEAM_ID
        POL_Wine_WaitExit "$TITLE"

elif [ "$INSTALL_METHOD" == "DVD" ]; then
        POL_SetupWindow_cdrom

        POL_SetupWindow_check_cdrom "SETUP.EXE"
        POL_Wine start /unix "$CDROM/SETUP.EXE"
        POL_Wine_WaitExit "SETUP.EXE"

        # Restore screen resolution (game's default is 800x600 ?)
        # POL_Shortcut_InsertBeforeWine "$SHORTCUT" "trap 'xrandr -s 0' EXIT"

        POL_Shortcut "$SHORTCUT_FILENAME" "$TITLE" "" "" "$SOFTWARE_CATEGORIES"

elif [ "$INSTALL_METHOD" == "DOWNLOAD" ]; then
        # POL_SetupWindow_message "$(eval_gettext '\nNote: this script will download the DEMO.')" "$TITLE"

        cd "$WINEPREFIX/drive_c"
        POL_Download "https://archive.org/download/Golden_Axe_Myth/Golden_Axe_Myth.7z"

        POL_SetupWindow_wait_next_signal "$(eval_gettext 'Extracting the archive...')" "$TITLE"
        # POL_System_unzip "Golden_Axe_Myth.7z" -d "$WINEPREFIX/drive_c/"

        if [ $(which 7z) ];then
                # REF: https://tecadmin.net/how-to-extract-filename-extension-in-shell-script/
                # filename=$(basename "$APP_ANSWER")
                # cp "$APP_ANSWER" "$WINEPREFIX/drive_c/"
                # cd "$WINEPREFIX/drive_c/"

                # NOK 7z x "$WINEPREFIX/drive_c/Golden_Axe_Myth.7z" -o"$WINEPREFIX/drive_c/game/"
                7z x "Golden_Axe_Myth.7z" -o"$WINEPREFIX/drive_c/game/"
        else
                POL_SetupWindow_message "$(eval_gettext 'Game Extract failed (Missing 7z).')" "$TITLE"
        fi

        # POL_Wine start /unix "setup.exe" # "/SILENT"
        # POL_Wine_WaitExit "$TITLE"

        # POL_Wine "motoracing2_setup.exe" # "/SILENT"
        # POL_Wine_WaitBefore "$TITLE"

        # Restore screen resolution (game's default is 1024x768)
        # POL_Shortcut_InsertBeforeWine "$SHORTCUT" "trap 'xrandr -s 0' EXIT"

        POL_Shortcut "$SHORTCUT_FILENAME" "$TITLE" "" "" "$SOFTWARE_CATEGORIES"
        POL_Shortcut_Document "$TITLE" "GAMythUser_Manual.pdf"

elif [ "$INSTALL_METHOD" == "LOCAL" ]; then
        # POL_SetupWindow_menu "$(eval_gettext 'What is the type of the file?.')" "$TITLE" "$(eval_gettext '.EXE')~$(eval_gettext '.ZIP')~$(eval_gettext '.RAR')" "~"
        POL_SetupWindow_menu "$(eval_gettext 'What is the type of the file?.')" "$TITLE" "$(eval_gettext '.7Z')~$(eval_gettext '.ZIP')" "~"
        # APP_ANSWER=".EXE"

if [ "$APP_ANSWER" == ".EXE" ]; then
        # Asking then installing local files of the game
        cd "$HOME"
        POL_SetupWindow_browse "$(eval_gettext 'Please select the setup file to run')" "$TITLE"
        # POL_SetupWindow_message "$(eval_gettext 'Note: we recommend you to uncheck all the checkboxes:\n[x] -> [ ]')" "$TITLE"
        SETUP_EXE="$APP_ANSWER"
        POL_Wine start /unix "$SETUP_EXE" # "/SILENT"
        POL_Wine_WaitExit "$TITLE"

        # Restore screen resolution (game's default is 1024x768)
        # POL_Shortcut_InsertBeforeWine "$SHORTCUT" "trap 'xrandr -s 0' EXIT"

        POL_Shortcut "$SHORTCUT_FILENAME" "$TITLE" "" "" "$SOFTWARE_CATEGORIES"
        # POL_Shortcut_Document "$TITLE" "Readme.txt"

elif [ "$APP_ANSWER" == "$(eval_gettext '.ZIP')" ]; then
        cd "$HOME"
        POL_SetupWindow_browse "$(eval_gettext 'Please select the .ZIP file')" "$TITLE"
        SETUP_EXE="$APP_ANSWER"
        cd "$POL_System_TmpDir"
        POL_SetupWindow_wait_next_signal "$(eval_gettext 'Extracting the archive...')" "$TITLE"
        POL_System_unzip "$APP_ANSWER" -d "$WINEPREFIX/drive_c/"
        POL_Shortcut "$SHORTCUT_FILENAME" "$TITLE" "" "" "$SOFTWARE_CATEGORIES"

        POL_Shortcut_Document "$TITLE" "GAMythUser_Manual.pdf"

elif [ "$APP_ANSWER" == "$(eval_gettext '.7Z')" ]; then
        cd "$HOME"
        POL_SetupWindow_browse "$(eval_gettext 'Please select the .7Z file')" "$TITLE"
        SETUP_EXE="$APP_ANSWER"
        cd "$POL_System_TmpDir"
        POL_SetupWindow_wait_next_signal "$(eval_gettext 'Extracting the archive...')" "$TITLE"
        # POL_System_unzip "$APP_ANSWER" -d "$WINEPREFIX/drive_c/"

#        7z x "$APP_ANSWER"  "$WINEPREFIX/drive_c/1/"


if [ $(which 7z) ];then

        # POL_SetupWindow_message "$(eval_gettext '7z for:  $APP_ANSWER')" "$TITLE"

        # REF: https://tecadmin.net/how-to-extract-filename-extension-in-shell-script/
        filename=$(basename "$APP_ANSWER")
        cp "$APP_ANSWER" "$WINEPREFIX/drive_c/"
        cd "$WINEPREFIX/drive_c/"
        7z x "$WINEPREFIX/drive_c/$filename" -o"$WINEPREFIX/drive_c/game/"

        # debug
        # POL_SetupWindow_message "$(eval_gettext 'output: "$WINEPREFIX/drive_c/OpenBOR/" ')" "$TITLE"

        # 7z x "APP_ANSWER" &> /dev/null
        # 7z x "$WINEPREFIX/drive_c/"$filename -o"$WINEPREFIX/drive_c/OpenBOR/"

else

        POL_SetupWindow_message "$(eval_gettext 'Game Extract failed (Missing 7z).')" "$TITLE"
fi


        # POL_System_unrar x "$APP_ANSWER" "$WINEPREFIX/drive_c/" || POL_Debug_Fatal "unrar is required to unarchive $TITLE (unrar package is not installed on the OS)."
        POL_Shortcut "$SHORTCUT_FILENAME" "$TITLE" "" "" "$SOFTWARE_CATEGORIES"
        POL_Shortcut_Document "$TITLE" "GAMythUser_Manual.pdf"
fi
fi

################
#      GPU     #
################

# Set Graphic Card information keys for wine
# POL_Wine_SetVideoDriver

# Asking about memory size of graphic card
# POL_SetupWindow_VMS $GAME_VMS

# Useful for Nvidia GPUs
# POL_Call POL_Install_physx

# Empty (no text):
# POL_SetupWindow_message "$LNG_SUCCES" "$TITLE"

# POL_SetupWindow_message "$(eval_gettext '\nInstallation is finished ! :)')" "$TITLE"

# POL_SetupWindow_message "$(eval_gettext 'WARNING: to avoid to have a big useless POL/POM log file, you should type \ninto Debug flags : fixme-all')" "$TITLE"

POL_SetupWindow_Close
exit 0
