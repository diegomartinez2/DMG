#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  untitled.py
#
#  Copyright 2024 Diego Martinez Gutierrez <diego.martinez@ehu.eus>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
#
# ---------------------------
# Importación de los módulos
# ---------------------------
#from tkinter import
#from tkinter.ttk import
from time import strftime
from numpy import floor #the greatest integer that does not exceed x.
#int(a//1) #the value of x with the fractional portion removed.


# -------
# Clases
# -------

# ----------
# Funciones
# ----------
def Earth2Dni(arg):
    delta = strftime('%H:%M:%S %p')
    har = floor(Millisec / MillisecPerHar)
    pass

def get_time():
    Day = strftime('%d')
    Month = strftime('%m')
    Year = strftime('%Y')
    Hour = strftime('%H')
    Minute = strftime('%M')
    Second = strftime('%S')
    return Day, Month, Year, Hour, Minute, Second

def FIX(arg):
    return int(arg//1)

def Gregorian_to_Julian(Day, Month, Year, Hour=0, Minute=0, Second=0):
    if Year<0:
        Year = -(Year - 1)
    if Month < 3:
        Month = Month + 12
        Year = Year - 1
    #whole days
        WD = Day + FIX(((153 * Month) - 457) / 5) + floor(365.25 * Year) - floor(0.01 * Year) + floor(0.0025 * Year)
    #fractional day
        FD = ((Hour * 3600) + (Minute * 60) + Second) / 86400
    # Julian Day
        JD = WD + FD
    return JD

def Julian_to_Gregorian(JD):
    Z = floor(JD)
    G = Z - 0.25
    A = floor(G / 36524.25)
    B = A - (0.25 * A)
    Year = floor((G + B) / 365.25)
    C = Z + B - floor(365.25 * Year)
    Month = FIX(((5 * C) + 456) / 153)
    Day = C - FIX(((153 * Month) - 457) / 5)
    if Month > 12:
        Year = Year + 1
        Month = Month - 12
    if Year < 1:
        Year = 1 - Year
    Z = (JD - floor(JD)) * 86400
    Hour = FIX(Z / 3600)
    R = Z - (Hour * 3600)
    Minute = FIX(R / 60)
    Second = R - (Minute * 60)
    return Day, Month, Year, Hour, Minute, Second

def Cavernian_to_AYN(Yahr, Vailee, Hahr,Gahrtahvo=0, Tahvo=0, Gorahn=0, Prorahn=0):
    #whole yahrtee
    WY = Yahr + ((Vailee - 1) * 29) + ((Hahr - 9647) * 290)
    #fractional yahr
    FY = ((Gahrtahvo * 15625) + (Tahvo * 625) + (Gorahn * 25) + Prorahn) / 78125
    #Atrian Yahr
    AY = WY + FY
    return AY

def AYN_to_Cavernian(AY):
    #Extract the number of whole yahrtee and calculate the date
    Z = floor(AY)
    G = Z - 0.25
    A = floor(G / 290)
    Hahr = 9647 + A
    C = Z - (A * 290)
    Vailee = floor((C - 0.25) / 29) + 1
    Yahr = C - ((Vailee - 1) * 29)
    #If AY includes a fractional (time of yahr) part, extract the fraction and calculate the time
    Z = (AY - floor(AY)) * 78125
    Gahrtahvo = FIX(Z / 15625)
    R = Z - (Gahrtahvo * 15625)
    Tahvo = FIX(R / 625)
    R = R - (Tahvo * 625)
    Gorahn = FIX(R / 25)
    Prorahn = R - (Gorahn * 25)
    return Yahr, Vailee, Hahr, Gahrtahvo, Tahvo, Gorahn, Prorahn

def Cavernian_to_Gregorian(AY):
    #Atrian Yahr
    AYD = AY - 1.0
    #to Julian
    JDD = AYD * 1.25945582758621
    JD = JDD + 727249.704166666
    return Julian_to_Gregorian(JD)

def Gregorian_to_Cavernian(Day, Month, Year, Hour=0, Minute=0, Second=0):
    JD = Gregorian_to_Julian(arg)
    JDD = JD - 727249.704166666
    AYD = JDD * 0.793993705929756
    AY = AYD + 1.0
    #Convert the calculated Atrian Yahr to a Cavernian date
    return AYN_to_Cavernian(AY)

# def Vailee_to_VaileeTEXT(arg):
#     VaileeDictionary = {
#     1:Leefo,
#     2:Leebro,
#     3:Leesahn,
#     4:Leetar,
#     5:Leevot,
#     6:Leevofo,
#     7:Leevobro,
#     8:Leevosahn,
#     9:Leevotar,
#     10:Leenovoo
#     }
#     pass

#-----------------------------------------------------

def main(args):
    VaileeDictionary = {1:"Leefo",2:"Leebro",3:"Leesahn",4:"Leetar",5:"Leevot",6:"Leevofo",7:"Leevobro",8:"Leevosahn",9:"Leevotar",10:"Leenovoo"}
    Day, Month, Year, Hour, Minute, Second = get_time()
    print("Year:", Year, "Month=", Month, "Day:", Day, "Hour=",Hour,":",Minute,":",Second)
    Yahr, Vailee, Hahr, Gahrtahvo, Tahvo, Gorahn, Prorahn = Gregorian_to_Cavernian(Day, Month, Year, Hour, Minute, Second)
    print("Yahr:",Yahr, "Vailee:",Vailee,"_",VaileeDictionary[Vailee], "Hahr:", Hahr, "Gahrtahvo:", Gahrtahvo, "Tahvo:", Tahvo, "Gorahn:", Gorahn, "Prorahn:", Prorahn)
    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
