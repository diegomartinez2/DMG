from tkinter import
from tkinter.ttk import
from time import strftime
from numpy import floor #the greatest integer that does not exceed x.
#int(a//1) #the value of x with the fractional portion removed.
def Earth2Dni(arg):
    delta = strftime('%H:%M:%S %p')
    har = floor(Millisec / MillisecPerHar)
   pass

def FIX(arg):
    return nt(arg//1)

def Gregorian_to_Julian(arg):
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
    Z = INT(JD)
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
    Z = (JD - INT(JD)) * 86400
    Hour = FIX(Z / 3600)
    R = Z - (Hour * 3600)
    Minute = FIX(R / 60)
    Second = R - (Minute * 60)
    return Year,Month,Hour,Minute,Second

def Cavernian_to_AYN(arg):
    #whole yahrtee
    WY = Yahr + ((Vailee - 1) * 29) + ((Hahr - 9647) * 290)
    #fractional yahr
    FY = ((Gahrtahvo * 15625) + (Tahvo * 625) + (Gorahn * 25) + Prorahn) / 78125
    #Atrian Yahr
    AY = WY + FY
    return AY

def Cavernian_to_Gregorian(arg):
    #Atrian Yahr
    AYD = AY - 1.0
    #to Julian
    JDD = AYD * 1.25945582758621
    JD = JDD + 727249.704166666
    return Julian_to_Gregorian(JD)

def Gregorian_to_Cavernian(arg):
    JD = Gregorian_to_Julian(arg)
    JDD = JD - 727249.704166666
    AYD = JDD * 0.793993705929756
    AY = AYD + 1.0
    #Convert the calculated Atrian Yahr to a Cavernian date
    pass
