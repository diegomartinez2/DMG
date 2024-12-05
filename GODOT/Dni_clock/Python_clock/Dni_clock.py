from tkinter import
from tkinter.ttk import
from time import strftime
from numpy import floor #the greatest integer that does not exceed x.
#int(a//1) #the value of x with the fractional portion removed.
def Earth2Dni(arg):
    delta = strftime('%H:%M:%S %p')
    har = floor(Millisec / MillisecPerHar)
   pass

def get_time(arg):
    Day = strftime('%d')
    Month = strftime('%m')
    Year = strftime('%Y')
    Hour = strftime('%H')
    Minute = strftime('%M')
    Second = strftime('%S')
    return Day, Month, Year, Hour, Minute, Second

def FIX(arg):
    return nt(arg//1)

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
    Z = (JD - INT(JD)) * 86400
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
    return Yahr, Vailee, Hahr

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
