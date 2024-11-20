from tkinter import
from tkinter.ttk import
from time import strftime
from numpy import floor #the greatest integer that does not exceed x.
#int(a//1) #the value of x with the fractional portion removed.
def Earth2Dni(arg):
    delta = strftime('%H:%M:%S %p')
    har = floor(Millisec / MillisecPerHar)
    if Year<0:
        Year = -(Year - 1)
    if Month < 3 then
        Month = Month + 12
        Year = Year - 1
#whole days
    WD = Day + FIX(((153 * Month) - 457) / 5) + floor(365.25 * Year) - floor(0.01 * Year) + floor(0.0025 * Year)
#fractional day
    FD = ((Hour * 3600) + (Minute * 60) + Second) / 86400
# Julian Day
    JD = WD + FD
   pass
def FIX(arg):
    return nt(arg//1)
def Gregorian_to_Julian(arg):
    pass
