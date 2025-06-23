#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  Dni_clock.py
#
#  Copyright 2025 Diego Martinez Gutierrez <diego.martinez@ehu.eus>
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
# D'ni clock converter from Gregorian to D'ni
# ---------------------------
# Importación de los módulos
# ---------------------------
#from tkinter import
#from tkinter.ttk import
import tkinter as tk
import datetime
from time import strftime
from numpy import floor #the greatest integer that does not exceed x.
#int(a//1) #the value of x with the fractional portion removed.


# -------
# Clases
# -------

# ----------
# Funciones
# ----------

#Conversion numerical base:
def toDigits(n, b):
    """Convert a positive number n to its digit representation in base b."""
    digits = []
    while n > 0:
        digits.insert(0, n % b)
        n  = n // b
    return digits

def fromDigits(digits, b):
    """Compute the number given by digits in base b."""
    n = 0
    for d in digits:
        n = b * n + d
    return n

#Conversion to CAVE timezone
def Earth2Dni(arg):
    delta = strftime('%H:%M:%S %p')
    har = floor(Millisec / MillisecPerHar)
    pass

def get_time():
    Day = int(strftime('%d'))
    Month = int(strftime('%m'))
    Year = int(strftime('%Y'))
    Hour = int(strftime('%H'))
    Minute = int(strftime('%M'))
    Second = int(strftime('%S'))
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
    return Hahr, Vailee, Yahr, Gahrtahvo, Tahvo, Gorahn, Prorahn

def Cavernian_to_Gregorian(AY):
    #Atrian Yahr
    AYD = AY - 1.0
    #to Julian
    JDD = AYD * 1.25945582758621
    JD = JDD + 727249.704166666
    return Julian_to_Gregorian(JD)

def Gregorian_to_Cavernian(Day, Month, Year, Hour=0, Minute=0, Second=0):
    JD = Gregorian_to_Julian(Day, Month, Year, Hour, Minute, Second)
    JDD = JD - 727249.704166666
    AYD = JDD * 0.793993705929756
    AY = AYD + 1.0
    #Convert the calculated Atrian Yahr to a Cavernian date
    return AYN_to_Cavernian(AY)

"""import math
from PIL import Image
import matplotlib.pyplot as plt
import os

def create_image_row(inputs, image_dir="images", output_file="output.png"):
"""    """
    Create a single image by arranging 9 input images (0.png to 24.png) in a row with separators.

    Args:
        inputs (list): List of 9 numbers (0 to 24) to be rounded down.
        image_dir (str): Directory containing images named 0.png to 24.png.
        output_file (str): Path to save the output image.
"""    """
    # Validate input length
    if len(inputs) != 9:
        raise ValueError("Exactly 9 inputs are required.")

    # Round down inputs and ensure they are within 0-24
    rounded_inputs = []
    for num in inputs:
        if not isinstance(num, (int, float)) or num < 0 or num > 24:
            raise ValueError("All inputs must be numbers between 0 and 24.")
        rounded_inputs.append(math.floor(num))

    # Load images
    image_files = [f"{i}.png" for i in rounded_inputs]
    images = []
    for file in image_files:
        file_path = os.path.join(image_dir, file)
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"Image {file} not found in {image_dir}.")
        img = Image.open(file_path)
        images.append(img)

    # Assume all images have the same size
    img_width, img_height = images[0].size

    # Calculate total width: 9 images + 8 separators (10 pixels each)
    separator_width = 10
    total_width = 9 * img_width + 8 * separator_width
    total_height = img_height

    # Create a new blank image
    result = Image.new("RGB", (total_width, total_height), color="white")

    # Paste images and add separators
    x_offset = 0
    for i, img in enumerate(images):
        result.paste(img, (x_offset, 0))
        x_offset += img_width
        if i < 8:  # Add separator after each image except the last
            # Draw a black vertical line
            for x in range(x_offset, x_offset + separator_width):
                for y in range(total_height):
                    result.putpixel((x, y), (0, 0, 0))
            x_offset += separator_width

    ## Save the image using matplotlib to ensure proper rendering
    #plt.figure(figsize=(total_width / 100, total_height / 100))
    #plt.imshow(result)
    #plt.axis("off")
    #plt.savefig(output_file, bbox_inches="tight", pad_inches=0)
    #plt.close()

    # Display the image using matplotlib
    plt.figure(figsize=(total_width / 100, total_height / 100))
    plt.imshow(result)
    plt.axis("off")
    plt.show()

## Example usage:
#if __name__ == "__main__":
#    # Example inputs
#    sample_inputs = [0.7, 5.2, 10.8, 15.3, 20.9, 2.1, 7.6, 12.4, 24.0]
#    create_image_row(sample_inputs, image_dir="images", output_file="output.png")
"""
#-----------------------------------------------------

def main(args):
    base = 25
    VaileeDictionary = {1:"Leefo",2:"Leebro",3:"Leesahn",4:"Leetar",5:"Leevot",6:"Leevofo",7:"Leevobro",8:"Leevosahn",9:"Leevotar",10:"Leenovoo"}
    Day, Month, Year, Hour, Minute, Second = get_time()
    print("Year:", Year, "Month=", Month, "Day:", Day, "Hour=",Hour,":",Minute,":",Second)
    Hahr, Vailee, Yahr, Gahrtahvo, Tahvo, Gorahn, Prorahn = Gregorian_to_Cavernian(Day, Month, Year, Hour, Minute, Second)
    if base==10:
        print("Hahr:", Hahr, "Yahr:",Yahr, "Vailee:",Vailee,"_",VaileeDictionary[Vailee], "Gahrtahvo:", Gahrtahvo, "Tahvo:", Tahvo, "Gorahn:", Gorahn, "Prorahn:", Prorahn)
    else:
        print("Hahr:", toDigits(Hahr,25), "Yahr:",Yahr, "Vailee:",Vailee,"_",VaileeDictionary[Vailee], "Gahrtahvo:", Gahrtahvo, "Tahvo:", Tahvo, "Gorahn:", Gorahn, "Prorahn:", Prorahn)
        #print("Hahr:", toDigits(Hahr,25), "Yahr:",toDigits(Yahr), "Vailee:",toDigits(Vailee),"_",VaileeDictionary[Vailee], "Gahrtahvo:", toDigits(Gahrtahvo), "Tahvo:", toDigits(Tahvo), "Gorahn:", toDigits(Gorahn), "Prorahn:", toDigits(Prorahn))
    return 0

# Add a funtion to print instruction of use.
def tutorial(arg):
    print("Instructions of use:")
    print("")
    pass

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))

"""
import tkinter as tk
import datetime

def actualizar_hora():
    ahora = datetime.datetime.now()
    Year = ahora.year
    Month = ahora.month
    Day = ahora.day
    Hour = ahora.hour
    Minute = ahora.minute
    Second = ahora.second
    Hahr, Vailee, Yahr, Gahrtahvo, Tahvo, Gorahn, Prorahn = Gregorian_to_Cavernian(Day, Month, Year, Hour, Minute, Second)

    etiqueta_fecha.config(text=f"Fecha: {Year}-{Month:02}-{Day:02}")
    etiqueta_hora.config(text=f"Hora: {Hour:02}:{Minute:02}:{Second:02}")
    etiqueta_fecha_Dni.config(text=f"Fecha Dni: DHahr: {toDigits(Hahr,25)}-Yahr:{Yahr}-Vailee:{Vailee}_{VaileeDictionary[Vailee]}Gahrtahvo:{Gahrtahvo}-Tahvo:{Tahvo}-Gorahn:{Gorahn}-Prorahn:{Prorahn}
)

    ventana.after(1000, actualizar_hora) # Actualizar cada segundo

ventana = tk.Tk()
ventana.title("Fecha y Hora")

# Etiquetas para la fecha y la hora
etiqueta_fecha = tk.Label(ventana, text="")
etiqueta_fecha.pack(pady=10)

etiqueta_hora = tk.Label(ventana, text="")
etiqueta_hora.pack(pady=10)

# Llamar a la función para actualizar la hora cada segundo
actualizar_hora()

ventana.mainloop()
"""
