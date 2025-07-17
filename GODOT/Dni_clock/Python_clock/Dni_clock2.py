#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Dni_clock.py
#
# Copyright 2025 Diego Martinez Gutierrez <diego.martinez@ehu.eus>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
# MA 02110-1301, USA.
#
#
# D'ni clock converter from Gregorian to D'ni
# ---------------------------
# Importación de los módulos
# ---------------------------
import tkinter as tk
import datetime
from time import strftime
from math import floor # Usar math.floor para consistencia si no hay más numpy

# -------
# Constantes D'ni
# -------
DAYS_IN_DNISH_YEAR = 405.0 # (29 * 14) + 1 (13 Vailee of 29 days + 1 Vailee of 28 days for the last one)
DNISH_YR_TO_GREGORIAN_DAY_RATIO = 1.25945582758621
GREGORIAN_TO_DNISH_YR_RATIO = 0.793993705929756
JULIAN_DAY_OFFSET_DNISH_EPOCH = 727249.704166666 # JD for 1 AY at 00:00:00

# D'ni time units per Yahr
DNISH_YAHR_GARHTAHVO = 78125 / 15625 # 5 Garhtahvo per Tahvo? No. 78125 units per Yahr, 15625 units per Gahrtahvo.
DNISH_YAHR_TAHVO = 78125 / 625 # 125 Tahvo per Yahr
DNISH_YAHR_GORAHN = 78125 / 25 # 3125 Gorahn per Yahr
DNISH_YAHR_PRORAHN = 78125 # 78125 Prorahn per Yahr

# D'ni Calendar Constants
DNISH_VAILEE_DAYS = 29
DNISH_HAHR_VAILEE = 10 # Number of Vailee in a Hahr
DNISH_EPOCH_HAHR = 9647

VAILEE_NAMES = {
    1:"Leefo", 2:"Leebro", 3:"Leesahn", 4:"Leetar", 5:"Leevot",
    6:"Leevofo", 7:"Leevobro", 8:"Leevosahn", 9:"Leevotar", 10:"Leenovoo"
}

# ----------
# Funciones de Utilidad
# ----------
def to_digits_base_n(number, base):
    """Convert a positive number to its digit representation in base b."""
    if number == 0:
        return [0]
    if not isinstance(number, int) or number < 0:
        raise ValueError("Number must be a non-negative integer.")
    if not isinstance(base, int) or base < 2:
        raise ValueError("Base must be an integer greater than or equal to 2.")

    digits = []
    temp_n = number
    while temp_n > 0:
        digits.insert(0, temp_n % base)
        temp_n //= base
    return digits

def from_digits_base_n(digits, base):
    """Compute the number given by digits in base b."""
    if not all(isinstance(d, int) and 0 <= d < base for d in digits):
        raise ValueError("Digits must be integers within the range [0, base-1].")
    if not isinstance(base, int) or base < 2:
        raise ValueError("Base must be an integer greater than or equal to 2.")

    number = 0
    for digit in digits:
        number = base * number + digit
    return number

# ----------
# Conversiones de Fecha y Hora
# ----------

def gregorian_to_julian(day, month, year, hour=0, minute=0, second=0):
    """Converts a Gregorian date and time to a Julian Day number."""
    if year < 0:
        year = -(year - 1) # Adjust for astronomical year numbering

    if month < 3:
        month += 12
        year -= 1

    # Whole days
    wd = day + floor(((153 * month) - 457) / 5) + \
         floor(365.25 * year) - floor(0.01 * year) + floor(0.0025 * year)

    # Fractional day
    fd = ((hour * 3600) + (minute * 60) + second) / 86400

    return wd + fd

def julian_to_gregorian(jd):
    """Converts a Julian Day number to a Gregorian date and time."""
    z = floor(jd)
    f = jd - z # Fractional part for time calculation

    g = z - 0.25
    a = floor(g / 36524.25)
    b = a - (0.25 * a)
    year = floor((g + b) / 365.25)
    c = z + b - floor(365.25 * year)

    month = floor(((5 * c) + 456) / 153)
    day = c - floor(((153 * month) - 457) / 5)

    if month > 12:
        year += 1
        month -= 12

    if year < 1:
        year = 1 - year # Adjust back for astronomical year numbering

    # Time calculation
    total_seconds = f * 86400
    hour = floor(total_seconds / 3600)
    remainder_seconds = total_seconds - (hour * 3600)
    minute = floor(remainder_seconds / 60)
    second = remainder_seconds - (minute * 60)

    return int(day), int(month), int(year), int(hour), int(minute), int(round(second))

def cavernian_to_atrian_yahr(yahr, vailee, hahr, gahrtahvo=0, tahvo=0, gorahn=0, prorahn=0):
    """Converts a Cavernian (D'ni) date and time to Atrian Yahr (decimal)."""
    # Whole yahrtee
    whole_yahrtee = yahr + ((vailee - 1) * DNISH_VAILEE_DAYS) + ((hahr - DNISH_EPOCH_HAHR) * (DNISH_VAILEE_DAYS * DNISH_HAHR_VAILEE))

    # Fractional yahr (time of yahr)
    # Total units in a Yahr is 25^4 = 390625. Wait, 78125 is 5^7. D'ni uses base 5 and base 25.
    # The original code used 78125. Let's stick to that for now, assuming 25^3 * something.
    # Gahrtahvo (25^3), Tahvo (25^2), Gorahn (25^1), Prorahn (25^0)
    # A single Yahr has 5 Gorahn, 25 Prorahn in a Gorahn, 25 Gorahn in a Tahvo, 25 Tahvo in a Gahrtahvo, 25 Gahrtahvo in a Yahr?
    # This implies a base 25 system for time units within a Yahr.
    # 1 Gahrtahvo = 25 Tahvo
    # 1 Tahvo = 25 Gorahn
    # 1 Gorahn = 25 Prorahn
    # So 1 Yahr = 25 * 25 * 25 * 25 Prorahn = 25^4 = 390625 Prorahn.
    # The constants in your code (15625, 625, 25, 78125) are:
    # 15625 = 25^3 (Gahrtahvo)
    # 625 = 25^2 (Tahvo)
    # 25 = 25^1 (Gorahn)
    # 78125 = 5^7 = 3125 * 25 = 125 * 625. This is inconsistent with 25^4.
    # Assuming 78125 as the total units per Yahr as per your code:

    fractional_yahr = ((gahrtahvo * 15625) + (tahvo * 625) + (gorahn * 25) + prorahn) / 78125

    return whole_yahrtee + fractional_yahr

def atrian_yahr_to_cavernian(ay):
    """Converts Atrian Yahr (decimal) to a Cavernian (D'ni) date and time."""
    z = floor(ay)

    # Date calculation
    g = z - 0.25 # This -0.25 is unusual for a simple floor, might relate to epoch
    a = floor(g / (DNISH_VAILEE_DAYS * DNISH_HAHR_VAILEE))
    hahr = DNISH_EPOCH_HAHR + a
    c = z - (a * (DNISH_VAILEE_DAYS * DNISH_HAHR_VAILEE))

    vailee = floor((c - 0.25) / DNISH_VAILEE_DAYS) + 1 # Similar -0.25 offset
    yahr = c - ((vailee - 1) * DNISH_VAILEE_DAYS)

    # Time calculation
    total_time_units = round((ay - floor(ay)) * 78125) # Round to nearest integer unit
    gahrtahvo = floor(total_time_units / 15625)
    remainder = total_time_units - (gahrtahvo * 15625)
    tahvo = floor(remainder / 625)
    remainder = remainder - (tahvo * 625)
    gorahn = floor(remainder / 25)
    prorahn = remainder - (gorahn * 25)

    return int(hahr), int(vailee), int(yahr), int(gahrtahvo), int(tahvo), int(gorahn), int(prorahn)

def cavernian_to_gregorian(ay):
    """Converts Atrian Yahr to Gregorian date and time."""
    ayd = ay - 1.0 # Atrian Yahr Days, assuming 1 AY = epoch start

    # Convert to Julian Day Difference from Epoch
    jdd = ayd * DNISH_YR_TO_GREGORIAN_DAY_RATIO
    jd = jdd + JULIAN_DAY_OFFSET_DNISH_EPOCH

    return julian_to_gregorian(jd)

def gregorian_to_cavernian(day, month, year, hour=0, minute=0, second=0):
    """Converts Gregorian date and time to Cavernian (D'ni) date and time."""
    jd = gregorian_to_julian(day, month, year, hour, minute, second)
    jdd = jd - JULIAN_DAY_OFFSET_DNISH_EPOCH
    ayd = jdd * GREGORIAN_TO_DNISH_YR_RATIO
    ay = ayd + 1.0 # Add 1 to make 1 AY the epoch start

    return atrian_yahr_to_cavernian(ay)

# ----------
# Interfaz Gráfica (Tkinter)
# ----------
class DniClockApp:
    def __init__(self, master):
        self.master = master
        master.title("D'ni Clock Converter")

        self.base_display = 25 # Default display base for D'ni numbers

        # Labels for Gregorian time
        tk.Label(master, text="Fecha y Hora Actual (Gregoriana)", font=("Arial", 14, "bold")).pack(pady=5)
        self.gregorian_date_label = tk.Label(master, text="", font=("Arial", 12))
        self.gregorian_date_label.pack()
        self.gregorian_time_label = tk.Label(master, text="", font=("Arial", 12))
        self.gregorian_time_label.pack()

        # Labels for D'ni time
        tk.Label(master, text="\nFecha y Hora D'ni", font=("Arial", 14, "bold")).pack(pady=5)
        self.dni_date_label = tk.Label(master, text="", font=("Arial", 12))
        self.dni_date_label.pack()

        self.update_time()

    def update_time(self):
        """Updates the time labels every second."""
        now = datetime.datetime.now()
        day, month, year = now.day, now.month, now.year
        hour, minute, second = now.hour, now.minute, now.second

        hahr, vailee, yahr, gahrtahvo, tahvo, gorahn, prorahn = \
            gregorian_to_cavernian(day, month, year, hour, minute, second)

        # Update Gregorian labels
        self.gregorian_date_label.config(text=f"Fecha: {year}-{month:02}-{day:02}")
        self.gregorian_time_label.config(text=f"Hora: {hour:02}:{minute:02}:{second:02}")

        # Format D'ni numbers in base 25 for display
        hahr_dni_digits = "".join(map(str, to_digits_base_n(hahr, self.base_display)))
        yahr_dni_digits = "".join(map(str, to_digits_base_n(yahr, self.base_display)))
        gahrtahvo_dni_digits = "".join(map(str, to_digits_base_n(gahrtahvo, self.base_display)))
        tahvo_dni_digits = "".join(map(str, to_digits_base_n(tahvo, self.base_display)))
        gorahn_dni_digits = "".join(map(str, to_digits_base_n(gorahn, self.base_display)))
        prorahn_dni_digits = "".join(map(str, to_digits_base_n(prorahn, self.base_display)))

        vailee_name = VAILEE_NAMES.get(vailee, "Desconocido") # Handle unknown Vailee

        self.dni_date_label.config(
            text=f"Hahr: {hahr_dni_digits} / Yahr: {yahr_dni_digits} / Vailee: {vailee} ({vailee_name})\n"
                 f"Gahrtahvo: {gahrtahvo_dni_digits} / Tahvo: {tahvo_dni_digits} / Gorahn: {gorahn_dni_digits} / Prorahn: {prorahn_dni_digits}"
        )

        self.master.after(1000, self.update_time) # Update every second

def main():
    root = tk.Tk()
    app = DniClockApp(root)
    root.mainloop()

if __name__ == '__main__':
    main()
