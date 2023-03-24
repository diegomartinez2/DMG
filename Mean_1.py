#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  Mean_1.py
#
#  Copyright 2023 Diego Martinez Gutierrez <diego.martinez@ehu.eus>
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
import numpy as np

def mse_1(actual, predicted):
    actual = np.array(actual)
    predicted = np.array(predicted)
    differences = np.subtract(actual, predicted)
    squared_differences = np.square(differences)
    return squared_differences.mean()

print(mse(df['y'], df['y_predicted']))

# A shorter version of the code above
#import numpy as np

def mse_2(actual, predicted):
    return np.square(np.subtract(np.array(actual), np.array(predicted))).mean()

print(mse(df['y'], df['y_predicted']))

#another
def mse_3(arg):
    y = [11,20,19,17,10]
    y_bar = [12,18,19.5,18,9]
    summation = 0  #variable to store the summation of differences
    n = len(y) #finding total number of items in list
    for i in range (0,n):  #looping through each element of the list
      difference = y[i] - y_bar[i]  #finding the difference between observed and predicted value
      squared_difference = difference**2  #taking square of the differene
      summation = summation + squared_difference  #taking a sum of all the differences
    MSE = summation/n  #dividing summation by total values to obtain average
    print("The Mean Square Error is: " , MSE)
    return MSE

#another
#import numpy as np


def main(args):
    # Given values
    Y_true = [1,1,2,2,4]  # Y_true = Y (original values)

    # Calculated values
    Y_pred = [0.6,1.29,1.99,2.69,3.4]  # Y_pred = Y'

    # Mean Squared Error

    MSE = np.square(np.subtract(Y_true,Y_pred)).mean()
    print("The Mean Square Error is: " , MSE)
    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
