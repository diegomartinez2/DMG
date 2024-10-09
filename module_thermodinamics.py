#!/usr/bin/env python
import numpy as np
from scipy.special import tanh, sinh, cosh

def w_to_a(w, T, n):
    a = np.zeros(n)
    if T == 0.0:
        a[:] = np.sqrt(1.0 / (2.0 * w))
    else:
        a[:] = np.sqrt((1.0 / tanh(0.5 * w * 315774.65221921849 / T)) / (2.0 * w))
    return a

def w_to_da(w, T, n):
    da = np.zeros(n)
    if T == 0.0:
        da[:] = -np.sqrt(1.0 / (8.0 * w**3.0))
    else:
        beta = 315774.65221921849 / T
        a = w * beta + sinh(w * beta)
        b = np.sqrt(1.0 / (32.0 * w**3.0 * (sinh(0.5 * w * beta)**3.0) * cosh(0.5 * w * beta)))
        da[:] = -a * b
    return da

def dW_f0_u0(w, T):
    if T == 0.0:
        return 0.25
    else:
        return 0.25 * (2.0 * nb(w, T) + 1.0 + 2.0 * (315774.65221921849 / T) * w * np.exp(w * 315774.65221921849 / T) * nb(w, T)**2.0)

def nb(w, T):
    if T == 0.0:
        return 0.0
    else:
        return 1.0 / (np.exp(w * 315774.65221921849 / T ) - 1.0)
