#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  untitled.py
#
#  Copyright 2022 Diego Martinez Gutierrez <diego.martinez@ehu.eus>
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

# -------
# Clases
# -------
class NombredeClase(object):
    """docstring for NombredeClase."""

    def __init__(self, arg):
        super(NombredeClase, self).__init__()
        self.arg = arg
class Dielectric_function_RPA_Jellium(object):
    """Dielectric Function Random Path Approximation with Jellium model"""
    def __init__(self,m):
        self.e = 1.602176462e-19 #C (SI)
        self.k_F = 4.5437957e14 #J-1 (SI) 1.1663787(6)e-5 GeV-2
        self.hbar = 1.054571817e-34 #Js (SI) 6.582119569e-16 eVs
        self.k2_FT = (4.0*m*(e**2)*k_F)/(np.pi*(hbar**2))  #Fermi-Thomas wavenumber
        self.v_0 = (hbar*k_F)/m #Fermi velocity
        return 0
    def real_part(self, q, omega, m):
        part_1 = 1-(((omega-((self.hbar*(q**2))/(2*m)))**2)/((q**2)*(self.v_0**2)))
        part_2 = 1-(omega+((self.hbar*(q**2))/(2*m)))**2/()(q**2)*(self.v_0**2))
        part_3 = (omega-q*v_0-(self.hbar*q/(2*m)))/(omega+q*self.v_0-(slef.hbar*q/(2*m)))
        part_4 = (omega+q*v_0+(self.hbar*q/(2*m)))/(omega-q*self.v_0+(slef.hbar*q/(2*m)))
        epsilon = 1 + self.k2_FT/(q**2) * (
        0.5*self.k_F/(4*q) * part_1*np.log(part_3) + part_2 *np.log(part_4)
        )
        return epsilon
    def imaginary_part(self, q, omega, m):
        part_1 = 1-(((omega-((self.hbar*(q**2))/(2*m)))**2)/((q**2)*(self.v_0**2)))
        if (q <= (2*self.k_F)):
            if (0.0 <= omega <= (q*self.v_0-(self.hbar*q**2)/(2*m))):
                epsilon = (np.pi/2) * (self.k2_FT/self.v_0) * (omega/(q**3))
                return epsilon
            elif ((q*self.v_0-(self.hbar*q**2)/(2*m)) <= omega <= (q*self.v_0+(self.hbar*q**2)/(2*m))):
                epsilon = (np.pi/4) * ((self.k2_FT*self.k_F)/q**3) * part_1
                return epsilon
            else:
                return 0
        if  (q >= (2*self.k_F)):
            if ((-q*self.v_0+(self.hbar*q**2)/(2*m)) <= omega <= (q*self.v_0+(self.hbar*q**2)/(2*m))):
                epsilon = (np.pi/4) * ((self.k2_FT*self.k_F)/q**3) * part_1
                return epsilon
            else:
                return 0
        else:
            raise RuntimeError("Something bad happened")
        return 0

# ----------
# Funciones
# ----------
def Plasma_frequency_SI(n_e):
    """
    Plasma frequency for 'cold' electrons approximation.
    e = 1.602176462e-19  C             (Electronic Charge)
    epsilon_0 = 8.854187813e-12  F⋅m−1 (Permittivity of Vaccum)
    m_e = 9.10938188e-31  Kg           (Mass of Electron)
    """
    omega_p_SI = np.sqrt((n_e*1.602176462e-19)/(9.10938188e-31*8.854187813e-12))
    return omega_p_SI

def Resonant_frequency_SI(n_e):
    """
    Look if this is OK
    """
    omega_0 = 1/(2*np.pi*np.sqrt((4*np.pi*(n_e*1.602176462e-19)**2)/(3*9.10938188e-31)))
    return omega_0

def dielectric_as1(omega,w_vector):
    epsilon_as = w_vector*w_vector*8.988e16     #c in SI m/s c**2 in J/Kg
    epsilon_as = epsilon_as/omega
    return epsilon_as

def dielectric_as2(omega,epsilon_inf,epsilon_0,omega_T):
    epsilon_as = epsilon_inf + (epsilon_inf-epsilon_0)/((omega**2/omega_T**2)-1)
    return epsilon_as2

def lydanne_sachs_teller(omega_T,epsilon_0,epsilon_inf):
    Omega_l = np.sqrt(omega_T*omega_T*(epsilon_0/epsilon_inf))
    return Omega_l

def Displacent_polarizability_omega_T(omega,epsilon_0,epsilon_inf):
    omega_T=np.sqrt(omega*omega*(epsilon_inf+2)/(epsilon_0+2)) # omega_T depends on omega??
    return omega_T

def dielectric_harmonic(omega, epsilon_inf, N, atom_a, atom_b, nu):
    """
    Input data:
     omega = Frequency
     epsilon_inf = dielctric constant of vacuum
     N =
     a = atom a -> M(a) mass of atom a
     b = atom b -> Z(b) atomic number of atom b
     nu = damping constant
     ---------
     Z() = Born effective charge
     M() = Atomic masses
     e() =
     omega_nu = resonant frequency
     Big_omega =
    """
    #electric_charge = 4.803e-10 #Fr (CGS)
    electric_charge = 1.602176462e-19 #C (SI)

    response1 = -(N/Big_omega) * electric_charge**2
    for a in range(atom_a):
        for b in range(atom_b):
            temp = ((Z(a)*Z(b))/np.sqrt(M(a)*M(b)))*G(a,b,omega,nu,mu)
            response2 += temp
    response_function = response1*response2

    epsilon=epsilon_inf+4*np.pi*response_function

def G(a,b,omega,nu,mu):
    for mu_index in range(mu):
        G += e(a,mu_index)*e(b,mu_index)*(1/(2*omega_nu(mu_index)))*((1/(omega-omega_nu(mu_index)+1.j*nu))-(1/(omega+omega_nu(mu_index)+1.j*nu)))
    return G

def Dielectric_function_SI(kappa,omega):
    epsilon_0 = 8.854187813e-12
    return epsilon_0*np.identity()+1.j*electric_conductivity(kappa,omega)/omega
    #return epsilon_0*np.identity()+epsilon_0*electric_susceptibility(kappa,omega)
def Dielectric_function_G_CGS():
    return np.identity()+4*np.pi()*electric_conductivity(kappa,omega)/omega
    #return np.identity()+4*np.pi()*electric_susceptibility(kappa,omega)

def electric_conductivity(kappa,omega):
    return (1.j*omega/(4*np.pi))*(1-epsilon) #this epsilon is the dielectric function?!!
def electric_susceptibility(kappa,omega):
    for i in range(atoms):
        sum_polarizability+=polarizability(atom_type(i))
    return sum_polarizability/(epsilon_0*Delta_V)*E_loc/E #S.I. <- note that the local field E_loc is defined by th dielectric function!!

def main(args):
    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
