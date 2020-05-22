#!/usr/bin/env python
"""
A simple ionic model

Pythonised stuff from Remco's spreadsheet
Based on Young et al. 2015
"""

import numpy as np

def melt_bond_length(p, coeffs=None):
    """
    Return the bond length (in m) given
    a polynomial expansion in pressure.
    See Remco's
    spreadsheet. If coeffs
    this uses Remco's parameters including
    the bond length fudge
    """
    if coeffs is None:
        coeffs = [2.0909, 0.00165, 0.0000019]

    r = 0.0000000001 * (coeffs[0] + p * coeffs[1] + p**2 * coeffs[2])
    return r

def ionic_model_force_constant(r):
    """
    Ionic model force constant following equation 31 of Young (2015).
    All parameters other than r are designed to follow Remco's spreadsheet.
    """
    zi = 2.0 # Cation valence
    zj = -2.0 # Anion valence
    n = 12 # Born-Mayer constant. Ultimatly LJ I think
    eps0 = 8.854187817E-12 # Vaccum permittivity (F/m)
    e = 1.60217662E-19 # electron charge (C)
    
    kf = (zi * zj * e**2 * (1-n)) / (4.0 * np.pi * eps0 * r**3)
    return kf

def ionic_model_beta(kf, T):
   """
   Calculate beta as per equation 27 of Young et al. (2015)
   """
   h = 6.62607004E-34 # Plank constant (m^2 kg s^-1)
   kb = 1.38064852E-23 # Boltzman constant (m^2 kg s^-2 K^-1)
   m1 = 24 # amu
   m2 = 26 # amu
   u = 1.6605402E-27 # Unified atomic mass (kg) - an amu in kg

   beta = 1 + (1/24) * (h / (kb * T))**2 * ((1 / (m1 * u)) - (1 / (m2 * u))) * (kf / (4 * np.pi**2))

   # Return lmn(beta) as that's what we do with everything else
   return np.log(beta)
