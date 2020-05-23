#!/usr/bin/env python
"""
A simple ionic model

Pythonised stuff from Remco's spreadsheet
Based on Young et al. 2015
"""

import numpy as np

def melt_bond_length(p, coeffs):
    """
    Return the bond length (in m) given
    a polynomial expansion in pressure.
    """

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

   # Return 1000.ln(beta) as that's what we do with everything else
   frac_factor = 1000.0 * np.log(beta)
   return frac_factor


def plot_force_constants(pressures, coeff_sets, names=None, styles=None,
                         colors=None, filename=None):

    import matplotlib
    if filename is not None:
        matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, axs = plt.subplots(ncols=2)
    if type(coeff_sets) is list or type(coeff_sets) is tuple:
        if (names is None) and (styles is None):
            for coeffs in coeff_sets:
                r = melt_bond_length(pressures, coeffs)
                k = ionic_model_force_constant(r)
                axs[0].plot(pressures, r/1-E10)
                axs[1].plot(pressures, k)
        elif styles is None:
            for coeffs, name in zip(coeff_sets, names):
                r = melt_bond_length(pressures, coeffs)
                k = ionic_model_force_constant(r)
                axs[0].plot(pressures, r/1E-10, label=name)
                axs[1].plot(pressures, k, label=name)
        else:
            for coeffs, name, style, color in zip(coeff_sets, names,
                                                styles, colors):
                r = melt_bond_length(pressures, coeffs)
                k = ionic_model_force_constant(r)
                axs[0].plot(pressures, r/1E-10, label=name, linestyle=style,
                         color=color)
                axs[1].plot(pressures, k, label=name, linestyle=style,
                         color=color)
    else:
        r = melt_bond_length(pressures, coeff_sets)
        k = ionic_model_force_constant(r)
        axs[0].plot(pressures, r/1E-10, "b-")
        axs[1].plot(pressures, k, "b-")

    axs[0].set_ylabel("Bond length (Angstroms)")
    axs[0].set_xlabel("P (GPa)")
    axs[1].set_ylabel("Force constant (N/m ???)")
    axs[1].set_xlabel("P (GPa)")
    plt.tight_layout()

    if names is not None:
        plt.legend(loc=2)

    if filename is not None:
        plt.savefig(filename)
    else:
        plt.show()


if __name__ == "__main__":

    # If run on the command line, just make some plots

    import castep_isotope_sub

    # From Remco's spreadsheet. r(p) = c0 + c1*p + c2*p^2
    # where r is the bond length (in angstroms here) and
    # p is pressure (in GPa)
    r_coefs_melt_dekoker = [1.9613, 0.00165, 0.0000019]
    r_corfs_melt_fudge =   [2.0909, 0.00165, 0.0000019]


    # Plot r and force constant as function of P
    plot_force_constants(np.linspace(0, 150, num=100),
                         [r_coefs_melt_dekoker, r_corfs_melt_fudge],
                         names=['de Koker', 'fudge'], styles=['-','-'],
                         colors=['k', 'b'], filename="ionic_model_rs.pdf")


    temps = np.concatenate((np.linspace(300.0, 500.0, num=40), 
                            np.linspace(501.0, 4000.0, num=40)))

    betas = []
    names = []
    colors = []
    styles = []

    names.append('de Koker 0 GPa')
    betas.append(ionic_model_beta(ionic_model_force_constant(
                 melt_bond_length(0, r_coefs_melt_dekoker)),temps))
    colors.append('k')
    styles.append('-')

    names.append('de Koker 50 GPa')
    betas.append(ionic_model_beta(ionic_model_force_constant(
                 melt_bond_length(50, r_coefs_melt_dekoker)),temps))
    colors.append('k')
    styles.append('--')

    names.append('de Koker 100 GPa')
    betas.append(ionic_model_beta(ionic_model_force_constant(
                 melt_bond_length(100, r_coefs_melt_dekoker)),temps))
    colors.append('k')
    styles.append(':')

    castep_isotope_sub.plot_beta(temps, betas , names,
            colors=colors, styles=styles, filename='test.pdf') 
