#!/usr/bin/env python
"""
A simple ionic model

Pythonised stuff from Remco's spreadsheet
Based on Young et al. 2015
"""

import numpy as np
import scipy.optimize as spopt

def melt_bond_length(p, coeffs):
    """
    Return the bond length (in m) given
    a polynomial expansion in pressure.
    """

    r = 0.0000000001 * (coeffs[0] + p * coeffs[1] + p**2 * coeffs[2])
    return r

def ionic_model_force_constant(r, correction=1.0, offset=0.0, n=12):
    """
    Ionic model force constant following equation 31 of Young (2015).
    All parameters other than r are designed to follow Remco's spreadsheet.

    The optional correction allows a (r independent) change to the force
    constant to be made (should be equivlent to changing the constant term
    in the r expansion but easer to justify).
    """
    zi = 2.0 # Cation valence
    zj = -2.0 # Anion valence
    # n is Born-Mayer constant. Ultimatly LJ I think. Default to 12
    eps0 = 8.854187817E-12 # Vaccum permittivity (F/m)
    e = 1.60217662E-19 # electron charge (C)
    
    kf = (zi * zj * e**2 * (1-n)) / (4.0 * np.pi * eps0 * r**3)
    kf = (kf * correction) + offset
    return kf


def calculate_force_constant_correction(target_beta, c, t, p=0, mode='correction'):
   """
   Find a correction term for the ionic model force constants such
   that beta for these r coefficents and at this temperature is equal
   to target_beta.
   """
   def get_my_beta(corec):
       r = melt_bond_length(p, c)
       if mode=='correction':
           kf = ionic_model_force_constant(r, correction=corec)
       elif mode=='offset':
           kf = ionic_model_force_constant(r, offset=corec)
       beta = ionic_model_beta(kf, t)
       return beta

   def beta_error(corec):
       return get_my_beta(corec) - target_beta
   
   if mode=='correction':
       corec_needed = spopt.brentq(beta_error, 0.1, 10)
   elif mode=='offset':
      ks_uncor = ionic_model_force_constant(melt_bond_length(p, c))
      corec_needed = spopt.brentq(beta_error, -0.999 * ks_uncor, 100*ks_uncor ) 
   return corec_needed

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
                         colors=None, filename=None, kcorrs=None, offsets=None):

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
            for coeffs, name, style, color, cor, off in zip(coeff_sets, names,
                                                styles, colors, kcorrs, offsets):
                r = melt_bond_length(pressures, coeffs)
                k = ionic_model_force_constant(r, correction=cor, offset=off)
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
    import alpha_plotting

    # From Remco's spreadsheet. r(p) = c0 + c1*p + c2*p^2
    # where r is the bond length (in angstroms here) and
    # p is pressure (in GPa)
    r_coefs_melt_dekoker = [1.9613, -0.00165, 0.0000019]
    r_coefs_melt_fudge =   [2.0909, -0.00165, 0.0000019]

    forsterite_r_0 = 2.1170E-10
    forsterite_r_50 = 1.9520E-10
    forsterite_r_100 = 1.870E-10

    beta_dek_lq = ionic_model_beta(ionic_model_force_constant(
                         melt_bond_length(0, r_coefs_melt_dekoker)),1573.0)
    beta_fudge_lq = ionic_model_beta(ionic_model_force_constant(
                         melt_bond_length(0, r_coefs_melt_fudge)),1573.0)
    beta_fo = ionic_model_beta(ionic_model_force_constant(
                         forsterite_r_0),1573.0)

    print("Reduced frac factor liquid de Koker:", beta_dek_lq, "per mil")
    print("Reduced frac factor liquid r fudge:", beta_fudge_lq, "per mil")
    print("Reduced frac factor forsterite r:", beta_fo, "per mil")
    print("Fractionation factor de Koker - fo:", beta_fo - beta_dek_lq, "per mil")
    print("Fractionation factor fudge - fo:", beta_fo - beta_fudge_lq, "per mil")

    meaured_alpha = 0.080
    meaured_alpha = 0.07677 # Test to check the bond length fudge matces
    target_beta = meaured_alpha + beta_fo
    k_corr = calculate_force_constant_correction(target_beta,
                                                 r_coefs_melt_dekoker, 1573.0) 
    print("We need a correction to the force constant of:", k_corr, "UNITS?")

    print("Corrected beta:", ionic_model_beta(ionic_model_force_constant(
                 melt_bond_length(0, r_coefs_melt_dekoker),
                 correction=k_corr), 1573.0))
    print("Corrected alpha:", beta_fo - ionic_model_beta(ionic_model_force_constant(
                 melt_bond_length(0, r_coefs_melt_dekoker),
                 correction=k_corr), 1573.0))


    # Plot r and force constant as function of P
    plot_force_constants(np.linspace(0, 150, num=100),
                         [r_coefs_melt_dekoker, r_coefs_melt_fudge, r_coefs_melt_dekoker],
                         names=['de Koker', 'fudge', 'k cor'], 
                         kcorrs=[1,1,k_corr],
                         styles=['-','-', ':'],
                         colors=['k', 'b', 'g'], filename="ionic_model_rs.pdf")


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

    names.append('de Koker 0 GPa kcor')
    betas.append(ionic_model_beta(ionic_model_force_constant(
                 melt_bond_length(0, r_coefs_melt_dekoker),
                 correction=k_corr),temps))
    colors.append('y')
    styles.append('-')

    names.append('de Koker 50 GPa kcor')
    betas.append(ionic_model_beta(ionic_model_force_constant(
                 melt_bond_length(50, r_coefs_melt_dekoker),
                 correction=k_corr),temps))
    colors.append('y')
    styles.append('--')

    names.append('de Koker 100 GPa kcor')
    betas.append(ionic_model_beta(ionic_model_force_constant(
                 melt_bond_length(100, r_coefs_melt_dekoker),
                 correction=k_corr),temps))
    colors.append('y')
    styles.append(':')

    names.append('fudge 0 GPa')
    betas.append(ionic_model_beta(ionic_model_force_constant(
                 melt_bond_length(0, r_coefs_melt_fudge)),temps))
    colors.append('b')
    styles.append('-')

    names.append('fudge 50 GPa')
    betas.append(ionic_model_beta(ionic_model_force_constant(
                 melt_bond_length(50, r_coefs_melt_fudge)),temps))
    colors.append('b')
    styles.append('--')

    names.append('fudge 100 GPa')
    betas.append(ionic_model_beta(ionic_model_force_constant(
                 melt_bond_length(100, r_coefs_melt_fudge)),temps))
    colors.append('b')
    styles.append(':')

    names.append('forsterite ionic model 0 GPa')
    betas.append(ionic_model_beta(ionic_model_force_constant(
                 forsterite_r_0),temps))
    colors.append('g')
    styles.append('-')

    names.append('forsterite ionic model 50 GPa')
    betas.append(ionic_model_beta(ionic_model_force_constant(
                 forsterite_r_50),temps))
    colors.append('g')
    styles.append('--')

    names.append('forsterite ionic model 100 GPa')
    betas.append(ionic_model_beta(ionic_model_force_constant(
                 forsterite_r_100),temps))
    colors.append('g')
    styles.append(':')

    castep_isotope_sub.plot_beta(temps, betas , names,
            colors=colors, styles=styles, filename='ionic_model_beta.pdf') 


    # Plot alpha for Fo # NB: this is liq - Fo not Fo - liq.
    alpha_t_0gpa = ((ionic_model_beta(ionic_model_force_constant(
                         melt_bond_length(0, r_coefs_melt_fudge)),temps)) -
                    (ionic_model_beta(ionic_model_force_constant(
                         forsterite_r_0),temps)))
    alpha_t_50gpa = ((ionic_model_beta(ionic_model_force_constant(
                         melt_bond_length(50, r_coefs_melt_fudge)),temps)) -
                    (ionic_model_beta(ionic_model_force_constant(
                         forsterite_r_50),temps)))
    alpha_t_100gpa = ((ionic_model_beta(ionic_model_force_constant(
                         melt_bond_length(100, r_coefs_melt_fudge)),temps)) -
                    (ionic_model_beta(ionic_model_force_constant(
                         forsterite_r_100),temps)))

    dkalpha_t_0gpa = ((ionic_model_beta(ionic_model_force_constant(
                         melt_bond_length(0, r_coefs_melt_dekoker)),temps)) -
                    (ionic_model_beta(ionic_model_force_constant(
                         forsterite_r_0),temps)))
    dkalpha_t_50gpa = ((ionic_model_beta(ionic_model_force_constant(
                         melt_bond_length(50, r_coefs_melt_dekoker)),temps)) -
                    (ionic_model_beta(ionic_model_force_constant(
                         forsterite_r_50),temps)))
    dkalpha_t_100gpa = ((ionic_model_beta(ionic_model_force_constant(
                         melt_bond_length(100, r_coefs_melt_dekoker)),temps)) -
                    (ionic_model_beta(ionic_model_force_constant(
                         forsterite_r_100),temps)))

    import matplotlib.pyplot as plt
    fix, ax1 = plt.subplots()
    ax1.plot(temps, alpha_t_0gpa, 'k-', label='[fudge lq - fo] 0 GPa')
    ax1.plot(temps, alpha_t_50gpa, 'b-', label='[fudge lq - fo] 50 GPa')
    ax1.plot(temps, alpha_t_100gpa, 'g-', label='[fudge lq - fo] 100 GPa')
    ax1.plot(temps, dkalpha_t_0gpa, 'k:', label='[dk lq - fo] 0 GPa')
    ax1.plot(temps, dkalpha_t_50gpa, 'b:', label='[dk lq - fo] 50 GPa')
    ax1.plot(temps, dkalpha_t_100gpa, 'g:', label='[dk lq - fo] 100 GPa')
    ax1.set_ylabel(r"$\Delta^{}$Mg (per mill) relative to {}".format('{26}', "liquid"))
    ax1.set_xlabel("T (K)")
    plt.legend()
    plt.savefig("ionic_model_alpha.pdf")
