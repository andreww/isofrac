#!/usr/bin/env python

import numpy as np

import bulk_run_phonons
import fit_beta_V
import process_PVT_castep
import bm3_eos as eos
import earthref

earth_model = earthref.EarthModel(earthref.ak135)

def depth_PT(depth):
    """Retrun liquidus P and T at a given depth in a magma ocean

    Liquidus data Andrault et at. 2011 (EPSL doi:10.1016/j.epsl.2011.02.006)
    who fit a modified Simmon and Glatzel equation:

         T = T0 (P/a+1_^(1/c) 

    (see section 3.4) with parameters listed below. This replaces a
    previous linear fit to data at 0 and 60 GPa.
    """
    
    P = earth_model(6371-depth) # Interpolating AK135...
    # We now have P, T is from TP plot
    T_0 = 1940.0 # virtual liqidus temperature at 0 GPa
    a = 26.0 # GPa
    c = 1.9
    T = T_0 * ((P / a) + 1)**(1/c)
    return T, P

def fit_beta(files, supercell=False):

    paths_and_seeds = bulk_run_phonons.process_paths(files)
    data = fit_beta_V.get_data(paths_and_seeds, supercell=supercell)
    A1, A2, A3, B1, B2, B3, C1, C2, C3 = fit_beta_V.fit_beta_T_V(data, plot=False)

    def get_beta_T_V(T, V):
        return fit_beta_V.ln_beta_V_function_wrap(T, V, A1, A2, A3, B1, 
                                                 B2, B3, C1, C2, C3)

    return np.vectorize(get_beta_T_V)

def fit_PVT_EOS_params(files):

    data = []
    for f in files:
        print(f)
        data = process_PVT_castep.parse_castep_file(f, data)

    Ts = [0, 500, 1000, 1500, 2000, 2500, 3000, 3500]
    Vs = []
    Fs = []
    K0s = []
    Kp0s = []
    E0s = []
    V0s = []
    for T in Ts:
        V, F = process_PVT_castep.get_VF(data, T)
        V0, E0, K0, Kp0 =  eos.fit_BM3_EOS(V, F, verbose=True)
        Vs.append(V)
        Fs.append(F)
        K0s.append(K0)
        Kp0s.append(Kp0)
        E0s.append(E0)
        V0s.append(V0)
    fV0, fE0, fK0, fKp0 = eos.fit_parameters_quad(Ts, V0s, E0s, K0s, Kp0s,
        plot=False)

    def get_volume(P, T):
        return eos.get_V(P, T, fV0, fK0, fKp0)

    return np.vectorize(get_volume)


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


if __name__ == "__main__":

    import glob

    import matplotlib
    import matplotlib.pyplot as plt

    # Depth range of interest
    depths = np.linspace(0.0, 2800.0, num=200)

    # Get our list of Ps and Ts
    Ts, Ps = depth_PT(depths)

    # Volume of MgO
    MgO_eos = fit_PVT_EOS_params(
        glob.glob('../free_energy/MgO/MgO_*GPa/MgO.castep'))
    MgO_Vs = MgO_eos(Ps, Ts)
    MgO_Vs_athermal = MgO_eos(Ps, np.zeros_like(Ts))

    # Volume of MgSiO3 Pv
    MgSiO3_eos = fit_PVT_EOS_params(
        glob.glob('../free_energy/MgSiO3/MgSiO3_*GPa/MgSiO3.castep'))
    MgSiO3_Vs = MgSiO3_eos(Ps, Ts)
    MgSiO3_Vs_athermal = MgSiO3_eos(Ps, np.zeros_like(Ts))
    
    # Volume of MgSiO3 Pv
    Mg2SiO4_eos = fit_PVT_EOS_params(
        glob.glob('../free_energy/Mg2SiO4/Mg2SiO4_*GPa/Mg2SiO4.castep'))
    Mg2SiO4_Vs = Mg2SiO4_eos(Ps, Ts)
    Mg2SiO4_Vs_athermal = Mg2SiO4_eos(Ps, np.zeros_like(Ts))

    # 1000.ln(beta) for MgO
    MgO_beta_fun = fit_beta(glob.glob('../free_energy/MgO/MgO_*GPa/MgO.castep')) 
    MgO_betas = MgO_beta_fun(Ts, MgO_Vs)
    MgO_betas_athermal = MgO_beta_fun(Ts, MgO_Vs_athermal)

    # 1000.ln(beta) for MgSiO3
    MgSiO3_beta_fun = fit_beta(glob.glob('../free_energy/MgSiO3/MgSiO3_*GPa/MgSiO3.castep')) 
    MgSiO3_betas = MgSiO3_beta_fun(Ts, MgSiO3_Vs)
    MgSiO3_betas_athermal = MgSiO3_beta_fun(Ts, MgSiO3_Vs_athermal)

    # 1000.ln(beta) for Mg2SiO3
    Mg2SiO4_beta_fun = fit_beta(glob.glob('../free_energy/Mg2SiO4/Mg2SiO4_*GPa/Mg2SiO4.castep')) 
    Mg2SiO4_betas = Mg2SiO4_beta_fun(Ts, Mg2SiO4_Vs)
    Mg2SiO4_betas_athermal = Mg2SiO4_beta_fun(Ts, Mg2SiO4_Vs_athermal)


    print("Done fitting... now some key data" )
    print("P(GPa) T(K) Depth(km), 1000.ln(alpha(Fo, MgO)), 1000.ln(alpha(Fo,MgPv)")
    for P, T, D, B_Fo, B_MgO, B_MgPv in zip(Ps, Ts, depths, Mg2SiO4_betas, MgO_betas, MgSiO3_betas):
        print(P, T, D, B_Fo-B_MgO, B_Fo-B_MgPv)

    print("Done fitting... now plotting")

    f, ax1 = plt.subplots()

    fs = 14
    fs_l = fs

    ax_depths = np.array([0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 
                          2400, 2600, 2800])
    ax_Ts, ax_Ps = depth_PT(ax_depths)
    ax1.invert_yaxis()
    ax1.plot(ax_Ts, ax_Ps, ':g')
    ax1.set_xlim(left=2000, right=5000)
    ax1.set_xlabel("T (K)", fontsize=fs)
    ax1.set_ylabel("P (GPa)", fontsize=fs)
    ax1.tick_params(axis='both', which='both', labelsize=fs_l)

    ax2 = ax1.twinx()
    ax2.invert_yaxis()
    ax2.plot(ax_Ts, ax_depths, alpha=0) # invisable
    ax2.set_ylabel("Depth (km)", fontsize=fs)
    ax2.tick_params(axis='both', which='both', labelsize=fs_l)
    ax3 = ax2.twiny()
    # For Tim's latest we want Mg25 - aparantly half the fractionation
    mg_25 = False
    if mg_25:
        ax3.set_xlabel(r"$\Delta^{}$Mg (per mill) relative to forsterite".format('{25}'), fontsize=fs)
        ax3.set_xlim(left=0.0, right=0.12/2.0)
    else:
        ax3.set_xlabel(r"$\Delta^{}$Mg (per mill) relative to forsterite".format('{26}'), fontsize=fs)
        ax3.set_xlim(left=0.0, right=0.12)

    ax3.tick_params(axis='both', which='both', labelsize=fs_l)
    if mg_25:
        ax3.plot((Mg2SiO4_betas - MgSiO3_betas)/2.0, depths, 'r-')
    else: 
        ax3.plot((Mg2SiO4_betas_athermal - MgO_betas_athermal), depths, 'b--')
        ax3.plot((Mg2SiO4_betas - MgO_betas), depths, 'b-')
        ax3.plot((Mg2SiO4_betas - MgSiO3_betas), depths, 'r-')
        ax3.plot((Mg2SiO4_betas_athermal - MgSiO3_betas_athermal), depths, 'r--')

    
        #ax2 = ax1.twinx()
        #ax2_tick_ds = np.array([200, 400, 600, 800, 1000])
        #ax2_tick_Ps, ax2_tick_Ts = depth_PT(ax2_tick_ds)
        #ax2_tick_labs = ["200", "400", "600", "800", "1000"]
        #ax2.set_ylabel("P (GPa)")
        #ax2.set_yticks(ax2_tick_ds)
        #ax2.set_yticks(ax2_tick_labs)

    f.tight_layout()

    f.savefig("alpha_geotherm_4_liqidus.pdf")

    #plt.show()

    # New plot for melt...

    melt_ln_betas = ionic_model_beta(ionic_model_force_constant(melt_bond_length(Ps)), Ts)

    f, ax1 = plt.subplots()

    fs = 14
    fs_l = fs

    ax_depths = np.array([0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 
                          2400, 2600, 2800])
    ax_Ts, ax_Ps = depth_PT(ax_depths)
    ax1.invert_yaxis()
    ax1.plot(ax_Ts, ax_Ps, ':g')
    ax1.set_xlim(left=2000, right=5000)
    ax1.set_xlabel("T (K)", fontsize=fs)
    ax1.set_ylabel("P (GPa)", fontsize=fs)
    ax1.tick_params(axis='both', which='both', labelsize=fs_l)

    ax2 = ax1.twinx()
    ax2.invert_yaxis()
    ax2.plot(ax_Ts, ax_depths, alpha=0) # invisable
    ax2.set_ylabel("Depth (km)", fontsize=fs)
    ax2.tick_params(axis='both', which='both', labelsize=fs_l)
    ax3 = ax2.twiny()
    # For Tim's latest we want Mg25 - aparantly half the fractionation
    mg_25 = False
    if mg_25:
        ax3.set_xlabel(r"$\Delta^{}$Mg (per mill) relative to liquid".format('{25}'), fontsize=fs)
        ax3.set_xlim(left=0.0, right=0.12/2.0)
    else:
        ax3.set_xlabel(r"$\Delta^{}$Mg (per mill) relative to liquid".format('{26}'), fontsize=fs)
        #ax3.set_xlim(left=0.0, right=0.12)

    ax3.tick_params(axis='both', which='both', labelsize=fs_l)
    if mg_25:
        ax3.plot((melt_ln_betas - MgSiO3_betas)/2.0, depths, 'r-')
    else: 
        ax3.plot((melt_ln_betas - Mg2SiO4_betas_athermal), depths, 'k--')
        ax3.plot((melt_ln_betas - Mg2SiO4_betas), depths, 'k-')
        ax3.plot((melt_ln_betas - MgO_betas_athermal), depths, 'b--')
        ax3.plot((melt_ln_betas - MgO_betas), depths, 'b-')
        ax3.plot((melt_ln_betas - MgSiO3_betas), depths, 'r-')
        ax3.plot((melt_ln_betas - MgSiO3_betas_athermal), depths, 'r--')

    
        #ax2 = ax1.twinx()
        #ax2_tick_ds = np.array([200, 400, 600, 800, 1000])
        #ax2_tick_Ps, ax2_tick_Ts = depth_PT(ax2_tick_ds)
        #ax2_tick_labs = ["200", "400", "600", "800", "1000"]
        #ax2.set_ylabel("P (GPa)")
        #ax2.set_yticks(ax2_tick_ds)
        #ax2.set_yticks(ax2_tick_labs)

    f.tight_layout()

    f.savefig("alpha_geotherm_4_liqidus_melt.pdf")

    plt.show()
