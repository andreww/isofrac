#!/usr/bin/env python

import numpy as np

import bulk_run_phonons
import fit_beta_V
import process_PVT_castep
import bm3_eos as eos

def depth_PT(depth):
    """Retrun liquidus P and T at a given depth in a magma ocean

    Liquidus data is taken from figure 3 of Andrault et at. 2011
    (EPSL doi:10.1016/j.epsl.2011.02.006) and is for a chondritic 
    material. We assime linear behaviour to 60 GPa.
    """
    dPdDepth = 60.0/1400.0 # GPa/km
    P = depth * dPdDepth # GPa
    # We now have P, T is from TP plot
    T_0 = 2000.0 # T_liquid is 2000 K at 0 GPa
    T_60 = 3500.0 # at 60 GPa
    dTdP = (T_60 - T_0) / 60.0
    T = T_0 + P*dTdP
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
        print f
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
    

if __name__ == "__main__":

    import glob

    import matplotlib
    import matplotlib.pyplot as plt

    # Depth range of interest
    depths = np.linspace(0.0, 1400.0, num=50)

    # Get our list of Ps and Ts
    Ts, Ps = depth_PT(depths)

    # Volume of MgO
    MgO_eos = fit_PVT_EOS_params(
        glob.glob('/Volumes/earawa/Castep-isotopes-work/free_energy/MgO/MgO_*GPa/MgO.castep'))
    MgO_Vs = MgO_eos(Ps, Ts)
    MgO_Vs_athermal = MgO_eos(Ps, np.zeros_like(Ts))

    # Volume of MgSiO3 Pv
    MgSiO3_eos = fit_PVT_EOS_params(
        glob.glob('/Volumes/earawa/Castep-isotopes-work/free_energy/MgSiO3/MgSiO3_*GPa/MgSiO3.castep'))
    MgSiO3_Vs = MgSiO3_eos(Ps, Ts)
    MgSiO3_Vs_athermal = MgSiO3_eos(Ps, np.zeros_like(Ts))
    
    # Volume of MgSiO3 Pv
    Mg2SiO4_eos = fit_PVT_EOS_params(
        glob.glob('/Volumes/earawa/Castep-isotopes-work/free_energy/Mg2SiO4/Mg2SiO4_*GPa/Mg2SiO4.castep'))
    Mg2SiO4_Vs = Mg2SiO4_eos(Ps, Ts)
    Mg2SiO4_Vs_athermal = Mg2SiO4_eos(Ps, np.zeros_like(Ts))

    # 1000.ln(beta) for MgO
    MgO_beta_fun = fit_beta(glob.glob('/Volumes/earawa/Castep-isotopes-work/MgO_DFPT/MgO_prim222_900eV_k333_q444_DFPT_*GPa/MgO.castep'), supercell=True) 
    MgO_betas = MgO_beta_fun(Ts, MgO_Vs)
    MgO_betas_athermal = MgO_beta_fun(Ts, MgO_Vs_athermal)

    # 1000.ln(beta) for MgSiO3
    MgSiO3_beta_fun = fit_beta(glob.glob('/Volumes/earawa/Castep-isotopes-work/MgSiO3_DFPT/MgSiO3_900eV_k333_q333_DFPT_*GPa/mgsio3.castep')) 
    MgSiO3_betas = MgSiO3_beta_fun(Ts, MgSiO3_Vs)
    MgSiO3_betas_athermal = MgSiO3_beta_fun(Ts, MgSiO3_Vs_athermal)

    # 1000.ln(beta) for Mg2SiO3
    Mg2SiO4_beta_fun = fit_beta(glob.glob('/Volumes/earawa/Castep-isotopes-work/Forsterite_DFPT/Fo_900ev_k424_q313_DFPT_*GPa/forsterite.castep')) 
    Mg2SiO4_betas = Mg2SiO4_beta_fun(Ts, Mg2SiO4_Vs)
    Mg2SiO4_betas_athermal = Mg2SiO4_beta_fun(Ts, Mg2SiO4_Vs_athermal)



    print "Done fitting... now plotting" 

    f, ax1 = plt.subplots()
    ax1.set_xlabel("T (K)")
    ax1.set_ylabel("P (GPa)")
    ax1.invert_yaxis()
    ax1.plot(Ts, Ps, ':g')
    ax2 = ax1.twinx()
    ax2.set_ylabel("Depth (km)")
    ax2.invert_yaxis()
    ax2.plot(Ts, depths, alpha=0) # invisable

    ax3 = ax1.twiny()
    # For checking volumes are sane:
    #ax3.plot(MgO_Vs, Ps)
    #ax3.plot(MgSiO3_Vs, Ps)
    #ax3.plot(Mg2SiO4_Vs, Ps)

    # For checking betas
    #ax3.plot(MgO_betas, Ps)
    #ax3.plot(MgO_betas_athermal, Ps)
    #ax3.plot(MgSiO3_betas, Ps)
    #ax3.plot(MgSiO3_betas_athermal, Ps)
    #ax3.plot(Mg2SiO4_betas, Ps)
    #ax3.plot(Mg2SiO4_betas_athermal, Ps)

    # Alpha plots
    ax3.plot((Mg2SiO4_betas_athermal - MgO_betas_athermal), Ps, 'b--')
    ax3.plot((Mg2SiO4_betas - MgO_betas), Ps, 'b-')
    ax3.plot((Mg2SiO4_betas_athermal - MgSiO3_betas_athermal), Ps, 'r--')
    ax3.plot((Mg2SiO4_betas - MgSiO3_betas), Ps, 'r-')
    ax3.set_xlabel(r'$1000.\ln(\alpha)$ (per mill) relative to Mg$_2$SiO$_4$')


    plt.show()
