#!/usr/bin/env python
"""
fit_beta_V

Script to fit 1000.ln(beta(T,V). Assumes phonons
tool has already been run in bulk.
"""
import os
import re
import glob
import subprocess
import numpy as np
import scipy.optimize as spopt

import calc_beta
import castep_isotope_sub

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

get_volumes_RE = re.compile( r"^\s+Current cell volume =\s+(\d+\.\d+)\s+A\*\*3",
        re.MULTILINE)

def get_volume(seedname):
    """Read the final cell volume from a .castep file

    This can be used to plot e.g. beta as a function of
    cell volume in cases where the volume in the .phonon
    file does not represent the optimisesd structure. 
    """
    fh = open(seedname+".castep", "r")
    lines = fh.read()
    fh.close()
    volume_grps = get_volumes_RE.findall(lines)
    volume = float(volume_grps[-1]) # Last one is probably the optimized volume
    return volume

def fit_V_beta_func(data):

    Ts = np.linspace(300.0, 4000.0, num=50)
#        betas = calc_beta.beta_T(Ts, 1, v, vs, w, ws)
    # Get betas for T = 300 - 4000

    (v, w, vs, ws) = get_freqs(seedname, fineqpoints)
    Ts = np.linspace(300.0, 4000.0, num=4000)
    betas = calc_beta.beta_T(Ts, 1, v, vs, w, ws)
    lnbetas = 1000.0 * np.log(betas)


    # Check error so that caller can check value
    calc_betas = ln_beta_function(Ts, popt[0], popt[1], popt[2])
    error = np.abs(lnbetas-calc_betas)
    max_error = np.max(error)

    return (popt, pconv, max_error)

def fit_beta_T_V(data):

    print "Fitting beta to T,V data"
    Ts = np.linspace(500.0, 3500.0, num=50)
    allTs = []
    allVs = []
    betas = []
    for V in data.keys():
        (v, w, vs, ws) = data[V]
        betas.extend(calc_beta.beta_T(Ts, 1, v, vs, w, ws))
        allTs.extend(Ts)
        allVs.extend(np.ones(np.size(Ts))*V)

    TVs = np.array([allTs, allVs])
    betas = np.array(betas)
    lnbetas = 1000.0 * np.log(betas)
    # Fit to functional form for ln(beta)
    popt, pconv = spopt.curve_fit(ln_beta_V_function, TVs, 
        lnbetas, p0=[1E14, 1.7E16, 0, -1E10, -1E10, 1E10, 1E6, 5E6, 1E6])

    # Check results...
    calc_betas = ln_beta_V_function(TVs, popt[0], popt[1], popt[2],
         popt[3], popt[4], popt[5], popt[6], popt[7], popt[8])
    max_error = np.max(np.abs(lnbetas - calc_betas))

    print "For function:\n  1000 ln(beta) = " + \
        "(A1+A2.V^-1+A3.V^-2)/T^6 + (B1+B2.V^-1+B3.V^-2)/T^4 + (C1+C2.V-1+C3.V-2)/T^2"
    print "parameters are: \n  A1 = {:7g} \t   B1 = {:7g} \t   C1 = {:7g}".format(
                popt[0], popt[3], popt[6])
    print "  A2 = {:7g} \t   B2 = {:7g} \t   C2 = {:7g}".format(
                popt[1], popt[4], popt[7])
    print "  A3 = {:7g} \t   B3 = {:7g} \t   C3 = {:7g}".format(
                popt[2], popt[5], popt[8])
    print "maximum error is: {:7g}".format(max_error)

    fig = plt.figure(figsize=(12.0,10.0), dpi=600)
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(TVs[0], TVs[1], lnbetas, c='r')
    ax.plot_trisurf(TVs[0], TVs[1], calc_betas, cmap=plt.get_cmap('Greens'))
    ax.set_xlabel('Temperature (K)')
    ax.set_ylabel('Volume (A$^3$)')
    ax.set_zlabel(r'$1000.\ln(\beta)$ (per mill)')
    ax.set_xlim(500, 3500)
    plt.savefig('beta_V_T.png')

    # Convergence is to ~0.05 per mil, so worse than this is a problem
    assert max_error < 0.06, ValueError

    return popt[0], popt[1], popt[2], popt[3], popt[4], \
           popt[5], popt[6], popt[7], popt[8]


def ln_beta_V_function_wrap(T, V, A1, A2, A3, B1, B2, B3, C1, C2, C3):
    return ln_beta_V_function([T, V], A1, A2, A3, B1, B2, B3, C1, C2, C3)
    

def ln_beta_V_function(TV, A1, A2, A3, B1, B2, B3, C1, C2, C3):
    "A function used to parameterise 1000.ln(beta(T,V))"
    b = (A1+A2/TV[1]+A3/TV[1]**2)/TV[0]**6 + \
        (B1+B2/TV[1]+B3/TV[1]**2)/TV[0]**4 + \
        (C1+C2/TV[1]+C3/TV[1]**2)/TV[0]**2
    return b


def run_and_report(seedname, fineqpoints=None):

    print "Using 'Phonons' for frequency calculation "
    print "of 24Mg and 26 Mg substitution into {} ...\n".format(seedname)
    popt, pconv, max_error = fit_beta_func(seedname, fineqpoints)
    print "For function:\n  1000 ln(beta) = A/T^6 + B/T^4 + C/T^2"
    print "parameters are: \n  A = {:7g} \n   B = {:7g} \n   C = {:7g}".format(
                popt[0], popt[1], popt[2])
    print "maximum error is: {:7g}".format(max_error)
    # Convergence is to ~0.01 per mil, so worse than this is a problem
    assert max_error < 0.01, ValueError

    return (popt, pconv)


def get_data(paths_and_seeds):
    """Return data from phonons runs

       Data is returned as a dictionary indexed 
       by cell volume, values are 4-tuples of arrays 
       of frequencies and wights of the normal and 
       heavy isotopes
    """
    old_dir = os.getcwd()
    data = {}
    for path, seedname in paths_and_seeds:
        os.chdir(path)
        print "Extracting data from {} in {}".format(seedname, path)
        vol = get_volume(seedname)
        # For MgO we have a 2*2*2 primitive cell so 
        #vol = vol / 2.0
	(v, w, vs, ws) = castep_isotope_sub.get_freqs(seedname)
        data[vol] = (v, w, vs, ws)
        os.chdir(old_dir)
    return data


if __name__ == "__main__":
    import argparse
    import bulk_run_phonons
    parser = argparse.ArgumentParser(description=
        "Fit 1000.ln(beta(T,V)")
    parser.add_argument('castep_files', nargs='+', action='store',
                    metavar='Castep file', help='A list of .castep file paths')
    args = parser.parse_args()

    paths_and_seeds = bulk_run_phonons.process_paths(args.castep_files)

    data = get_data(paths_and_seeds)
        
    A1, A2, A3, B1, B2, B3, C1, C2, C3 = fit_beta_T_V(data)
    print 2500, 25, ln_beta_V_function_wrap(2600, 272.331, A1, A2, A3, B1, B2, B3, C1, C2, C3)
    print 3200, 25, ln_beta_V_function_wrap(3200, 280.110, A1, A2, A3, B1, B2, B3, C1, C2, C3)
    print 3000, 60, ln_beta_V_function_wrap(3000, 237.137, A1, A2, A3, B1, B2, B3, C1, C2, C3)
    print 4000, 60, ln_beta_V_function_wrap(4000, 241.042, A1, A2, A3, B1, B2, B3, C1, C2, C3)
    print 
    print 2500, 25, ln_beta_V_function_wrap(2600, 259.405506654, A1, A2, A3, B1, B2, B3, C1, C2, C3)
    print 3200, 25, ln_beta_V_function_wrap(3200, 259.405506654, A1, A2, A3, B1, B2, B3, C1, C2, C3)
    print 3000, 60, ln_beta_V_function_wrap(3000, 228.45872272, A1, A2, A3, B1, B2, B3, C1, C2, C3)
    print 4000, 60, ln_beta_V_function_wrap(4000, 228.45872272, A1, A2, A3, B1, B2, B3, C1, C2, C3)
