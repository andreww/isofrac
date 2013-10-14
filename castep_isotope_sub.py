#!/usr/bin/env python
"""
castep_isotope_sub

Trivial script to do castep isotope calculation on top of phonons
and fit the results
"""
import os
import subprocess
import numpy as np
import scipy.optimize as spopt

import calc_beta

def produce_dotcell(seedname, mass):
    """
    produce_dotcell: reads <seedname>.cell (CASTEP cell file)
    and writes a new .cell file to <filename> appending
    %block SPECIES_MASS block with one designed for heavy 
    or light Mg isotope substitution (mass must be 24Mg or
    26Mg).
    """

    if mass == "24Mg":
        isotope_line = "    Mg:iso  24.0\n"
        filename = seedname+"__isotope_l.cell"
    elif mass == "26Mg":
        isotope_line = "    Mg:iso  26.0\n"
        filename = seedname+"__isotope_h.cell"
    else:
        raise NotImplementedError("Only 24Mg and 26Mg can be used")

    inputfile = open(seedname+".cell", "r")
    outputfile = open(filename, "w")
    for line in inputfile:
        outputfile.write(line)
    inputfile.close
    outputfile.write("%block SPECIES_MASS\n")
    outputfile.write("    O   15.9994000\n")
    outputfile.write("    Si   28.0855000\n")
    outputfile.write("    Mg   24.3050000\n")
    outputfile.write(isotope_line)
    outputfile.write("%endblock  SPECIES_MASS\n")
    outputfile.close
    return()

def run_phonons(seedname):
    """Does phonon calculation for seedname with 24Mg and 26Mg
       and returns the frequencies and weights"""

    phonons_path = "/share/apps/atomistic/CASTEP-6.11-serial/" + \
        "linux_x86_64_gfortran--serial/phonons"

    # Setup and run for light isotope
    produce_dotcell(seedname, "24Mg")
    os.symlink(seedname+".param", seedname+"__isotope_l.param")
    subprocess.check_call([phonons_path, seedname+"__isotope_l"])

    # Setup and run for heavy isotope
    produce_dotcell(seedname, "26Mg")
    os.symlink(seedname+".param", seedname+"__isotope_h.param")
    subprocess.check_call([phonons_path, seedname+"__isotope_h"])

    # Pull out the frequencies and weights
    v, w = calc_beta.read_frequences(seedname+"__isotope_l.phonon")
    vs, ws = calc_beta.read_frequences(seedname+"__isotope_h.phonon")

    return(v, w, vs, ws)

def fit_beta_func(seedname):

    # Get betas for T = 300 - 4000
    (v, w, vs, ws) = run_phonons(seedname)
    Ts = np.linspace(300.0, 4000.0, num=4000)
    betas = calc_beta.beta_T(Ts, 1, v, vs, w, ws)
    lnbetas = 1000.0 * np.log(betas)

    popt, pconv = spopt.curve_fit(beta_function, Ts, lnbetas, p0=[1E14, -1E10, 1E6])
    return (popt, pconv)
    

def beta_function(T, A, B, C):
    b = A/T**6 + B/T**4 + C/T**2
    return b

if __name__ == "__main__":
    import sys
    seedname = sys.argv[1]
    print fit_beta_func(seedname)


