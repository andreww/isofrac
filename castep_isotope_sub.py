#!/usr/bin/env python
"""
castep_isotope_sub

Trivial script to do castep isotope calculation on top of phonons
and fit the results
"""
import os
import glob
import subprocess
import numpy as np
import scipy.optimize as spopt

import calc_beta

def produce_dotcell(seedname, mass, fineqpoints=None):
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
    if fineqpoints is not None:
        outputfile.write(
            "phonon_fine_kpoint_mp_grid {0[0]:d} {0[1]:d} {0[2]:d}\n".format(
            fineqpoints))
        # FIXME - should deal with this
        # phonon_fine_kpoint_mp_offset x y z
    outputfile.close()
    return()

def run_phonons(seedname, fineqpoints=None):
    """Does phonon calculation for seedname with 24Mg and 26Mg"""
       
    phonons_path = "/nfs/see-fs-02_users/earawa/lvs/Code/castep/CASTEP-6.11/bin/linux_x86_64_ifort12/phonons_AMW"

    # Setup and run for light isotope
    produce_dotcell(seedname, "24Mg", fineqpoints)
    os.symlink(seedname+".param", seedname+"__isotope_l.param")
    subprocess.check_call([phonons_path, seedname+"__isotope_l"])

    # Setup and run for heavy isotope
    produce_dotcell(seedname, "26Mg", fineqpoints)
    os.symlink(seedname+".param", seedname+"__isotope_h.param")
    subprocess.check_call([phonons_path, seedname+"__isotope_h"])


def get_freqs(seedname, fineqpoints=None):
    """After phonon calculation for seedname with 24Mg and 26Mg
       returns the frequencies and weights"""

    if (os.path.isfile(seedname+"__isotope_l.phonon") and
        os.path.isfile(seedname+"__isotope_h.phonon")):
        print "Reusing existing .phonon files"
        if fineqpoints is not None:
            "phonon_fine_kpoint_mp_grid ignored"
    else:
        # First, run phonons
        if fineqpoints is not None:
            print "phonon_fine_kpoint_mp_grid: {0[0]:d} {0[1]:d} {0[2]}".format(
                fineqpoints)
        run_phonons(seedname, fineqpoints)

    # Pull out the frequencies and weights
    v, w, vol = calc_beta.read_frequences(seedname+"__isotope_l.phonon")
    vs, ws, vol = calc_beta.read_frequences(seedname+"__isotope_h.phonon")

    print "Total q-points: {0:d}\n".format(len(w))

    return(v, w, vs, ws)

def fit_beta_func(seedname, fineqpoints=None):

    # Get betas for T = 300 - 4000
    (v, w, vs, ws) = get_freqs(seedname, fineqpoints)
    Ts = np.linspace(300.0, 4000.0, num=4000)
    betas = calc_beta.beta_T(Ts, 1, v, vs, w, ws)
    lnbetas = 1000.0 * np.log(betas)

    # Fit to functional form for ln(beta)
    popt, pconv = spopt.curve_fit(ln_beta_function, Ts, 
        lnbetas, p0=[1E14, -1E10, 1E6])

    # Check error so that caller can check value
    calc_betas = ln_beta_function(Ts, popt[0], popt[1], popt[2])
    error = np.abs(lnbetas-calc_betas)
    max_error = np.max(error)

    return (popt, pconv, max_error)
    

def ln_beta_function(T, A, B, C):
    "A function used to parameterise 1000.ln(beta(T))"
    b = A/T**6 + B/T**4 + C/T**2
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


def plot_beta(Ts, betas, names=None, filename=None):

    import matplotlib
    if filename is not None:
        matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    Tsm1 = 1E6/(Ts**2.0)
    fix, ax1 = plt.subplots()
    if type(betas) is list or type(betas) is tuple:
        if names is None:
            for beta in betas:
                ax1.plot(Tsm1, beta)
        else:
            for beta, name in zip(betas, names):
                ax1.plot(Tsm1, beta, label=name)
    else:
        ax1.plot(Tsm1, betas, "b-")
    ax1.set_ylabel("1000 * ln beta (per mill)")
    ax1.set_xlabel("1000000 / T^2 (1^6 K^-2)")
    ax1.set_xlim(right=Tsm1.max())
    x1locs, x1labels = plt.xticks()

    if names is not None:
        plt.legend(loc=2)

    ax2 = ax1.twiny()
    ax2.set_xlabel("T (K)")
    x2vals = []
    x2locs = []
    max_pos = 1E6/(Ts.min()**2.0)
    min_pos = 1E6/(Ts.max()**2.0)
    for xloc in np.linspace(min_pos, max_pos, 6):
        thisval = np.sqrt((1.0/(xloc/1E6)))
        if thisval != float("inf"):
            x2vals.append("{:4.0f}".format(thisval))
            x2locs.append(xloc)
    plt.xticks(x2locs, x2vals)

    if filename is not None:
        plt.savefig(filename)
    else:
        plt.show()

def cleanup(seedname):

    print "Cleaning up..."
    for file in glob.glob(seedname + "__isotope_h.*") + \
            glob.glob(seedname + "__isotope_l.*"):
        os.unlink(file)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=
              'Reduced partition function from Castep.')
    parser.add_argument('seedname', metavar='seedname', 
                   action='store', help='Castep seedname')
    parser.add_argument('--plot', help='generate a matplotlib graph',
                   action="store_true")
    parser.add_argument('-c', '--cleanup', help='remove phonons output',
                   action="store_true")
    parser.add_argument('-f', '--fineqpoints', nargs=3, type=int,
                   help='Add fine q-points')
    args = parser.parse_args()

    seedname = args.seedname
    print "\nCalculation of reduced partition function"
    print "=========================================\n"
    (popt, pconv) = run_and_report(seedname, fineqpoints=args.fineqpoints)
    # Tabulate output
    print "\n T(K)        1000 ln beta"
    print " -------------------------"
    for t in [300, 500, 1000, 2600, 3700]:
        print " {:5f}     {:5f}".format(t, ln_beta_function(t, 
                   popt[0], popt[1], popt[2]))

    # Plot output
    Ts = np.linspace(300.0, 4000.0, num=40)
    betas = ln_beta_function(Ts, popt[0], popt[1], popt[2])
    if args.plot:
        plot_beta(Ts, betas)

    if args.cleanup:
        cleanup(seedname)


