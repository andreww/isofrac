#!/usr/bin/env python

import numpy as np
#import matplotlib
#matplotlib.use('agg')

import matplotlib.pyplot as plt

import castep_isotope_sub

def plot_beta(Ts, betas):

    Tsm1 = 1E6/(Ts**2.0)
    print Tsm1
    fix, ax1 = plt.subplots()
    ax1.plot(Tsm1, betas, "b-")
    ax1.set_ylabel("1000 * ln beta (per mill)")
    ax1.set_xlabel("1000000 / T^2 (1^6 K^-2)")
    ax1.set_xlim(right=Tsm1.max())
    x1locs, x1labels = plt.xticks()

    ax2 = ax1.twiny()
    #ax2.plot(Ts, betas, "b-")
    ax2.set_xlabel("T (K)")
    x2vals = []
    x2locs = []
    max_pos = 1E6/(Ts.min()**2.0)
    print max_pos
    min_pos = 1E6/(Ts.max()**2.0)
    for xloc in np.linspace(min_pos, max_pos, 6):
        thisval = np.sqrt((1.0/(xloc/1E6)))
        if thisval != float("inf"):
            x2vals.append("{:4.0f}".format(thisval))
            x2locs.append(xloc)
    print x2locs
    plt.xticks(x2locs, x2vals)
    print plt.xticks()
    plt.show()

if __name__ == "__main__":
    import sys
    seedname = sys.argv[1]
    (popt, pconv) = castep_isotope_sub.fit_beta_func(seedname)
    Ts = np.linspace(300.0, 4000.0, num=40)
    betas = castep_isotope_sub.beta_function(Ts, popt[0], popt[1], popt[2])
    print Ts
    print betas
    plot_beta(Ts, betas)

