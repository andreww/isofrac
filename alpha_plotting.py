#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import castep_isotope_sub

def plot_alpha(Ts, betas_ref, betas_plot, name_ref, names_plot):

    fix, ax1 = plt.subplots()

    for bet_plot, name_plot in zip(betas_plot, names_plot):
        alpha = betas_ref - bet_plot
        ax1.plot(Ts, alpha, label=name_plot)

    ax1.set_ylabel("1000 * ln alpha (per mill) relative to {}".format(name_ref))
    ax1.set_xlabel("T (K)")
    plt.legend()
    plt.show()


if __name__ == "__main__":
    import readline

    names = []
    Ts = np.linspace(1000.0, 4000.0, num=500)
    betas = []

    name_ref = raw_input("Enter name of reference phase: ")
    A = float(raw_input("A parameter of reference: "))
    B = float(raw_input("B parameter of reference: "))
    C = float(raw_input("C parameter of reference: "))
    betas_ref = castep_isotope_sub.beta_function(Ts, A, B, C)

    done = False
    while not done:
       name = raw_input("Enter name (empty string when done): ")
       if name == "":
           done = True
           next
       else:
           A = float(raw_input("A parameter: "))
           B = float(raw_input("B parameter: "))
           C = float(raw_input("C parameter: "))
           names.append(name)
           betas.append(castep_isotope_sub.beta_function(Ts, A, B, C))

    plot_alpha(Ts, betas_ref, betas, name_ref, names) 
