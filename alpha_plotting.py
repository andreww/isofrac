#!/usr/bin/env python

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import castep_isotope_sub

def plot_alpha(Ts, betas_ref, betas_plot, name_ref, names_plot):

    fix, ax1 = plt.subplots()

    styles = ['k--', 'k-', 'b--', 'b-', 'g--', 'g-', 'r--', 'r-']
    i = 0

    print()
    for pressure, beta_ref in betas_ref.items(): 
        for bet_plot, name_plot in zip(betas_plot[pressure], names_plot[pressure]):
            alpha = beta_ref - bet_plot
            print("pressure", pressure, "name", name_plot, "style", styles[i]) 
            ax1.plot(Ts, alpha, styles[i], label=name_plot + str(pressure))
            i = i + 1

    ax1.set_ylabel(r"$\Delta^{}$Mg (per mill) relative to {}".format('{26}', name_ref))
    ax1.set_xlabel("T (K)")
    #plt.legend()
    #plt.show()
    plt.savefig('alpha_4.eps')


if __name__ == "__main__":
    import readline

    Ts = np.linspace(1000.0, 4000.0, num=500)
    betas_p = {}
    betas_p_names = {}
    betas_ref = {}

    name_ref = input("Enter name of reference phase: ")
    done = False
    while not done:
       pressure = input("Pressure: ")
       if pressure == "":
           done = True
           print("Done reference")
           continue
       A = float(input("A parameter of reference: "))
       B = float(input("B parameter of reference: "))
       C = float(input("C parameter of reference: "))
       betas_ref[pressure] = castep_isotope_sub.ln_beta_function(Ts, A, B, C)

    Pdone = False
    while not Pdone:
        pressure = input("Pressure: ")
        if pressure == "":
           Pdone = True
           continue

        names = []
        betas = []
        done = False
        while not done:
            name = input("Enter name (empty string when done): ")
            if name == "":
                done = True
                continue
            A = float(input("A parameter: "))
            B = float(input("B parameter: "))
            C = float(input("C parameter: "))
            names.append(name)
            betas.append(castep_isotope_sub.ln_beta_function(Ts, A, B, C))
        betas_p_names[pressure] = names
        betas_p[pressure] = betas

    plot_alpha(Ts, betas_ref, betas_p, name_ref, betas_p_names) 
