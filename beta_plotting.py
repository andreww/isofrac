#!/usr/bin/env python

import numpy as np
import castep_isotope_sub

if __name__ == "__main__":
    import readline
    import argparse

    parser = argparse.ArgumentParser(description=
              'Plot reduced partition function from parameters.')
    parser.add_argument('--filename', help='Ouput plot to a file.')
    args = parser.parse_args()


    names = []
    Ts = np.linspace(300.0, 4000.0, num=40)
    betas = []

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

    castep_isotope_sub.plot_beta(Ts, betas, names, filename=args.filename) 
