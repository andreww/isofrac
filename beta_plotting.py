#!/usr/bin/env python

import numpy as np
import castep_isotope_sub

if __name__ == "__main__":
    import readline
    import argparse

    parser = argparse.ArgumentParser(description=
              'Plot reduced partition function from parameters.')
    parser.add_argument('--filename', help='Ouput plot to a file.')
    parser.add_argument('--linestyle', help='Include style and color.',
                        action="store_true")
    args = parser.parse_args()


    names = []
    colors = []
    styles = []
    Ts = np.concatenate((np.linspace(300.0, 500.0, num=40), np.linspace(501.0, 4000.0, num=40)))
    betas = []

    done = False
    while not done:
       name = raw_input("Enter name (empty string when done): ")
       if name == "":
           done = True
           next
       elif name[0] =="*":
           name = name[1:]
           if args.linestyle:
               colors.append(raw_input("Color:"))
               styles.append(raw_input("Style:"))
           A = float(raw_input("Wu A parameter: "))
           B = float(raw_input("Wu B parameter: "))
           C = float(raw_input("Wu C parameter: "))
           names.append(name)
           betas.append(castep_isotope_sub.wu_ln_beta_function(Ts, A, B, C))
       else:
           if args.linestyle:
               colors.append(raw_input("Color:"))
               styles.append(raw_input("Style:"))
           A = float(raw_input("A parameter: "))
           B = float(raw_input("B parameter: "))
           C = float(raw_input("C parameter: "))
           names.append(name)
           betas.append(castep_isotope_sub.ln_beta_function(Ts, A, B, C))

    if args.linestyle:
        castep_isotope_sub.plot_beta(Ts, betas, names, 
            colors=colors, styles=styles, filename=args.filename) 
    else:
        castep_isotope_sub.plot_beta(Ts, betas, names, 
            filename=args.filename) 
