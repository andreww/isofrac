#!/usr/bin/env python
"""
bulk_run_phonons

Trivial script run pairs of phonon calculations 
from castep output with different isotopic masses
and leave the results around to allow isotopic 
fractionation to be calculated
"""
import os

import castep_isotope_sub

def process_paths(castep_files):
    """Given paths to a list of .castep files, 
       try to return the path and seedname used
    """ 
    paths_and_seeds = []
    for file in castep_files:
        (path, head) = os.path.split(file)
        (seedname, ext) = os.path.splitext(head)
        assert (seedname != ""), ValueError
        if (path == ""): path = "."
        paths_and_seeds.append((path, seedname))
    return paths_and_seeds


def walk_cleanup(paths_and_seeds):
    old_dir = os.getcwd()
    for path, seedname in paths_and_seeds:
        os.chdir(path)
        print("Cleaning up {} isotope calculations in {}".format(seedname, path))
        castep_isotope_sub.cleanup(seedname)
        os.chdir(old_dir)


def walk_phonons(paths_and_seeds):
    old_dir = os.getcwd()
    for path, seedname in paths_and_seeds:
        os.chdir(path)
        print("Running phonons for {} in {}".format(seedname, path))
        castep_isotope_sub.run_phonons(seedname)
        os.chdir(old_dir)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=
        "Phonon bulk runner for reduced partition functions from Castep.")
    parser.add_argument('castep_files', nargs='+', action='store',
                    metavar='Castep file', help='A list of .castep file paths')
    parser.add_argument('-c', '--cleanup', help='remove phonons output',
                   action="store_true")
    args = parser.parse_args()

    paths_and_seeds  = process_paths(args.castep_files)

    if args.cleanup:
        walk_cleanup(paths_and_seeds)
    else:
        walk_phonons(paths_and_seeds)



