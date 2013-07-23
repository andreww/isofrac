#!/usr/bin/env python

import re
import math as m

def read_frequences(filename):
    """Return data from <seedname>.phonon

    The CASTEP and related PHONON codes both
    generate a file containing phonon frequences and 
    related information. This function returns a 
    list of frequencies from a file (assumed to be
    in .phonon file format"""

    # phonon frequences are on lines by themselves 
    # as integers followed by reals. Only case in
    # the file that is like this
    # But we may also (for gamma point) end up with 
    # flippin IR activities too. We ignore these if 
    # we find em.
    get_freq_RE = re.compile(r"^\s+\d+\s+([\+\-]?\d+\.\d+)(?:\s+[\+\-]?\d+\.\d+)?\s*$",
                     re.MULTILINE)

    fh = open(filename, 'r')
    freq_grps = get_freq_RE.findall(fh.read())
    fh.close
    freqs = []
    for freq in freq_grps:
        freqs.append(float(freq))
    return freqs

def beta(T, N, Nq, freq, freqstar):

    h = 4.135667516E-15 # eV.s
    k = 8.6173324E-3 # eV/K
    cm2ev = 1.23984E-4 # *cm^-1 to give eV
    cm2Hz = 0.03E12 # *cm^-1 to give Hz (1/s)

    assert len(freq)==len(freqstar)

    beta = 1.0
    for vs, v in zip(freqstar, freq):
        if v <= 0.0:
            continue
        if vs <= 0.0:
            continue
        vs = vs*cm2Hz
        v = v*cm2Hz
        evs = m.exp((-1.0*h*vs)/(2.0*k*T))
        ev = m.exp((-1.0*h*v)/(2.0*k*T))
        beta = beta * (vs/v) * (evs / (1.0-evs)) * ((1.0-ev)/ev)

    beta = beta**(1.0/(Nq*N))

    return beta
    

if __name__ == "__main__":
    import sys
    v = read_frequences(sys.argv[1])
    vs = read_frequences(sys.argv[2])

    for T in [10, 100, 300, 500, 1000, 1500, 2000, 2500]:
        b = beta(T, 4, 4, v, vs)
        print T, b, m.log(b)*1E3
        

