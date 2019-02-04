#!/usr/bin/env python

import re
import math as m
import numpy as np

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

    get_freq_weights_RE = re.compile(r"^\s+q-pt=\s+\d+\s+[\+\-]?\d+\.\d+\s+[\+\-]?\d+\.\d+\s+[\+\-]?\d+\.\d+\s+(\d+\.\d+)\s*$", re.MULTILINE)

    get_lattice_vecs_RE = re.compile(r"^\s+Unit cell vectors \(A\)\n\s*(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s*\n\s*(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s*\n\s*(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)", re.MULTILINE)

    fh = open(filename, 'r')
    filelines = fh.read()
    freq_grps = get_freq_RE.findall(filelines)
    wgt_grps = get_freq_weights_RE.findall(filelines)
    lvec_grps = get_lattice_vecs_RE.findall(filelines)[0]
    fh.close
    wgts = []
    for wgt in wgt_grps:
        wgts.append(float(wgt))
    freqs = []
    for freq in freq_grps:
        freqs.append(float(freq))

    lvec = np.array([[float(lvec_grps[0]), float(lvec_grps[1]), float(lvec_grps[2])],
                     [float(lvec_grps[3]), float(lvec_grps[4]), float(lvec_grps[5])],
                     [float(lvec_grps[6]), float(lvec_grps[7]), float(lvec_grps[8])]]
                    )
    vol = np.linalg.det(lvec)
    return freqs, wgts, vol

def beta(T, N, freq, freqstar, wgt, wgtstar):

    h = 4.135667516E-15 # eV.s
    k = 8.6173324E-5 # eV/K
    cm2ev = 1.23984E-4 # *cm^-1 to give eV
    cm2Hz = 0.03E12 # *cm^-1 to give Hz (1/s)

    assert len(freq)==len(freqstar)
    assert wgt == wgtstar
 
    N_qpt = len(wgt) 
    N_fr = len(freq)//N_qpt
    i = 0
    beta = 1.0
    for Nqwt in wgt:
        this_bt = 1.0
        for vs, v in zip(freqstar[i*N_fr:i*N_fr+N_fr], freq[i*N_fr:i*N_fr+N_fr]):
            if v <= 0.0:
               continue
            if vs <= 0.0:
                continue
            vs = vs*cm2Hz
            v = v*cm2Hz
            evs = m.exp((-1.0*h*vs)/(2.0*k*T))
            evsb = m.exp((-1.0*h*vs)/(k*T))
            ev = m.exp((-1.0*h*v)/(2.0*k*T))
            evb = m.exp((-1.0*h*v)/(k*T))
            this_bt = this_bt * (vs/v) * (evs / (1.0-evsb)) * ((1.0-evb)/ev)
        beta = beta*this_bt**Nqwt
        i = i + 1

    beta = beta**(1.0/N)

    return beta

def beta_T(Ts, N, freq, freqstar, wgt, wgtstar):

    betas = np.zeros_like(Ts)
    i = 0
    for T in Ts:
        betas[i] = beta(T, N, freq, freqstar, wgt, wgtstar)
        i = i + 1

    return betas

if __name__ == "__main__":
    import sys
    v, w, vol = read_frequences(sys.argv[1])
    vs, ws, vol = read_frequences(sys.argv[2])

    for T in [15, 30, 60, 120, 240, 300, 500, 670, 1000, 1500, 2000, 2500, 2600, 3000, 3500, 3700, 4000]:
        b = beta(T, 1, v, vs, w, ws)
        print(T, b, m.log(b)*1E3)
        

    # Or, the 'vectorised' version...
    Ts = np.array([15.0, 30.0, 60.0, 120.0, 240.0, 300.0, 500.0, 670.0, 1000.0, 1500.0, 2000.0, 2500.0, 2600.0, 3000.0, 3500.0, 3700.0, 4000.0])
    betas = beta_T(Ts, 1, v, vs, w, ws)
    print(Ts)
    print(betas)
    print(np.log(betas)*1E3)
