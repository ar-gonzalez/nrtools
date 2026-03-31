#!/usr/bin/python
# ConvertEOSTable_SGRID2LORENE.py

"""
Convert EOS tables from SGRID to LORENE format
"""

__author__ = 'AG'
__date__   = '03/2026'

import sys
import numpy as np


def unit_conv_fact():
    """
    Return conversion factors (same as original script)
    """
    return dict(
        Length_cgs = 1.476625038500000e+05,
        Time_cgs   = 4.925490949141889e-06,
        Mass_cgs   = 1.988920000000000e+33,
        Energy_cgs = 1.787552150093231e+54,
        Surface_cgs= 2.180421504325126e+10,
        Volume_cgs = 3.219664987770316e+15,
        Press_cgs  = 5.551981826938919e+38,
        Mdens_cgs  = 6.177412890952259e+17,
        Edens_cgs  = 5.551981826938919e+38,
        Amom_cgs   = 8.804571936403334e+48,
        Length_fm  = 1.476999442301651e+18,
        Length_km  = 1.476625038500000e+00,
        Density_kgm3 = 6.177412890952258e+20,
        Energy_MeV = 1.115702344637173e+60,
        Volume_fm3 = 3.219664987770317e+54,
        fm_to_cm   = 1e-13,
        barmass_g  = 1.66e-24
    )


def read_sgrid_file(fname):
    """
    Read SGRID EOS file (skip comment lines)
    """
    return np.loadtxt(fname, comments='#')


def write_lorene_eosfile(fname, data):
    """
    Write EOS file in LORENE format
    """
    n = len(data)

    with open(fname, 'w') as fh:
        fh.write('#\n#\n#\n')
        fh.write(f'{n} <-- Number of lines\n')
        fh.write('#\n')
        fh.write('#        n_B [fm^{-3}]  rho [g/cm^3]   p [dyn/cm^2]\n')
        fh.write('#\n')

    with open(fname, 'a') as fh:
        np.savetxt(fh, data, fmt='%d %.16e %.16e %.16e')


def convert_sgrid2lorene(datain):
    """
    Convert SGRID data to LORENE
    Input columns:
        rho0, epsl, P
    Output columns:
        i, nB, rho, p
    """
    units = unit_conv_fact()

    npts = len(datain)
    dataout = np.zeros((npts, 4))

    rho0 = datain[:, 0]
    epsl = datain[:, 1]
    P    = datain[:, 2]

    # ---- n_B from rho0 ----
    # rho0 = (mB * nB / fm_to_cm^3) / Mdens_cgs
    nB = rho0 * units["Mdens_cgs"] * (units["fm_to_cm"]**3) / units["barmass_g"]

    # ---- rho from epsl ----
    # epsl = (rho/Mdens_cgs)/rho0 - 1
    rho = (epsl + 1.0) * rho0 * units["Mdens_cgs"]

    # ---- pressure ----
    p = P * units["Press_cgs"]

    # fill array
    dataout[:, 0] = np.arange(npts)  # index i
    dataout[:, 1] = nB
    dataout[:, 2] = rho
    dataout[:, 3] = p

    return dataout


if __name__ == '__main__':

    nf = len(sys.argv) - 1

    if nf < 1:
        print("Usage:\n " + sys.argv[0] + " <input SGRID files>")
        sys.exit()

    print("Files to process:", sys.argv[1:])

    for f in sys.argv[1:]:

        data = read_sgrid_file(f)
        convdata = convert_sgrid2lorene(data)

        name = f.rsplit(".", 1)[0] + "_lorene.txt"

        write_lorene_eosfile(name, convdata)

        print(f"{f}  =>  {name}")
        print("n={} col={}".format(convdata.shape[0], convdata.shape[1]))