from __future__ import division
import numpy as np
from scipy.integrate import quad, dblquad
from scipy.constants import parsec
from scipy.special import erfinv
import matplotlib.pyplot as plt
import healpy as hp
import sys
import os
import threading, subprocess
from math import cos, sin, erf, atan
import time
import os
import ctypes as C
import scipy.interpolate as interpolate
import subprocess
from scipy.stats import lognorm

AngRes = 0.5

def ymw16(DM, l, b):
    """
    Calcultes the distance to a source using the provided
    despersion measure. This function assumes the ymw16 model
    for electron density distribution in the galaxy.

    Requires an installed version of the ymw16 model in the local
    directory.
    
    A base error on the measured distance is 20% as well as the mean
    variation in distance calculated from the neighbouring pixels.

    Arguments
    ---------
    DM [float] : Dispersion measure # in pp/cm^3
    l [float] : longitude (galactic coordinates)
    b [float] : latitude (galactic coordinates)
    
    Returns
    -------
    Dist [float] : Distance to MSP # in kpc
    """
    try:
        p = subprocess.check_output(["./ymw16", "Gal", str(l), str(b),str(DM), "1"])
    except subprocess.CalledProcessError as e:
        if e.returncode == 53 or e.returncode == 137:
            p = e.output
        else:
            raise e
    Dist = float(p[105:114])/1e3
    return Dist

def dist_pdf(DM, l, b, n_dbins, plots=False):
    """
    Calcultes the distance to a source using the provided
    despersion measure. This function assumes the ymw16 model
    for electron density distribution in the galaxy.

    Requires an installed version of the ymw16 model in the local
    directory.
    
    The pdf is calculated by varying the parameters in the ymw16 model
    and recalculating the distance to the MSP which is the biggest source
    of systematics in the distance calculation

    Arguments
    ---------
    DM [array] : Dispersion measure # in pp/cm^3
    l [array] : longitude (galactic coordinates)
    b [array] : latitude (galactic coordinates)
    
    Returns
    -------
    Dist [array, array] : PDFs of the distance to MSPs # in kpc
    """
    filename = 'ymw16par.txt'
    names = np.loadtxt(filename, usecols=(0,), dtype='str')
    vals = np.loadtxt(filename, usecols=(1,))

    n1_list = np.linspace(0.008, 0.016, 18)
    H1_list = np.linspace(1200, 2000, 18)
    dist_pdfs = []

    for DM, l, b in zip(DM, l, b):
        D_list = []
        l = float(l)
        b = float(b)
        DM = float(DM)
        for n1 in n1_list:
            vals[names=='n1'] = n1
            for H1 in H1_list:
                vals[names=='H1'] = H1
                np.savetxt(filename, zip(names, vals), 
                    delimiter=" ", fmt="%s") 
                D_list.append(float(Distance(DM=DM, 
                    l=l, b=b, ne2001=False)))
        d_hist = np.histogram(D_list, bins=n_dbins, density=True)
        if plots:
            d_bins = np.linspace(0,max(D_list), num=n_dbins)
            plt.figure()
            plt.hist(D_list, bins=d_bins, normed=True)
            plt.ylabel('N')
            plt.xlabel('D [kpc]')
            plt.savefig('%s_test.pdf'%(str(DM)))
            plt.close()
        dist_pdfs.append(d_hist[0])

    return dist_pdfs

# def Distance(DM, l, b, ne2001=True):
#     """
#     Short wrapper for the above to functions
#     """
#     D = ymw16(DM, l, b)
#     return D


# DM = np.array([4., 2., 3.])
# l = np.array([10., 20., 30.])
# b = np.array([15., 25., 35.])
# newthing = dist_pdf(DM, l, b, 20, plots=True)
# print newthing