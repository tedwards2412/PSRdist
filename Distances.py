from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import ctypes as C
import subprocess
import py_ymw16 as ymw

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
    DM [array] : Dispersion measure # in pp/cm^3
    l [array] : longitude (galactic coordinates)
    b [array] : latitude (galactic coordinates)
    
    Returns
    -------
    Dist [array] : Distance to MSP # in pc
    """
    # try:
    #     p = subprocess.check_output(["./ymw16", "Gal", str(l), str(b),str(DM), "1"])
    # except subprocess.CalledProcessError as e:
    #     if e.returncode == 53 or e.returncode == 137:
    #         p = e.output
    #     else:
    #         raise e
    # Dist = float(p[105:114])/1e3
    Dist = ymw.distances(l=l, b=b, DM=DM)
    return Dist

def dist_pdf(name, DM, l, b, n_dbins, MC=True, save_files=False, plots=False):
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
    Dist [array, array] : PDFs of the distance to MSPs # in pc
    """
    filename = 'ymw16_v1.3/ymw16par.txt'
    names = np.loadtxt(filename, usecols=(0,), dtype='str')
    vals = np.loadtxt(filename, usecols=(1,))

    n1_list = np.linspace(0.008, 0.016, 18)
    H1_list = np.linspace(1200, 2000, 18)
    dist_pdfs = []

    D_list = []
    for n1 in n1_list:
        vals[names=='n1'] = n1
        for H1 in H1_list:
            vals[names=='H1'] = H1
            np.savetxt(filename, zip(names, vals), 
                delimiter=" ", fmt="%s") 
            D_list.append(ymw16(DM=DM, l=l, b=b))
    D_list = np.array(D_list)
    for i in range(DM.size):
        d_hist = np.histogram(D_list[:,i], bins=n_dbins, density=True)
        dist_pdfs.append(d_hist)

    if plots:
        for i in range(DM.size):
            print i
            plt.figure()
            plt.hist(D_list[:,i], bins=n_dbins, normed=True)
            plt.ylabel('N')
            plt.title('%s_test.pdf'%(str(name[i])))
            plt.xlabel('D [kpc]')
            plt.savefig('%s_test.pdf'%(str(name[i])))
            plt.close()

    if save_files:
        print "hello"
        # Do some saving

    return dist_pdfs
names = np.array(["B0628-28", "B1237+25", "B1237+201"])
DM = np.array([4., 2., 3.])
l = np.array([10., 20., 30.])
b = np.array([15., 25., 35.])
newthing = dist_pdf(names, DM, l, b, 20, plots=True, save_files=True)
print newthing