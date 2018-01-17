from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import ctypes as C
import subprocess
import py_ymw16 as ymw

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
    Dist [array] : Distance to MSP # in kpc
    """
    Dist = ymw.distances(l=l, b=b, DM=DM)
    return Dist

def dist_pdf(name, DM, l, b, n_dbins,
    MC=True, save_files=False, plots=False,
    n_grid=10, n_MC=1000):
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
    name [array] : List of names of sources for labelling
    DM [array] : Dispersion measure # in pp/cm^3
    l [array] : longitude (galactic coordinates)
    b [array] : latitude (galactic coordinates)
    n_grid [int] : number of grid points for each variable
    n_MC [float] : number of monte carlo samples

    Keyword arguments:
    MC -- Monte Carlo vs grid scan (default True)
    save_files -- saves pdfs as txt files (default False)
    plots -- pdf plots for each pulsar (default False)

    Returns
    -------
    Dist [array, array] : PDFs of the distance to MSPs # in kpc
    """
    # original best fit values
    filename_bf = 'ymw16_v1.3/ymw16par_bestfit.txt'
    names_bf = np.loadtxt(filename_bf, usecols=(0,), dtype='str')
    vals_bf = np.loadtxt(filename_bf, usecols=(1,))

    # Load the input file for the ymw16_v1.3 code
    filename_input = 'ymw16_v1.3/ymw16par.txt'
    names = np.loadtxt(filename_input, usecols=(0,), dtype='str')
    vals = np.loadtxt(filename_input, usecols=(1,))
    if (names==names_bf).all():
        vals = vals_bf
        print "Resetting input file... Done!"
    else:
        print "Cannot reset input file, quitting"
        quit()

    # Load the file containing parameter ranges and best fit values
    filename_ranges = 'ymw16_v1.3/ymw16par_ranges.txt'
    names_ranges = np.loadtxt(filename_ranges, usecols=(0,), dtype='str')
    val_lo, val_hi, val_bf, val_unc = np.loadtxt(
        filename_ranges, usecols=(1,2,3,4)).T

    # Create list
    D_list = []
    dist_pdfs = []
    if MC:
        tag="_MC"

        # Mask that says which parameters to vary
        mask1 = (val_lo != val_hi)
        print "Number of free parameters: %i"%mask1.sum()

        # Draw from a uniform distribution
        val_MC = np.random.uniform(low=val_lo[mask1], high=val_hi[mask1],
            size=(n_MC, mask1.sum()))

        # Match the MC values to input file values
        ind_order = np.array([(n==names).argmax() for n in names_ranges[mask1]])

        # Run MC
        for i in range(n_MC):
            # update values
            vals[ind_order] = val_MC[i] # Set values
            np.savetxt(filename_input, zip(names, vals),
                delimiter=" ", fmt="%s")
            D_list.append(ymw16(DM=DM, l=l, b=b))
        D_list = np.array(D_list)


    else: # Perform a grid scan
        tag="_grid"
        n1_list = np.linspace(0.008, 0.016, n_grid)
        H1_list = np.linspace(1200, 2000, n_grid)
        # ngn_list = np.linspace(1, 3, 3)
        # Wgn_list = np.linspace(10, 20, 3)
        # Agn_list = np.linspace(120, 130, 3)
        # nLI_list = np.linspace(0, 3, 4)
        nlb1_list = np.linspace(0.5, 1.5, 3)
        nlb2_list = np.linspace(1, 3, 3)

        for n1 in n1_list:
            vals[names=='n1'] = n1
            for H1 in H1_list:
                vals[names=='H1'] = H1
                # for ngn in ngn_list:
                    # vals[names=='ngn'] = ngn
                # for Wgn in Wgn_list:
                #     vals[names=='Wgn'] = Wgn
                # for Agn in Agn_list:
                #     vals[names=='Agn'] = Agn
                # for nLI in nLI_list:
                #     vals[names=='nLI'] = nLI
                for nlb1 in nlb1_list:
                    vals[names=='nlb1'] = nlb1
                    for nlb2 in nlb2_list:
                        vals[names=='nlb2'] = nlb2
                        np.savetxt(filename_input, zip(names, vals),
                            delimiter=" ", fmt="%s")
                        D_list.append(ymw16(DM=DM, l=l, b=b))
        D_list = np.array(D_list)

    # --> Make histograms
    for i in range(DM.size):
        d_hist = np.histogram(D_list[:,i], bins=n_dbins, density=True)
        dist_pdfs.append(d_hist)

    # --> make plots
    if plots:
        for i in range(DM.size):
            plt.figure()
            plt.hist(D_list[:,i], bins=n_dbins, normed=True)
            plt.ylabel('N')
            plt.title('%s'%(str(name[i])))
            plt.xlabel('D [kpc]')
            plt.savefig('plots/%s%s.pdf'%(name[i], tag))
            plt.close()

    # --> save files
    if save_files:
        for i in range(DM.size):
            np.savetxt("output/%s%s.dat"%(name[i], tag),
                zip(dist_pdfs[i][0], dist_pdfs[i][1]),
                    delimiter=" ")

    return dist_pdfs

names = np.array(["B0628-28", "B1237+25"])
DM = np.array([34.36, 9.2755])
l = np.array([236.95, 252.44])
b = np.array([-16.75, 86.54])
# newthing = dist_pdf(names, DM, l, b,
#     n_dbins=40, n_grid=181, n_MC=10000,
#     plots=True, save_files=True)
newthing = dist_pdf(names, DM, l, b,
    n_dbins=40, n_grid=40,
    plots=True, save_files=True, MC=False)
