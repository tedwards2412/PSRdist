from __future__ import division
import ctypes
import numpy as np
import matplotlib.pyplot as plt
import py_ymw16 as ymw
import os

# Get directory name
dirname = os.path.dirname(os.path.realpath(__file__))

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
    Dist = ymw.distances(l=l, b=b, DM=DM, dirname=dirname+"/ymw16_v1.3/")
    # Dist = ymw.distances(l=l, b=b, DM=DM)
    return np.array(Dist)


def dist_bf(name, DM, l, b):
    """
    Calcultes the distance to a source using the provided
    despersion measure. This function assumes the ymw16 model
    for electron density distribution in the galaxy with their best-fit
    parameters.

    Arguments
    ---------
    name [array] : List of names of sources for labelling
    DM [array] : Dispersion measure # in pp/cm^3
    l [array] : longitude (galactic coordinates)
    b [array] : latitude (galactic coordinates)

    Returns
    -------
    Dist [array] : Distance to MSPs # in kpc
    """
    # original best fit values
    filename_bf = dirname + '/ymw16_v1.3/ymw16par_bestfit.txt'
    names_bf = np.loadtxt(filename_bf, usecols=(0,), dtype='str')
    vals_bf = np.loadtxt(filename_bf, usecols=(1,))
    # --> remove incomplete Ele_arm and Wid_arm
    vals_bf = vals_bf[(names_bf!='Wid_arm') & (names_bf!='Ele_arm')]
    names_bf = names_bf[(names_bf!='Wid_arm') & (names_bf!='Ele_arm')]
    # --> have to load Ele_arm separately (since it has one column for each narm)
    with open(filename_bf) as fp:
        for line in fp.readlines():
          if "Ele_arm" in line:
              narm = np.array(line.split(" ")[1:], dtype='float')
          if "Wid_arm" in line:
              warm = np.array(line.split(" ")[1:], dtype='float')


    # reset the input file for the ymw16_v1.3 code and load the values
    filename_input = dirname + '/ymw16_v1.3/ymw16par.txt'
    np.savetxt(filename_input, zip(names_bf, vals_bf),
            delimiter=" ", fmt="%s")
    with open(filename_input, 'a') as f:
        f.write('Ele_arm %.6f %.6f %.6f %.6f %.6f\n'%tuple(narm))
        f.write('Wid_arm %i %i %i %i %i\n'%tuple(warm))
    print "Resetting input file: Done!"

    return ymw16(DM=DM, l=l, b=b)

# def dist_pdf(name, DM, l, b, n_dbins,
def dist_pdf(name, DM, l, b, d_edges=np.linspace(0, 10, 101),
    MC=True, MC_mode = "bestfit",
    save_files=False, plots=False,
    n_grid=10, n_MC=1000):
    """
    Calcultes the distance to a source using the provided
    despersion measure. This function assumes the ymw16 model
    for electron density distribution in the galaxy.

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
    MC_mode -- bestfit, uniform, gaussian# (# = error in %)
    MC -- Monte Carlo vs grid scan (default True)
    save_files -- saves pdfs as txt files (default False)
    plots -- pdf plots for each pulsar (default False)

    Returns
    -------
    dist_pdfs [array] : PDFs of the distance to MSPs
    dist_edges [array] : Edges of the pdf histogram # [kpc]
    """
    # original best fit values
    filename_bf = dirname + '/ymw16_v1.3/ymw16par_bestfit.txt'
    names_bf = np.loadtxt(filename_bf, usecols=(0,), dtype='str')
    vals_bf = np.loadtxt(filename_bf, usecols=(1,))
    # --> remove incomplete Ele_arm and Wid_arm
    vals_bf = vals_bf[(names_bf!='Wid_arm') & (names_bf!='Ele_arm')]
    names_bf = names_bf[(names_bf!='Wid_arm') & (names_bf!='Ele_arm')]
    # --> have to load Ele_arm separately (since it has one column for each narm)
    with open(filename_bf) as fp:
        for line in fp.readlines():
          if "Ele_arm" in line:
              narm = np.array(line.split(" ")[1:], dtype='float')
          if "Wid_arm" in line:
              warm = np.array(line.split(" ")[1:], dtype='float')


    # reset the input file for the ymw16_v1.3 code and load the values
    filename_input = dirname + '/ymw16_v1.3/ymw16par.txt'
    np.savetxt(filename_input, zip(names_bf, vals_bf),
            delimiter=" ", fmt="%s")
    with open(filename_input, 'a') as f:
        f.write('Ele_arm %.6f %.6f %.6f %.6f %.6f\n'%tuple(narm))
        f.write('Wid_arm %i %i %i %i %i\n'%tuple(warm))
    print "Resetting input file: Done!"

    names = np.loadtxt(filename_input, usecols=(0,), dtype='str')
    vals = np.loadtxt(filename_input, usecols=(1,))
    vals = vals[(names!='Wid_arm') & (names!='Ele_arm')]
    names = names[(names!='Wid_arm') & (names!='Ele_arm')]
    with open(filename_input) as fp:
        for line in fp.readlines():
          if "Ele_arm" in line:
              narm = np.array(line.split(" ")[1:], dtype='float')
          if "Wid_arm" in line:
              warm = np.array(line.split(" ")[1:], dtype='float')

    # Load the file containing parameter ranges and best fit values
    filename_ranges = dirname + '/ymw16_v1.3/ymw16par_ranges.txt'
    names_ranges = np.loadtxt(filename_ranges, usecols=(0,), dtype='str')
    val_lo, val_hi, val_bf, val_unc = np.loadtxt(
        filename_ranges, usecols=(1,2,3,4)).T
    narm_unc = np.loadtxt('%s/ymw16_v1.3/ymw16par_Ele_arm_unc.txt'%dirname)
    narm_lo, narm_hi = np.loadtxt('%s/ymw16_v1.3/ymw16par_Ele_arm_ranges.txt'%dirname).T

    # Create list
    D_list = []
    dist_pdfs = []
    dist_edges = []
    if MC:
        tag="_MC_%s"%MC_mode

        # Mask that says which parameters to vary
        mask1 = (val_lo != val_hi)
        print "Number of free parameters: %i"%mask1.sum()

        # Match the MC values to input file values
        ind_order = np.array([(n==names).argmax() for n in names_ranges[mask1]])

        if MC_mode == 'bestfit':
            # Draw from a gaussian distribution with an uncertainty
            # corresponding to the uncertainty in table 2 of YMW16
            val_MC = np.random.normal(loc=val_bf[mask1], scale=val_unc[mask1],
                size=(n_MC, mask1.sum()))
            val_narm_MC = np.random.normal(loc=narm, scale=narm_unc,
                size=(n_MC, len(narm)))

        elif MC_mode[:8]=='gaussian':
            # Draw from a gaussian distribution with an uncertainty
            # of a given percentage (appended after gaussian)
            sig = float(MC_mode[8:])/100.
            val_MC = np.random.normal(loc=val_bf[mask1],
                scale=val_bf[mask1] * sig,
                size=(n_MC, mask1.sum()))
            val_narm_MC = np.random.normal(loc=narm, scale=narm * sig,
                size=(n_MC, len(narm)))

        elif MC_mode == 'uniform':
            # Draw from a uniform distribution
            val_MC = np.random.uniform(low=val_lo[mask1], high=val_hi[mask1],
                size=(n_MC, mask1.sum()))
            val_narm_MC = np.random.uniform(low=narm_lo, high=narm_hi,
                size=(n_MC, len(narm)))


        # Run MC
        for i in range(n_MC):
            # update values
            vals[ind_order] = val_MC[i] # Set values
            np.savetxt(filename_input, zip(names, vals),
                delimiter=" ", fmt="%s")
            with open(filename_input, 'a') as f:
                f.write('Ele_arm %.6f %.6f %.6f %.6f %.6f\n'%tuple(val_narm_MC[i]))
                f.write('Wid_arm %i %i %i %i %i\n'%tuple(warm))
            D_list.append(ymw16(DM=DM, l=l, b=b))
        D_list = np.array(D_list)


    else: # Perform a grid scan
        tag="_grid"
        n1_list = np.linspace(0.008, 0.016, n_grid)
        H1_list = np.linspace(1200, 2000, n_grid)
        detlb2_list = np.linspace(10, 60, 11)

        for n1 in n1_list:
            vals[names=='n1'] = n1
            for H1 in H1_list:
                vals[names=='H1'] = H1

                for detlb2 in detlb2_list:
                    vals[names=='detlb2'] = detlb2
                    np.savetxt(filename_input, zip(names, vals),
                        delimiter=" ", fmt="%s")
                    D_list.append(ymw16(DM=DM, l=l, b=b))
        D_list = np.array(D_list)

    # --> Make histograms
    for i in range(DM.size):
        d_hist, d_edges = np.histogram(D_list[:,i],
            bins=d_edges, density=True)
            # bins=n_dbins, density=True)
        dist_pdfs.append(d_hist)
        dist_edges.append(d_edges)
    dist_pdfs = np.array(dist_pdfs)
    dist_edges = np.array(dist_edges)

    # --> make plots
    if plots:
        for i in range(DM.size):
            print tag
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
            # np.savetxt("output/%s_%s.dat"%(name[i], tag),
            #     zip(dist_pdfs[i][0], dist_pdfs[i][1]),
            #         delimiter=" ")
            np.savetxt("output/%s_pdf_%s.dat"%(name[i], tag),
                dist_pdfs[i], delimiter=" ")
            np.savetxt("output/%s_Dedges_%s.dat"%(name[i], tag),
                dist_edges[i], delimiter=" ")

    return dist_pdfs, dist_edges

# names = np.array(["B0628-28", "B1237+25"])
# DM = np.array([34.36, 9.2755])
# l = np.array([236.95, 252.44])
# b = np.array([-16.75, 86.54])
#
# newthing = dist_pdf(names, DM, l, b,
#     n_dbins=40, n_grid=181, n_MC=10000,
#     plots=True, save_files=True)
#
# newthing = dist_pdf(names, DM, l, b,
#     n_dbins=40, n_grid=90,
#     plots=True, save_files=True, MC=False)
