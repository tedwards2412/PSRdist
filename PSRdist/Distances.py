from __future__ import division
import ctypes
import numpy as np
import matplotlib.pyplot as plt
import py_ymw16 as ymw
from scipy import stats
from scipy.integrate import simps
import os
import scipy.optimize as so
from shutil import copyfile
import tempfile
import time
# from sklearn.neighbors.kde import KernelDensity

# Get directory name
dirname = os.path.dirname(os.path.realpath(__file__))
workingdir = os.getcwd()

# ------------------------------------------------------------------------------
# H E L P E R  F U N C T I O N S
# ------------------------------------------------------------------------------
def histedges_equalN(x, nbin):
    """
    Get histogram edges s.t. there are equal number or events in each bin
    From : https://stackoverflow.com/questions/39418380/histogram-with-equal-number-of-points-in-each-bin
    """
    npt = len(x)
    return np.interp(np.linspace(0, npt, nbin + 1),
                     np.arange(npt),
                     np.sort(x))

# ------------------------------------------------------------------------------
# M A I N
# ------------------------------------------------------------------------------
def ymw16(DM, l, b, ndir=1, best_fit = False):
    """
    Calcultes the distance (DM) to a source using the provided
    dispersion measure (distance). This function assumes the ymw16 model
    for electron density distribution in the galaxy.

    Requires an installed version of the ymw16 model in the local
    directory.

    A base error on the measured distance is 20% as well as the mean
    variation in distance calculated from the neighbouring pixels.

    Arguments
    ---------
    * DM [array] : Dispersion measure [pp/cm^3] or Distance [kpc]
    * l [array] : longitude (galactic coordinates)
    * b [array] : latitude (galactic coordinates)
    * ndir : (1) DM to distance, (2) distance to DM

    Returns
    -------
    * Dist [array] : Distance to PSR # in kpc
    """
    if best_fit:
        # reset to best fit value
        copyfile("%s/ymw16_v1.3/ymw16par_bestfit.txt"%dirname, "%s/ymw16_v1.3/ymw16par.txt"%dirname)

    Dist = ymw.distances(l=l, b=b, DM=DM, dirname=dirname+"/ymw16_v1.3/", ndir=ndir)
    # Dist = ymw.distances(l=l, b=b, DM=DM)
    return np.array(Dist)

def dist_bf(name, DM, l, b, ndir=1):
    """
    Calcultes the distance to a source using the provided
    dispersion measure. This function assumes the ymw16 model
    for electron density distribution in the galaxy with their best-fit
    parameters.

    Arguments
    ---------
    * name [array] : List of names of sources for labelling
    * DM [array] : Dispersion measure # in pp/cm^3
    * l [array] : longitude (galactic coordinates)
    * b [array] : latitude (galactic coordinates)
    * ndir : (1) DM to distance, (2) distance to DM

    Returns
    -------
    Dist [array] : Distance to PSRs # in kpc
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

    return ymw16(DM=DM, l=l, b=b, ndir=ndir)

def dist_pdf(name, DM, l, b,
    mode="kde", nd=100,
    MC_mode = "bestfit", n_MC=1000,
    save_files=False, plots=False, error=False,
    outdir=dirname):
    """
    Calcultes the distance to a source using the provided
    dispersion measure. This function assumes the ymw16 model
    for electron density distribution in the galaxy.

    The pdf is calculated by varying the parameters in the ymw16 model
    and recalculating the distance to the PSR which is the biggest source
    of systematics in the distance calculation

    Arguments
    ---------
    * name [array] : List of names of sources for labelling
    * DM [array] : Dispersion measure # in pp/cm^3
    * l [array] : longitude (galactic coordinates)
    * b [array] : latitude (galactic coordinates)
    * ndir : (1) DM to distance, (2) distance to DM
    * n_MC [float] : number of monte carlo samples
    * mode [str] : return PDF from kernel density estimator (kde) or
                 as a hisgogram (hist).
    * nd [int] : number of distance values (kde) or bins (hist)

    Keyword arguments:
    * MC_mode -- bestfit, uniform, gaussian# (# = error in %)
    * save_files -- saves pdfs as txt files (default False)
    * plots -- pdf plots for each pulsar (default False)

    Returns
    -------
    * dist_pdfs [array] : PDFs of the distance to PSRs
    * dist [array] : distances for each pdf point (if kde) or
                   bin edges (if hist) # [kpc]
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

    # --- Run the Monte Carlo --- #
    # Create list
    D_list = []
    dist_pdfs = []
    dist = []

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
        D_list.append(ymw16(DM=DM, l=l, b=b, ndir=1))
    D_list = np.array(D_list)


    # Get PDF
    if mode =='kde':
        for i in range(DM.size):
            _dist = np.linspace(max(D_list[:,i].min()/1.2, 0),
                D_list[:,i].max()*1.2, nd)
            kde_func = stats.gaussian_kde(D_list[:,i])

            dist.append(_dist)
            dist_pdfs.append(kde_func(_dist))

    elif mode=='hist':
        for i in range(DM.size):
            d_hist, d_edges = np.histogram(D_list[:,i],
                bins=nd, density=True)
                # bins=n_dbins, density=True)
            dist_pdfs.append(d_hist)
            dist.append(d_edges)
    dist_pdfs = np.array(dist_pdfs)
    dist = np.array(dist)

    if error:
        errors = []
        def confidence_interval(x, pdf, xbins, pdf_peak, confidence_level):
            binwidth = np.diff(xbins)[0]
            return (pdf[np.logical_and(xbins < pdf_peak + x, xbins > pdf_peak - x)]*binwidth).sum() - confidence_level
        for i in range(DM.size):
            pdf_peak = dist[i,np.argmax(dist_pdfs[i,:])]
            one_sigma = so.brentq(confidence_interval, 0., 10., args=(dist_pdfs[i,:], dist[i,:], pdf_peak, 0.68))
            errors.append([pdf_peak,one_sigma])
        errors = np.array(errors)

    # --> make plots
    if plots:
        for i in range(DM.size):
            kde_func = stats.gaussian_kde(D_list[:,i])
            _dist = np.linspace(max(D_list[:,i].min()/1.2, 0),
                D_list[:,i].max()*1.2, nd)

            plt.figure()
            plt.hist(D_list[:,i], bins=nd, density=True, alpha=0.8)
            if error:
                plt.axvline(x=errors[i,0]+errors[i,1], color='r')
                plt.axvline(x=errors[i,0]-errors[i,1], color='r')
            plt.plot(_dist, kde_func(_dist))
            plt.ylabel('N')
            plt.title('%s'%(str(name[i])))
            plt.xlabel('D [kpc]')
            plt.savefig(dirname + '/../plots/%s_%s.pdf'%(name[i], MC_mode))
            plt.close()

    # --> save files
    if save_files:
        for i in range(DM.size):
            if mode=="kde":
                np.savetxt("%s/%s_pdf_%s_%s.dat"%(outdir, name[i], mode, MC_mode),
                    zip(dist[i], dist_pdfs[i]), delimiter=" ")

            elif mode=="hist":
                np.savetxt("%s/%s_pdf_%s_%s.dat"%(outdir, name[i], mode, MC_mode),
                    dist_pdfs[i], delimiter=" ")
                np.savetxt("%s/%s_Dedges_%s_%s.dat"%(outdir, name[i], mode, MC_mode),
                    dist[i], delimiter=" ")


    if error:
        return dist_pdfs, dist, errors
    else:
        return dist_pdfs, dist

def DM_pdf(name,
    l, b,
    d_bf, sigma_d=None,
    # mode="hist",
    MC_mode = "bestfit", n_MC=1000,
    outdir=dirname,
    nd=100, ndm=50,
    equal_d_spacing=True):
    """
    Calcultes the distance to a source using the provided
    dispersion measure. This function assumes the ymw16 model
    for electron density distribution in the galaxy.

    The pdf is calculated by varying the parameters in the ymw16 model
    and recalculating the distance to the PSR which is the biggest source
    of systematics in the distance calculation

    Arguments
    ---------
    * name [array] : List of names of sources for labelling
    * l [array] : longitude (galactic coordinates)
    * b [array] : latitude (galactic coordinates)
    * d_bf [array]: best fit distance 
    * n_MC [float] : number of monte carlo samples

    Keyword arguments:
    * MC_mode -- bestfit, uniform, gaussian# (# = error in %)
    * save_files -- saves pdfs as txt files (default False)
    * plots -- pdf plots for each pulsar (default False)

    Returns
    -------
    * P(DM|D) [array] : 2-D Histogram where $\int dDM P(DM|D) = 1$
    * dist_edges [array] : bin edges # [kpc]
    * DM_edges [array] : bin edges # [pp/cm^3]
    """
    mode = "hist"
    # Set up distance: we draw distance assuming a lognormal around the best-fit
    if sigma_d == None:
        sigma_d = np.log(2)/2. # Meaning that a value that is at 2xd_bf is 2 sigma away
    dist = np.random.lognormal(np.log(d_bf), sigma_d, size = int(n_MC))
    _mask = (dist>20)
    while _mask.sum() > 1: # Make sure we don't have very many
        dist[_mask] = np.random.lognormal(np.log(d_bf), sigma_d, size = int(_mask.sum()))
        _mask = (dist>20)

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

    # --- Run the Monte Carlo --- #
    # Create list
    DM_list = []
    dm_pdfs = []
    dm = []

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
    print "Starting MC"
    t1 = time.time()
    # DM_list = np.zeros((n_MC, len(dist)))
    DM_list = np.zeros(n_MC)
    print "MC values for YMW16-shape: ", val_MC.shape
    for i in range(n_MC):
        # update values
        vals[ind_order] = val_MC[i] # Set values
        np.savetxt(filename_input, zip(names, vals),
            delimiter=" ", fmt="%s")
        with open(filename_input, 'a') as f:
            f.write('Ele_arm %.6f %.6f %.6f %.6f %.6f\n'%tuple(val_narm_MC[i]))
            f.write('Wid_arm %i %i %i %i %i\n'%tuple(warm))
        DM_list[i] = ymw16(DM=np.array([dist[i]*1e3]), l=l, b=b, ndir=2) # convert kpc to pc

    # print DM_list
    t2 = time.time()
    print "MC took: %.2fs"%(t2-t1)
    print "This is: %.2es per DM"%((t2-t1) / n_MC)

    # dist = np.array([dist for i in range(n_MC)]) # reshpae dist

    # Currently not functioning
    # Get PDF
    # if mode =='kde':
    #     print "\nWARNING: KDE mode is currently not working properly. Please use mode=hist instead.\n"
    #     print "Starting on gaussian_kde"
    #     t1 = time.time()
    #     dist = histedges_equalN(dist, nd)
    #     DM_vals = histedges_equalN(DM_list, ndm) # Equal number of events between values
    #     values =  np.array(zip(dist.flatten(), DM_list.flatten())).T # Get 2d values into correct shape for gaussian_kde
    #     kde_func = stats.gaussian_kde(values)
    #     t2 = time.time()
    #     print "Performing kde took: %.2fs"%(t2-t1)
    #     return kde_func, dist, DM_vals

    if mode=='hist':
        # print "Mode hist currently not defined for DM_pdf(). Quitting"
        # quit()
        print "Starting on histogram"
        t1 = time.time()
        # print DM_list
        # dist_edges = histedges_equalN(dist, np.min(nd, int(n_MC / 100.)))
        # DMedges = histedges_equalN(DM_list, np.min(ndm, int(n_MC / 100.)))
        dist_edges = histedges_equalN(dist, nd)
        if equal_d_spacing:
            dist_edges = np.linspace(dist_edges[1], dist_edges[-1], nd)

        DM_edges = histedges_equalN(DM_list, ndm)
        H, dist_edges, DMedges = np.histogram2d(dist.flatten(), DM_list.flatten(),
            bins=np.array([dist_edges, DM_edges]), normed=True)
        t2 = time.time()
        print "Getting histogram took: %.2fs"%(t2-t1)
        return H, dist_edges, DMedges


def dist_parallax(name, P, P_errors, nd=100, save_files=False, plots=False,
    outdir=dirname):
    """
    Arguments
    ---------
    * name [array] : List of names of sources for labelling
    * P [array] : List of names of sources for labelling
    * P_error [ndarray] : List of names of sources for labelling
        Shape[i] = (# of sources,) for symmetric
        Shape[i] = (# of sources, 2) for asymmetric
        - [:,0] is the lower error, [:,1] is the upper error
        Example np.array([[0.1,0.2],0.1,0.1])
        First element represents asymmetric errors, second
        and third are symmetric

    Returns
    -------
    * dist_pdfs [array] : PDFs of the true distance to PSRs given a measurement of w
    * dist [array] : distances for each pdf point (if kde) or
                   bin edges (if hist) # [kpc]
    """
    # Check for symmetric or asymmetric
    sigma_P = np.zeros([P.size,2])
    for i in range(P.size):
        if isinstance(P_errors[i], list):
            sigma_P[i,:] = P_errors[i]
        else:
            sigma_P[i,:] = P_errors[i]

    dist_pdfs = []
    dist_list = []
    # d_list = np.linspace(0.,15.,num=nd)
    # d_width = d_list[1] - d_list[0]
    # d_cen = d_list + d_width/2

    for i in range(P.size):
        # print np.where(P[i]-3*sigma_P[i,0]>0, P[i]+3*sigma_P[i,0], 0.001)
        # d_list = np.linspace(np.where(P[i]-3*sigma_P[i,0]>0, P[i]-3*sigma_P[i,0], 0.001),
        #     1./(P[i]+3*sigma_P[i,1]), num=nd)
        d_list = np.linspace(1./(P[i]+3*sigma_P[i,1]),
            np.where(P[i]-3*sigma_P[i,0]>0, 1./(P[i]-3*sigma_P[i,0]), 15.),
            num=nd)


        d_pdf = (np.heaviside(1./d_list - P[i],0.5) *
            np.exp(-(((P[i] - 1/d_list)/sigma_P[i,1])**2)/2.) / (d_list**2.) +
            np.heaviside(P[i] - 1./d_list , 0.5) *
            np.exp(-(((P[i] - 1/d_list)/sigma_P[i,0])**2)/2.) / (d_list**2.))

        # A = sum(d_hist*d_width) # Normalising factor
        d_pdf /= simps(d_pdf, d_list)
        # print name[i], d_list.min(), d_list.max(), d_pdf.max()
        dist_list.append(d_list)
        dist_pdfs.append(d_pdf)

    dist_pdfs = np.array(dist_pdfs)
    dist_list = np.array(dist_list)

    # print dist_pdfs.shape
    if save_files:
        for i in range(P.size):
            np.savetxt("%s/%s_pdf_parallax.dat"%(outdir, name[i]),
                zip(dist_list[i], dist_pdfs[i]), delimiter=" ")

    # --> make plots
    if plots:
        for i in range(P.size):
            plt.figure()
            plt.plot(dist_list[i], dist_pdfs[i])
            plt.ylabel('N')
            plt.title('%s'%(str(name[i])))
            plt.xlabel('D [kpc]')
            # plt.savefig(dirname + '/../plots/%s_parallax.pdf'%(name[i]))
            plt.savefig(dirname + '/../plots/%s_parallax.pdf'%(name[i]))
            plt.close()

    return dist_pdfs, dist_list

def px_pdf(name, px, dist, sigma_px):
    """
    Arguments
    ---------
    * name [array] : List of names of sources for labelling
    * px [float] : observed parallax [mas]
    * dist [array] : True distance to a source [kpc]
    * sigma_px [array] : Error on parallax [mas]

    Returns
    -------
    * px_pdf [array] : PDFs of the measured parallax to PSRs given a true distance Dt (not correctly normalized!)
    * dist [array] : distances for each pdf point (if kde) or
                   bin edges (if hist) # [kpc]
    """
    # Check for symmetric or asymmetric

    px_pdf = (np.heaviside(px - 1./dist,0.5) *
        np.exp(-(((px - 1/dist)/sigma_px[1])**2)/2.) +
        np.heaviside(1./dist - px, 0.5) *
        np.exp(-(((px - 1/dist)/sigma_px[0])**2)/2.)) # P(w_meas|w(D))

    return px_pdf, dist
