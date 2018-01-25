import numpy as np
from scipy.integrate import simps
from MSPDist import Distances as D

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font', family='serif', serif='cmr10', size=10)
rc('text', usetex=True)


names = np.array(["B0628-28", "B1237+25", "B1237+10"])
DM = np.array([34.36, 9.2755, 10.2755])
l = np.array([236.95, 252.44, 200.44])
b = np.array([-16.75, 86.54, 19.54])
mode='kde'
for MC_mode in [
    "uniform", "bestfit"
    ]:

    tag = "_%s_%s"%(mode, MC_mode)
    dist_pdfs, dist = D.dist_pdf(name=names, DM=DM, l=l, b=b,
        nd=100, n_MC=1000, MC_mode=MC_mode,
        mode=mode, plots=True, error=True)

    D_bf = D.dist_bf(name=names, DM=DM, l=l, b=b)

