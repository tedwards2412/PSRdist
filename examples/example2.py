import numpy as np
from scipy.integrate import simps
from MSPDist import Distances as D

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font', family='serif', serif='cmr10', size=10)
rc('text', usetex=True)


names = np.array(["B0628-28", "B1237+25"])
DM = np.array([34.36, 9.2755])
l = np.array([236.95, 252.44])
b = np.array([-16.75, 86.54])
for MC_mode in ["uniform", "bestfit", "gaussian5", "gaussian10",
    "gaussian20"]:

    tag = "_MC_%s"%MC_mode
    dist_pdfs, dist_edges = D.dist_pdf(name=names, DM=DM, l=l, b=b,
        d_edges=np.linspace(0, 10, 201),
        MC=True, n_MC=50000, MC_mode=MC_mode)

    D_bf = D.dist_bf(name=names, DM=DM, l=l, b=b)

    d_mean = ((dist_edges.T[1:]+dist_edges.T[:-1])/2).T

    for i in range(len(names)):
        print "%s: integral PDF: %.2f"%(names[i], simps(dist_pdfs[i], d_mean[i]))
        plt.figure(figsize=(4,3))
        plt.plot(d_mean[i], dist_pdfs[i], "r-")
        plt.plot([D_bf[i], D_bf[i]], [-10000, 10000], "k:",
            label='ymw16: best-fit', zorder=-10)
        plt.title('%s'%names[i])
        plt.xlabel(r'$D$ [kpc]')
        plt.ylabel(r'pdf')
        plt.ylim(0, dist_pdfs[i].max() * 1.2)
        plt.legend(loc=1)
        plt.tight_layout()
        plt.savefig('plots/%s%s.pdf'%(names[i], tag))
