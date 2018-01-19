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
mode='kde'
for MC_mode in [
    "uniform", "bestfit"
        # "gaussian5", "gaussian10",
    # "gaussian20",
    # "gaussian15"
    ]:

    tag = "_%s_%s"%(mode, MC_mode)
    dist_pdfs, dist = D.dist_pdf(name=names, DM=DM, l=l, b=b,
        nd=100, n_MC=50000, MC_mode=MC_mode,
        mode=mode, plots=True)

    D_bf = D.dist_bf(name=names, DM=DM, l=l, b=b)

    for i in range(len(names)):
        print "%s: integral PDF: %.2f"%(names[i], simps(dist_pdfs[i], dist[i]))
        plt.figure(figsize=(4,3))
        plt.plot(dist[i], dist_pdfs[i], "r-")
        plt.plot([D_bf[i], D_bf[i]], [-10000, 10000], "k:",
            label='ymw16: best-fit', zorder=-10)
        plt.title('%s'%names[i])
        plt.xlabel(r'$D$ [kpc]')
        plt.ylabel(r'pdf')
        plt.ylim(0, dist_pdfs[i].max() * 1.2)
        plt.legend(loc=1)
        plt.tight_layout()
        plt.savefig('plots/%s%s.pdf'%(names[i], tag))
