import numpy as np
from scipy.integrate import simps
from PSRdist import Distances as D

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font', family='serif', serif='cmr10', size=10)
rc('text', usetex=True)


names = np.array(["B0628-28", "B1237+25", "B1237+10"])
# DM = np.array([34.36, 9.2755, 10.2755])
# l = np.array([236.95, 252.44, 200.44])
# b = np.array([-16.75, 86.54, 19.54])
P = np.array([0.4, 0.4, 0.4])
Perr = np.array([[3.,0.1], 0.1, 0.5])
dist_pdfs, dist = D.dist_parallax(names, P, Perr,
    nd=100, plots=True, save_files=True)

