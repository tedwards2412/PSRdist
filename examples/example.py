import numpy as np
import py_ymw16 as ymw
import time

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font', family='serif', serif='cmr10', size=10)
rc('text', usetex=True)

# Load parameter file
filename = 'ymw16_v1.3/ymw16par.txt'
names = np.loadtxt(filename, usecols=(0,), dtype='str')
vals = np.loadtxt(filename, usecols=(1,))
names = np.loadtxt(filename, usecols=(0,), dtype='str')
vals = np.loadtxt(filename, usecols=(1,))

# Setup grid parameters
ngrid=91
n1_list = np.linspace(0.008, 0.016, ngrid)
H1_list = np.linspace(1200, 2000, ngrid)

# Setup sources to evaluate
sources = np.array([
            ["B0628-28", 34.36, 236.95, -16.75],
            ["B1237+25", 9.2755, 252.44, 86.54]]) # Name, DM, glon [deg], glat [deg]

for source, DM, l, b in sources:
    t1 = time.time()
    print " --- %s --- "%source
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
            D_list.append(ymw.distance(l=l, b=b, DM=DM, vbs=0))
    
    t2 = time.time()
    print "Evaluating %i grid points took %.2fs"%(ngrid**2, t2-t1)
    D_list = np.array(D_list)
    D_list*=1
    plt.figure()
    plt.hist(D_list, bins=np.linspace(0, 10, 161))
    plt.title(r"%s"%source)
    plt.ylabel('N')
    plt.xlabel('D [kpc]')
    plt.savefig('examples/%s_test.pdf'%(source))
    plt.close()