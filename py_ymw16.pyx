# ----------------------------------------------------------------------------------------
# py_ymw16.pyx
# ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
#
# Cython wrapper for py_ymw16.pyx
# Author: Richard Bartels (r.t.bartels [at] uva.nl)
#
# ----------------------------------------------------------------------------------------

import numpy as np
cimport numpy as np
cimport cython

cimport py_ymw16

cpdef double distance(double l, double b, double DM, 
    double DM_Host=0,
    int ndir=1, int np=1, int vbs=0) nogil:
    """
    returns distance in kpc
    """
    

    # Fix names of directories
    cdef char* dirname = "./ymw16_v1.3/"
    cdef char* text = ""
    
    return dmdtau(gl=l, gb=b, dordm=DM, DM_Host=DM_Host, 
        ndir=ndir, np=np, vbs=vbs, dirname=dirname, text=text) * 1e-3 # convert pc to kpc