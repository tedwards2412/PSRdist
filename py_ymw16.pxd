# ----------------------------------------------------------------------------------------
# py_ymw16.pyx
# ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
#
# Cython wrapper for py_ymw16.pyx
# Author: Richard Bartels (r.t.bartels [at] uva.nl)
#
# ----------------------------------------------------------------------------------------

cdef extern from "ymw16_v1.3/cn.h":
    double dmdtau(double gl, double gb, double dordm, double DM_Host, int ndir, int np, int vbs, char *dirname, char *text) nogil;
#     void dmdtau(double gl, double gb, double dordm, double DM_Host, int ndir, int np, int vbs, char *dirname, char *text) nogil;
