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

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.initializedcheck(False)

cpdef double distance(double l, double b, double DM,
    double DM_Host=0,
    int ndir=1, int _np=1, int vbs=0) nogil:
    """
    returns distance in kpc
    """


    # Fix names of directories
    cdef char* dirname = "./ymw16_v1.3/"
    cdef char* text = ""

    return dmdtau(gl=l, gb=b, dordm=DM, DM_Host=DM_Host,
        ndir=ndir, np=_np, vbs=vbs, dirname=dirname, text=text) * 1e-3 # convert pc to kpc

cpdef double[:] distances(double[:] l, double[:] b, double[:] DM,
  double DM_Host=0, int ndir=1, int _np=1, int vbs=0):
  """
  Return distances for multiple sources
  """

  cdef int N = len(l) ; # Number of sources
  cdef double[:] D = np.zeros(N) ;
  cdef Py_ssize_t iN ;

  with nogil:
    for iN in range(N):
      D[iN] = distance(l=l[iN], b=b[iN], DM=DM[iN],
              DM_Host=DM_Host, ndir=ndir, _np=_np, vbs=vbs) ;

  return D
