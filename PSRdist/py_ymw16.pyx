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
    int ndir=1, int _np=1, int vbs=0,
    char* dirname="./ymw16_v1.3/") nogil:
    """

    Arguments
    ---------
    DM [double/float] : Dispersion measure [pp/cm^3] or Distance in [pc]
    l [double/float] : longitude (galactic coordinates) [deg]
    b [double/float] : latitude (galactic coordinates) [deg]
    DM_Host : Contribution of Fast-radio burst host galaxy to the observed DM (only used for IGM mode)
    ndir : (1) DM to distance, (2) distance to DM
    _np : (-1) IGM, (0) Magellanic clouds, (1) Galaxy
    vbs : verbose (default is 0)
    dirname : directory of the ymw16 code input files

    Returns
    -------
    Distance in kpc
    """
    # Fix names of directories
    cdef char* text = ""

    if ndir==1:
      return dmdtau(gl=l, gb=b, dordm=DM, DM_Host=DM_Host,
          ndir=ndir, np=_np, vbs=vbs, dirname=dirname, text=text) * 1e-3 # convert pc to kpc

    elif ndir==2:
      return dmdtau(gl=l, gb=b, dordm=DM, DM_Host=DM_Host,
          ndir=ndir, np=_np, vbs=vbs, dirname=dirname, text=text) # cm-3 pc

cpdef double[:] distances(double[:] l, double[:] b, double[:] DM,
  double DM_Host=0, int ndir=1, int _np=1, int vbs=0,
  char* dirname="./ymw16_v1.3/"):
  """

  Arguments
  ---------
  DM [array] : Dispersion measure [pp/cm^3] or Distances [pc]
  l [array] : longitude (galactic coordinates) [deg]
  b [array] : latitude (galactic coordinates) [deg]
  DM_Host : Contribution of Fast-radio burst host galaxy to the observed DM (only used for IGM mode)
  ndir : (1) DM to distance, (2) distance to DM
  _np : (-1) IGM, (0) Magellanic clouds, (1) Galaxy
  vbs : verbose (default is 0)
  dirname : directory to the ymw16 code input files

  Returns
  -------
  Distances for multiple sources in kpc
  """

  cdef int N = len(l) ; # Number of sources
  cdef double[:] D = np.zeros(N) ;
  cdef Py_ssize_t iN ;

  with nogil:
    for iN in range(N):
      D[iN] = distance(l=l[iN], b=b[iN], DM=DM[iN],
              DM_Host=DM_Host, ndir=ndir, _np=_np, vbs=vbs,
              dirname=dirname) ;

  return D
