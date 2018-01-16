import os
from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy

extensions = [
    Extension("*", ["*.pyx", 
        "ymw16_v1.3/dmdtau.c", "ymw16_v1.3/dora.c", "ymw16_v1.3/fermibubble.c",
        "ymw16_v1.3/frb_d.c", "ymw16_v1.3/galcen.c", "ymw16_v1.3/gum.c", 
        "ymw16_v1.3/lmc.c", "ymw16_v1.3/localbubble.c", "ymw16_v1.3/ne_crd.c",
        "ymw16_v1.3/nps.c", "ymw16_v1.3/smc.c", "ymw16_v1.3/spiral.c", 
        "ymw16_v1.3/thick.c", "ymw16_v1.3/thin.c", "ymw16_v1.3/ymw16_ne.c",
        "ymw16_v1.3/ymw16par.c"],
    include_dirs=[numpy.get_include()], 
)]

print(extensions)

setup(
    ext_modules = cythonize(extensions),
)
