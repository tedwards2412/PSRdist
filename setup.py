#!/usr/bin/env python
import os
from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy

extensions = [
    Extension("MSPDist.*", ["MSPDist/*.pyx", 
        "MSPDist/ymw16_v1.3/dmdtau.c", "MSPDist/ymw16_v1.3/dora.c", "MSPDist/ymw16_v1.3/fermibubble.c",
        "MSPDist/ymw16_v1.3/frb_d.c", "MSPDist/ymw16_v1.3/galcen.c", "MSPDist/ymw16_v1.3/gum.c", 
        "MSPDist/ymw16_v1.3/lmc.c", "MSPDist/ymw16_v1.3/localbubble.c", "MSPDist/ymw16_v1.3/ne_crd.c",
        "MSPDist/ymw16_v1.3/nps.c", "MSPDist/ymw16_v1.3/smc.c", "MSPDist/ymw16_v1.3/spiral.c", 
        "MSPDist/ymw16_v1.3/thick.c", "MSPDist/ymw16_v1.3/thin.c", "MSPDist/ymw16_v1.3/ymw16_ne.c",
        "MSPDist/ymw16_v1.3/ymw16par.c"],
    include_dirs=[numpy.get_include()], 
)]

print(extensions)

setup(name='MSPDist',
      version='1.0',
      description='Python Distribution Utilities',
      author='Richard Bartels and Thomas Edwards',
      author_email='tedwards2412@gmail.com',
      url='https://github.com/tedwards2412/MSPDist',
      packages=['MSPDist'],
      data_files=[('lib/python2.7/site-packages/MSPDist/ymw16_v1.3/', ['MSPDist/ymw16_v1.3/ymw16par_ranges.txt', 'MSPDist/ymw16_v1.3/ymw16par_Ele_arm_unc.txt','MSPDist/ymw16_v1.3/ymw16par_bestfit.txt','MSPDist/ymw16_v1.3/ymw16par_Ele_arm_ranges.txt','MSPDist/ymw16_v1.3/spiral.txt','MSPDist/ymw16_v1.3/ymw16par.txt'])],
    ext_modules = cythonize(extensions),
    long_description=open('README.md').read()
)
