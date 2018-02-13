#!/usr/bin/env python
import os
from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy

extensions = [
    Extension("PSRDist.*", ["PSRDist/*.pyx", 
        "PSRDist/ymw16_v1.3/dmdtau.c", "PSRDist/ymw16_v1.3/dora.c", "PSRDist/ymw16_v1.3/fermibubble.c",
        "PSRDist/ymw16_v1.3/frb_d.c", "PSRDist/ymw16_v1.3/galcen.c", "PSRDist/ymw16_v1.3/gum.c", 
        "PSRDist/ymw16_v1.3/lmc.c", "PSRDist/ymw16_v1.3/localbubble.c", "PSRDist/ymw16_v1.3/ne_crd.c",
        "PSRDist/ymw16_v1.3/nps.c", "PSRDist/ymw16_v1.3/smc.c", "PSRDist/ymw16_v1.3/spiral.c", 
        "PSRDist/ymw16_v1.3/thick.c", "PSRDist/ymw16_v1.3/thin.c", "PSRDist/ymw16_v1.3/ymw16_ne.c",
        "PSRDist/ymw16_v1.3/ymw16par.c"],
    include_dirs=[numpy.get_include()], 
)]

print(extensions)

setup(name='PSRDist',
      version='1.0',
      description='Python Distribution Utilities',
      author='Richard Bartels and Thomas Edwards',
      author_email='tedwards2412@gmail.com',
      url='https://github.com/tedwards2412/PSRDist',
      packages=['PSRDist'],
      data_files=[('lib/python2.7/site-packages/PSRDist/ymw16_v1.3/', ['PSRDist/ymw16_v1.3/ymw16par_ranges.txt', 'PSRDist/ymw16_v1.3/ymw16par_Ele_arm_unc.txt','PSRDist/ymw16_v1.3/ymw16par_bestfit.txt','PSRDist/ymw16_v1.3/ymw16par_Ele_arm_ranges.txt','PSRDist/ymw16_v1.3/spiral.txt','PSRDist/ymw16_v1.3/ymw16par.txt'])],
    ext_modules = cythonize(extensions),
    long_description=open('README.md').read()
)
