"""# MSPDist
Calculates probability density functions for the distances to Millisecond Pulsars.
Code consists of a wrapper for the ymw16 electron density model (Yao, Manchester and Wang, 2017, Astrophys. J., 835, 29; arXiv:1610.09448). It has been adopted to scan over the parameters of the electron density model.

## Requirements and installation
- numpy
- cython

Compile the cython code by running ./make.sh. Then py_ymw16 can be imported by python.
Currently there is an issue with reading and writing to files in installation directories
For now simply add MSPDist directory to your PYTHONPATH

## Errors
We assume the ymw16 model for the electron density distribution of our galaxy
To calculate the errors on the distances to the pulsars we assume that the
parameters of the ymw16 are uncertain and perform a grid scan to produce an array
of distances to each MSP. We then create a pdf through a normalised histogram of this list.

## Changes in the original ymw16 code.
The version of the ymw16 code provided contains some changes compared to the original in the following files:
- cn.h
- dmdtau.c: now returns a double and print statements have been turned off
- ymw16par.txt

The original best fit parameters from the ymw16 paper can be found in ymw16par_bestfit.txt.
"""

from PSRDist.Distances import *
import py_ymw16 as ymw
