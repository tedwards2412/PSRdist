# PSRDist
Calculates a variety of properties relating to the distances to Pulsars.
Code consists of a wrapper for the ymw16 electron density model (Yao, Manchester and Wang, 2017, Astrophys. J., 835, 29; arXiv:1610.09448).
It has been adopted to scan over the parameters of the electron density model.
The documentaiton can be found [here](https://tedwards2412.github.io/PSRdist/).
If you use this package for research purposes please include a reference to the original ymw16 paper as well as Bartels, Edwards & Weniger (in prep). 

## Requirements and installation
- numpy
- cython

Compile the cython code by running ./make.sh.

```
./make.sh
```

Currently working on getting it pip installable.
For now simply add the module location to your PYTHONPATH by adding 

```
PYTHONPATH="/PATH/TO/MODULE/PSRdist:$PYTHONPATH"
```

to your bash_profile or bashrc scripts. The module should then be imported easily as a wrapper:

```
import py_ymw16 as ymw
```
or with the use of the additional MC functions:

```
from PSRDist import Distances as D
```

Examples on how to use the code are provided and found [here](https://github.com/tedwards2412/PSRdist/tree/master/examples).

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
