# MSPDist
Calculates probability density functions for the distances to Millisecond Pulsars

## Requirements
- numpy
- cython

## Errors
We assume the ymw16 model for the electron density distribution of our galaxy
To calculate the errors on the distances to the pulsars we assume that the 
parameters of the ymw16 are uncertain and performa grid scan to produce an array
of distances to each MSP. We then create a pdf through a normalised histogram of this list.
