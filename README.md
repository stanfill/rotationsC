## *intervals*
========================================================

This folder contains the code necessary to run the simulations for the intervals research

## *rotations*
========================================================
The *rotations* package is also available in this folder.  It uses C++ to greatly increase computing time, especially for the estimators.  Because of differences in compilers, only a Mac version is available at this time.  Sorry.  To install, both a Fortran and C++ compiler are required.


### Installation instructions:


#### For Mac: 
```
library(devtools)
install_github('rotationsC','stanfill',subdir='rotations')
library(rotations)
```