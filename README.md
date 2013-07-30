## *intervals*
========================================================

This folder contains the code necessary to run the simulations for the intervals research

## *rotations2*
========================================================
The *rotations2* package is also available in this folder.  It uses C++ to greatly increase computing time, especially for the estimators.  Because of differences in compilers, a Windows and a Mac version need to be kept separate to be widely accessible.  To install, both a Fortran and C++ compiler are required.

**To contributors**: Please *do not* make changes to the files in the *'rotations2Mac'* folder!  Or, if you do be sure to make the same changes to the files in the  *'rotations2'* folder.

### Installation instructions:

#### For Windows:
```
library(devtools)
install_github('rotationsC','stanfill',subdir='rotations2')
library(rotations2)
```

#### For Mac: 
```
library(devtools)
install_github('rotationsC','stanfill',subdir='rotations2Mac')
library(rotations2)
```