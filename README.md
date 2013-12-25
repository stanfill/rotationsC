*rotations*
========================================================
A stable version (0.2) of the *rotations* package is available for download from CRAN.  For MAC users, the current version (0.2.1) of *rotations* can be downloaded from GitHub.  Windows users will need to download the repo, compile the C++ code and install locally.

### Installation instructions for Mac: 
```
library(devtools)
install_github('rotationsC','stanfill',subdir='rotations')
library(rotations)
```

### Change log:
The differences between the CRAN (0.2) and GitHub (1.0) versions of the *rotations* package are as follows:

* Fixes in documentation for Bayes point estimate

* The functions `SO3` and `Q4` no longer exist.  All of their functaionality has been moved to `as.SO3` and `as.Q4`

* Create a `plot.Q4` function that uses `plot.SO3` after casting the object to class `SO3`

*intervals*
--------------------------------------------------------
This folder contains the code necessary to run the simulations for the confidence region research

