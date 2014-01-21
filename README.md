*rotations*
========================================================
A stable version (0.2) of the *rotations* package is available for download from CRAN.  For Mac users, the current version (1.0) of *rotations* can be downloaded from GitHub.  Windows users will need to download the repo, compile the C++ code and install locally.

### Installation instructions for Mac: 
```
library(devtools)
install_github('rotationsC','stanfill',subdir='rotations')
library(rotations)
```

### Change log:
The differences between the CRAN (0.2) and GitHub (1.0) versions of the *rotations* package are as follows:

* Fixes in documentation for Bayes point estimate

* The functions `SO3` and `Q4` no longer exist.  All of their functionality has been moved to `as.SO3` and `as.Q4`

* Create a `plot.Q4` function that uses `plot.SO3` after casting the object to class `SO3`

* `print.Q4` and `print.SO3` no longer print the object class

* `print.SO3` now names the columns `R11` through `R33` to signify which element in the matrix each row corresponds to

* Made `is.SO3` more rigorous in testing for conditions of `SO3`

* `angle` and `axis` have been renamed `mis.angle` and `mis.axis`, respectively, to avoid naming clashes with the `graphics` package

* Adopt `period.sep` naming convention for all functions in package.  The affected functions were formerly known as `sum_dist`, `exp_skew`, `cayley_kappa`, `fisher_kappa` and `vmises_kappa`.  New names are the same but with `.` in place of `_`.

* `dist` renamed to `rot.dist` to avoid clashes with `stats` package

*intervals*
--------------------------------------------------------
This folder contains the code necessary to run the simulations for the confidence region research

