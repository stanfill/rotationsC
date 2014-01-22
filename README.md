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
#### New Features - 

* Vignette added that introduces the package

* Create a `plot.Q4` function that uses `plot.SO3` after casting the object to class `SO3`

#### Bug Fixes -

* Updated Bayes sampling method to avoid seg faults

#### Major Changes -

* The functions `SO3` and `Q4` no longer exist.  All of their functionality has been moved to `as.SO3` and `as.Q4`

* `angle` and `axis` have been renamed `mis.angle` and `mis.axis`, respectively, to avoid naming clashes with the `graphics` package

* Adopt `period.sep` naming convention for all functions in package.  The affected functions were formerly known as `sum_dist`, `exp_skew`, `cayley_kappa`, `fisher_kappa` and `vmises_kappa`.  New names are the same but with `.` in place of `_`.

* `dist` renamed to `rot.dist` to avoid clashes with `stats` package

#### Minor Changes -

* Fixes in documentation for Bayes point estimate

* `print.Q4` and `print.SO3` no longer print the object class

* `print.SO3` now names the columns `R11` through `R33` to signify which element in the matrix each row corresponds to

* `is.SO3` more rigorously tests for conditions of `SO3`

* Arithmetic for `SO3` objects now is possible for samples of rotations

*intervals*
--------------------------------------------------------
This folder contains the code necessary to run the simulations for the confidence region research

