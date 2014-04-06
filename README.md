*rotations*
========================================================
A stable version (1.1) of the *rotations* package is available for download from CRAN.  For Mac users, the current version (1.2) of *rotations* can be downloaded from GitHub.  Windows users will need to download the repo, compile the C++ code and install locally.

### Installation instructions for Mac: 
```
library(devtools)
install_github('rotationsC','stanfill',subdir='rotations')
library(rotations)
```

### Change log:

#### Minor Changes -

* Sampling from the von Mises Fisher and matrix Fisher distributions is now done in C++ 
* `method` argument in `regions` function was changed from `trans` to `transformation`, a call to `match.arg()` still allows for `trans`
* `type` argument in `regions` function has been changed from `theory` to `asymptotic`

#### Bug Fixes -

* A typo in the (not run) Bayes credible region example has been fixed
* `print` method for `Q4` objects now respects the `digits` option