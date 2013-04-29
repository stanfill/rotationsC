# How to Build Package with Rcpp
========================================================
This document is meant for me to remember how to build a package based on Rcpp with my system.  Other users should use with caution.

## Using `Rcpp.package.skeleton`

This is a simple example of how to create a package from *Rcpp.package.skeleton*.  Make sure your working directory is where you want the package folder to be created then run

```
library(Rcpp)
Rcpp.package.skeleton("rotations2",example_code=F,cpp_files=c("estimators.cpp","FisherMethod.cpp"))
```
This will not allow for RcppArmadillo so make the following two changes to the *DESCRIPTION* file:

* *Depends: Rcpp (>= 0.10.3)* to *Depends: Rcpp (>= 0.10.3), RcppArmadillo*
* *LinkingTo: Rcpp* to *LinkingTo: Rcpp, RcppArmadillo*

In the *src* folder add the following line to the Makevars file:
```
PKG_LIBS = $(shell $(R_HOME)/bin/Rscript -e "Rcpp:::LdFlags()" ) $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)
```

In the same folder add the following line to Makevars.win:
```
PKG_LIBS = $(shell $(R_HOME)/bin${R_ARCH_BIN}/Rscript.exe -e "Rcpp:::LdFlags()") $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)
```

Now run
```
compileAttributes()
```
In the "/inst/include/rotations2_RcppExports.h" file change "../inst/include/rotations2.h" to "rotations2.h"

Finally, run the following code
```
library(roxygen2)
library(devtools)
roxygenize(rotations2)
install()
```

## Using `RcppArmadillo.package.skeleton`

The basic idea is use the command
```
library(RcppArmadillo)
RcppArmadillo.package.skeleton('rotations2')
compileAttributes(verbose=T)
```
This will make the package compatible with RcppArmadillo but other changes will need to be made.  Haven't looked into this much