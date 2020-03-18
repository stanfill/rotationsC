# rotations

<!-- badges: start -->
<!-- badges: end -->

The goal of rotations is to provide tools for working with rotation data. A stable version (1.4) of the *rotations* package is available for download from CRAN.  For Windows users, the current version (1.5) of *rotations* can be downloaded from GitHub.  Mac users will need to download the repo, compile the C++ code and install.

## Installation

```
remotes::install_github("stanfill/rotationsC", subdir = "rotations")
```

## Example

This is a basic example which shows you how to solve a common problem:

```{r}
library(rotations)
## basic example code
```

## Change log

* The Maxwell-Boltzmann distribution has been added: `dmaxwell`, `pmaxwell` and `rmaxwell`.
* Due to issue with the bessel functions, the Cayley distribution was previously used to approximate the matrix-Fisher distribution for large kappa.  Now the Maxwell-Boltzmann distribution is used for `kappa>200` and the approximation is far superior.
