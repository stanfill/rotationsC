*rotations*
========================================================
A stable version (1.4) of the *rotations* package is available for download from CRAN.  For Windows users, the current version (1.5) of *rotations* can be downloaded from GitHub.  Mac users will need to download the repo, compile the C++ code and install.

### Installation instructions for Mac: 
```
devtools::install_github('stanfill/rotationsC',subdir='rotations')
```

### Change log:

* The Maxwell-Boltzmann distribution has been added: `dmaxwell`, `pmaxwell` and `rmaxwell`
* Due to issue with the bessel functions, the Cayley distribution was previously used to approximate the matrix-Fisher distribution for large kappa.  Now the Maxwell-Boltzmann distribution is used for `kappa>200` and the approximation is far superior.
