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

* Samling from the von Mises Fisher distribution is now done in C++ to increase speed

