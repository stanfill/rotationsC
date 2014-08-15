*rotations*
========================================================
A stable version (1.2) of the *rotations* package is available for download from CRAN.  For Mac users, the current version (1.3) of *rotations* can be downloaded from GitHub.  Windows users will need to download the repo, compile the C++ code and install locally.

### Installation instructions for Mac: 
```
library(devtools)
install_github('rotationsC','stanfill',subdir='rotations')
library(rotations)
```

### Change log:


#### Major Changes -
* The function `discord` has been added, which computes a measure of discord for random rotations
* The previously proposed extrinsic measure of discord is included as well as the novel intrinsic measure

#### Minor Changes -

* Vignette has been updated
* `tail.Q4` and `tail.SO3` have been added
