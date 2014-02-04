*rotations*
========================================================
A stable version (1.0) of the *rotations* package is available for download from CRAN.  For Mac users, the current version (1.1) of *rotations* can be downloaded from GitHub.  Windows users will need to download the repo, compile the C++ code and install locally.

### Installation instructions for Mac: 
```
library(devtools)
install_github('rotationsC','stanfill',subdir='rotations')
library(rotations)
```

### Change log:

#### Minor Changes -

* The options for the `method` argument of the `region` functions have been changed to `"trans"` and `"direct"` to align with the current names in my dissertation.  `"trans"` is short for "transformation," which is used to access the methods based on a transformation of a directional statistics result while `"direct"`  is used for the methods that rely on theory for SO(3) directly.

#### Major Changes -

* The `plot.SO3` function now has an `interactive` option that renders 3D plots of the supplied data using the `rgl` and `sphereplot` packages.  Currently, all options available for the static `ggplot2` plots are availble for the interactive plots, but only one column can be displayed at a time.  As an example see

```
library(rotations)
Rs<-ruars(20,rcayley,kappa=10)
plot(Rs,show_estimates='proj.mean', mean_region='moment theory',alp=.1,intractive=FALSE)
plot(Rs,show_estimates='proj.mean', mean_region='moment theory',alp=.1,intractive=TRUE)
```
