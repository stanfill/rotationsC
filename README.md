*rotations*
========================================================
A staple version of the *rotations* package is available for download from CRAN.  For MAC users, the next version of *rotations* can be downloaded as shown below.

### Installation instructions for Mac: 
```
library(devtools)
install_github('rotationsC','stanfill',subdir='rotations')
library(rotations)
```

### Change log:
The differences between the CRAN and GitHub versions of the *rotations* package are as follows:

* The `]` operator has been redefined to maintain the `SO3` or `Q4` class of the object.

* Addition `+` and subtraction `-` have been redefined for the multiplicative group $SO(3)$.  That is, for $R_1$ and $R_2$ in $SO(3)$, $R_1+R_2=R_2R_1$, $R_1-R_2=R_2^\top R_1$ and $-R_1=R_1^\top$.

* The `plot` function can now print multiple columns of a sample of rotations simultaneously through the argument `col`.  For example `col=c(1,3)` will print two labelled eyeballs, one for the $x$- and one for the $z$-axis.

* Bayesian credible regions have been added to the `region` function by setting `method='Bayes'`.  Additional arguments `tuneS`, `tuneK`, `burnin` and `B` will set the tuning parameters for the central orientation, concentration parameter, burnin for the MCMC and total number of draws from the posterior, respectively.  See the help file for more details.

*intervals*
--------------------------------------------------------
This folder contains the code necessary to run the simulations for the intervals research

