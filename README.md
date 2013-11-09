*rotations*
========================================================
A stable version of the *rotations* package is available for download from CRAN.  For MAC users, the next version of *rotations* can be downloaded from GitHub.  Windows users will need to download the repo, compile the C++ code and install locally.

### Installation instructions for Mac: 
```
library(devtools)
install_github('rotationsC','stanfill',subdir='rotations')
library(rotations)
```

### Change log:
The differences between the CRAN and GitHub versions of the *rotations* package are as follows:

* The `'['` operator has been redefined to maintain the `SO3` or `Q4` class of the object.

* Addition `+` and subtraction `-` have been redefined for the multiplicative group $SO(3)$.  That is, for $R_1$ and $R_2$ in $SO(3)$, $R_1+R_2=R_2R_1$, $R_1-R_2=R_2^\top R_1$ and $-R_1=R_1^\top$.

* The `plot` function can now print multiple columns of a sample of rotations simultaneously through the argument `col`.  For example `col=c(1,3)` will print two labelled eyeballs, one for the $x$- and one for the $z$-axis.

* The `head` and `str` functions have been modified to properly handle objects of class `'SO3'` and `'Q4'`.

* Quaternions are unique up to their sign, that is if q is a quaternion the q=-q.  So `==` has be redifined such that `q==-q` will return `TRUE`.

* Functions for Bayesian analysis
  * Each of the following functions requires the user to supply `type`, `tuneS`, `tuneK`, `burnin` and `B`.  These arguments  determine the transition probability form, tuning parameters for the central orientation, concentration parameter, burnin for the MCMC and total number of draws from the posterior, respectively.  Currently only `'Cayley'`, `'Fisher'` and `'Mises'` are valid options for `type`.  See the respective help files for more details.
  * The function `MCMCSO3` implements a Gibbs-within-MCMC algorithm to get draws from the posterior distributions for the concentration parameter $S$ and concentration parameter $\kappa$.  A list of four items is returned: `S` is a B-by-p matrix where each row corresponds to a draw from the posterior for the central orientation S, `kappa` is a vector of B draws from the posterior for the concentration parameter $\kappa$ and the transition probabilities for the central orientation and concentration are given by `Saccept` and `Kaccept`,respectively.
  * Bayesian point estimates for the central orientation and concentration parameter are given by `bayes.mean`.
  * Bayesian credible regions have been added to the `region` by setting `method='Bayes'`.  

*intervals*
--------------------------------------------------------
This folder contains the code necessary to run the simulations for the intervals research

