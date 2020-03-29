# BWK

## Code from Boys, Wilkinson and Kirkwood (2008)

Repository to share code from:

Boys, R. J., Wilkinson, D. J., Kirkwood, T. B. L. (2008) Bayesian inference for a discretely observed stochastic kinetic model, Statistics and Computing, 18(2):125-135. http://dx.doi.org/10.1007/s11222-007-9043-x

Initially just Wilkinson's C code for block updating MCMC algorithm.

Building this requires the [GSL](https://www.gnu.org/software/gsl/). On Debian-derived Linux distributions, this should just require the packages `gsl-bin` and `libgsl-dev`. There will be similar packages on other Linux and Unix-like systems. Then, just typing `make` in the relevant directory should build and run the code. See the `Makefile` for further details.


I hope to find Boys' Fortran code for RJ-MCMC algorithm at some point.


