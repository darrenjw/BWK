# BWK

## Code from Boys, Wilkinson and Kirkwood (2008)

Repository to share code from:

Boys, R. J., Wilkinson, D. J., Kirkwood, T. B. L. (2008) Bayesian inference for a discretely observed stochastic kinetic model, *Statistics and Computing*, **18**(2):125-135. http://dx.doi.org/10.1007/s11222-007-9043-x

Wilkinson's C code for the exact block updating MCMC algorithm is in the directories `exact` and `partial`. Building this requires the [GSL](https://www.gnu.org/software/gsl/). On Debian-derived Linux distributions, this should just require the packages `gsl-bin` and `libgsl-dev`. There will be similar packages on other Linux and Unix-like systems. Then, just typing `make` in the relevant directory should build and run the code. See the `Makefile` and `Readme` in each directory for further details.

Wilkinson's approximate code can be found in the directories `approx` and `partial-approx`. Again, see the `Readmes` and `Makefiles` for further details.

I have now found Boys' main Fortran 77 code for a RJ-MCMC approach. It is in the directory `partial-rj`. I've managed to tweak it to get it to build and run. Seems to be running correctly. It should be easy to build and run on any Linux-like system with `f77` installed. Type `make`. Again, see the `Readme` in that directory for further details.




