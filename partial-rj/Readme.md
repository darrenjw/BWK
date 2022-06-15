# Boys' RJ-MCMC algorithm for exact inference

## Fortran code

I (DJW) have added a `Makefile` for slightly simplifying the building and running of the code.

Just typing `make` should do a test build and run of a toy problem.

After that, you need to carefully inspect the Makefile and the input files to figure out how to run the particular scenarios of interest. eg. `make rjb40-epi.tab` will use `in6` to run with 40 intervals, doing the partially observed case with the initial state fixed.

