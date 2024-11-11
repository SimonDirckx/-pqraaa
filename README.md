# PQR-AAA
## Description
PQR-AAA is a package to construct Set-valued (aka vector-valued) rational approximations using a greedy approach known as 'AAA' (adaptive Antoulas-Anderson).
##Dependencies
### Eigen
PQR-AAA depends on the [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) library for its linear algebra, as well as a custom variant of a subset of the Eigen Library (can be found in 'pqraaa-support').
## Getting started with PQR-AAA
In the 'examples' folder, 4 reference problems from the [NLEVP](https://github.com/ftisseur/nlevp) toolbox are included. After initial build, these can be compiled by running e.g.
```
make examples/randomDelays
./examples/randomDelays
```
from command line. These 4 examples run both SV-AAA and QR-AAA on parameter dependent systems, and compare the timings.
