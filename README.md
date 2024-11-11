# PQR-AAA
## Description
PQR-AAA is a package to construct Set-valued (aka vector-valued) rational approximations using a greedy approach known as 'AAA' (adaptive Antoulas-Anderson).
##Dependencies
### Eigen
PQR-AAA depends on the [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) library for its linear algebra, as well as a custom variant of a subset of the Eigen Library (can be found in 'pqraaa-support').
## Getting started with PQR-AAA
# pre-loaded examples
In the 'examples' folder, 4 reference problems from the [NLEVP](https://github.com/ftisseur/nlevp) toolbox are included. After initial build, these can be compiled by running e.g.
```
make examples/randomDelays
./examples/randomDelays
```
from command line. These 4 examples run both SV-AAA and QR-AAA on parameter dependent systems, and compare the timings.
# computing your own QR-AAA approximation
Creating your own QR-AAA example can be done in 3 steps:
* Include the correct headers and , for simplicity, define a number of aliases
```
#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <PQRAAA/sv-aaa/qraaa.hh>

using namespace std;
using namespace Eigen;
using namespace boost;

template< typename Tval >
using Mat = Eigen::Matrix<Tval,Eigen::Dynamic,Eigen::Dynamic>;

template< typename Tval >
using Vec = Eigen::Vector<Tval,Eigen::Dynamic>;
```
* Create the function matrix F: an nZxN matrix containing, as columns, N functions discretized on on a set of nZ points in the complex plane
```
//example: 100000 functions discretized on equispaced grid of 1000 points
int N = 100000;
int nZ = 1000;
Eigen::ArrayXd Z = Eigen::ArrayXd::LinSpaced(nZ,-1.,1.);
Mat<CTval> F = Mat<CTval>::Zero(nZ,N);
//fill F in some way...
```
