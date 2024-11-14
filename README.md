# PQR-AAA
## Description
PQR-AAA is a package to construct Set-valued (aka vector-valued) rational approximations using a greedy approach known as 'AAA' (adaptive Antoulas-Anderson).
## Dependencies
### Eigen
PQR-AAA depends on the [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) library for its linear algebra, as well as a custom variant of a subset of the Eigen Library (can be found in 'pqraaa-support').
## Getting started with PQR-AAA
### pre-loaded examples
In the 'examples' folder, 4 reference problems from the [NLEVP](https://github.com/ftisseur/nlevp) toolbox are included. After initial build, these can be compiled by running e.g.
```
make examples/randomDelays
./examples/randomDelays 1 4 1e-4
```
from command line.  The boolean (here '1') determines whether qr-aaa is used. If '0', regular sv-aaa is used. The input '4' determines the number of cores used (1 corresponding to serial qr-aaa). The input '1e-4' determines the tolerance. NOTE: the current implementation still exhibits some false sharing.
### computing your own QR-AAA approximation
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
Mat<complex<double>> F = Mat<complex<double>>::Zero(nZ,N);
//fill F in some way...
```

* After possibly rescaling the columns of F (application dependent, but unit norm for each column is usually a good idea), set tolerance, set max_degree, number of cores, whether to use qr (always true when n_cores>1) and compute the qr-aaa approximation
```
QRAAA::infoType info;
QRAAA::AAAopts opts;
opts.tol = tol;
opts.max_degree = 30;
opts.qr = true;
opts.n_cores = 1;
auto repr_f=QRAAA::sv_aaa(F,Z,opts,info);
```
* Optionally, output the info of the qr-aaa approximation
```
QRAAA::summarize(info);
```
* The nodes, weights, and coefficients (as glas2 objects) of the qr-aaa approximation can be accessed by
```
auto nodes = repr_f.nodes()
auto weights = repr_f.weights()
auto coeffs = repr_f.coefficients()
```
