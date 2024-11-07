# pqraaa-support

## Description

This repository contains supporting packages on which [pqraaa](https://github.com/SimonDirckx/-pqraaa) (partly) depends. This is purely for convenience purposes to aggregate all these dependencies in one place and guarantee availability of package versions that are compatible with BEACHpack.

## Packages

We briefly list all the packages.

### Eigen

For linear algebra, BEACHpack uses the C++ template library [Eigen](https://eigen.tuxfamily.org). The [*eigen-3.4.0* folder](eigen-3.4.0/) contains the 3.4.0 version. 

### CORK++

The PQR-AAA algorithm relies on the SV-AAA algorithm implemented in the CORK++ package, located in the [*cork* folder](cork/). There is no real indication which version of CORK++ this is, as this was lost to time. A seemingly more up-to-date version is available at [scm.cs.kuleuven.be/scm/git/cork](https://scm.cs.kuleuven.be/scm/git/cork) under directory *cork* and should also work just fine.

### GLAS2

CORK++ uses the GLAS2-package for linear algebra instead of Eigen. This package is located in [*glas2*](glas2/) and the current version originates from the one accompanying the [CORK repository](https://scm.cs.kuleuven.be/scm/git/cork) previously mentioned, under directory *glas2*.

### Boost.Numeric.Bindings

To be able to use BLAS and LAPACK, CORK++ and GLAS2 need bindings to them. This is supported by Boost.Numeric.Bindings in the [*boost/numeric/bindings* folder](boost/boost/numeric/bindings/). The version in this folder is a modification of the version accompanying the [CORK repository](https://scm.cs.kuleuven.be/scm/git/cork) (under *boost/numeric/bindings*). Details on the modifications are found at [gitlab.kuleuven.be/u0156790/lapack-c-bindings](https://gitlab.kuleuven.be/u0156790/lapack-c-bindings) where the same version is found under *include/boost-new*.

### MUMPS

While normally not necessary in PQR-AAA, MUMPS allows efficient operations with sparse matrices. CORK++ relies on MUMPS to support sparse matrices. We provide it in this repository under [*MUMPS_5.0.1*](MUMPS_5.0.1/) for completenesss sake.
