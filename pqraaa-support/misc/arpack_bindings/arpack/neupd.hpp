//
// Copyright (c) 2015
// Karl Meerbergen
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
//

#ifndef ARPACK_NEUPD_HPP
#define ARPACK_NEUPD_HPP

#include <arpack/data.hpp>
#include <boost/numeric/bindings/detail/config/fortran.hpp>
#include <boost/numeric/bindings/value_type.hpp>
#include <boost/numeric/bindings/size.hpp>
#include <boost/numeric/bindings/begin.hpp>
#include <cassert>

//
// ARPACK computational routines
//

/*
c\Description: 
c
c  This subroutine returns the converged approximations to eigenvalues
c  of A*z = lambda*B*z and (optionally):
c
c      (1) The corresponding approximate eigenvectors;
c
c      (2) An orthonormal basis for the associated approximate
c          invariant subspace;
c
c      (3) Both.
c
c  There is negligible additional cost to obtain eigenvectors.  An orthonormal
c  basis is always computed.  There is an additional storage cost of n*nev
c  if both are requested (in this case a separate array Z must be supplied).
c
c  The approximate eigenvalues and eigenvectors of  A*z = lambda*B*z
c  are derived from approximate eigenvalues and eigenvectors of
c  of the linear operator OP prescribed by the MODE selection in the
c  call to DNAUPD.  DNAUPD must be called before this routine is called.
c  These approximate eigenvalues and vectors are commonly called Ritz
c  values and Ritz vectors respectively.  They are referred to as such
c  in the comments that follow.  The computed orthonormal basis for the
c  invariant subspace corresponding to these Ritz values is referred to as a
c  Schur basis.
c
c  See documentation in the header of the subroutine DNAUPD for 
c  definition of OP as well as other terms and the relation of computed
c  Ritz values and Ritz vectors of OP with respect to the given problem
c  A*z = lambda*B*z.  For a brief description, see definitions of 
c  IPARAM(7), MODE and WHICH in the documentation of DNAUPD.
c
c\Usage:
c  call dneupd 
c     ( RVEC, HOWMNY, SELECT, DR, DI, Z, LDZ, SIGMAR, SIGMAI, WORKEV, BMAT, 
c       N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM, IPNTR, WORKD, WORKL, 
c       LWORKL, INFO )
c
c\Arguments:
c  RVEC    LOGICAL  (INPUT) 
c          Specifies whether a basis for the invariant subspace corresponding 
c          to the converged Ritz value approximations for the eigenproblem 
c          A*z = lambda*B*z is computed.
c
c             RVEC = .FALSE.     Compute Ritz values only.
c
c             RVEC = .TRUE.      Compute the Ritz vectors or Schur vectors.
c                                See Remarks below. 
c 
c  HOWMNY  Character*1  (INPUT) 
c          Specifies the form of the basis for the invariant subspace 
c          corresponding to the converged Ritz values that is to be computed.
c
c          = 'A': Compute NEV Ritz vectors; 
c          = 'P': Compute NEV Schur vectors;
c          = 'S': compute some of the Ritz vectors, specified
c                 by the logical array SELECT.
c
c  SELECT  Logical array of dimension NCV.  (INPUT)
c          If HOWMNY = 'S', SELECT specifies the Ritz vectors to be
c          computed. To select the Ritz vector corresponding to a
c          Ritz value (DR(j), DI(j)), SELECT(j) must be set to .TRUE.. 
c          If HOWMNY = 'A' or 'P', SELECT is used as internal workspace.
c
c  DR      Double precision array of dimension NEV+1.  (OUTPUT)
c          If IPARAM(7) = 1,2 or 3 and SIGMAI=0.0  then on exit: DR contains 
c          the real part of the Ritz  approximations to the eigenvalues of 
c          A*z = lambda*B*z. 
c          If IPARAM(7) = 3, 4 and SIGMAI is not equal to zero, then on exit:
c          DR contains the real part of the Ritz values of OP computed by 
c          DNAUPD. A further computation must be performed by the user
c          to transform the Ritz values computed for OP by DNAUPD to those
c          of the original system A*z = lambda*B*z. See remark 3 below.
c
c  DI      Double precision array of dimension NEV+1.  (OUTPUT)
c          On exit, DI contains the imaginary part of the Ritz value 
c          approximations to the eigenvalues of A*z = lambda*B*z associated
c          with DR.
c
c          NOTE: When Ritz values are complex, they will come in complex 
c                conjugate pairs.  If eigenvectors are requested, the 
c                corresponding Ritz vectors will also come in conjugate 
c                pairs and the real and imaginary parts of these are 
c                represented in two consecutive columns of the array Z 
c                (see below).
c
c  Z       Double precision N by NEV+1 array if RVEC = .TRUE. and HOWMNY = 'A'. (OUTPUT)
c          On exit, if RVEC = .TRUE. and HOWMNY = 'A', then the columns of 
c          Z represent approximate eigenvectors (Ritz vectors) corresponding 
c          to the NCONV=IPARAM(5) Ritz values for eigensystem 
c          A*z = lambda*B*z. 
c 
c          The complex Ritz vector associated with the Ritz value 
c          with positive imaginary part is stored in two consecutive 
c          columns.  The first column holds the real part of the Ritz 
c          vector and the second column holds the imaginary part.  The 
c          Ritz vector associated with the Ritz value with negative 
c          imaginary part is simply the complex conjugate of the Ritz vector 
c          associated with the positive imaginary part.
c
c          If  RVEC = .FALSE. or HOWMNY = 'P', then Z is not referenced.
c
c          NOTE: If if RVEC = .TRUE. and a Schur basis is not required,
c          the array Z may be set equal to first NEV+1 columns of the Arnoldi
c          basis array V computed by DNAUPD.  In this case the Arnoldi basis
c          will be destroyed and overwritten with the eigenvector basis.
c
c  LDZ     Integer.  (INPUT)
c          The leading dimension of the array Z.  If Ritz vectors are
c          desired, then  LDZ >= max( 1, N ).  In any case,  LDZ >= 1.
c
c  SIGMAR  Double precision  (INPUT)
c          If IPARAM(7) = 3 or 4, represents the real part of the shift. 
c          Not referenced if IPARAM(7) = 1 or 2.
c
c  SIGMAI  Double precision  (INPUT)
c          If IPARAM(7) = 3 or 4, represents the imaginary part of the shift. 
c          Not referenced if IPARAM(7) = 1 or 2. See remark 3 below.
c
c  WORKEV  Double precision work array of dimension 3*NCV.  (WORKSPACE)
c
c  **** The remaining arguments MUST be the same as for the   ****
c  **** call to DNAUPD that was just completed.               ****
c
c  NOTE: The remaining arguments
c
c           BMAT, N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM, IPNTR,
c           WORKD, WORKL, LWORKL, INFO
c
c         must be passed directly to DNEUPD following the last call
c         to DNAUPD.  These arguments MUST NOT BE MODIFIED between
c         the the last call to DNAUPD and the call to DNEUPD.
c
c  Three of these parameters (V, WORKL, INFO) are also output parameters:
c
c  V       Double precision N by NCV array.  (INPUT/OUTPUT)
c
c          Upon INPUT: the NCV columns of V contain the Arnoldi basis
c                      vectors for OP as constructed by DNAUPD .
c
c          Upon OUTPUT: If RVEC = .TRUE. the first NCONV=IPARAM(5) columns
c                       contain approximate Schur vectors that span the
c                       desired invariant subspace.  See Remark 2 below.
c
c          NOTE: If the array Z has been set equal to first NEV+1 columns
c          of the array V and RVEC=.TRUE. and HOWMNY= 'A', then the
c          Arnoldi basis held by V has been overwritten by the desired
c          Ritz vectors.  If a separate array Z has been passed then
c          the first NCONV=IPARAM(5) columns of V will contain approximate
c          Schur vectors that span the desired invariant subspace.
c
c  WORKL   Double precision work array of length LWORKL.  (OUTPUT/WORKSPACE)
c          WORKL(1:ncv*ncv+3*ncv) contains information obtained in
c          dnaupd.  They are not changed by dneupd.
c          WORKL(ncv*ncv+3*ncv+1:3*ncv*ncv+6*ncv) holds the
c          real and imaginary part of the untransformed Ritz values,
c          the upper quasi-triangular matrix for H, and the
c          associated matrix representation of the invariant subspace for H.
c
c          Note: IPNTR(9:13) contains the pointer into WORKL for addresses
c          of the above information computed by dneupd.
c          -------------------------------------------------------------
c          IPNTR(9):  pointer to the real part of the NCV RITZ values of the
c                     original system.
c          IPNTR(10): pointer to the imaginary part of the NCV RITZ values of
c                     the original system.
c          IPNTR(11): pointer to the NCV corresponding error bounds.
c          IPNTR(12): pointer to the NCV by NCV upper quasi-triangular
c                     Schur matrix for H.
c          IPNTR(13): pointer to the NCV by NCV matrix of eigenvectors
c                     of the upper Hessenberg matrix H. Only referenced by
c                     dneupd if RVEC = .TRUE. See Remark 2 below.
c          -------------------------------------------------------------
c
c  INFO    Integer.  (OUTPUT)
c          Error flag on output.
c
c          =  0: Normal exit.
c
c          =  1: The Schur form computed by LAPACK routine dlahqr
c                could not be reordered by LAPACK routine dtrsen.
c                Re-enter subroutine dneupd with IPARAM(5)=NCV and 
c                increase the size of the arrays DR and DI to have 
c                dimension at least dimension NCV and allocate at least NCV 
c                columns for Z. NOTE: Not necessary if Z and V share 
c                the same space. Please notify the authors if this error
c                occurs.
c
c          = -1: N must be positive.
c          = -2: NEV must be positive.
c          = -3: NCV-NEV >= 2 and less than or equal to N.
c          = -5: WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'
c          = -6: BMAT must be one of 'I' or 'G'.
c          = -7: Length of private work WORKL array is not sufficient.
c          = -8: Error return from calculation of a real Schur form.
c                Informational error from LAPACK routine dlahqr.
c          = -9: Error return from calculation of eigenvectors.
c                Informational error from LAPACK routine dtrevc.
c          = -10: IPARAM(7) must be 1,2,3,4.
c          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatible.
c          = -12: HOWMNY = 'S' not yet implemented
c          = -13: HOWMNY must be one of 'A' or 'P' if RVEC = .true.
c          = -14: DNAUPD did not find any eigenvalues to sufficient
c                 accuracy.
*/

// Value-type variants of saupd
#define ARPACK_DNEUPD FORTRAN_ID( dneupd )
#define ARPACK_SNEUPD FORTRAN_ID( sneupd )
#define ARPACK_DNEUPD_C FORTRAN_ID( dneupd_c )
#define ARPACK_SNEUPD_C FORTRAN_ID( sneupd_c )

extern "C" {

void ARPACK_SNEUPD( fortran_bool_t* RVEC, char* HOWMNY, fortran_bool_t* SELECT, float* DR, float* DI, float* Z, int* LDZ
                  , float* SIGMAR, float* SIGMAI, float* WORKEV, char const* BMAT, int const* N, char const* WHICH, const int* NEV, float const* TOL
                  , float* RESID, int const* NCV, float* V, const int* LDV, int* IPARAM, int* IPNTR, float* WORKD
                  , float* WORKL, int const* LWORKL, int* INFO
                  ) ;

void ARPACK_SNEUPD_C( int* RVEC, char* HOWMNY, int* SELECT_C, fortran_bool_t* SELECT, float* DR, float* DI, float* Z, int* LDZ
                  , float* SIGMAR, float* SIGMAI, float* WORKEV, char const* BMAT, int const* N, char const* WHICH, const int* NEV, float const* TOL
                  , float* RESID, int const* NCV, float* V, const int* LDV, int* IPARAM, int* IPNTR, float* WORKD
                  , float* WORKL, int const* LWORKL, int* INFO
                  ) ;

void ARPACK_DNEUPD( fortran_bool_t* RVEC, char* HOWMNY, fortran_bool_t* SELECT, double* DR, double* DI, double* Z, int* LDZ
                  , double* SIGMAR, double* SIGMAI, double* WORKEV, char const* BMAT, int const* N, char const* WHICH, const int* NEV, double const* TOL
                  , double* RESID, int const* NCV, double* V, const int* LDV, int* IPARAM, int* IPNTR, double* WORKD
                  , double* WORKL, int const* LWORKL, int* INFO
                  ) ;

void ARPACK_DNEUPD_C( int* RVEC, char* HOWMNY, int* SELECT_C, fortran_bool_t* SELECT, double* DR, double* DI, double* Z, int* LDZ
                  , double* SIGMAR, double* SIGMAI, double* WORKEV, char const* BMAT, int const* N, char const* WHICH, const int* NEV, double const* TOL
                  , double* RESID, int const* NCV, double* V, const int* LDV, int* IPARAM, int* IPNTR, double* WORKD
                  , double* WORKL, int const* LWORKL, int* INFO
                  ) ;

}

// C++ bindings

namespace ARPACK {

  namespace detail {

    template <typename T>
    struct neupd_driver {
    } ;

    template <>
    struct neupd_driver< float > {
      static void apply( int RVEC, char HOWMNY, int* SELECT_C, fortran_bool_t* SELECT, float* DR, float* DI, float* WORKEV, float* Z, int LDZ
                       , float SIGMAR, float SIGMAI
                       , char BMAT, int N, char const* WHICH, int NEV, float TOL
                       , float* RESID, int NCV, float* V, int LDV, int* IPARAM, int* IPNTR, float* WORKD
                       , float* WORKL, int LWORKL, int& INFO )
      {
        ARPACK_SNEUPD_C( &RVEC, &HOWMNY, SELECT_C, SELECT, DR, DI, Z, &LDZ, &SIGMAR, &SIGMAI, WORKEV, &BMAT, &N, WHICH
                     , &NEV, &TOL, RESID, &NCV, V, &LDV, IPARAM, IPNTR, WORKD, WORKL, &LWORKL, &INFO
                     ) ;
      }
    } ;

    template <>
    struct neupd_driver< double > {
      static void apply( int RVEC, char HOWMNY, int* SELECT_C, fortran_bool_t* SELECT, double* DR, double* DI, double* WORKEV, double* Z, int LDZ
                       , double SIGMAR, double SIGMAI
                       , char BMAT, int N, char const* WHICH, int NEV, double TOL
                       , double* RESID, int NCV, double* V, int LDV, int* IPARAM, int* IPNTR, double* WORKD
                       , double* WORKL, int LWORKL, int& INFO )
      {
        ARPACK_DNEUPD_C( &RVEC, &HOWMNY, SELECT_C, SELECT, DR, DI, Z, &LDZ, &SIGMAR, &SIGMAI, WORKEV, &BMAT, &N, WHICH
                     , &NEV, &TOL, RESID, &NCV, V, &LDV, IPARAM, IPNTR, WORKD, WORKL, &LWORKL, &INFO
                     ) ;
      }
    } ;

  } // detail

  // Binding for eigenvalues of a single matrix, no shift-and-invert.
  template <typename Select, typename SelectC, typename DR, typename DI, typename Z, typename Problem, typename Data, typename Work>
  int neupd( bool RVEC, SelectC& select_c, Select& select, DR& dr, DI& di, Z& z, Problem const& problem, Data& data, Work& workev ) {
    typedef typename Data::value_type value_type ;

    assert( boost::numeric::bindings::size(select)==data.ncv_ ) ;
    assert( boost::numeric::bindings::size(select_c)==data.ncv_ || boost::numeric::bindings::size(select_c)==0 ) ;

    assert( !RVEC || boost::numeric::bindings::size_column(z)==data.nev_+1 ) ;
    assert( boost::numeric::bindings::size(dr)==data.nev_+1 ) ;
    assert( boost::numeric::bindings::size(di)==data.nev_+1 ) ;

    assert( boost::numeric::bindings::size(workev)>=3*data.ncv_ ) ;

    char how_many ;
    if (boost::numeric::bindings::size(select_c)==0) {
      how_many = 'A' ;
    } else {
      how_many = 'S' ;
    }

    detail::neupd_driver<value_type>::apply( int(RVEC), how_many
                                           , boost::numeric::bindings::begin_value(select_c)
                                           , boost::numeric::bindings::begin_value(select)
                                           , boost::numeric::bindings::begin_value(dr)
                                           , boost::numeric::bindings::begin_value(di)
                                           , boost::numeric::bindings::begin_value(workev)
                                           , boost::numeric::bindings::begin_value(z)
                                           , boost::numeric::bindings::stride_major(z), problem.sigma_r(), problem.sigma_i()
                                           , data.bmat, data.n_, data.which_, data.nev_, data.tol_, data.resid_, data.ncv_, data.v_, data.ldv_
                                           , data.iparam_, data.ipntr_, data.workd_, data.workl_, data.lworkl_, data.info_
                                           ) ;
    return data.info_ ;
  } // seupd()

} // namespace ARPACK

#endif

