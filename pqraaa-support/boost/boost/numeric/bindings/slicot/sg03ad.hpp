//
// Copyright (c) 2011-
// Toon Knapen, Karl Meerbergen, Kresimir Fresl,
// Thomas Klimpel and Rutger ter Borg
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_NUMERIC_BINDINGS_SLICOT_SG03AD_HPP
#define BOOST_NUMERIC_BINDINGS_SLICOT_SG03AD_HPP

#include <boost/numeric/bindings/detail/config/fortran.hpp>
#include <boost/numeric/bindings/detail/array.hpp>
#include <boost/numeric/bindings/lapack/workspace.hpp>

//
// SLICOT computational routines
//

// Value-type variants of hbgst
#define SLICOT_SG03AD FORTRAN_ID( sg03ad )


extern "C" {

//
// LAPACK computational routines
//

// Value-type variants of hbgst
void SLICOT_SG03AD( const char* dico, const char* job, const char* fact, const char* trana, const char* uplo,
        const fortran_int_t* n,
        double* a, const fortran_int_t* lda,
        double* e, const fortran_int_t* lde,
        double* q, const fortran_int_t* ldq,
        double* z, const fortran_int_t* ldz,
        double* x, const fortran_int_t* ldx, double* scale, double* sep, double* ferr,
        double* wr, double* wi, double *beta,
        fortran_int_t* iwork, double* dwork, const fortran_int_t* ldwork,
        fortran_int_t* info ) ;

} // extern "C"

namespace boost {
namespace numeric {
namespace bindings {
namespace slicot {

/*
C     PURPOSE
C
C     To solve for X either the generalized continuous-time Lyapunov
C     equation
C
C             T                T
C        op(A)  X op(E) + op(E)  X op(A) = SCALE * Y,                (1)
C
C     or the generalized discrete-time Lyapunov equation
C
C             T                T
C        op(A)  X op(A) - op(E)  X op(E) = SCALE * Y,                (2)
C
C     where op(M) is either M or M**T for M = A, E and the right hand
C     side Y is symmetric. A, E, Y, and the solution X are N-by-N
C     matrices. SCALE is an output scale factor, set to avoid overflow
C     in X.
C
C     Estimates of the separation and the relative forward error norm
C     are provided.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     DICO    CHARACTER*1
C             Specifies which type of the equation is considered:
C             = 'C':  Continuous-time equation (1);
C             = 'D':  Discrete-time equation (2).
C
C     JOB     CHARACTER*1
C             Specifies if the solution is to be computed and if the
C             separation is to be estimated:
C             = 'X':  Compute the solution only;
C             = 'S':  Estimate the separation only;
C             = 'B':  Compute the solution and estimate the separation.
C
C     FACT    CHARACTER*1
C             Specifies whether the generalized real Schur
C             factorization of the pencil A - lambda * E is supplied
C             on entry or not:
C             = 'N':  Factorization is not supplied;
C             = 'F':  Factorization is supplied.
C
C     TRANS   CHARACTER*1
C             Specifies whether the transposed equation is to be solved
C             or not:
C             = 'N':  op(A) = A,    op(E) = E;
C             = 'T':  op(A) = A**T, op(E) = E**T.
C
C     UPLO    CHARACTER*1
C             Specifies whether the lower or the upper triangle of the
C             array X is needed on input:
C             = 'L':  Only the lower triangle is needed on input;
C             = 'U':  Only the upper triangle is needed on input.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the matrix A.  N >= 0.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, if FACT = 'F', then the leading N-by-N upper
C             Hessenberg part of this array must contain the
C             generalized Schur factor A_s of the matrix A (see
C             definition (3) in section METHOD). A_s must be an upper
C             quasitriangular matrix. The elements below the upper
C             Hessenberg part of the array A are not referenced.
C             If FACT = 'N', then the leading N-by-N part of this
C             array must contain the matrix A.
C             On exit, the leading N-by-N part of this array contains
C             the generalized Schur factor A_s of the matrix A. (A_s is
C             an upper quasitriangular matrix.)
C
C     LDA     INTEGER
C             The leading dimension of the array A.  LDA >= MAX(1,N).
C
C     E       (input/output) DOUBLE PRECISION array, dimension (LDE,N)
C             On entry, if FACT = 'F', then the leading N-by-N upper
C             triangular part of this array must contain the
C             generalized Schur factor E_s of the matrix E (see
C             definition (4) in section METHOD). The elements below the
C             upper triangular part of the array E are not referenced.
C             If FACT = 'N', then the leading N-by-N part of this
C             array must contain the coefficient matrix E of the
C             equation.
C             On exit, the leading N-by-N part of this array contains
C             the generalized Schur factor E_s of the matrix E. (E_s is
C             an upper triangular matrix.)
C
C     LDE     INTEGER
C             The leading dimension of the array E.  LDE >= MAX(1,N).
C
C     Q       (input/output) DOUBLE PRECISION array, dimension (LDQ,N)
C             On entry, if FACT = 'F', then the leading N-by-N part of
C             this array must contain the orthogonal matrix Q from
C             the generalized Schur factorization (see definitions (3)
C             and (4) in section METHOD).
C             If FACT = 'N', Q need not be set on entry.
C             On exit, the leading N-by-N part of this array contains
C             the orthogonal matrix Q from the generalized Schur
C             factorization.
C
C     LDQ     INTEGER
C             The leading dimension of the array Q.  LDQ >= MAX(1,N).
C
C     Z       (input/output) DOUBLE PRECISION array, dimension (LDZ,N)
C             On entry, if FACT = 'F', then the leading N-by-N part of
C             this array must contain the orthogonal matrix Z from
C             the generalized Schur factorization (see definitions (3)
C             and (4) in section METHOD).
C             If FACT = 'N', Z need not be set on entry.
C             On exit, the leading N-by-N part of this array contains
C             the orthogonal matrix Z from the generalized Schur
C             factorization.
C
C     LDZ     INTEGER
C             The leading dimension of the array Z.  LDZ >= MAX(1,N).
C
C     X       (input/output) DOUBLE PRECISION array, dimension (LDX,N)
C             On entry, if JOB = 'B' or 'X', then the leading N-by-N
C             part of this array must contain the right hand side matrix
C             Y of the equation. Either the lower or the upper
C             triangular part of this array is needed (see mode
C             parameter UPLO).
C             If JOB = 'S', X is not referenced.
C             On exit, if JOB = 'B' or 'X', and INFO = 0, 3, or 4, then
C             the leading N-by-N part of this array contains the
C             solution matrix X of the equation.
C             If JOB = 'S', X is not referenced.
C
C     LDX     INTEGER
C             The leading dimension of the array X.  LDX >= MAX(1,N).
C
C     SCALE   (output) DOUBLE PRECISION
C             The scale factor set to avoid overflow in X.
C             (0 < SCALE <= 1)
C
C     SEP     (output) DOUBLE PRECISION
C             If JOB = 'S' or JOB = 'B', and INFO = 0, 3, or 4, then
C             SEP contains an estimate of the separation of the
C             Lyapunov operator.
C
C     FERR    (output) DOUBLE PRECISION
C             If JOB = 'B', and INFO = 0, 3, or 4, then FERR contains an
C             estimated forward error bound for the solution X. If XTRUE
C             is the true solution, FERR estimates the relative error
C             in the computed solution, measured in the Frobenius norm:
C             norm(X - XTRUE) / norm(XTRUE)
C
C     ALPHAR  (output) DOUBLE PRECISION array, dimension (N)
C     ALPHAI  (output) DOUBLE PRECISION array, dimension (N)
C     BETA    (output) DOUBLE PRECISION array, dimension (N)
C             If FACT = 'N' and INFO = 0, 3, or 4, then
C             (ALPHAR(j) + ALPHAI(j)*i)/BETA(j), j=1,...,N, are the
C             eigenvalues of the matrix pencil A - lambda * E.
C             If FACT = 'F', ALPHAR, ALPHAI, and BETA are not
C             referenced.
*/


  struct sg03ad_info {
    fortran_int_t info ;
    double scale ;
    double sep ;
    double ferr ;
  } ;

  template <typename A, typename E, typename Q, typename Z, typename X, typename WR, typename WI, typename Beta>
  sb03md_info sg03ad( char dico, char job, char fact, char trana, A& a, E& e, Q& q, Z& z, X& x, WR& alphar, WI& alphai, Beta& beta, lapack::minimal_workspace ) {
    BOOST_STATIC_ASSERT( (bindings::is_column_major< A >::value) );
    BOOST_STATIC_ASSERT( (bindings::is_column_major< U >::value) );
    BOOST_STATIC_ASSERT( (bindings::is_column_major< C >::value) );
    assert( bindings::size_column(a)==bindings::size_row(a) ) ;
    assert( bindings::size_column(u)==bindings::size_row(u) ) ;
    assert( bindings::size_column(c)==bindings::size_row(c) ) ;
    assert( bindings::size_column(c)==bindings::size_column(a) ) ;
    assert( bindings::size_column(u)==bindings::size_column(a) ) ;
    assert( bindings::size(wr)==bindings::size(wi) ) ;
    assert( bindings::size(wr)==bindings::size_row(a) ) ;
    assert( job=='X' || job=='S' || job=='B' ) ;
    assert( fact=='F' || fact=='N' ) ;

    sb03md_info info ;
    int n = bindings::size_row(a) ;

    int liwork = 1 ;
    if (job!='X')
      liwork = n*n ;
    bindings::detail::array< fortran_int_t > iwork( liwork ) ;
    int lwork = 1 ;
    if (job=='X') {
      if (fact=='F') {
        if (dico=='C') lwork = n*n ;
        else lwork = std::max(n*n,2*n) ;
      } else/*if (fact=='N')*/ {
        lwork = std::max(n*n,3*n) ;
      }
    } else {
      if (fact=='F') {
        if (dico=='C') lwork = n*n ;
        else lwork = 2*n*n+2*n ;
      } else/*if (fact=='N')*/ {
        if (dico=='C') lwork = std::max(2*n*n,3*n) ;
        else lwork = 2*n*n+2*n ;
      }
    }
    bindings::detail::array< double > dwork( lwork ) ;
    fortran_int_t lda( bindings::stride_major(a) ) ;
    fortran_int_t ldu( bindings::stride_major(u) ) ;
    fortran_int_t ldc( bindings::stride_major(c) ) ;
    SLICOT_SG03AD( &dico, &job, &fact, &trana, &n, bindings::begin_value(a), &lda
                 , bindings::begin_value(u), &ldu
                 , bindings::begin_value(c), &ldc
                 , &info.scale, &info.sep, &info.ferr, bindings::begin_value(wr), bindings::begin_value(wi)
                 , bindings::begin_value(iwork)
                 , bindings::begin_value(dwork), &lwork, &info.info ) ;
    return info ;
  } // sb03od()

} // namespace slicot
} // namespace bindings
} // namespace numeric
} // namespace boost

#endif
