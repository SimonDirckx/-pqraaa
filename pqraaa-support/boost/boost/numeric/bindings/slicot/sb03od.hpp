//
// Copyright (c) 2012-
// Toon Knapen, Karl Meerbergen, Kresimir Fresl,
// Thomas Klimpel and Rutger ter Borg
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_NUMERIC_BINDINGS_SLICOT_SB03OD_HPP
#define BOOST_NUMERIC_BINDINGS_SLICOT_SB03OD_HPP

#include <boost/numeric/bindings/slicot/tools.hpp>
#include <boost/numeric/bindings/detail/config/fortran.hpp>
#include <boost/numeric/bindings/detail/array.hpp>
#include <boost/numeric/bindings/lapack/workspace.hpp>

//
// SLICOT computational routines
//

// Value-type variants of hbgst
#define SLICOT_SB03OD FORTRAN_ID( sb03od )


extern "C" {

//
// LAPACK computational routines
//

// Value-type variants of hbgst
void SLICOT_SB03OD( const char* dico, const char* fact, const char* trans,
        const fortran_int_t* n, const fortran_int_t* m,
        double* a, const fortran_int_t* lda,
        double* q, const fortran_int_t* ldq,
        double* b, const fortran_int_t* ldb, double* scale,
        double* wr, double* wi,
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
C     To solve for X = op(U)'*op(U) either the stable non-negative
C     definite continuous-time Lyapunov equation
C                                   2
C        op(A)'*X + X*op(A) = -scale *op(B)'*op(B)                   (1)
C
C     or the convergent non-negative definite discrete-time Lyapunov
C     equation
C                                   2
C        op(A)'*X*op(A) - X = -scale *op(B)'*op(B)                   (2)
C
C     where op(K) = K or K' (i.e., the transpose of the matrix K), A is
C     an N-by-N matrix, op(B) is an M-by-N matrix, U is an upper
C     triangular matrix containing the Cholesky factor of the solution
C     matrix X, X = op(U)'*op(U), and scale is an output scale factor,
C     set less than or equal to 1 to avoid overflow in X. If matrix B
C     has full rank then the solution matrix X will be positive-definite
C     and hence the Cholesky factor U will be nonsingular, but if B is
C     rank deficient then X may be only positive semi-definite and U
C     will be singular.
C
C     In the case of equation (1) the matrix A must be stable (that
C     is, all the eigenvalues of A must have negative real parts),
C     and for equation (2) the matrix A must be convergent (that is,
C     all the eigenvalues of A must lie inside the unit circle).
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     DICO    CHARACTER*1
C             Specifies the type of Lyapunov equation to be solved as
C             follows:
C             = 'C':  Equation (1), continuous-time case;
C             = 'D':  Equation (2), discrete-time case.
C
C     FACT    CHARACTER*1
C             Specifies whether or not the real Schur factorization
C             of the matrix A is supplied on entry, as follows:
C             = 'F':  On entry, A and Q contain the factors from the
C                     real Schur factorization of the matrix A;
C             = 'N':  The Schur factorization of A will be computed
C                     and the factors will be stored in A and Q.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the matrix A and the number of columns in
C             matrix op(B).  N >= 0.
C
C     M       (input) INTEGER
C             The number of rows in matrix op(B).  M >= 0.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading N-by-N part of this array must
C             contain the matrix A. If FACT = 'F', then A contains
C             an upper quasi-triangular matrix S in Schur canonical
C             form; the elements below the upper Hessenberg part of the
C             array A are not referenced.
C             On exit, the leading N-by-N upper Hessenberg part of this
C             array contains the upper quasi-triangular matrix S in
C             Schur canonical form from the Shur factorization of A.
C             The contents of array A is not modified if FACT = 'F'.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,N).
C
C     Q       (input or output) DOUBLE PRECISION array, dimension
C             (LDQ,N)
C             On entry, if FACT = 'F', then the leading N-by-N part of
C             this array must contain the orthogonal matrix Q of the
C             Schur factorization of A.
C             Otherwise, Q need not be set on entry.
C             On exit, the leading N-by-N part of this array contains
C             the orthogonal matrix Q of the Schur factorization of A.
C             The contents of array Q is not modified if FACT = 'F'.
C
C     LDQ     INTEGER
C             The leading dimension of array Q.  LDQ >= MAX(1,N).
C
C     B       (input/output) DOUBLE PRECISION array, dimension (LDB,N)
C             if TRANS = 'N', and dimension (LDB,max(M,N)), if
C             TRANS = 'T'.
C             On entry, if TRANS = 'N', the leading M-by-N part of this
C             array must contain the coefficient matrix B of the
C             equation.
C             On entry, if TRANS = 'T', the leading N-by-M part of this
C             array must contain the coefficient matrix B of the
C             equation.
C             On exit, the leading N-by-N part of this array contains
C             the upper triangular Cholesky factor U of the solution
C             matrix X of the problem, X = op(U)'*op(U).
C             If M = 0 and N > 0, then U is set to zero.
C
C     LDB     INTEGER
C             The leading dimension of array B.
C             LDB >= MAX(1,N,M), if TRANS = 'N';
C             LDB >= MAX(1,N),   if TRANS = 'T'.
C
C     SCALE   (output) DOUBLE PRECISION
C             The scale factor, scale, set less than or equal to 1 to
C             prevent the solution overflowing.
C
C     WR      (output) DOUBLE PRECISION array, dimension (N)
C     WI      (output) DOUBLE PRECISION array, dimension (N)
C             If FACT = 'N', and INFO >= 0 and INFO <= 2, WR and WI
C             contain the real and imaginary parts, respectively, of
C             the eigenvalues of A.
C             If FACT = 'F', WR and WI are not referenced.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, or INFO = 1, DWORK(1) returns the
C             optimal value of LDWORK.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             If M > 0, LDWORK >= MAX(1,4*N + MIN(M,N));
C             If M = 0, LDWORK >= 1.
C             For optimum performance LDWORK should sometimes be larger.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 1:  if the Lyapunov equation is (nearly) singular
C                   (warning indicator);
C                   if DICO = 'C' this means that while the matrix A
C                   (or the factor S) has computed eigenvalues with
C                   negative real parts, it is only just stable in the
C                   sense that small perturbations in A can make one or
C                   more of the eigenvalues have a non-negative real
C                   part;
C                   if DICO = 'D' this means that while the matrix A
C                   (or the factor S) has computed eigenvalues inside
C                   the unit circle, it is nevertheless only just
C                   convergent, in the sense that small perturbations
C                   in A can make one or more of the eigenvalues lie
C                   outside the unit circle;
C                   perturbed values were used to solve the equation;
C             = 2:  if FACT = 'N' and DICO = 'C', but the matrix A is
C                   not stable (that is, one or more of the eigenvalues
C                   of A has a non-negative real part), or DICO = 'D',
C                   but the matrix A is not convergent (that is, one or
C                   more of the eigenvalues of A lies outside the unit
C                   circle); however, A will still have been factored
C                   and the eigenvalues of A returned in WR and WI.
C             = 3:  if FACT = 'F' and DICO = 'C', but the Schur factor S
C                   supplied in the array A is not stable (that is, one
C                   or more of the eigenvalues of S has a non-negative
C                   real part), or DICO = 'D', but the Schur factor S
C                   supplied in the array A is not convergent (that is,
C                   one or more of the eigenvalues of S lies outside the
C                   unit circle);
C             = 4:  if FACT = 'F' and the Schur factor S supplied in
C                   the array A has two or more consecutive non-zero
C                   elements on the first sub-diagonal, so that there is
C                   a block larger than 2-by-2 on the diagonal;
C             = 5:  if FACT = 'F' and the Schur factor S supplied in
C                   the array A has a 2-by-2 diagonal block with real
C                   eigenvalues instead of a complex conjugate pair;
C             = 6:  if FACT = 'N' and the LAPACK Library routine DGEES
C                   has failed to converge. This failure is not likely
C                   to occur. The matrix B will be unaltered but A will
C                   be destroyed.
C
*/


  struct sb03od_info {
    fortran_int_t info ;
    double scale ;
    double sep ;
    double ferr ;
  } ;

  template <typename A, typename Q, typename B, typename WR, typename WI>
  sb03md_info sb03od( char dico, char fact, A& a, Q& q, double const& scale, B& b, WR& wr, WI& wi, lapack::minimal_workspace ) {
    BOOST_STATIC_ASSERT( (bindings::is_column_major< A >::value) );
    BOOST_STATIC_ASSERT( (bindings::is_column_major< Q >::value) );
    BOOST_STATIC_ASSERT( (bindings::is_column_major< B >::value) );
    BOOST_STATIC_CONSTANT( char, trans = detail::trans<A>::value ) ;
    BOOST_STATIC_ASSERT( (trans==detail::trans<B>::value) );
    assert( bindings::size_column(a)==bindings::size_row(a) ) ;
    assert( bindings::size_column(q)==bindings::size_row(q) ) ;
    assert( bindings::size_row(b)==bindings::size_column(a) ) ;
    assert( bindings::size_column(q)==bindings::size_column(a) ) ;
    assert( bindings::size(wr)==bindings::size(wi) ) ;
    assert( bindings::size(wr)==bindings::size_row(a) ) ;
    assert( dico=='C' || job=='D' ) ;
    assert( fact=='F' || fact=='N' ) ;

    sb03md_info info ;
    int n = bindings::size_row(a) ;
    int m = bindings::size_row(b) ;

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
    std::ptr_diff_t lwork = std::max*(1, 4*n + std::min(m,n));
    bindings::detail::array< double > dwork( lwork ) ;
    fortran_int_t lda( bindings::stride_major(a) ) ;
    fortran_int_t ldq( bindings::stride_major(q) ) ;
    fortran_int_t ldb( bindings::stride_major(b) ) ;
    SLICOT_SB03MD( &dico, &fact, &trans, &n, bindings::begin_value(a), &lda
                 , bindings::begin_value(q), &ldq
                 , bindings::begin_value(b), &ldb
                 , &info.scale, bindings::begin_value(wr), bindings::begin_value(wi)
                 , bindings::begin_value(dwork), &lwork, &info.info ) ;
    return info ;
  } // sb03od()

} // namespace slicot
} // namespace bindings
} // namespace numeric
} // namespace boost

#endif
