//
// Copyright (c) 2012-
// Toon Knapen, Karl Meerbergen, Kresimir Fresl,
// Thomas Klimpel and Rutger ter Borg
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_NUMERIC_BINDINGS_SLICOT_SB04MD_HPP
#define BOOST_NUMERIC_BINDINGS_SLICOT_SB04MD_HPP

#include <boost/numeric/bindings/detail/config/fortran.hpp>
#include <boost/numeric/bindings/detail/array.hpp>
#include <boost/numeric/bindings/lapack/workspace.hpp>

//
// SLICOT computational routines
//

// Value-type variants of hbgst
#define SLICOT_SB04MD FORTRAN_ID( sb04md )


extern "C" {

//
// LAPACK computational routines
//

// Value-type variants of hbgst
void SLICOT_SB04PD( const char* dico, const char* facta, const char* factb, const char* trana, const char* trana, const fortran_int_t const* isgn,
        const fortran_int_t* m, const fortran_int_t* n,
        double* a, const fortran_int_t* lda,
        double* u, const fortran_int_t* ldu,
        double* b, const fortran_int_t* ldb,
        double* v, const fortran_int_t* ldv,
        double* c, const fortran_int_t* ldc,
        double* scale,
        double* dwork, fortran_int_t* ldwork,
        fortran_int_t* info ) ;

} // extern "C"

namespace boost {
namespace numeric {
namespace bindings {
namespace slicot {

/*
C     PURPOSE
C
C     To solve for X either the real continuous-time Sylvester equation
C
C        op(A)*X + ISGN*X*op(B) = scale*C,                           (1)
C
C     or the real discrete-time Sylvester equation
C
C        op(A)*X*op(B) + ISGN*X = scale*C,                           (2)
C
C     where op(M) = M or M**T, and ISGN = 1 or -1. A is M-by-M and
C     B is N-by-N; the right hand side C and the solution X are M-by-N;
C     and scale is an output scale factor, set less than or equal to 1
C     to avoid overflow in X. The solution matrix X is overwritten
C     onto C.
C
C     If A and/or B are not (upper) quasi-triangular, that is, block
C     upper triangular with 1-by-1 and 2-by-2 diagonal blocks, they are
C     reduced to Schur canonical form, that is, quasi-triangular with
C     each 2-by-2 diagonal block having its diagonal elements equal and
C     its off-diagonal elements of opposite sign.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     DICO    CHARACTER*1
C             Specifies the equation from which X is to be determined
C             as follows:
C             = 'C':  Equation (1), continuous-time case;
C             = 'D':  Equation (2), discrete-time case.
C
C     FACTA   CHARACTER*1
C             Specifies whether or not the real Schur factorization
C             of the matrix A is supplied on entry, as follows:
C             = 'F':  On entry, A and U contain the factors from the
C                     real Schur factorization of the matrix A;
C             = 'N':  The Schur factorization of A will be computed
C                     and the factors will be stored in A and U;
C             = 'S':  The matrix A is quasi-triangular (or Schur).
C
C     FACTB   CHARACTER*1
C             Specifies whether or not the real Schur factorization
C             of the matrix B is supplied on entry, as follows:
C             = 'F':  On entry, B and V contain the factors from the
C                     real Schur factorization of the matrix B;
C             = 'N':  The Schur factorization of B will be computed
C                     and the factors will be stored in B and V;
C             = 'S':  The matrix B is quasi-triangular (or Schur).
C
C     TRANA   CHARACTER*1
C             Specifies the form of op(A) to be used, as follows:
C             = 'N':  op(A) = A    (No transpose);
C             = 'T':  op(A) = A**T (Transpose);
C             = 'C':  op(A) = A**T (Conjugate transpose = Transpose).
C
C     TRANB   CHARACTER*1
C             Specifies the form of op(B) to be used, as follows:
C             = 'N':  op(B) = B    (No transpose);
C             = 'T':  op(B) = B**T (Transpose);
C             = 'C':  op(B) = B**T (Conjugate transpose = Transpose).
C
C     ISGN    INTEGER
C             Specifies the sign of the equation as described before.
C             ISGN may only be 1 or -1.
C
C     Input/Output Parameters
C
C     M       (input) INTEGER
C             The order of the matrix A, and the number of rows in the
C             matrices X and C.  M >= 0.
C
C     N       (input) INTEGER
C             The order of the matrix B, and the number of columns in
C             the matrices X and C.  N >= 0.
C
C     A       (input or input/output) DOUBLE PRECISION array,
C             dimension (LDA,M)
C             On entry, the leading M-by-M part of this array must
C             contain the matrix A. If FACTA = 'S', then A contains
C             a quasi-triangular matrix, and if FACTA = 'F', then A
C             is in Schur canonical form; the elements below the upper
C             Hessenberg part of the array A are not referenced.
C             On exit, if FACTA = 'N', and INFO = 0 or INFO >= M+1, the
C             leading M-by-M upper Hessenberg part of this array
C             contains the upper quasi-triangular matrix in Schur
C             canonical form from the Schur factorization of A. The
C             contents of array A is not modified if FACTA = 'F' or 'S'.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,M).
C
C     U       (input or output) DOUBLE PRECISION array, dimension
C             (LDU,M)
C             If FACTA = 'F', then U is an input argument and on entry
C             the leading M-by-M part of this array must contain the
C             orthogonal matrix U of the real Schur factorization of A.
C             If FACTA = 'N', then U is an output argument and on exit,
C             if INFO = 0 or INFO >= M+1, it contains the orthogonal
C             M-by-M matrix from the real Schur factorization of A.
C             If FACTA = 'S', the array U is not referenced.
C
C     LDU     INTEGER
C             The leading dimension of array U.
C             LDU >= MAX(1,M), if FACTA = 'F' or 'N';
C             LDU >= 1,        if FACTA = 'S'.
C
C     B       (input or input/output) DOUBLE PRECISION array,
C             dimension (LDB,N)
C             On entry, the leading N-by-N part of this array must
C             contain the matrix B. If FACTB = 'S', then B contains
C             a quasi-triangular matrix, and if FACTB = 'F', then B
C             is in Schur canonical form; the elements below the upper
C             Hessenberg part of the array B are not referenced.
C             On exit, if FACTB = 'N', and INFO = 0 or INFO = M+N+1,
C             the leading N-by-N upper Hessenberg part of this array
C             contains the upper quasi-triangular matrix in Schur
C             canonical form from the Schur factorization of B. The
C             contents of array B is not modified if FACTB = 'F' or 'S'.
C
C     LDB     (input) INTEGER
C             The leading dimension of the array B.  LDB >= max(1,N).
C
C     V       (input or output) DOUBLE PRECISION array, dimension
C             (LDV,N)
C             If FACTB = 'F', then V is an input argument and on entry
C             the leading N-by-N part of this array must contain the
C             orthogonal matrix V of the real Schur factorization of B.
C             If FACTB = 'N', then V is an output argument and on exit,
C             if INFO = 0 or INFO = M+N+1, it contains the orthogonal
C             N-by-N matrix from the real Schur factorization of B.
C             If FACTB = 'S', the array V is not referenced.
C
C     LDV     INTEGER
C             The leading dimension of array V.
C             LDV >= MAX(1,N), if FACTB = 'F' or 'N';
C             LDV >= 1,        if FACTB = 'S'.
C
C     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
C             On entry, the leading M-by-N part of this array must
C             contain the right hand side matrix C.
C             On exit, if INFO = 0 or INFO = M+N+1, the leading M-by-N
C             part of this array contains the solution matrix X.
C
C     LDC     INTEGER
C             The leading dimension of array C.  LDC >= MAX(1,M).
C
C     SCALE   (output) DOUBLE PRECISION
C             The scale factor, scale, set less than or equal to 1 to
C             prevent the solution overflowing.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0 or M+N+1, then: DWORK(1) returns the
C             optimal value of LDWORK; if FACTA = 'N', DWORK(1+i) and
C             DWORK(1+M+i), i = 1,...,M, contain the real and imaginary
C             parts, respectively, of the eigenvalues of A; and, if
C             FACTB = 'N', DWORK(1+f+j) and DWORK(1+f+N+j), j = 1,...,N,
C             with f = 2*M if FACTA = 'N', and f = 0, otherwise, contain
C             the real and imaginary parts, respectively, of the
C             eigenvalues of B.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             LDWORK >= MAX( 1, a+MAX( c, b+d, b+e ) ),
C             where a = 1+2*M, if FACTA =  'N',
C                   a = 0,     if FACTA <> 'N',
C                   b = 2*N,   if FACTB =  'N', FACTA =  'N',
C                   b = 1+2*N, if FACTB =  'N', FACTA <> 'N',
C                   b = 0,     if FACTB <> 'N',
C                   c = 3*M,   if FACTA =  'N',
C                   c = M,     if FACTA =  'F',
C                   c = 0,     if FACTA =  'S',
C                   d = 3*N,   if FACTB =  'N',
C                   d = N,     if FACTB =  'F',
C                   d = 0,     if FACTB =  'S',
C                   e = M,     if DICO  =  'C', FACTA <> 'S',
C                   e = 0,     if DICO  =  'C', FACTA =  'S',
C                   e = 2*M,   if DICO  =  'D'.
C             An upper bound is
C             LDWORK = 1+2*M+MAX( 3*M, 5*N, 2*N+2*M ).
C             For good performance, LDWORK should be larger, e.g.,
C             LDWORK = 1+2*M+MAX( 3*M, 5*N, 2*N+2*M, 2*N+M*N ).
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = i:  if INFO = i, i = 1,...,M, the QR algorithm failed
C                   to compute all the eigenvalues of the matrix A
C                   (see LAPACK Library routine DGEES); the elements
C                   2+i:1+M and 2+i+M:1+2*M of DWORK contain the real
C                   and imaginary parts, respectively, of the
C                   eigenvalues of A which have converged, and the
C                   array A contains the partially converged Schur form;
C             = M+j:  if INFO = M+j, j = 1,...,N, the QR algorithm
C                   failed to compute all the eigenvalues of the matrix
C                   B (see LAPACK Library routine DGEES); the elements
C                   2+f+j:1+f+N and 2+f+j+N:1+f+2*N of DWORK contain the
C                   real and imaginary parts, respectively, of the
C                   eigenvalues of B which have converged, and the
C                   array B contains the partially converged Schur form;
C                   as defined for the parameter DWORK,
C                   f = 2*M, if FACTA =  'N',
C                   f = 0,   if FACTA <> 'N';
C             = M+N+1:  if DICO = 'C', and the matrices A and -ISGN*B
C                   have common or very close eigenvalues, or
C                   if DICO = 'D', and the matrices A and -ISGN*B have
C                   almost reciprocal eigenvalues (that is, if lambda(i)
C                   and mu(j) are eigenvalues of A and -ISGN*B, then
C                   lambda(i) = 1/mu(j) for some i and j);
C                   perturbed values were used to solve the equation
C                   (but the matrices A and B are unchanged).
C
C     METHOD
C
C     An extension and refinement of the algorithms in [1,2] is used.
C     If the matrices A and/or B are not quasi-triangular (see PURPOSE),
C     they are reduced to Schur canonical form
C
C        A = U*S*U',  B = V*T*V',
C
C     where U, V are orthogonal, and S, T are block upper triangular
C     with 1-by-1 and 2-by-2 blocks on their diagonal. The right hand
C     side matrix C is updated accordingly,
C
C        C = U'*C*V;
C
C     then, the solution matrix X of the "reduced" Sylvester equation
C     (with A and B in (1) or (2) replaced by S and T, respectively),
C     is computed column-wise via a back substitution scheme. A set of
C     equivalent linear algebraic systems of equations of order at most
C     four are formed and solved using Gaussian elimination with
C     complete pivoting. Finally, the solution X of the original
C     equation is obtained from the updating formula
C
C        X = U*X*V'.
C
C     If A and/or B are already quasi-triangular (or in Schur form), the
C     initial factorizations and the corresponding updating steps are
C     omitted.
C
C     REFERENCES
C
C     [1] Bartels, R.H. and Stewart, G.W.  T
C         Solution of the matrix equation A X + XB = C.
C         Comm. A.C.M., 15, pp. 820-826, 1972.
C
C     [2] Anderson, E., Bai, Z., Bischof, C., Demmel, J., Dongarra, J.,
C         Du Croz, J., Greenbaum, A., Hammarling, S., McKenney, A.,
C         Ostrouchov, S., and Sorensen, D.
C         LAPACK Users' Guide: Second Edition.
C         SIAM, Philadelphia, 1995.
C
C     NUMERICAL ASPECTS
C
C     The algorithm is stable and reliable, since orthogonal
C     transformations and Gaussian elimination with complete pivoting
C     are used. If INFO = M+N+1, the Sylvester equation is numerically
C     singular.
C
C
C
*/


  struct sb04pd_info {
    fortran_int_t info ;
  } ;

  template <typename A, typename B, typename C, typename Z>
  sb03md_info sb04pd( char dico, char job, char facta, char factb, char trana, char tranb, int isgn, A& a, U& u, B& b, V& v, C& c, lapack::minimal_workspace ) {
    BOOST_STATIC_ASSERT( (bindings::is_column_major< A >::value) );
    BOOST_STATIC_ASSERT( (bindings::is_column_major< B >::value) );
    BOOST_STATIC_ASSERT( (bindings::is_column_major< C >::value) );
    BOOST_STATIC_ASSERT( (bindings::is_column_major< U >::value) );
    BOOST_STATIC_ASSERT( (bindings::is_column_major< V >::value) );

    sb04pd_info info ;
    int n = bindings::size_row(a) ;
    int m = bindings::size_row(b) ;
    assert( bindings::size_column(a)==n ) ;
    assert( bindings::size_column(b)==m ) ;
    assert( bindings::size_row(c)==n ) ;
    assert( bindings::size_column(c)==m ) ;
    assert( bindings::size_column(z)==m ) ;
    assert( bindings::size_row(z)==m ) ;

    fortran_int_t lda( bindings::stride_major(a) ) ;
    fortran_int_t ldb( bindings::stride_major(b) ) ;
    fortran_int_t ldc( bindings::stride_major(c) ) ;
    fortran_int_t ldz( bindings::stride_major(z) ) ;

    bindings::detail::array< fortran_int_t > iwork( 4*n ) ;
    fortran_int_t ldwork = std::max( std::max(n+m, 2*n*n+8*n), std::max(1, 5*m) ) ;
    bindings::detail::array< double > dwork( ldwork ) ;

    SLICOT_SB04MD( &n, &m, bindings::begin_value(a), &lda
                 , bindings::begin_value(b), &ldb
                 , bindings::begin_value(c), &ldc
                 , bindings::begin_value(z), &ldz
                 , bindings::begin_value(iwork)
                 , bindings::begin_value(dwork), &ldwork, &info.info ) ;
    return info ;
  } // sb04md()

} // namespace slicot
} // namespace bindings
} // namespace numeric
} // namespace boost

#endif
