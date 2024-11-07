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
void SLICOT_SB04MD( const fortran_int_t* n, const fortran_int_t* m,
        double* a, const fortran_int_t* lda,
        double* b, const fortran_int_t* ldb,
        double* c, const fortran_int_t* ldc,
        double* z, const fortran_int_t* ldz,
        fortran_int_t* iwork, double* dwork, fortran_int_t* ldwork,
        fortran_int_t* info ) ;

} // extern "C"

namespace boost {
namespace numeric {
namespace bindings {
namespace slicot {

/*
C     PURPOSE
C
C     To solve for X the continuous-time Sylvester equation
C
C        AX + XB = C
C
C     where A, B, C and X are general N-by-N, M-by-M, N-by-M and
C     N-by-M matrices respectively.
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     A       (input/output) DOUBLE PRECISION array, dimension (N,N)
C             On entry, the leading N-by-N part of this array must
C             contain the coefficient matrix A of the equation.
C             On exit, the leading N-by-N upper Hessenberg part of this
C             array contains the matrix H, and the remainder of the
C             leading N-by-N part, together with the elements 2,3,...,N
C             of array DWORK, contain the orthogonal transformation
C             matrix U (stored in factored form).
C
C     B       (input/output) DOUBLE PRECISION array, dimension (M,M)
C             On entry, the leading M-by-M part of this array must
C             contain the coefficient matrix B of the equation.
C             On exit, the leading M-by-M part of this array contains
C             the quasi-triangular Schur factor S of the matrix B'.
C
C     C       (input/output) DOUBLE PRECISION array, dimension (N,M)
C             On entry, the leading N-by-M part of this array must
C             contain the coefficient matrix C of the equation.
C             On exit, the leading N-by-M part of this array contains
C             the solution matrix X of the problem.
C
C     Z       (output) DOUBLE PRECISION array, dimension (M,M)
C             The leading M-by-M part of this array contains the
C             orthogonal matrix Z used to transform B' to real upper
C             Schur form.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             > 0:  if INFO = i, 1 <= i <= M, the QR algorithm failed to
C                   compute all the eigenvalues (see LAPACK Library
C                   routine DGEES);
C             > M:  if a singular matrix was encountered whilst solving
C                   for the (INFO-M)-th column of matrix X.
C
C     METHOD
C
C     The matrix A is transformed to upper Hessenberg form H = U'AU by
C     the orthogonal transformation matrix U; matrix B' is transformed
C     to real upper Schur form S = Z'B'Z using the orthogonal
C     transformation matrix Z. The matrix C is also multiplied by the
C     transformations, F = U'CZ, and the solution matrix Y of the
C     transformed system
C
C        HY + YS' = F
C
C     is computed by back substitution. Finally, the matrix Y is then
C     multiplied by the orthogonal transformation matrices, X = UYZ', in
C     order to obtain the solution matrix X to the original problem.
C
C     REFERENCES
C
C     [1] Golub, G.H., Nash, S. and Van Loan, C.F.
C         A Hessenberg-Schur method for the problem AX + XB = C.
C         IEEE Trans. Auto. Contr., AC-24, pp. 909-913, 1979.
C
C     NUMERICAL ASPECTS
C                                         3       3      2         2
C     The algorithm requires about (5/3) N  + 10 M  + 5 N M + 2.5 M N
C     operations and is backward stable.
C
*/


  struct sb04md_info {
    fortran_int_t info ;
  } ;

  template <typename A, typename B, typename C, typename Z>
  sb04md_info sb04md( A& a, B& b, C& c, Z& z ) {
    BOOST_STATIC_ASSERT( (bindings::is_column_major< A >::value) );
    BOOST_STATIC_ASSERT( (bindings::is_column_major< B >::value) );
    BOOST_STATIC_ASSERT( (bindings::is_column_major< C >::value) );
    BOOST_STATIC_ASSERT( (bindings::is_column_major< Z >::value) );

    sb04md_info info ;
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
