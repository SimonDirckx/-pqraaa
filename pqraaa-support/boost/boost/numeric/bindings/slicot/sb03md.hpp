//
// Copyright (c) 2011-
// Toon Knapen, Karl Meerbergen, Kresimir Fresl,
// Thomas Klimpel and Rutger ter Borg
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_NUMERIC_BINDINGS_SLICOT_SB03MD_HPP
#define BOOST_NUMERIC_BINDINGS_SLICOT_SB03MD_HPP

#include <boost/numeric/bindings/detail/config/fortran.hpp>
#include <boost/numeric/bindings/detail/array.hpp>
#include <boost/numeric/bindings/lapack/workspace.hpp>

//
// SLICOT computational routines
//

// Value-type variants of hbgst
#define SLICOT_SB03MD FORTRAN_ID( sb03md )


extern "C" {

//
// LAPACK computational routines
//

// Value-type variants of hbgst
void SLICOT_SB03MD( const char* dico, const char* job, const char* fact, const char* trana,
        const fortran_int_t* n,
        double* a, const fortran_int_t* lda,
        double* u, const fortran_int_t* ldu,
        double* c, const fortran_int_t* ldc, double* scale, double* sep, double* ferr,
        double* wr, double* wi,
        fortran_int_t* iwork, double* dwork, const fortran_int_t* ldwork,
        fortran_int_t* info ) ;

} // extern "C"

namespace boost {
namespace numeric {
namespace bindings {
namespace slicot {

/*
  To solve for X either the real continuous-time Lyapunov equation

     op(A)'*X + X*op(A) = scale*C                             (1)

  or the real discrete-time Lyapunov equation

     op(A)'*X*op(A) - X = scale*C                             (2)

  and/or estimate an associated condition number, called separation,
  where op(A) = A or A' (A**T) and C is symmetric (C = C').
  (A' denotes the transpose of the matrix A.) A is N-by-N, the right
  hand side C and the solution X are N-by-N, and scale is an output
  scale factor, set less than or equal to 1 to avoid overflow in X.

  DICO    CHARACTER*1
          Specifies the equation from which X is to be determined
          as follows:
          = 'C':  Equation (1), continuous-time case;
          = 'D':  Equation (2), discrete-time case.

  JOB     CHARACTER*1
          Specifies the computation to be performed, as follows:
          = 'X':  Compute the solution only;
          = 'S':  Compute the separation only;
          = 'B':  Compute both the solution and the separation.

  FACT    CHARACTER*1
          Specifies whether or not the real Schur factorization
          of the matrix A is supplied on entry, as follows:
          = 'F':  On entry, A and U contain the factors from the
                  real Schur factorization of the matrix A;
          = 'N':  The Schur factorization of A will be computed
                  and the factors will be stored in A and U.

  TRANA   CHARACTER*1
          Specifies the form of op(A) to be used, as follows:
          = 'N':  op(A) = A    (No transpose);
          = 'T':  op(A) = A**T (Transpose);
          = 'C':  op(A) = A**T (Conjugate transpose = Transpose).

*/


  struct sb03md_info {
    fortran_int_t info ;
    double scale ;
    double sep ;
    double ferr ;
  } ;

  template <typename A, typename U, typename C, typename WR, typename WI>
  sb03md_info sb03md( char dico, char job, char fact, char trana, A& a, U& u, C& c, WR& wr, WI& wi, lapack::minimal_workspace ) {
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
    SLICOT_SB03MD( &dico, &job, &fact, &trana, &n, bindings::begin_value(a), &lda
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
