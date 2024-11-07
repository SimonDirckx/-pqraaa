//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2015.
//
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_lapack_schur_hpp
#define cork_lapack_schur_hpp

#include <glas2/vector.hpp>
#include <glas2/matrix.hpp>
#include <glas2/bindings/vector.hpp>
#include <glas2/bindings/matrix.hpp>
#include <boost/numeric/bindings/lapack/driver/gees.hpp>
#include <boost/numeric/bindings/lapack/driver/gges.hpp>
#include <type_traits>
#include <cassert>

namespace CORK { namespace lapack {

  namespace detail {

    template <typename T, typename EnableIf=void>
    struct schur {} ;

    // For real types
    template <typename T>
    struct schur< T, typename std::enable_if< std::is_floating_point<T>::value >::type > {
      template <typename A, typename B, typename X, typename Y, typename E>
      int operator() ( A& a, B& b, X& x, Y& y, E& e ) const {
        int info ;
        int sdim ;
        int n = a.num_rows() ;
        glas2::vector<T> alpha_r( n ) ;
        glas2::vector<T> alpha_i( n ) ;
        glas2::vector<T> beta( n ) ;

        info = boost::numeric::bindings::lapack::gges( 'V', 'V', 'N', 0, a, b, sdim, alpha_r, alpha_i, beta, y, x ) ;
        if (info) return info ;

        glas2::real(e) = alpha_r / beta ;
        glas2::imag(e) = alpha_i / beta ;

        return info ;
      }

      template <typename A, typename X, typename E>
      int operator() ( A& a, X& x, E& e ) const {
        int info ;
        int sdim ;
        glas2::vector<T> wr( e.size() ) ;
        glas2::vector<T> wi( e.size() ) ;
        assert( e.size()==a.num_rows() ) ;
        assert( e.size()==a.num_columns() ) ;
        assert( e.size()==x.num_rows() ) ;
        assert( e.size()==x.num_columns() ) ;

        info = boost::numeric::bindings::lapack::gees( 'V', 'N', 0, a, sdim, wr, wi, x ) ;
        if (info) return info ;

        glas2::real(e) = wr ; glas2::imag(e) = wi ;

        return info ;
      }
    } ;

    template <>
    struct schur< std::complex<double> > {
      template <typename A, typename B, typename X, typename Y, typename E>
      int operator() ( A& a, B& b, X& x, Y& y, E& e ) const {
        int info ;
        int sdim ;
        int n = a.num_rows() ;
        glas2::vector< std::complex<double> > alpha( n ) ;
        glas2::vector< std::complex<double> > beta( n ) ;
        info = boost::numeric::bindings::lapack::gges( 'V', 'V', 'N', 0, a, b, sdim, alpha, beta, y, x ) ;

        e = alpha / beta ;

        return info ;
      }

      template <typename A, typename X, typename E>
      int operator() ( A& a, X& x, E& e ) const {
        int info ;
        int sdim ;
        info = boost::numeric::bindings::lapack::gees( 'V', 'N', 0, a, sdim, e, x ) ;

        return info ;
      }
    } ;
  } // namespace detail

  //
  // Generic function schur (as in Matlab)
  // Compute the eigendecomposition
  //  A X = X S
  // where S is upper triangular for complex A
  // and quasi upper triangular for real A
  //
  // A is overwritten by S.
  //
  // Concept requirements:
  // StridedDenseMatrix(A)
  // StridedDenseMatrix(X)
  // StridedDenseVector(E) (E must be complex valued)
  //
  template <typename A, typename X, typename E>
  int schur( A& a, X& x, E& e ) {
    return detail::schur< typename A::value_type >() ( a, x, e ) ;
  } // schur()



  //
  // Generic function eig (as in Matlab)
  // Compute the eigendecomposition
  //  A * X = Y * S
  //  B * X = Y * T
  // where S and T are upper triangular and complex if A and B are complex
  // and real quasi upper triangular when A and B are real.
  //
  // A is overwritten by S and B is overwritten by T.
  //
  // Concept requirements:
  // StridedDenseMatrix(A)
  // StridedDenseMatrix(B)
  // StridedDenseMatrix(X)
  // StridedDenseMatrix(Y)
  // StridedDenseVector(E) (E must be complex valued)
  //
  template <typename A, typename B, typename X, typename Y, typename E>
  int schur( A& a, B& b, X& x, Y& y, E& e ) {
    static_assert( std::is_same< typename A::value_type, typename B::value_type >::value, "CORK::lapack::schur: A and B must have the same value_type" ) ;
    return detail::schur< typename A::value_type >() ( a, b, x, y, e ) ;
  } // schur()

} } // namespace CORK::lapack

#endif
