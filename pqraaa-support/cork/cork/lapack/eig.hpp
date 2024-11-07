//  (C) Copyright Karl Meerbergen 2007. 
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_lapack_eig_hpp
#define cork_lapack_eig_hpp

#include <glas2/vector.hpp>
#include <glas2/matrix.hpp>
#include <glas2/bindings/vector.hpp>
#include <glas2/bindings/matrix.hpp>
#include <boost/numeric/bindings/lapack/driver/geev.hpp>
#include <boost/numeric/bindings/lapack/driver/ggev.hpp>
#include <cassert>

#undef GLAS_INCLUDE_LEVEL
#define GLAS_INCLUDE_LEVEL GLAS_INCLUDE_TOOLBOX

namespace CORK { namespace lapack {

  namespace detail {

    // For real types.
    template <typename A, typename X, typename E>
    typename std::enable_if< std::is_floating_point<typename A::value_type>::value
                           && !std::is_floating_point<typename X::value_type>::value
                           , int >::type eig( A& a, X& x, E& e ) {
      typedef typename A::value_type value_type ;
      glas2::shared_vector< value_type > beta( e.size() ) ;
      glas2::shared_vector< value_type > wr( e.size() ) ;
      glas2::shared_vector< value_type > wi( e.size() ) ;
      glas2::shared_matrix< value_type > xr( e.size(), e.size() ) ;

      int info = ::boost::numeric::bindings::lapack::geev( 'N', 'V', a, wr, wi, xr, xr ) ;
      if (info) return info ;
      real(e) = wr ;
      imag(e) = wi ;
      for ( int i=0; i<e.size(); ++i ) {
        glas2::real( x( glas2::all(), i ) ) = xr( glas2::all(), i) ;
        if ( wi(i)!=0.0 ) {
          glas2::real( x( glas2::all(), i+1 ) ) = xr( glas2::all(), i) ;
          glas2::imag( x( glas2::all(), i ) ) = xr( glas2::all(), i+1) ;
          glas2::imag( x( glas2::all(), i+1 ) ) = -xr( glas2::all(), i+1) ;
          ++i ;
        } else {
          fill( glas2::imag( x( glas2::all(), i ) ),  0.0 ) ;
        }
      }
      return info ;
    }

    template <typename A, typename X, typename E>
    typename std::enable_if< !std::is_floating_point<typename A::value_type>::value, int >::type eig( A& a, X& x, E& e ) {
      return ::boost::numeric::bindings::lapack::geev( 'N', 'V', a, e, x, x ) ;
    }


    template <typename A, typename B, typename X, typename E>
    typename std::enable_if< std::is_floating_point< typename X::value_type >::value
                           && std::is_floating_point< typename A::value_type >::value
                           , int >::type eig( A& a, B& b, X& x, E& e )
    {
      typedef typename A::value_type value_type ;
      glas2::shared_vector< value_type > beta( e.size() ) ;
      glas2::shared_vector< value_type > alphar( e.size() ) ;
      glas2::shared_vector< value_type > alphai( e.size() ) ;

      int info = ::boost::numeric::bindings::lapack::ggev( 'N', 'V', a, b, alphar, alphai, beta, x, x ) ;
      if (info) return info ;
      real(e) = alphar / beta ;
      imag(e) = alphai / beta ;
      return info ;
    }

    template <typename A, typename B, typename X, typename E>
    typename std::enable_if< !std::is_floating_point< typename X::value_type >::value
                           && std::is_floating_point<typename A::value_type >::value
                           && std::is_floating_point<typename B::value_type >::value
                           , int >::type eig( A& a, B& b, X& x, E& e ) {
      typedef typename A::value_type value_type ;
      glas2::shared_vector< value_type > beta( e.size() ) ;
      glas2::shared_vector< value_type > alphar( e.size() ) ;
      glas2::shared_vector< value_type > alphai( e.size() ) ;
      glas2::shared_matrix< value_type > xr( e.size(), e.size() ) ;

      int info = ::boost::numeric::bindings::lapack::ggev( 'N', 'V', a, b, alphar, alphai, beta, xr, xr ) ;
      if (info) return info ;
      real(e) = alphar / beta ;
      imag(e) = alphai / beta ;
      for ( int i=0; i<e.size(); ++i ) {
        real( x( glas2::all(), i ) ) = xr( glas2::all(), i ) ;
        if ( alphai(i)!=0.0 ) {
          real( x( glas2::all(), i+1 ) ) = xr( glas2::all(), i ) ;
          imag( x( glas2::all(), i ) ) = xr( glas2::all(), i+1 ) ;
          imag( x( glas2::all(), i+1 ) ) = -xr( glas2::all(), i+1 ) ;
          e(i+1) = conj(e(i));
          ++i ;
        } else {
          fill( imag( x( glas2::all(), i ) ), 0.0 ) ;
        }
      }
      return info ;
    }

    template <typename A, typename B, typename X, typename E>
    typename std::enable_if< !std::is_floating_point<typename A::value_type >::value, int >::type eig( A& a, B& b, X& x, E& e ) {
      glas2::vector< typename A::value_type > beta( e.size() ) ;

      int info = ::boost::numeric::bindings::lapack::ggev( 'N', 'V', a, b, e, beta, x, x ) ;
      if (info) return info ;
      std::cout << e << std::endl ;
      for (typename E::size_type i=0; i<e.size(); ++i) {
  //      if (e(i).imag()!=0) {
  //        assert( e(i)==conj(e(i+1)) ) ;
  //        e(i) /= beta(i) ;
  //        e(i+1) = conj(e(i));
  //      } else {
          e(i) /= beta(i) ;
  //      }
      }
      return info ;
    }

    template <typename A, typename B, typename X, typename Y, typename E>
    typename std::enable_if< !std::is_floating_point<typename A::value_type >::value, int >::type eig( A& a, B& b, X& x, Y& y, E& e ) {
      glas2::vector< typename A::value_type > beta( e.size() ) ;

      int info = ::boost::numeric::bindings::lapack::ggev( 'V', 'V', a, b, e, beta, y, x ) ;
      if (info) return info ;
      for (typename E::size_type i=0; i<e.size(); ++i) {
        e(i) /= beta(i) ;
      }
      return info ;
    }
  } // namespace detail

  //
  // Generic function eig (as in Matlab)
  // Compute the eigendecomposition
  //  A X = X E
  // where E is diagonal and complex (even when A and B are real)
  //
  // Concept requirements:
  // StridedDenseMatrix(A)
  // StridedDenseMatrix(X)
  // StridedDenseVector(E) (E and X must be complex valued)
  //
  // This function dispatches depending on the structures of A and B
  //
  template <typename A, typename X, typename E>
  int eig( A& a, X& x, E& e ) {
    assert( a.num_rows()==a.num_columns() ) ;

    assert( a.num_rows()==x.num_rows() ) ;
    assert( a.num_rows()==x.num_columns() ) ;

    assert( a.num_rows()==e.size() ) ;

    return detail::eig( a, x, e ) ;
  } // eig()


  //
  // Generic function eig (as in Matlab)
  // Compute the eigendecomposition
  //  A X = B X E
  // where E is diagonal and complex (even when A and B are real)
  //
  // Concept requirements:
  // StridedDenseMatrix(A)
  // StridedDenseMatrix(B)
  // StridedDenseMatrix(X)
  // StridedDenseMatrix(Y)
  // StridedDenseVector(E) (E, X and Y must be complex valued)
  //
  // This function dispatches depending on the structures of A and B
  //
  template <typename A, typename B, typename X, typename E>
  int eig( A& a, B& b, X& x, E& e ) {
    assert( a.num_rows()==a.num_columns() ) ;

    assert( a.num_rows()==b.num_rows() ) ;
    assert( a.num_rows()==b.num_columns() ) ;

    assert( a.num_rows()==x.num_rows() ) ;
    assert( a.num_rows()==x.num_columns() ) ;

    assert( a.num_rows()==e.size() ) ;

    return detail::eig( a, b, x, e ) ;
  } // eig()

  template <typename A, typename B, typename X, typename Y, typename E>
  int eig( A& a, B& b, X& x, Y& y, E& e ) {
    assert( a.num_rows()==a.num_columns() ) ;

    assert( a.num_rows()==b.num_rows() ) ;
    assert( a.num_rows()==b.num_columns() ) ;

    assert( a.num_rows()==x.num_rows() ) ;
    assert( a.num_rows()==x.num_columns() ) ;

    assert( a.num_rows()==y.num_rows() ) ;
    assert( a.num_rows()==y.num_columns() ) ;

    assert( a.num_rows()==e.size() ) ;

    return detail::eig( a, b, x, y, e ) ;
  } // eig()

} } // namespace CORK::lapack

#endif
