//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2015.
//
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_lapack_qr_hpp
#define cork_lapack_qr_hpp

#include <glas2/vector.hpp>
#include <glas2/matrix.hpp>
#include <glas2/bindings/vector.hpp>
#include <glas2/bindings/matrix.hpp>
#include <boost/numeric/bindings/lapack/computational/geqrf.hpp>
#include <boost/numeric/bindings/lapack/computational/unmqr.hpp>
#include <boost/numeric/bindings/lapack/computational/ormqr.hpp>
#include <boost/numeric/bindings/tag.hpp>
#include <type_traits>
#include <cassert>

namespace CORK { namespace lapack {

  namespace detail {

    template <typename T, typename EnableIf=void>
    struct extract_q_from_geqrf {} ;

    // For real types
    template <typename T>
    struct extract_q_from_geqrf< T, typename std::enable_if< std::is_floating_point<T>::value >::type > {
      template <typename A, typename Tau, typename M>
      int operator() ( A const& a, Tau const& tau, M& m ) const {
        m = glas2::identity_matrix<T>( m.num_rows(), m.num_columns() ) ;
        int info = boost::numeric::bindings::lapack::ormqr( boost::numeric::bindings::tag::right(), a,tau, m ) ;
        return info;
      }
    } ;

    template <typename T>
    struct extract_q_from_geqrf< std::complex<T>, typename std::enable_if< std::is_floating_point<T>::value >::type > {
      template <typename A, typename Tau, typename M>
      int operator() ( A const& a, Tau const& tau, M& m ) const {
        m = glas2::identity_matrix<T>( m.num_rows(), m.num_columns() ) ;
        return boost::numeric::bindings::lapack::unmqr( boost::numeric::bindings::tag::right(), a,tau, m ) ;
      }
    } ;
  } // namespace detail

  //
  // Input: a
  //
  // Output:
  //  a: contains R
  //  q: contains Q
  template <typename A, typename Q>
  int qr( A& a, Q& q ) {
    typedef typename A::value_type value_type ;
    assert( q.num_rows()==a.num_rows() || q.num_rows()==std::max(a.num_rows(),a.num_columns()) ) ;
    assert( q.num_columns()==a.num_columns() || q.num_columns()==std::max(a.num_rows(),a.num_columns()) ) ;

    glas2::vector< value_type > tau( std::min(a.num_rows(), a.num_columns() ) ) ;

    int info = boost::numeric::bindings::lapack::geqrf( a, tau ) ;
    if (info!=0) return info ;

    info = detail::extract_q_from_geqrf< typename A::value_type >() ( a, tau, q ) ;
    if (info) return info ;

    typedef typename A::size_type size_type ;
    size_type n = std::min(a.num_rows(), a.num_columns() ) ;
    for (size_type i=0; i<n; ++i) {
      fill( a( glas2::range_from_end(i+1,0), i ), 0.0 ) ;
    }

    return 0 ;
  } // qr()

} } // namespace CORK::lapack

#endif
