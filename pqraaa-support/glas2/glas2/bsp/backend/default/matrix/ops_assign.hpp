//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_bsp_backend_default_matrix_ops_assign_hpp
#define glas2_bsp_backend_default_matrix_ops_assign_hpp

#include <glas2/bsp/backend/default/default_backend.hpp>
#include <glas2/bsp/matrix/concept/bsp_matrix.hpp>
#include <glas2/backend/current_backend.hpp>
#include <glas2/concept/is.hpp>
#include <glas2/scalar/concept/scalar.hpp>
#include <type_traits>
#include <cassert>

namespace glas2 {

  template <typename To, typename From>
  typename std::enable_if< ::glas2::is<glas2::bsp::BSPMatrix,To>::value && ::glas2::is<bsp::BSPMatrix,From>::value, To >::type plus_assign( current_backend, To to, From const& from ) {
    // We have to test whether there is no BSP needed in evaluating from
    assert( to.distribution()==from.distribution() ) ;
    to.local() += from.local() ;
    return to ;
  }

  template <typename To, typename From>
  typename std::enable_if< ::glas2::is<bsp::BSPMatrix,To>::value && ::glas2::is<bsp::BSPMatrix,From>::value, To >::type minus_assign( current_backend, To to, From const& from ) {
    // We have to test whether there is no BSP needed in evaluating from
    assert( to.distribution()==from.distribution() ) ;
    to.local() -= from.local() ;
    return to ;
  }

  template <typename To, typename From>
  typename std::enable_if< ::glas2::is<bsp::BSPMatrix,To>::value && ::glas2::is<bsp::BSPMatrix,From>::value, To >::type multiplies_assign( current_backend, To to, From const& from ) {
    // We have to test whether there is no BSP needed in evaluating from
    assert( to.distribution()==from.distribution() ) ;
    to.local() *= from.local() ;
    return to ;
  }

  template <typename To, typename From>
  typename std::enable_if< ::glas2::is<bsp::BSPMatrix,To>::value && ::glas2::is<Scalar,From>::value, To >::type multiplies_assign( current_backend, To to, From const& from ) {
    // We have to test whether there is no BSP needed in evaluating from
    to.local() *= from ;
    return to ;
  }

  template <typename To, typename From>
  typename std::enable_if< ::glas2::is<bsp::BSPMatrix,To>::value && ::glas2::is<bsp::BSPMatrix,From>::value, To >::type divides_assign( current_backend, To to, From const& from ) {
    // We have to test whether there is no BSP needed in evaluating from
    assert( to.distribution()==from.distribution() ) ;
    to.local() /= from.local() ;
    return to ;
  }

  template <typename To, typename From>
  typename std::enable_if< ::glas2::is<bsp::BSPMatrix,To>::value && ::glas2::is<Scalar,From>::value, To >::type divides_assign( current_backend, To to, From const& from ) {
    // We have to test whether there is no BSP needed in evaluating from
    to.local() /= from ;
    return to ;
  }

} // namespace glas2

#endif
