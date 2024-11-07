//  (C) Copyright Karl Meerbergen 2020.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_approximation_unique_test_set_hpp
#define cork_approximation_unique_test_set_hpp

#include <glas2/vector.hpp>
#include <algorithm>

namespace CORK { namespace approximation {

  template <typename V>
  decltype(auto) unique_test_set( V& vec ) {
    auto it = std::unique( vec.begin(), vec.end() ) ;
    int new_size = it - vec.begin() ;
    glas2::shared_vector< typename V::value_type > v_new( new_size ) ;
    v_new = vec( glas2::range(0,new_size) ) ;
    return v_new ;
  } // unique_test_set()

  template <typename V, typename Predicate>
  decltype(auto) unique_test_set( V& vec, Predicate const& pred ) {
    auto it = std::unique( vec.begin(), vec.end(), pred ) ;
    int new_size = it - vec.begin() ;
    glas2::shared_vector< typename V::value_type > v_new( new_size ) ;
    v_new = vec( glas2::range(0,new_size) ) ;
    return v_new ;
  } // unique_test_set()

} } // namespace CORK::approximation

#endif
