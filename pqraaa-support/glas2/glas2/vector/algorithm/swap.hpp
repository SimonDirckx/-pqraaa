//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_vector_algorithm_swap_hpp
#define glas2_vector_algorithm_swap_hpp

#include <glas2/concept/is.hpp>
#include <glas2/vector/concept/dense_vector.hpp>
#include <type_traits>
#include <cassert>

namespace glas2 {

  template <typename V, typename W>
  void swap( V&& v, W&& w ) {
    assert( v.size()==w.size() ) ;
    for (int i=0; i<v.size(); ++i) {
      typename std::decay<V>::type::value_type elt = v(i) ;
      v(i) = w(i) ;
      w(i) = elt ;
    }
  }

} // namespace glas2

#endif
