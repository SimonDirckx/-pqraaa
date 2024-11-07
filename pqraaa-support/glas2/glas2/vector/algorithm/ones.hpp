#ifndef glas2_vector_algorithm_ones_hpp
#define glas2_vector_algorithm_ones_hpp

//#include <glas2/vector/type/static_const_vector.hpp>
#include <glas2/vector/type/constant_vector.hpp>

namespace glas2 {

  template <typename T, typename I>
  decltype (auto) ones( I n ) {
    //return static_const_vector<I,T,1>( n ) ;
    return constant_vector<I,T>( n, 1 ) ;
  }

} // namespace glas2

#endif
