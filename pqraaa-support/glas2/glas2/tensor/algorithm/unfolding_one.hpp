//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_tensor_algorithm_unfolding_one_hpp
#define glas2_tensor_algorithm_unfolding_one_hpp

#include <glas2/tensor/concept/dense_tensor.hpp>
#include <glas2/tensor/concept/tensor_transform.hpp>
#include <glas2/concept/is.hpp>
#include <type_traits>

namespace glas2 {

  // Unfolding following 1 coordinate in vertical direction and the other coordinates in horizontal direction
  template <typename S>
  class unfolding_one_tag {
    public:
      unfolding_one_tag( S coordinate )
      : coordinate_( coordinate )
      {}

    public:
      S coordinate() const { return coordinate_ ; }

    private:
      S coordinate_ ;
  } ;

  template <typename X, typename S>
  typename std::enable_if< is<DenseTensor,X>::value, typename tensor_transform< unfolding_one_tag<S>, X >::result_type >::type unfolding_one( X x, S coordinate ) {
    return tensor_transform< unfolding_one_tag<S>, X>::apply(x,coordinate) ;
  }

} // namespace glas2

#endif
