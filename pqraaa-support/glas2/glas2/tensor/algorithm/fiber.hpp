//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_tensor_algorithm_fiber_hpp
#define glas2_tensor_algorithm_fiber_hpp

#include <glas2/tensor/concept/dense_tensor.hpp>
#include <glas2/tensor/concept/tensor_transform.hpp>
#include <glas2/concept/is.hpp>
#include <type_traits>

namespace glas2 {

  template <typename I, typename S>
  class fiber_tag {
    public:
      fiber_tag( I const& index, S coordinate )
      : index_( index )
      , coordinate_( coordinate )
      {}

    public:
      I const& index() const { return index_ ; }
      S coordinate() const { return coordinate_ ; }

    private:
      I index_ ;
      S coordinate_ ;
  } ;

  template <typename I, typename S>
  fiber_tag<I,S> make_fiber_tag( I const& index, S coordinate ) {
    return fiber_tag<I,S>( index, coordinate ) ;
  }

  template <typename X, typename I, typename S>
  typename std::enable_if< is<DenseTensor,X>::value, typename tensor_transform< fiber_tag<I,S>, X >::result_type >::type fiber( X x, I index, S coordinate ) {
    return tensor_transform< fiber_tag<I,S>, X>::apply(x,index,coordinate) ;
  }

} // namespace glas2

#endif
