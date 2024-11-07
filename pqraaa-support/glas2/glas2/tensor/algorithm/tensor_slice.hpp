//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_tensor_algorithm_tensor_slice_hpp
#define glas2_tensor_algorithm_tensor_slice_hpp

#include <glas2/tensor/concept/dense_tensor.hpp>
#include <glas2/tensor/concept/tensor_transform.hpp>
#include <glas2/concept/is.hpp>
#include <type_traits>

namespace glas2 {

  template <typename I, typename S>
  class tensor_slice_tag {
    public:
      tensor_slice_tag( I const& index, S coordinate1, S coordinate2 )
      : index_( index )
      , coordinate1_( coordinate1 )
      , coordinate2_( coordinate2 )
      {}

    public:
      I const& index() const { return index_ ; }
      S coordinate1() const { return coordinate1_ ; }
      S coordinate2() const { return coordinate2_ ; }

    private:
      I index_ ;
      S coordinate1_ ;
      S coordinate2_ ;
  } ;

  template <typename I, typename S>
  tensor_slice_tag<I,S> make_tensor_slice_tab( I const& index, S coordinate1, S coordinate2 ) {
    return tensor_slice_tag<I,S>( index, coordinate1, coordinate2 ) ;
  }

  template <typename X, typename I, typename S>
  typename std::enable_if< is<DenseTensor,X>::value, typename tensor_transform< tensor_slice_tag<I,S>, X >::result_type >::type tensor_slice( X x, I index, S coordinate1, S coordinate2 ) {
    return tensor_transform< tensor_slice_tag<I,S>, X>::apply(x,index,coordinate1, coordinate2) ;
  }

} // namespace glas2

#endif
