//  (C) Copyright Karl Meerbergen (2011).
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas_toolbox_tensor_concept_storage_size_hpp
#define glas_toolbox_tensor_concept_storage_size_hpp

#include <glas/toolbox/tensor/concept/shape.hpp>
#include <glas/algorithm/product.hpp>
#include <glas/concept/value_type.hpp>

namespace glas { 

template <typename X>
struct storage_size_type
{
  typedef typename shape_type<X>::type          shape_type ;
  typedef typename value_type<shape_type>::type type ;
};

template <typename X>
typename storage_size_type<X>::type storage_size( X const& x) {
  return product( shape(x) ) ;
} // storage_size()

} // namespace glas

#endif
