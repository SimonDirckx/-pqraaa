//  (C) Copyright Sam Corveleyn 2015.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas3_array_array_network_concept_Tucker_tensor_hpp
#define glas3_array_array_network_concept_Tucker_tensor_hpp

#include <glas3/array/array_network/concept/tensor_network.hpp>

namespace glas3 {

  struct TuckerTensor
  : TensorNetwork
  {
    typedef TuckerTensor type ;
  } ;

} // namespace glas3

#endif
