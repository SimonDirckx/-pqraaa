//  (C) Copyright Sam Corveleyn 2015.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas3_array_array_network_algorithm_assign_hpp
#define glas3_array_array_network_algorithm_assign_hpp

#include <glas3/concept/is.hpp>
#include <glas3/array/concept/array.hpp>
#include <glas3/array/dense_array/concept/dense_array.hpp>
#include <glas3/array/array_network/concept/tensor_network.hpp>

#include <glas3/array/array_network/algorithm/full.hpp>
#include <glas3/array/dense_array/algorithm/assign.hpp>

#include <type_traits>

namespace glas3 {

template <typename To, typename From>
typename std::enable_if< is< DenseArray, To >::value && is< TensorNetwork, From >::value >::type assign( To const& to, From const& from ) {
	auto from_full = full( from ) ;
	assign( to, from_full ) ;
}

} // namespace glas3

#endif
