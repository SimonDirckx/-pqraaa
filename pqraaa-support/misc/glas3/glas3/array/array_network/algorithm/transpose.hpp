//  (C) Copyright Sam Corveleyn 2015.
//  Use, modification and distribution are subject to the
//  GLAS Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt)

#ifndef glas3_array_array_network_algorithm_transpose_hpp
#define glas3_array_array_network_algorithm_transpose_hpp

#include <glas3/concept/is.hpp>
#include <glas3/array/concept/array.hpp>
#include <glas3/array/array_network/concept/tensor_network.hpp>

#include <glas3/array/dense_array/type/range.hpp>

#include <glas3/array/array_network/algorithm/permute.hpp>

#include <type_traits>
#include <initializer_list>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>

namespace glas3 {

template < typename T >
auto transpose( T const& t )
-> typename std::enable_if< is< TensorNetwork, T >::value, decltype( permute( t, glas3::range<typename T::ndims_type>( t.shape().size() - 1, -1, -1 ) ) ) >::type {
	return permute( t, glas3::range<typename T::ndims_type>( t.shape().size() - 1, -1, -1 ) ) ;
}

} // namespace glas3

#endif
