//  (C) Copyright Sam Corveleyn 2015.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas3_backend_default_backend_array_dense_array_assign_hpp
#define glas3_backend_default_backend_array_dense_array_assign_hpp

#include <glas3/backend/default_backend/default_backend.hpp>
#include <glas3/array/concept/array.hpp>
#include <glas3/array/dense_array/concept/dense_array.hpp>
#include <glas3/array/array_network/concept/tensor_network.hpp>
#include <glas3/concept/is.hpp>

#include <type_traits>
#include <cassert>

namespace glas3 {

template <typename To, typename From>
typename std::enable_if< is< DenseArray, To >::value && is< DenseArray, From >::value >::type assign( default_backend, To const& to, From const& from ) {
	assert( to.size() == from.size() ) ;
	for (typename To::size_type i = 0; i < to.size(); ++i) { to[i] = from[i] ; }
}

} // namespace glas3

#endif
