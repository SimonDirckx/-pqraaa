//  (C) Copyright Sam Corveleyn 2015.
//  Use, modification and distribution are subject to the
//  GLAS Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt)

#ifndef glas3_array_dense_array_algorithm_randomize_hpp
#define glas3_array_dense_array_algorithm_randomize_hpp

#include <glas3/backend/current_backend.hpp>
#include <glas3/backend/default_backend/array/dense_array/randomize.hpp>

#include <type_traits>
#include <utility>

namespace glas3 {

template <typename V, typename G >
typename std::enable_if< is< DenseArray, V >::value >::type randomize( V const& v, G& g ) {
	return randomize( current_backend(), v, g ) ;
}

} // namespace glas3

#endif
