//  (C) Copyright Sam Corveleyn 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas3_array_dense_array_algorithm_fill_hpp
#define glas3_array_dense_array_algorithm_fill_hpp

#include <glas3/backend/current_backend.hpp>
#include <glas3/backend/default_backend/array/dense_array/fill.hpp>
#include <glas3/array/dense_array/concept/dense_array.hpp>
#include <glas3/concept/is.hpp>
#include <type_traits>

namespace glas3 {

template <typename V, typename T>
typename std::enable_if< is< DenseArray, V >::value >::type fill( V const& v, T const& value ) {
	fill( current_backend(), v, value ) ;
}

} // namespace glas3

#endif
