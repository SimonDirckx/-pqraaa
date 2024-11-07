//  (C) Copyright Sam Corveleyn 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas3_array_dense_array_algorithm_assign_hpp
#define glas3_array_dense_array_algorithm_assign_hpp

#include <glas3/backend/default_backend/array/dense_array/assign.hpp>
#include <glas3/backend/current_backend.hpp>
#include <glas3/array/dense_array/concept/dense_array.hpp>
#include <glas3/concept/is.hpp>
#include <type_traits>

namespace glas3 {

template <typename To, typename From>
typename std::enable_if< is< DenseArray, To >::value && is< DenseArray, From >::value >::type assign( To const& to, From const& from ) {
	assign( current_backend(), to, from ) ;
}

} // namespace glas3

#endif
