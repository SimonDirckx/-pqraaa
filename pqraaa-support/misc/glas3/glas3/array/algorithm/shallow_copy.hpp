//  (C) Copyright Sam Corveleyn 2015.
//  Use, modification and distribution are subject to the
//  GLAS Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt)

#ifndef glas3_array_algorithm_shallow_copy_hpp
#define glas3_array_algorithm_shallow_copy_hpp

#include <glas3/array/concept/array.hpp>
#include <glas3/concept/is.hpp>
#include <type_traits>

namespace glas3 {

template <typename X>
typename std::enable_if< is<Array, X>::value, X >::type
shallow_copy( X const& x ) {
    return x.shallow_copy() ;
}

} // namespace glas3

#endif
