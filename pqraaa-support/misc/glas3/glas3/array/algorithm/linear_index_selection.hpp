//  (C) Copyright Sam Corveleyn 2015.
//  Use, modification and distribution are subject to the
//  GLAS Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt)

#ifndef glas3_array_algorithm_linear_index_selection_hpp
#define glas3_array_algorithm_linear_index_selection_hpp

#include <glas3/concept/is.hpp>
#include <glas3/concept/concept.hpp>
#include <glas3/array/concept/array.hpp>

#include <type_traits>

namespace glas3 {

template <typename T, typename S>
typename std::enable_if< !is< Array, T >::value, void >::type
linear_index_selection ( T const& t, S const& s)
{} ;

} // namespace glas3

#endif
