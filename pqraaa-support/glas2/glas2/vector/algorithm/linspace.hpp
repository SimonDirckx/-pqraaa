//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_vector_algorithm_linspace_hpp
#define glas2_vector_algorithm_linspace_hpp

#include <glas2/vector/type/linspace_type.hpp>
#include <glas2/concept/is.hpp>
#include <glas2/concept/is_arithmetic.hpp>
#include <type_traits>

namespace glas2 {

  template <typename T, typename S>
  typename std::enable_if< std::is_integral<S>::value && glas2::is_arithmetic<T>::value, linspace_type<T,S> >::type linspace( T const& begin, T const& end, S number ) {
    return linspace_type<T,S>( begin, end, number ) ;
  }

} // namespace glas2

#endif
