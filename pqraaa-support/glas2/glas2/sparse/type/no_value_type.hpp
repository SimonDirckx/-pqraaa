//  (C) Copyright Karl Meerbergen 2015.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_sparse_type_no_value_type
#define glas2_sparse_type_no_value_type

#include <glas2/scalar/concept/scalar.hpp>

namespace glas2 {

  struct no_value_type {
  } ;

  template <typename T>
  T operator*( no_value_type, T const& t ) { return t ; }

/*  template <typename T>
  T operator*( T const& t, no_value_type ) { return t ; }*/

  template <>
  struct glas_concept< no_value_type >
  : Scalar
  {};

} // namespace glas2


#endif
