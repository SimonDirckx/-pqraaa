//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_vector_type_all_hpp
#define glas2_vector_type_all_hpp

#include <glas2/vector/concept/dense_vector.hpp>
#include <glas2/vector/concept/is_no_range_or_slice.hpp>
#include <glas2/concept/concept.hpp>

namespace glas2 {

  class all {
    public:
      typedef int value_type ;
      typedef int size_type ;

    public:
      inline value_type operator()( size_type i ) const { return i ; }
  } ;

  template <>
  struct glas_concept< all >
  : DenseVector
  {} ;

  template <>
  struct is_no_range_or_slice< all >
  : std::false_type
  {} ;

} // namespace glas2

#endif
