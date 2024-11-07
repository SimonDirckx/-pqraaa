//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_vector_container_range_hpp
#define glas2_vector_container_range_hpp

#include <glas2/vector/concept/dense_vector.hpp>
#include <glas2/concept/concept.hpp>

namespace glas2 {

  class range {
    public:
      typedef int value_type ;
      typedef int size_type ;

    public:
      range( value_type begin, value_type end )
      : begin_( begin )
      , end_( end )
      {}

    public:
      size_type size() const { return end_-begin_ ; }

      value_type operator[]( size_type i ) const { return begin_ + i ; }
      value_type operator()( size_type i ) const { return begin_ + i ; }

    private:
      value_type begin_ ;
      value_type end_ ;
  } ;

  template <>
  struct glas_concept< range >
  : DenseVector
  {} ;

} // namespace glas2

#endif
