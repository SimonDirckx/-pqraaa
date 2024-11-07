//  (C) Copyright Karl Meerbergen 2015.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_vector_type_range_from_end_hpp
#define glas2_vector_type_range_from_end_hpp

#include <glas2/vector/concept/vector.hpp>
#include <glas2/vector/concept/is_no_range_or_slice.hpp>
#include <glas2/concept/concept.hpp>

namespace glas2 {

  class range_from_end {
    public:
      typedef int value_type ;
      typedef int size_type ;

    public:
      inline range_from_end( value_type begin, value_type end )
      : begin_( begin )
      , end_( end )
      {
        assert( begin_>=0 ) ;
        assert( end_>=0 ) ;
      }

    public:
      inline size_type begin() const { return begin_ ; }
      inline size_type from_end() const { return end_ ; }

      inline value_type operator()( size_type i ) const { return begin_ + i ; }

    private:
      value_type begin_ ;
      value_type end_ ;
  } ;

  template <>
  struct glas_concept< range_from_end >
  : DenseVector
  {} ;

  template <>
  struct is_no_range_or_slice< range_from_end >
  : std::false_type
  {} ;

} // namespace glas2

#endif
