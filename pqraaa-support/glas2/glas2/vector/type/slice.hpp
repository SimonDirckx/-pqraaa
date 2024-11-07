//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_vector_type_slice_hpp
#define glas2_vector_type_slice_hpp

#include <glas2/vector/concept/dense_vector.hpp>
#include <glas2/vector/concept/is_no_range_or_slice.hpp>
#include <glas2/concept/concept.hpp>
#include <cassert>

namespace glas2 {

  class slice {
    public:
      typedef int value_type ;
      typedef int size_type ;

    public:
      inline slice( value_type begin, value_type end, value_type step )
      : begin_( begin )
      , end_( end )
      , step_( step )
      {}

    public:
      inline size_type size() const { return ( end_-begin_+step_-1 ) / step_ ; }
      inline size_type begin() const { return begin_ ; }
      inline size_type step() const { return step_ ; }

      inline value_type operator()( size_type i ) const {
        assert( begin_ + i*step_>=0 && begin_ + i*step_<end_ ) ;
        return begin_ + i*step_ ;
      }

    private:
      value_type begin_ ;
      value_type end_ ;
      value_type step_ ;
  } ;

  template <>
  struct glas_concept< slice >
  : DenseVector
  {} ;

  template <>
  struct is_no_range_or_slice< slice >
  : std::false_type
  {} ;

} // namespace glas2

#endif
