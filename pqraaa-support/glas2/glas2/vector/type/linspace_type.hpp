//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_vector_type_linspace_type_hpp
#define glas2_vector_type_linspace_type_hpp

#include <glas2/vector/concept/dense_vector.hpp>
#include <glas2/concept/concept.hpp>
#include <cmath>

namespace glas2 {

  template <typename T, typename S=int>
  class linspace_type {
    public:
      typedef T                                value_type ;
      typedef decltype(std::abs(value_type())) real_type ;
      typedef S                                size_type ;

    public:
      inline linspace_type( value_type begin, value_type end, size_type number )
      : begin_( begin )
      , end_( end )
      , number_( number )
      , step_( (end_-begin_)/real_type(std::max<size_type>(number_,2)-1) )
      {
        assert( number_>0 ) ;
      }

    public:
      inline size_type size() const { return number_ ; }
      inline size_type begin() const { return begin_ ; }

      inline value_type operator()( size_type i ) const {
        assert( i>=0 && i<size() ) ;
        return begin_ + real_type(i)*step_ ;
      }

    private:
      value_type begin_ ;
      value_type end_ ;
      size_type  number_ ;
      value_type step_ ;
  } ;

  template <typename T, typename S>
  struct glas_concept< linspace_type<T,S> >
  : DenseVector
  {} ;

} // namespace glas2

#endif
