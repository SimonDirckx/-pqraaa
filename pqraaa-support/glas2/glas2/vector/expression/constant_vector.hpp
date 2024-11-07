//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_vector_expression_constant_vector_hpp
#define glas2_vector_expression_constant_vector_hpp

#include <glas2/concept/is.hpp>
#include <glas2/concept/concept.hpp>
#include <glas2/scalar/concept/scalar.hpp>
#include <glas2/vector/concept/dense_vector.hpp>
#include <type_traits>

namespace glas2 {

  template <typename I, typename T>
  class constant_vector {
    public:
      constant_vector( I const& n, T const& v )
      : size_( n )
      , value_( v )
      {}

    public:
      typedef I size_type ;
      typedef T value_type ;

      size_type size() const {
        return size_ ;
      }
      value_type operator() ( size_type ) const { return value_ ; }
      value_type operator[] ( size_type ) const { return value_ ; }

    private:
      size_type  size_ ;
      value_type value_ ;
  } ;

  template <typename I, typename T>
  struct glas_concept< constant_vector<I,T> >
  : DenseVector
  {} ;

} // namespace glas2

#endif
