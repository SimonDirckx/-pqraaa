//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_vector_type_flip_type_hpp
#define glas2_vector_type_flip_type_hpp

#include <glas2/concept/is.hpp>
#include <glas2/vector/concept/dense_vector.hpp>
#include <type_traits>

namespace glas2 {

  template <typename V>
  class flip_type
  {
    public:
      explicit flip_type( V v )
      : vector_( v )
      {}

    public:
      typedef typename V::size_type  size_type ;
      typedef typename V::value_type value_type ;

      size_type size() const {
        return vector_.size() ;
      }

      value_type operator() ( size_type i ) const { return vector_(vector_.size()-1-i) ; }
      value_type& operator() ( size_type i) { return vector_(vector_.size()-1-i) ; }
      
    public:
      flip_type& operator=( flip_type const& that ) {
        assign( *this, that ) ;
        return *this ;
      }

      template <typename E>
      flip_type& operator=( E const& that ) {
        assign( *this, that ) ;
        return *this ;
      }

    public:
      V vector() const { return vector_; }

    private:
      V  vector_ ;
  } ;

  template <typename V>
  struct glas_concept< flip_type<V>
                , typename std::enable_if< is<DenseVector,V>::value>::type
                >
  : DenseVector
  {} ;

} // namespace glas2

#endif
