//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_vector_type_numbers_vector_hpp
#define glas2_vector_type_numbers_vector_hpp

#include <glas2/concept/is.hpp>
#include <glas2/concept/concept.hpp>
#include <glas2/scalar/concept/scalar.hpp>
#include <glas2/vector/concept/dense_vector.hpp>
#include <type_traits>

namespace glas2 {

  template <typename V>
  class numbers_vector {
    public:
      numbers_vector( V const& v )
      : list_( v )
      {}

    public:
      typedef typename V::size_type  size_type ;
      typedef typename V::value_type value_type ;

      size_type size() const {
        return list_.size() ;
      }
      value_type operator() ( size_type i ) const {
        assert( i>=0 && i<size() ) ;
        return list_.begin()[i] ;
      }

    private:
      V const& list_ ;
  } ;

  template <typename V>
  struct glas_concept< numbers_vector<V> >
  : DenseVector
  {} ;

} // namespace glas2

#endif
