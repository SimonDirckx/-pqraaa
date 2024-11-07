//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_vector_expression_vector_vector_expression_hpp
#define glas2_vector_expression_vector_vector_expression_hpp

#include <glas2/concept/is.hpp>
#include <glas2/concept/concept.hpp>
#include <glas2/scalar/concept/scalar.hpp>
#include <glas2/vector/concept/dense_vector.hpp>
#include <type_traits>

namespace glas2 {

  template <typename V1, typename V2, typename Op>
  class vector_vector_expression {
    public:
      vector_vector_expression( V1 const& v1, V2 const& v2 )
      : v1_( v1 )
      , v2_( v2 )
      {}

      static_assert( is<DenseVector,V1>::value, "vector_vector_expression: V1 is not a Vector" );
      static_assert( is<DenseVector,V2>::value, "vector_vector_expression: V2 is not a Vector" );

    public:
      typedef typename V1::size_type                                                    size_type ;
      typedef decltype( Op() ( typename V2::value_type(), typename V2::value_type() ) ) value_type ;

      size_type size() const {
        assert( v1_.size()==v2_.size() ) ;
        return v1_.size() ;
      }
      value_type operator() ( size_type i) const { return Op()( v1_(i), v2_(i) ) ; }
      value_type operator[] ( size_type i) const { return Op()( v1_[i], v2_[i] ) ; }

    private:
      V1 v1_ ;
      V2 v2_ ;
  } ;

  template <typename V1, typename V2, typename Op>
  struct glas_concept< vector_vector_expression<V1,V2,Op>, typename std::enable_if< is<DenseVector,V1>::value && is<DenseVector,V2>::value >::type >
  : DenseVector
  {} ;

} // namespace glas2

#endif
