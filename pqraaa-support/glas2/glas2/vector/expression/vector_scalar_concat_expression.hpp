//  (C) Copyright Karl Meerbergen 2016.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_vector_scalar_expression_vector_concat_expression_hpp
#define glas2_vector_scalar_expression_vector_concat_expression_hpp

#include <glas2/concept/is.hpp>
#include <glas2/concept/concept.hpp>
#include <glas2/type/pass_reference.hpp>
#include <glas2/vector/concept/dense_vector.hpp>
#include <glas2/scalar/concept/scalar.hpp>
#include <type_traits>

namespace glas2 {

  template <typename V1, typename V2>
  class vector_scalar_concat_expression {
    public:
      vector_scalar_concat_expression( V1 const& v1, V2 const& v2 )
      : v1_( v1 )
      , v2_( v2 )
      {}

      static_assert( is<DenseVector,V1>::value, "vector_concat_expression: V1 is not a Vector" );
      static_assert( is<Scalar,V2>::value, "vector_concat_expression: V2 is not a Scalar" );

    public:
      typedef typename V1::size_type                        size_type ;
      typedef decltype( typename V1::value_type() + V2() )  value_type ;

      size_type size() const {
        return v1_.size()+1 ;
      }
      value_type operator() ( size_type i ) const {
        assert( i>=0 && i < size() ) ;
        if (i<v1_.size()) return v1_(i) ;
        else return v2_ ;
      }

    private:
      typedef typename pass_reference<V1>::type v1_ref ;
      v1_ref v1_ ;
      V2     v2_ ;
  } ;

  template <typename V1, typename V2>
  struct glas_concept< vector_scalar_concat_expression<V1,V2>, typename std::enable_if< is<DenseVector,V1>::value && is<Scalar,V2>::value >::type >
  : DenseVector
  {} ;

} // namespace glas2

#endif
