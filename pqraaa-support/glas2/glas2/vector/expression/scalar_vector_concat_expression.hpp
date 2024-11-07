//  (C) Copyright Karl Meerbergen 2016.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_scalar_vector_expression_vector_concat_expression_hpp
#define glas2_scalar_vector_expression_vector_concat_expression_hpp

#include <glas2/concept/is.hpp>
#include <glas2/concept/concept.hpp>
#include <glas2/type/pass_reference.hpp>
#include <glas2/scalar/concept/scalar.hpp>
#include <glas2/vector/concept/dense_vector.hpp>
#include <type_traits>

namespace glas2 {

  template <typename V1, typename V2>
  class scalar_vector_concat_expression {
    public:
      scalar_vector_concat_expression( V1 const& v1, V2 const& v2 )
      : v1_( v1 )
      , v2_( v2 )
      {}

      static_assert( is<Scalar,V1>::value, "vector_concat_expression: V1 is not a Scalar" );
      static_assert( is<DenseVector,V2>::value, "vector_concat_expression: V2 is not a Vector" );

    public:
      typedef typename V2::size_type                        size_type ;
      typedef decltype( V1() + typename V2::value_type() )  value_type ;

      size_type size() const {
        return 1+v2_.size() ;
      }
      value_type operator() ( size_type i ) const {
        assert( i>=0 && i < size() ) ;
        if (i==0) return v1_ ;
        else return v2_( i-1 ) ;
      }

    private:
      typedef typename pass_reference<V2>::type v2_ref ;
      V1     v1_ ;
      v2_ref v2_ ;
  } ;

  template <typename V1, typename V2>
  struct glas_concept< scalar_vector_concat_expression<V1,V2>, typename std::enable_if< is<Scalar,V1>::value && is<DenseVector,V2>::value >::type >
  : DenseVector
  {} ;

} // namespace glas2

#endif
