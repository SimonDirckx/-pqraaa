//  (C) Copyright Karl Meerbergen 2016.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_vector_expression_scalar_kron_expression_hpp
#define glas2_vector_expression_scalar_kron_expression_hpp

#include <glas2/concept/is.hpp>
#include <glas2/concept/concept.hpp>
#include <glas2/scalar/concept/scalar.hpp>
#include <glas2/vector/concept/dense_vector.hpp>
#include <type_traits>

namespace glas2 {

  template <typename V1, typename V2>
  class scalar_concat_expression {
    public:
      scalar_concat_expression( V1 const& v1, V2 const& v2 )
      : v1_( v1 )
      , v2_( v2 )
      {}

      static_assert( is<Scalar,V1>::value, "scalar_concat_expression: V1 is not a Scalar" );
      static_assert( is<Scalar,V2>::value, "scalar_concat_expression: V2 is not a Scalar" );

    public:
      typedef short                    size_type ;
      typedef decltype( V1() * V2() )  value_type ;

      size_type size() const { return 2 ; }

      value_type operator() ( size_type i ) const {
        assert( i>=0 && i<2 ) ;
        if (i==0) return v1_ ;
        else return v2_ ;
      }

    private:
      V1 v1_ ;
      V2 v2_ ;
  } ;

  template <typename V1, typename V2>
  struct glas_concept< scalar_concat_expression<V1,V2>, typename std::enable_if< is<Scalar,V1>::value && is<Scalar,V2>::value >::type >
  : DenseVector
  {} ;

} // namespace glas2

#endif
