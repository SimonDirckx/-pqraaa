//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_vector_expression_vector_scalar_expression_hpp
#define glas2_vector_expression_vector_scalar_expression_hpp

#include <glas2/concept/is.hpp>
#include <glas2/concept/concept.hpp>
#include <glas2/type/pass_reference.hpp>
#include <glas2/scalar/concept/scalar.hpp>
#include <glas2/vector/concept/dense_vector.hpp>
#include <type_traits>

namespace glas2 {

  template <typename V, typename S, typename Op>
  class vector_scalar_expression {
    public:
      vector_scalar_expression( V const& v, S const& s )
      : vector_( v )
      , scalar_(s)
      {}

      static_assert( is<Scalar,S>::value, "vector_scalar_expression: S is not a Scalar" );
      static_assert( is<DenseVector,V>::value, "vector_scalar_expression: V is not a DenseVector" );

    public:
      typedef typename V::size_type                              size_type ;
      typedef decltype( Op() ( typename V::value_type(), S() ) ) value_type ;

      size_type size() const {
        return vector_.size() ;
      }
      value_type operator() ( size_type i) const { return Op()( vector_(i), scalar_ ) ; }
      value_type operator[] ( size_type i) const { return Op()( vector_[i], scalar_ ) ; }

    private:
      typedef typename pass_reference<V>::type v_ref ;
      v_ref    vector_ ;
      S        scalar_ ;
  } ;

  template <typename V, typename S, typename Op>
  struct glas_concept< vector_scalar_expression<V,S,Op>, typename std::enable_if< is<DenseVector,V>::value >::type >
  : DenseVector
  {} ;

} // namespace glas2

#endif
