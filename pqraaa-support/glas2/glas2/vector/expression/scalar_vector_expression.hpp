//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_vector_expression_scalar_vector_expression_hpp
#define glas2_vector_expression_scalar_vector_expression_hpp

#include <glas2/concept/is.hpp>
#include <glas2/concept/concept.hpp>
#include <glas2/type/pass_reference.hpp>
#include <glas2/vector/concept/dense_vector.hpp>
#include <type_traits>

namespace glas2 {

  template <typename S, typename V, typename Op>
  class scalar_vector_expression {
    public:
      scalar_vector_expression( S const& s, V const& v )
      : scalar_(s)
      , vector_( v )
      {}

    public:
      typedef typename V::size_type                              size_type ;
      typedef decltype( Op() ( S(), typename V::value_type() ) ) value_type ;

      size_type size() const {
        return vector_.size() ;
      }
      value_type operator() ( size_type i) const { return Op()( scalar_, vector_(i) ) ; }
      value_type operator[] ( size_type i) const { return Op()( scalar_, vector_[i] ) ; }

    public:
      typedef typename pass_reference<V>::type v_ref ;
      S      scalar_ ;
      v_ref  vector_ ;
  } ;

  template <typename S, typename V, typename Op>
  struct glas_concept< scalar_vector_expression<S,V,Op>, typename std::enable_if< is<DenseVector,V>::value >::type >
  : DenseVector
  {} ;

} // namespace glas2

#endif
