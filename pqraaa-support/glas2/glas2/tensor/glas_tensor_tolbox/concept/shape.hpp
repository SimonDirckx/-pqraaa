//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas_toolbox_tensor_concept_shape_hpp
#define glas_toolbox_tensor_concept_shape_hpp

#include <glas/config.hpp>
#define GLAS_MY_INCLUDE_LEVEL GLAS_INCLUDE_ALGORITHM_FWD
#include <glas/concept/check_include_level.hpp>

#include <glas/toolbox/tensor/concept/tensor_expression.hpp>
#include <glas/concept/strip_const_reference.hpp>

#undef GLAS_INCLUDE_LEVEL
#define GLAS_INCLUDE_LEVEL GLAS_INCLUDE_ALGORITHM_FWD

namespace glas {

  template <typename X>
  struct shape_functor {
    typedef typename X::shape_type const& result_type ;

    result_type operator() (X const& x) const {
      return x.shape() ;
    }
  } ;

  template <typename E, typename Enable=void>
  struct shape_result_type
  : result_of< shape_functor< typename analyse_argument< E >::result_type > >
  {} ;

  template <typename E, typename Enable=void>
  struct shape_type
  : strip_const_reference< typename shape_result_type<E>::type >
  {} ;

  template <typename X>
  typename boost::lazy_enable_if< TensorExpression< X >, shape_result_type< X > >::type
  shape( X const& x ) {
    return shape_functor<X>()( x ) ;
  }

} // namespace glas::algorithm

#endif
