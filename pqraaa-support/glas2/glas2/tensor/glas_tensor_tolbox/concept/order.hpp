//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas_toolbox_tensor_concept_order_hpp
#define glas_toolbox_tensor_concept_order_hpp

#include <glas/config.hpp>
#define GLAS_MY_INCLUDE_LEVEL GLAS_INCLUDE_ALGORITHM_FWD
#include <glas/concept/check_include_level.hpp>

#include <glas/toolbox/tensor/concept/tensor_expression.hpp>
#include <glas/concept/strip_const_reference.hpp>

#undef GLAS_INCLUDE_LEVEL
#define GLAS_INCLUDE_LEVEL GLAS_INCLUDE_ALGORITHM_FWD

namespace glas {

  template <typename X>
  struct order_functor {
    typedef typename X::order_type result_type ;

    result_type const& operator() (X const& x) const {
      return x.order() ;
    }
  } ;

  template <typename E, typename Enable=void>
  struct order_type
  : result_of< order_functor< typename analyse_argument< E >::result_type > >
  {} ;

  template <typename X>
  typename boost::lazy_enable_if< TensorExpression< X >, order_type< X > >::type
  order( X const& x ) {
    return order_functor<X>()( x ) ;
  }

} // namespace glas::algorithm

#endif
