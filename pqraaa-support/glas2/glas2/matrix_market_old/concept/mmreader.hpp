//  (C) Copyright Karl Meerbergen (2009).
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_matrix_market_concept_mmreader_hpp
#define glas2_matrix_market_concept_mmreader_hpp

#include <type_traits>

namespace glas2 { namespace matrix_market { 

template <typename X, typename EnableIf=void>
struct mmreader_functor {
  typedef typename const_closure_type< X >::type result_type ;

  result_type operator()( X const& x) const {
    return x ;
  } // operator()()
} ; // struct mmreader_functor

template <typename X>
typename mmreader_result_type<X>::type mmreader( X x) {
  return mmreader_functor<X>() ( x ) ;
} // mmreader()

} } // namespace glas2::matrix_market

#endif
