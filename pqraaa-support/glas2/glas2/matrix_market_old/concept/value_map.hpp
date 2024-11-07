//  (C) Copyright Karl Meerbergen (2009).
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas_toolbox_matrix_market_concept_value_map_hpp
#define glas_toolbox_matrix_market_concept_value_map_hpp

#ifdef GLAS_COMPLEX
#include <glas/valuetype/complex.hpp>
#endif
#include <glas/concept/reference_argument.hpp>
#include <glas/concept/pass_by_value_argument.hpp>
#include <glas/concept/result_of.hpp>
#include <glas/concept/analyse_argument.hpp>
#include <glas/concept/strip_const_reference.hpp>
#include <glas/concept/closure.hpp>
#include <glas/concept/const_closure.hpp>
#include <boost/utility/enable_if.hpp>
#include <cassert>

namespace glas { namespace matrix_market { 

#ifdef GLAS_COMPLEX
  typedef std::complex<double> matrix_market_value_type ;
#else
  typedef double matrix_market_value_type ;
#endif

template <typename R, typename X, typename EnableIf=void>
struct value_map_functor {
  // R operator() ( X const& x, matrix_market_result_type value ) const ;
} ; // struct value_map_functor

template <typename R, typename X>
R value_map( X const& x, matrix_market_value_type const& value ) {
  return value_map_functor<R, X>() ( x, value ) ;
} // value_map()

template <typename R, typename X>
R value_map( X const& x, double const& real, double const& imag ) {
#ifdef GLAS_COMPLEX
  return value_map<R>( x, matrix_market_value_type( real, imag ) ) ;
#else
  assert( imag == 0.0 ) ;
  return value_map<R>( x, real ) ;
#endif
} // value_map()

} } // namespace glas::matrix_market

#endif
