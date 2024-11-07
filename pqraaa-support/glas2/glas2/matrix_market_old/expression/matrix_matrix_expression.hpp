//  (C) Copyright Karl Meerbergen 2009.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas_toolbox_matrix_market_expression_matrix_matrix_expression_hpp
#define glas_toolbox_matrix_market_expression_matrix_matrix_expression_hpp

#include <glas/config.hpp>
#define GLAS_MY_INCLUDE_LEVEL GLAS_INCLUDE_TOOLBOX
#include <glas/concept/check_include_level.hpp>

#ifdef GLAS_COMPLEX
#include <glas/valuetype/complex.hpp>
#endif

#include <glas/toolbox/matrix_market/concept/matrix_market_expression.hpp>
#include <glas/toolbox/matrix_market/concept/mmreader.hpp>
#include <glas/toolbox/matrix_market/concept/value_map.hpp>
#include <glas/expression/matrix_matrix_expression.hpp>
#include <glas/concept/num_rows.hpp>
#include <glas/concept/num_columns.hpp>
#include <glas/concept/pass_by_value_argument.hpp>
#include <glas/concept/value_type.hpp>
#include <iosfwd>

#undef GLAS_INCLUDE_LEVEL
#define GLAS_INCLUDE_LEVEL GLAS_INCLUDE_TOOLBOX

namespace glas { namespace matrix_market {

  template <typename M, typename EltOp>
  struct MatrixMarketExpression< glas::matrix_matrix_expression< M, EltOp > >
  : MatrixMarketExpression< M >
  {} ;


  template <typename M, typename EltOp>
  struct mmreader_functor< glas::matrix_matrix_expression< M, EltOp > >
  {
    typedef typename mmreader_result_type< M >::type result_type ;

    result_type operator() ( glas::matrix_matrix_expression< M, EltOp > const& e ) const {
      return mmreader( e.argument() ) ;
    }
  } ;
  template <typename M, typename EltOp>
  struct mmreader_functor< glas::matrix_matrix_expression< M, EltOp > const >
  : mmreader_functor< glas::matrix_matrix_expression< M, EltOp > >
  {} ;


  template <typename R, typename M, typename EltOp>
  struct value_map_functor< R, glas::matrix_matrix_expression< M, EltOp > >
  {
    typedef R result_type ;

    R operator() ( glas::matrix_matrix_expression< M, EltOp > const& e, matrix_market_value_type value ) const {
      return EltOp()( value_map< typename value_type< M >::type >( e.argument(), value ) ) ;
    }
  } ;

} } // namespace glas

#endif
