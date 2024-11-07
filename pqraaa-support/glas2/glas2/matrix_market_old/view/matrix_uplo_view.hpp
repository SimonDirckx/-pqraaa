//  (C) Copyright Karl Meerbergen 2009.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_matrix_market_view_matrix_uplo_view_hpp
#define glas2_matrix_market_view_matrix_uplo_view_hpp

/*#ifdef GLAS_COMPLEX
#include <glas2/valuetype/complex.hpp>
#endif*/

#include <glas2/matrix_market/concept/matrix_market_expression.hpp>
#include <glas2/matrix_market/concept/mmreader.hpp>
#include <glas2/matrix_market/concept/value_map.hpp>
#include <glas2/matrix_market/concept/part_of.hpp>
#include <glas2/concept/value_type.hpp>
#include <iosfwd>

#undef GLAS_INCLUDE_LEVEL
#define GLAS_INCLUDE_LEVEL GLAS_INCLUDE_TOOLBOX

namespace glas2 { namespace matrix_market {

  template <typename M, typename S>
  struct MatrixMarketExpression< glas2::matrix_uplo_view<M,S> >
  : MatrixMarketExpression< M >
  {} ;


  template <typename M, typename S>
  struct mmreader_functor< glas2::matrix_uplo_view<M,S> >
  {
    typedef typename mmreader_result_type< M >::type result_type ;

    result_type operator() ( glas2::matrix_uplo_view<M,S> const& e ) const {
      return mmreader( e.argument() ) ;
    }
  } ;
  template <typename M, typename S>
  struct mmreader_functor< glas2::matrix_uplo_view<M,S> const >
  : mmreader_functor< glas2::matrix_uplo_view<M,S> >
  {} ;


  template <typename R, typename M, typename S>
  struct value_map_functor< R, glas2::matrix_uplo_view<M,S> >
  {
    typedef R result_type ;

    R operator() ( glas2::matrix_uplo_view<M,S> const& e, matrix_market_value_type value ) const {
      return value_map< R >( e.argument(), value ) ;
    }
  } ;


  template <typename M>
  struct part_of_functor< glas2::matrix_uplo_view<M,glas2::upper_tag> >
  {
    typedef bool result_type ;

    bool operator() ( glas2::matrix_uplo_view<M,glas2::upper_tag> const& e, int row, int column ) const {
      return row<=column ;
    }
  } ;

  template <typename M>
  struct part_of_functor< glas2::matrix_uplo_view<M,glas2::lower_tag> >
  {
    typedef bool result_type ;

    bool operator() ( glas2::matrix_uplo_view<M,glas2::lower_tag> const& e, int row, int column ) const {
      return row>=column ;
    }
  } ;

  template <typename M>
  struct part_of_functor< glas2::matrix_uplo_view<M,glas2::strictly_lower_tag> >
  {
    typedef bool result_type ;

    bool operator() ( glas2::matrix_uplo_view<M,glas2::strictly_lower_tag> const& e, int row, int column ) const {
      return row>column ;
    }
  } ;

  template <typename M>
  struct part_of_functor< glas2::matrix_uplo_view<M,glas2::strictly_upper_tag> >
  {
    typedef bool result_type ;

    bool operator() ( glas2::matrix_uplo_view<M,glas2::strictly_upper_tag> const& e, int row, int column ) const {
      return row<column ;
    }
  } ;

} } // namespace glas2

#endif
