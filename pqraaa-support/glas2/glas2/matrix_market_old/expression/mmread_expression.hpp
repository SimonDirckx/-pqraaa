//  (C) Copyright Karl Meerbergen 2009.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_matrix_market_expression_mmread_expression_hpp
#define glas2_matrix_market_expression_mmread_expression_hpp

/*#ifdef GLAS_COMPLEX
#include <glas2/scalar/complex.hpp>
#endif*/

#include <glas2/matrix_market/concept/matrix_market_expression.hpp>
#include <glas2/matrix_market/concept/mmreader.hpp>
#include <glas2/matrix_market/concept/value_map.hpp>
#include <glas2/matrix_market/type/mminformation.hpp>
#include <glas2/concept/concept.hpp>
#include <iosfwd>

#undef GLAS_INCLUDE_LEVEL
#define GLAS_INCLUDE_LEVEL GLAS_INCLUDE_TOOLBOX

namespace glas2 { namespace matrix_market {

  class mmread_expression
  : public mminformation
  {
    public:
      typedef matrix_market_value_type        value_type ;
      typedef matrix_market_value_type&       reference ;
      typedef matrix_market_value_type const& const_reference ;
      typedef general_tag                     structure_type ;

    public:
      inline mmread_expression( std::istream& is )
      : mminformation( is )
      , is_( &is )
      {}

    public:
      mminformation const& info() const { return *this ; }

    public:
      inline std::istream& stream() const { return *is_ ; }

    private:
      std::istream* is_ ;
  } ;

  template <>
  struct MatrixMarketExpression< mmread_expression >
  : boost::mpl::true_
  {} ;

  template <>
  struct value_map_functor< double, mmread_expression >
  {
    typedef double result_type ;

    inline result_type operator() ( mmread_expression const& e, matrix_market_value_type v ) const {
      assert( glas2::imag(v)==0.0 ) ;
      return glas2::real(v) ;
    }
  } ;

#ifdef GLAS_COMPLEX
  template <>
  struct value_map_functor< std::complex< double >, mmread_expression >
  {
    typedef std::complex< double > result_type ;

    inline result_type operator() ( mmread_expression const& e, matrix_market_value_type v ) const {
      return v ;
    }
  } ;
#endif

} } // namespace glas2::matrix_market

namespace glas2 {

  template <>
  struct MatrixExpression< matrix_market::mmread_expression >
  : boost::mpl::true_
  {} ;

  template <>
  struct pass_by_value_argument< matrix_market::mmread_expression >
  : boost::mpl::true_
  {} ;

} // namespace glas2

#endif
