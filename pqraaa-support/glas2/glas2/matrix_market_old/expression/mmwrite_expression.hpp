//  (C) Copyright Karl Meerbergen 2009.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_matrix_market_expression_mmwrite_expression_hpp
#define glas2_matrix_market_expression_mmwrite_expression_hpp

#include <glas2/concept/is.hpp>
#include <glas2/dense/concept/dense_matrix.hpp>
#include <glas2/sparse/concept/coordinate_sparse_matrix.hpp>
#include <glas2/matrix_market/type/mminformation.hpp>
#include <iostream>
#include <iomanip>
#include <limits>
#include <string>
#include <sstream>
#include <type_traits>

namespace glas2 { namespace matrix_market {

  template <typename T>
  struct value_type_string {
    static std::string value() { return "real" ; }
  } ;

  template <typename T>
  struct print_value {
    std::string operator() ( T const& e ) const {
      std::ostringstream s ; s << std::setprecision(16) << e ;
      return s.str() ;
    }
  } ;

#ifdef GLAS_COMPLEX
  template <typename T>
  struct value_type_string< std::complex<T> > {
    static std::string value() { return "complex" ; }
  } ;

  template <typename T>
  struct print_value< std::complex<T> > {
    std::string operator()( std::complex<T> const& e ) const {
      std::ostringstream s ; s << std::setprecision(16) << e.real() << " " << e.imag() ;
      return s.str() ;
    }
  } ;
#endif

  template <typename E>
  class mmwrite_expression
  {
    public:
      inline mmwrite_expression( E const& e, std::string const& string )
      : e_( e )
      , string_( string )
      {}

    public:
      inline std::string const& string() const { return string_ ; }
      inline E const& expression() const { return e_ ; }

    private:
      E const&           e_ ;
      std::string const& string_ ;
  } ;


  template <typename E>
  typename std::enable_if< is< CoordinateSparseMatrix, E >::value, std::ostream& >::type operator<<( std::ostream& os, mmwrite_expression<E> const& e )
  {
    os << "%%MatrixMarket matrix coordinate " << value_type_string< typename glas2::value_type<E>::type >::value() << " general\n" ;
    os << "% " << e.string() << "\n" ;
    os << glas2::num_rows(e.expression()) << " " << glas2::num_columns(e.expression()) << " " << glas2::nnz(e.expression()) << "\n" ;
    for ( typename glas2::nnz_type<E>::type i=0; i<glas2::nnz(e.expression()); ++i ) {
      os << glas2::at(glas2::row_index_array(e.expression()),i)+1-glas2::index_base<E>::value << " "
                << glas2::at(glas2::column_index_array(e.expression()),i)+1-glas2::index_base<E>::value << " "
                << print_value<typename glas2::value_type<E>::type>() ( glas2::at(glas2::value_array(e.expression()),i) ) << "\n" ;
    }

    return os ;
  }


  template <typename E>
  typename std::enable_if< is< DenseMatrixExpression, E >::value, std::ostream& >::type operator<<( std::ostream& os, mmwrite_expression<E> const& e )
  {
    os << "%%MatrixMarket matrix array " << value_type_string< typename E::value_type >::value() << " general\n" ;
    os << "% " << e.string() << "\n" ;
    os << e.expression().num_rows() << " " << e.expression().num_columns() << "\n" ;
    for ( typename E::size_type i=0; i<e.expression().num_columns(); ++i ) {
      for ( typename E::size_type j=0; j<e.expression().num_rows(); ++j ) {
        os << print_value<typename E::value_type>() ( e.expression()( j, i ) ) << "\n" ;
      }
    }

    return os ;
  }


  template <typename E>
  typename std::enable_if< is< DenseVectorExpression, E >::value, std::ostream& >::type operator<<( std::ostream& os, mmwrite_expression<E> const& e )
  {
    os << "%%MatrixMarket matrix array " << value_type_string< typename E::value_type >::value() << " general\n" ;
    os << "% " << e.string() << "\n" ;
    os << e.expression().size() << " 1\n" ;
    for ( typename E::size_type i=0; i<e.expression().size(); ++i ) {
      os << print_value<typename E::value_type>() ( e.expression()( i ) ) << "\n" ;
    }

    return os ;
  }


} } // namespace glas2::matrix_market

#endif
