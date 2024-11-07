//  (C) Copyright Karl Meerbergen 2009.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas_toolbox_matrix_market_algorithm_mmread_hpp
#define glas_toolbox_matrix_market_algorithm_mmread_hpp

#include <glas/config.hpp>
#define GLAS_MY_INCLUDE_LEVEL GLAS_INCLUDE_TOOLBOX
#include <glas/concept/check_include_level.hpp>

#include <glas/toolbox/matrix_market/concept/matrix_market_expression.hpp>
#include <glas/toolbox/matrix_market/concept/mmreader.hpp>
#include <glas/toolbox/matrix_market/concept/value_map.hpp>
#include <glas/toolbox/matrix_market/concept/part_of.hpp>
#include <glas/toolbox/matrix_market/expression/mmread_expression.hpp>
#include <glas/toolbox/matrix_market/tools/read_value.hpp>
#include <glas/sparse/algorithm/push_back.hpp>
#include <glas/sparse/algorithm/sort_and_compress.hpp>
#include <glas/sparse/algorithm/reset.hpp>
#include <glas/algorithm/column.hpp>
#ifdef GLAS_BACKEND_GLAS
#include <glas/backend/glas/concept/cursor.hpp>
#endif
#include <glas/backend_fwd/glas/backend_glas_tag.hpp>
#include <glas/backend_fwd/backend_assign_functor.hpp>
#include <glas/backend_fwd/assign_tags.hpp>
#include <glas/concept/strided_dense_matrix_collection.hpp>
#include <glas/concept/coordinate_sparse_structure_collection.hpp>
#include <glas/concept/compressed_sparse_structure_expression.hpp>
#include <glas/concept/sparse_matrix_collection.hpp>
#include <glas/concept/integral_equal.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/utility/enable_if.hpp>

#undef GLAS_INCLUDE_LEVEL
#define GLAS_INCLUDE_LEVEL GLAS_INCLUDE_TOOLBOX

namespace glas { namespace matrix_market {

  inline mmread_expression mmread( std::istream& is ) {
    return mmread_expression( is ) ;
  }

} } // namespace glas::matrix_market


namespace glas {

  template <typename A, typename C, typename E>
  struct backend_assign_functor< backend_glas_tag, A, C, E
                               , typename boost::enable_if< boost::mpl::and_< matrix_market::MatrixMarketExpression<E>
                                                                            , StridedDenseMatrixCollection< C >
                                                                            > >::type
                               > {
    C& operator() ( C& c, E const& e ) const {
      typename matrix_market::mmreader_result_type<E>::type const& mmreader( matrix_market::mmreader( e ) ) ;
      assert( integral_equal( mmreader.num_rows(), c.num_rows() ) ) ;
      assert( integral_equal( mmreader.num_columns(), c.num_columns() ) ) ;
      assert( integral_equal( mmreader.num_nz(), mmreader.num_rows()*mmreader.num_columns() ) ) ;
      assert( mmreader.rep() == matrix_market::ARRAY ) ;
      assert( mmreader.field() != matrix_market::PATTERN ) ;
      assert( mmreader.symmetry() == matrix_market::GENERAL ) ; // To be done for other formats.

      typedef typename glas::value_type<C>::type value_type ;
      double r_val, i_val( 0.0 ) ;

      for ( int i=0; i<e.num_columns(); ++i ) {
        typename glas_backend::cursor_result_type< typename glas::column_result_type<C>::type >::type c_col( glas_backend::cursor( glas::column( c, i ) ) ) ;
        for ( int j=0; j<e.num_rows(); ++j ) {
          mmreader.stream() >> r_val ;
          if ( mmreader.field() == matrix_market::COMPLEX ) mmreader.stream() >> i_val ;
          if ( matrix_market::part_of( e, j, i ) ) c_col[j] = matrix_market::value_map<value_type>( e, r_val, i_val ) ;
        }
      }
      return c ;
    }
  } ;

  template <typename C, typename E>
  struct backend_assign_functor< backend_glas_tag, assign_tag, C, E
                               , typename boost::enable_if< boost::mpl::and_< matrix_market::MatrixMarketExpression<E>
                                                                            , CoordinateSparseStructureCollection< C >
                                                                            , SparseMatrixCollection< C >
                                                                            >
                                                   >::type
                               > {
    C& operator() ( C& c, E const& e ) const {
      typename matrix_market::mmreader_result_type<E>::type const& mmreader( matrix_market::mmreader( e ) ) ;
      assert( integral_equal( mmreader.num_rows(), c.num_rows() ) ) ;
      assert( integral_equal( mmreader.num_columns(), c.num_columns() ) ) ;
      assert( mmreader.rep() == matrix_market::COORDINATE ) ;
      assert( mmreader.field() != matrix_market::PATTERN ) ;

      reset( c ) ;

      typedef typename glas::value_type<C>::type value_type ;
      int row, col ;
      double r_val, i_val( 0.0 ) ;

      for ( int i=0; i<mmreader.num_nz(); ++i ) {
        mmreader.stream() >> row >> col >> r_val ; --row ; --col ;
        if ( mmreader.field() == matrix_market::COMPLEX ) mmreader.stream() >> i_val ;
        if ( matrix_market::part_of( e, row, col ) ) glas::push_back( c, row, col, matrix_market::value_map<value_type>( e, r_val, i_val ) ) ;
        if (mmreader.symmetry()==matrix_market::SYMMETRIC && row!=col) {
          if ( matrix_market::part_of( e, col, row ) ) glas::push_back( c, col, row, matrix_market::value_map<value_type>( e, r_val, i_val ) ) ;
        }
      }

      glas::sort_and_compress( c ) ;
      return c ;
    }
  } ;

  template <typename C, typename E>
  struct backend_assign_functor< backend_glas_tag, assign_tag, C, E
                               , typename boost::enable_if< boost::mpl::and_< CoordinateSparseStructureCollection< C >
                                                                            , boost::mpl::not_< SparseMatrixExpression< C > >
                                                                            >
                                                   >::type
                               > {
    C& operator() ( C& c, E const& e ) const {
      typename matrix_market::mmreader_result_type<E>::type const& mmreader( matrix_market::mmreader( e ) ) ;
      assert( integral_equal( mmreader.num_rows(), c.num_rows() ) ) ;
      assert( integral_equal( mmreader.num_columns(), c.num_columns() ) ) ;
      assert( mmreader.rep() == matrix_market::COORDINATE ) ;
      assert( mmreader.symmetry() == matrix_market::SYMMETRIC || mmreader.symmetry() == matrix_market::GENERAL ) ;

      reset( c ) ;
      int row, col ;
      std::string str_dump ;

      for ( int i=0; i<mmreader.num_nz(); ++i ) {
        mmreader.stream() >> row >> col ; --row ; --col ;
        std::getline( mmreader.stream(), str_dump ) ;
        if ( matrix_market::part_of( e, row, col) ) glas::push_back( c, row, col ) ;
        if (mmreader.symmetry()==matrix_market::SYMMETRIC && row!=col) {
          if ( matrix_market::part_of( e, col, row) ) glas::push_back( c, col, row ) ;
        }
      }

      glas::sort_and_compress( c ) ;
      return c ;
    }
  } ;

  //
  // Matrix whose structure is already created.
  //
  template <typename C, typename E>
  struct backend_assign_functor< backend_glas_tag, assign_tag, C, E
                               , typename boost::enable_if< boost::mpl::and_< CompressedSparseStructureExpression< C >
                                                                            , SparseMatrixCollection< C >
                                                                            , matrix_market::MatrixMarketExpression< E >
                                                                            >
                                                   >::type
                               > {
    C& operator() ( C& c, matrix_market::mmread_expression const& e ) const {
      typename matrix_market::mmreader_result_type<E>::type const& mmreader( matrix_market::mmreader( e ) ) ;
      assert( integral_equal( mmreader.num_rows(), c.num_rows() ) ) ;
      assert( integral_equal( mmreader.num_columns(), c.num_columns() ) ) ;
      assert( mmreader.rep() == matrix_market::COORDINATE ) ;
      assert( mmreader.field() != matrix_market::PATTERN ) ;
      assert( mmreader.symmetry() == matrix_market::SYMMETRIC || mmreader.symmetry() == matrix_market::GENERAL ) ;

      glas::backend_glas::assign( glas::value_array(c), 0.0 ) ; 

      typedef typename glas::value_type<C>::type value_type ;
      int row, col ;
      double r_val, i_val( 0.0 ) ;

      for ( int i=0; i<mmreader.num_nz(); ++i ) {
        mmreader.stream() >> row >> col >> r_val ; --row ; --col ;
        if ( mmreader.field() == matrix_market::COMPLEX ) mmreader.stream() >> i_val ;
        if ( matrix_market::part_of( e, row, col) ) c( row, col ) += matrix_market::value_map<value_type>( e, r_val, i_val ) ;
        if ( mmreader.symmetry()==matrix_market::SYMMETRIC && row!=col) {
          if ( matrix_market::part_of( e, col, row) ) c( col, row ) += matrix_market::value_map<value_type>( e, r_val, i_val ) ;
        }
      }

      return c ;
    }
  } ;

} // namespace glas

#endif
