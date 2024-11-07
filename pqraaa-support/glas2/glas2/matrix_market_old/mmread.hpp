//  (C) Copyright Karl Meerbergen 2009.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_matrix_market_mmread_hpp
#define glas2_matrix_market_mmread_hpp

#include <glas2/matrix_market/type/mminformation.hpp>
#include <glas2/sparse/concept/coordinate_sparse_matrix.hpp>
#include <glas2/concept/is.hpp>

namespace glas2 { namespace matrix_market {

  template <typename M>
  typename std::enable_if< is< CoordinateSparseMatrix, M>::value, std::istream& >::type mmread( M m, std::istream& is ) {
    mminformation info( is ) ;

    assert( mmreader.rep() == matrix_market::COORDINATE ) ;
    assert( mmreader.symmetry() == matrix_market::SYMMETRIC || mmreader.symmetry() == matrix_market::GENERAL ) ;

    m.reset( info.num_rows(), info.num_columns(), info.num_nz() ) ;

    typedef typename glas2::value_type<C>::type value_type ;
    int row, col ;
    double r_val, i_val( 0.0 ) ;

    for ( int i=0; i<mmreader.num_nz(); ++i ) {
      is >> row >> col >> r_val ; --row ; --col ;
      if ( mmreader.field() == matrix_market::COMPLEX ) mmreader.stream() >> i_val ;
      if ( matrix_market::part_of( e, row, col ) ) glas2::push_back( c, row, col, matrix_market::value_map<value_type>( e, r_val, i_val ) ) ;
      if (mmreader.symmetry()==matrix_market::SYMMETRIC && row!=col) {
        if ( matrix_market::part_of( e, col, row ) ) glas2::push_back( c, col, row, matrix_market::value_map<value_type>( e, r_val, i_val ) ) ;
      }
    }

    m.sort_and_compress() ;

    return is ;
  }

} } // namespace glas2::matrix_market


namespace glas2 {

  template <typename A, typename C, typename E>
  struct backend_assign_functor< backend_glas2_tag, A, C, E
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

      typedef typename glas2::value_type<C>::type value_type ;
      double r_val, i_val( 0.0 ) ;

      for ( int i=0; i<e.num_columns(); ++i ) {
        typename glas2_backend::cursor_result_type< typename glas2::column_result_type<C>::type >::type c_col( glas2_backend::cursor( glas2::column( c, i ) ) ) ;
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
  struct backend_assign_functor< backend_glas2_tag, assign_tag, C, E
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

      typedef typename glas2::value_type<C>::type value_type ;
      int row, col ;
      double r_val, i_val( 0.0 ) ;

      for ( int i=0; i<mmreader.num_nz(); ++i ) {
        mmreader.stream() >> row >> col >> r_val ; --row ; --col ;
        if ( mmreader.field() == matrix_market::COMPLEX ) mmreader.stream() >> i_val ;
        if ( matrix_market::part_of( e, row, col ) ) glas2::push_back( c, row, col, matrix_market::value_map<value_type>( e, r_val, i_val ) ) ;
        if (mmreader.symmetry()==matrix_market::SYMMETRIC && row!=col) {
          if ( matrix_market::part_of( e, col, row ) ) glas2::push_back( c, col, row, matrix_market::value_map<value_type>( e, r_val, i_val ) ) ;
        }
      }

      glas2::sort_and_compress( c ) ;
      return c ;
    }
  } ;

  template <typename C, typename E>
  struct backend_assign_functor< backend_glas2_tag, assign_tag, C, E
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
        if ( matrix_market::part_of( e, row, col) ) glas2::push_back( c, row, col ) ;
        if (mmreader.symmetry()==matrix_market::SYMMETRIC && row!=col) {
          if ( matrix_market::part_of( e, col, row) ) glas2::push_back( c, col, row ) ;
        }
      }

      glas2::sort_and_compress( c ) ;
      return c ;
    }
  } ;

  //
  // Matrix whose structure is already created.
  //
  template <typename C, typename E>
  struct backend_assign_functor< backend_glas2_tag, assign_tag, C, E
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

      glas2::backend_glas2::assign( glas2::value_array(c), 0.0 ) ; 

      typedef typename glas2::value_type<C>::type value_type ;
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

} // namespace glas2

#endif
