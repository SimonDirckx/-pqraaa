//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_matrix_type_indirect_matrix_hpp
#define glas2_matrix_type_indirect_matrix_hpp

#include <glas2/matrix/algorithm/assign.hpp>
#include <glas2/matrix/algorithm/select.hpp>
#include <glas2/vector/type/range.hpp>
#include <glas2/vector/type/slice.hpp>
#include <glas2/vector/type/all.hpp>
#include <glas2/vector/type/indirect_vector.hpp>
#include <glas2/matrix/view/matrix_selection.hpp>
#include <glas2/matrix/concept/dense_matrix.hpp>
#include <glas2/concept/concept.hpp>
#include <glas2/concept/is.hpp>
#include <type_traits>
#include <cassert>

namespace glas2 {

  template <typename M, typename R, typename C>
  class indirect_matrix {
    public:
      typedef typename std::decay<M>::type      matrix_type ;
      typedef typename std::decay<R>::type      rows_type ;
      typedef typename std::decay<C>::type      columns_type ;
      typedef typename matrix_type::value_type  value_type ;
      typedef typename matrix_type::size_type   size_type ;
      typedef typename matrix_type::orientation orientation ;

    public:
      indirect_matrix( M m, R r, C c )
      : matrix_( m )
      , rows_( r )
      , columns_( c )
      {}

      // Copy reference !!
      indirect_matrix( indirect_matrix const& that )
      : matrix_( that.matrix_ )
      , rows_( that.rows_ )
      , columns_( that.columns_ )
      {}

    public:
      matrix_type const& matrix() const { return matrix_ ; }
      rows_type const& rows() const { return rows_ ; }
      columns_type const& columns() const { return columns_ ; }

      typename std::add_lvalue_reference<M>::type matrix() { return matrix_ ; }
      typename std::add_lvalue_reference<R>::type rows() { return rows_ ; }
      typename std::add_lvalue_reference<C>::type columns() { return columns_ ; }

    public:
      // Matrix
      size_type num_rows() const { return rows_.size() ; }
      size_type num_columns() const { return columns_.size() ; }

      value_type const& operator() ( size_type i, size_type j ) const {
        assert( i>=0 && i<num_rows() ) ;
        assert( j>=0 && j<num_columns() ) ;
        return matrix_( rows_(i), columns_(j) ) ;
      }

      value_type& operator() ( size_type i, size_type j ) {
        assert( i>=0 && i<num_rows() ) ;
        assert( j>=0 && j<num_columns() ) ;
        return matrix_( rows_(i), columns_(j) ) ;
      }

      indirect_matrix& operator=( indirect_matrix const& that ) {
        assign( *this, that ) ;
        return *this ;
      }

      template <typename E>
      indirect_matrix operator=( E const& that ) {
        return assign( *this, that ) ;
      }

      template <typename I1, typename I2>
      typename std::enable_if< !(std::is_integral<I1>::value && std::is_integral<I2>::value)
                             , typename matrix_selection< indirect_matrix, I1, I2 >::result_type
                             >::type operator()( I1 const& s1, I2 const& s2 ) const {
        return matrix_selection< indirect_matrix, I1, I2 >::apply( *this, s1, s2 ) ;
      }

    private:
      M matrix_ ;
      R rows_ ;
      C columns_ ;
  } ;


  template <typename MM, typename RR, typename CC, typename R, typename C>
  struct matrix_selection< indirect_matrix<MM,RR,CC>, R, C,
        typename std::enable_if< is<DenseVector,R>::value && is<DenseVector,C>::value >::type > {
    typedef indirect_matrix<MM,RR,CC>                                                                                 matrix_type ;
    typedef typename matrix_selection< typename matrix_type::matrix_type
                                     , typename vector_selection< typename matrix_type::rows_type, R >::result_type
                                     , typename vector_selection< typename matrix_type::columns_type, C >::result_type
                                     >::result_type                                                                   result_type ;

    static result_type apply( matrix_type m, R const& r, C const& c ) {
      return m.matrix()( m.rows()(r), m.columns()(c) ) ;
    }
  } ;


  template <typename MM, typename RR, typename CC, typename R, typename I>
  struct matrix_selection< indirect_matrix<MM,RR,CC>, R, I,
        typename std::enable_if< is<DenseVector,R>::value && std::is_integral<I>::value >::type > {
    typedef indirect_matrix<MM,RR,CC>                                                                                 matrix_type ;
    typedef typename matrix_selection< typename matrix_type::matrix_type
                                     , typename vector_selection< typename matrix_type::rows_type, R >::result_type
                                     , typename matrix_type::columns_type::value_type
                                     >::result_type                                                                   result_type ;

    static result_type apply( matrix_type m, R const& r, I c ) {
      return m.matrix()( m.rows()(r), m.columns()(c) ) ;
    }
  } ;


  template <typename MM, typename RR, typename CC, typename C, typename I>
  struct matrix_selection< indirect_matrix<MM,RR,CC>, I, C,
        typename std::enable_if< is<DenseVector,C>::value && std::is_integral<I>::value >::type > {
    typedef indirect_matrix<MM,RR,CC>                                                                                 matrix_type ;
    typedef typename matrix_selection< typename matrix_type::matrix_type
                                     , typename matrix_type::rows_type::value_type
                                     , typename vector_selection< typename matrix_type::columns_type, C >::result_type
                                     >::result_type                                                                   result_type ;

    static result_type apply( matrix_type m, I r, C const& c ) {
      return m.matrix()( m.rows()(r), m.columns()(c) ) ;
    }
  } ;


  template <typename M, typename R, typename C>
  struct glas_concept< indirect_matrix<M,R,C> >
  : DenseMatrix
  {};

} // namespace glas2


#endif
