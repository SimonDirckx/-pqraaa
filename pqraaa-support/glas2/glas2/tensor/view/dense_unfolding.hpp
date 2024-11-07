//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_tensor_view_dense_unfolding_hpp
#define glas2_tensor_view_dense_unfolding_hpp

#include <glas2/matrix/algorithm/assign.hpp>
#include <glas2/matrix/algorithm/select.hpp>
#include <glas2/vector/algorithm/prod.hpp>
#include <glas2/vector/type/range.hpp>
#include <glas2/vector/type/slice.hpp>
#include <glas2/vector/type/all.hpp>
#include <glas2/vector/type/strided_vector.hpp>
#include <glas2/matrix/view/matrix_selection.hpp>
#include <glas2/matrix/concept/dense_matrix.hpp>
#include <glas2/concept/concept.hpp>
#include <glas2/concept/is.hpp>
#include <type_traits>
#include <cassert>

namespace glas2 {

  template <typename T>
  class dense_unfolding {
    public:
      typedef T                      tensor_type ;
      typedef typename T::value_type value_type ;
      typedef typename T::size_type  size_type ;

    public:
      dense_unfolding( T tensor, size_type coordinate )
      : tensor_( tensor )
      , coordinate_( coordinate )
      {
        assert( coordinate>0 && coordinate<tensor.order() ) ;
      }

      // Copy reference !!
      double_strided_matrix( double_strided_matrix const& that )
      : tensor_( that.tensor_ )
      , coordinate_( that.coordinate_ )
      {}

    public:
      tensor_type tensor() const { return tensor_ ; }

      // Matrix
      size_type num_rows() const { return tensor_.shape()(coordinate_) ; }
      size_type num_columns() const { return glas2::product( tensor_.shape() ) / num_rows() ; }

      size_type row_stride() const { 
        size_type stride = 1 ;
        for (size_type i=0; i<coordinate; ++i) stride_rows *= tensor_.shape()(i) ;
      }

      size_type column_stride( size_type row ) const { 
      }

      value_type const& operator() ( size_type i, size_type j ) const {
        assert( i>=0 && i<num_rows_ ) ;
        assert( j>=0 && j<num_columns_ ) ;
      }

      value_type& operator() ( size_type i, size_type j ) {
        assert( i>=0 && i<num_rows_ ) ;
        assert( j>=0 && j<num_columns_ ) ;
      }

      dense_unfolding& operator=( dense_unfolding const& that ) {
        assert( that.num_rows()==num_rows() ) ;
        assert( that.num_columns()==num_columns() ) ;
        return assign( *this, that ) ;
      }

      template <typename E>
      double_strided_matrix operator=( E const& that ) {
        assert( that.num_rows()==num_rows() ) ;
        assert( that.num_columns()==num_columns() ) ;
        return assign( *this, that ) ;
      }

      template <typename I1, typename I2>
      typename std::enable_if< !(std::is_integral<I1>::value && std::is_integral<I2>::value)
                             , typename matrix_selection< double_strided_matrix, I1, I2 >::result_type
                             >::type operator()( I1 const& s1, I2 const& s2 ) {
        return matrix_selection< double_strided_matrix, I1, I2 >::apply( *this, s1, s2 ) ;
      }

      template <typename I1, typename I2>
      typename std::enable_if< !(std::is_integral<I1>::value && std::is_integral<I2>::value)
                             , typename matrix_selection< double_strided_matrix, I1, I2 >::result_type
                             >::type operator()( I1 const& s1, I2 const& s2 ) const {
        return matrix_selection< double_strided_matrix, I1, I2 >::apply( *this, s1, s2 ) ;
      }

    private:
      tensor_type tensor_ ;
      size_type   coordinate_ ;
  } ;

  template <typename T>
  struct concept< dense_unfolding<T> >
  : DenseMatrix
  {};


  template <typename T, typename S, typename O>
  struct matrix_selection< double_strided_matrix<T,S,O>, all, all > {
    typedef double_strided_matrix<T,S,O> result_type ;

    static result_type apply( double_strided_matrix<T,S,O> m, all, all ) {
      return m ;
    }
  } ;


  template <typename T, typename S, typename O, typename R, typename I>
  struct matrix_selection< double_strided_matrix<T,S,O>, R, I,
        typename std::enable_if< is<DenseVector,R>::value && std::is_integral<I>::value >::type > {
    typedef typename vector_selection< strided_vector<T,S>, R>::result_type result_type ;

    static result_type apply( double_strided_matrix<T,S,O> m, R const& r, I col ) {
      assert( col<m.num_columns() && col>=0 ) ;
      return strided_vector<T,S>( &m(0,col), m.stride_columns(), m.num_rows() )( r ) ;
    }
  } ;


  template <typename T, typename S, typename O, typename I, typename R>
  struct matrix_selection< double_strided_matrix<T,S,O>, I, R,
        typename std::enable_if< is<DenseVector,R>::value && std::is_integral<I>::value >::type > {
    typedef typename vector_selection< strided_vector<T,S>, R>::result_type result_type ;

    static result_type apply( double_strided_matrix<T,S,O> m, I row, R const& r ) {
      assert( row<m.num_rows() && row>=0 ) ;
      return strided_vector<T,S>( &m(row,0), m.stride_rows(), m.num_columns() )( r ) ;
    }
  } ;


  template <typename T, typename S, typename O, typename R2>
  struct matrix_selection< double_strided_matrix<T,S,O>, range, R2 > {
    typedef double_strided_matrix<T,S,O>                                 temp_type ;
    typedef typename matrix_selection< temp_type, all, R2 >::result_type result_type ;

    static result_type apply( double_strided_matrix<T,S,O> m, range const& r, R2 const& r2 ) {
      assert( r.begin()+r.size()<=m.num_rows() ) ;
      return result_type( &m(r.begin(),0), m.stride_rows(), m.stride_columns(), r.size(), m.num_columns() )( all(), r2 ) ;
    }
  } ;



  template <typename T, typename S, typename O>
  struct matrix_selection< double_strided_matrix<T,S,O>, all, range > {
    typedef double_strided_matrix<T,S,O> result_type ;

    static result_type apply( double_strided_matrix<T,S,O> m, all, range const& r ) {
      assert( r.begin()+r.size()<=m.num_columns() ) ;
      return result_type( m.ptr()+r.begin()*m.stride_columns(), m.stride_rows(), m.stride_columns(), m.num_rows(), r.size() ) ;
    }
  } ;


/*
  template <typename T, typename S>
  struct matrix_transform< select_tag, double_strided_matrix<T,S,column_major> > {
    typedef double_strided_matrix<T,S,row_major> result_type ;
    static result_type apply( double_strided_matrix<T,S,column_major> m ) {
      return result_type( m.ptr(), m.num_columns(), m.num_rows() ) ;
    }
  } ;
*/
} // namespace glas2


#endif
