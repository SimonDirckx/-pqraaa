//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_tensor_type_dense_unfolding_hpp
#define glas2_tensor_type_dense_unfolding_hpp

#include <glas2/tensor/concept/tensor_unfolding_matrix.hpp>
#include <glas2/matrix/algorithm/assign.hpp>
#include <glas2/vector/algorithm/fill.hpp>
#include <glas2/vector/algorithm/inner_prod.hpp>
#include <glas2/vector/algorithm/prod.hpp>
#include <glas2/matrix/algorithm/select.hpp>
#include <glas2/vector/type/range.hpp>
#include <glas2/vector/type/slice.hpp>
#include <glas2/vector/type/all.hpp>
#include <glas2/vector/type/contiguous_vector.hpp>
#include <glas2/vector/container/shared_vector.hpp>
#include <glas2/matrix/view/matrix_selection.hpp>
#include <glas2/matrix/concept/dense_matrix.hpp>
#include <glas2/concept/concept.hpp>
#include <glas2/concept/is.hpp>
#include <type_traits>
#include <cassert>
#include <vector>

namespace glas2 {

  template <typename T, typename S>
  class dense_unfolding {
    public:
      typedef T            value_type ;
      typedef S            size_type ;
      typedef column_major orientation ;

    public:
      template <typename Shape, typename Stride>
      dense_unfolding( value_type* ptr, Shape const& shape, Stride const& stride, size_type coordinate )
      : ptr_(ptr)
      , num_rows_( shape( coordinate ) )
      , num_columns_( glas2::prod(shape) / num_rows_ )
//#ifdef NDEBUG
      , columns_( num_columns_ )
//#else
      //, columns_( num_columns_+1 )
//#endif
      {
        assert( coordinate>=0 && coordinate<shape.size() ) ;
        assert( stride.size()==shape.size() ) ;

        row_stride_ = 1 ;
        for (size_type i=0; i<coordinate; ++i) row_stride_ *= shape(i) ;

        columns_(0) = 0 ;
        glas2::shared_vector<size_type> index( shape.size()+1 ) ;
        fill( index, 0 ) ;

        //size_type nnz = num_rows_ * num_columns_ ;
        for (size_type col=0; col<num_columns_; ++col) {
          columns_(col) = inner_prod( stride, index(glas2::range(0,stride.size())) ) ;
          if (coordinate==0) index(0) = shape(0);
          else index(0)++ ;
          int k=0; while ( k<shape.size() && index(k)==shape(k) ) {
            index(k) = 0; ++k;
            if (coordinate==k) index(k) = shape(coordinate);
            else index(k)++ ;
          }
        }
//#ifndef NDEBUG
        //columns_(num_columns_) = nnz ;
//#endif

        //std::cout << "row stride " << row_stride_ << std::endl ;
        //std::cout << "columns " << columns_ << std::endl ;
      }

      template <typename Cols>
      dense_unfolding( value_type* ptr, size_type num_rows, size_type row_stride, Cols const& columns )
      : ptr_( ptr )
      , num_rows_( num_rows )
      , num_columns_( columns.size() )
      , row_stride_( row_stride )
//#ifdef NDEBUG
      , columns_( num_columns_ )
//#else
      //, columns_( num_columns_+1 )
//#endif
      {
        columns_ = columns ;
      }

      // Copy reference !!
      dense_unfolding( dense_unfolding const& that )
      : ptr_(that.ptr_)
      , num_rows_( that.num_rows_ )
      , num_columns_( that.num_columns_ )
      , row_stride_( that.row_stride_ )
      , columns_( that.columns_ )
      {}

    public:
      value_type* ptr() { return ptr_ ; }
      value_type const* const ptr() const { return ptr_ ; }
      size_type stride_rows() const{ return row_stride_ ;}
      size_type begin_column( size_type i ) const{
        assert( i>=0 && i<num_columns_ ) ;
        return columns_(i) ;
      }

      typename shared_vector< size_type >::base_type  columns() const { return columns_ ; }

      // Matrix
      size_type num_rows() const { return num_rows_ ; }
      size_type num_columns() const { return num_columns_ ; }

      value_type const& operator() ( size_type i, size_type j ) const {
        assert( i>=0 && i<num_rows_ ) ;
        assert( j>=0 && j<num_columns_ ) ;
        //assert( i*row_stride_ + columns_(j) < columns_(num_columns_) ) ;
        return ptr_[ i*row_stride_ + columns_(j) ] ;
      }

      value_type& operator() ( size_type i, size_type j ) {
        assert( i>=0 && i<num_rows_ ) ;
        assert( j>=0 && j<num_columns_ ) ;
        //assert( i*row_stride_ + columns_(j) < columns_(num_columns_) ) ;
        return ptr_[ i*row_stride_ + columns_(j) ] ;
      }

      dense_unfolding& operator=( dense_unfolding const& that ) {
        assert( that.num_rows()==num_rows() ) ;
        assert( that.num_columns()==num_columns() ) ;
        return assign( *this, that ) ;
      }

      template <typename E>
      dense_unfolding operator=( E const& that ) {
        assert( that.num_rows()==num_rows() ) ;
        assert( that.num_columns()==num_columns() ) ;
        return assign( *this, that ) ;
      }

      template <typename I1, typename I2>
      typename std::enable_if< !(std::is_integral<I1>::value && std::is_integral<I2>::value)
                             , typename matrix_selection< dense_unfolding, I1, I2 >::result_type
                             >::type operator()( I1 const& s1, I2 const& s2 ) {
        return matrix_selection< dense_unfolding, I1, I2 >::apply( *this, s1, s2 ) ;
      }

      template <typename I1, typename I2>
      typename std::enable_if< !(std::is_integral<I1>::value && std::is_integral<I2>::value)
                             , typename matrix_selection< dense_unfolding, I1, I2 >::result_type
                             >::type operator()( I1 const& s1, I2 const& s2 ) const {
        return matrix_selection< dense_unfolding, I1, I2 >::apply( *this, s1, s2 ) ;
      }

    private:
      value_type*                 ptr_ ;
      size_type                   num_rows_ ;
      size_type                   num_columns_ ;
      size_type                   row_stride_ ;
      shared_vector< size_type >  columns_ ;
  } ;

  template <typename T, typename S>
  struct concept< dense_unfolding<T,S> >
  : TensorUnfoldingMatrix
  {};


  template <typename T, typename S>
  struct matrix_selection< dense_unfolding<T,S>, all, all > {
    typedef dense_unfolding<T,S> result_type ;

    static result_type apply( dense_unfolding<T,S> m, all, all ) {
      return m ;
    }
  } ;


  template <typename T, typename S, typename R, typename I>
  struct matrix_selection< dense_unfolding<T,S>, R, I,
        typename std::enable_if< is<DenseVector,R>::value && std::is_integral<I>::value >::type > {
    typedef typename vector_selection< strided_vector<T,S>, R>::result_type result_type ;

    static result_type apply( dense_unfolding<T,S> m, R const& r, I col ) {
      assert( col<m.num_columns() && col>=0 ) ;
      return strided_vector<T,S>( m.ptr()+m.columns(col), m.stride_rows(), m.num_rows() )( r ) ;
    }
  } ;

/*
  template <typename T, typename S, typename O, typename I, typename R>
  struct matrix_selection< dense_unfolding<T,S,O>, I, R,
        typename std::enable_if< is<DenseVector,R>::value && std::is_integral<I>::value >::type > {
    typedef typename vector_selection< indirect_vector<T,S>, R>::result_type result_type ;

    static result_type apply( dense_unfolding<T,S,O> m, I row, R const& r ) {
      assert( row<m.num_rows() && row>=0 ) ;
      return indirect_vector<T,S>( &m(row,0), m.stride_rows(), m.num_columns() )( r ) ;
    }
  } ;
*/

  template <typename T, typename S, typename R2>
  struct matrix_selection< dense_unfolding<T,S>, range, R2 > {
    typedef dense_unfolding<T,S> result_type ;

    static result_type apply( dense_unfolding<T,S> m, range const& r, R2 const& r2 ) {
      assert( r.begin()+r.size()<=m.num_rows() ) ;
      return result_type( m.ptr()+r.begin(), r.size(), m.stride_rows(), m.columns()(r2) ) ;
    }
  } ;


  template <typename T, typename S>
  struct matrix_selection< dense_unfolding<T,S>, all, range > {
    typedef dense_unfolding<T,S> result_type ;

    static result_type apply( dense_unfolding<T,S> m, all, range const& r ) {
      assert( r.begin()+r.size()<=m.num_columns() ) ;
      return result_type( m.ptr(), m.num_rows(), m.stride_rows(), m.columns()(r) ) ;
    }
  } ;

} // namespace glas2


#endif
