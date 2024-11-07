//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_matrix_type_double_strided_matrix_hpp
#define glas2_matrix_type_double_strided_matrix_hpp

#include <glas2/matrix/algorithm/assign.hpp>
#include <glas2/matrix/algorithm/ops_assign.hpp>
#include <glas2/matrix/algorithm/select.hpp>
#include <glas2/vector/type/range.hpp>
#include <glas2/vector/type/slice.hpp>
#include <glas2/vector/type/all.hpp>
#include <glas2/vector/type/strided_vector.hpp>
#include <glas2/matrix/view/matrix_selection.hpp>
#include <glas2/matrix/concept/double_strided_dense_matrix.hpp>
#include <glas2/concept/concept.hpp>
#include <glas2/concept/is.hpp>
#include <type_traits>
#include <cassert>

namespace glas2 {

  template <typename T, typename S, typename O>
  class double_strided_matrix {
    public:
      typedef T              value_type ;
      typedef S              size_type ;
      typedef O              orientation ;
      typedef general_matrix structure ;

    public:
      double_strided_matrix()
      : ptr_(0)
      , stride_rows_(0)
      , stride_columns_(0)
      , num_rows_(0)
      , num_columns_(0)
      {}

      double_strided_matrix( T* ptr, S sr, S sc, S m, S n )
      : ptr_(ptr)
      , stride_rows_(sr)
      , stride_columns_(sc)
      , num_rows_( m )
      , num_columns_( n )
      {}

      // Copy reference !!
      double_strided_matrix( double_strided_matrix const& that )
      : ptr_(that.ptr_)
      , stride_rows_(that.stride_rows_)
      , stride_columns_(that.stride_columns_)
      , num_rows_( that.num_rows_ )
      , num_columns_( that.num_columns_ )
      {}

    public:
      void reset( T* ptr, S sr, S sc, S m, S n ) {
        ptr_ = ptr ;
        stride_rows_ = sr ;
        stride_columns_ = sc ;
        num_rows_ = m ;
        num_columns_ = n ;
      }

    public:
      value_type* ptr() { return ptr_ ; }
      value_type const* const ptr() const { return ptr_ ; }
      size_type stride_rows() const{ return stride_rows_ ;}
      size_type stride_columns() const{ return stride_columns_ ;}

      // Matrix
      size_type num_rows() const { return num_rows_ ; }
      size_type num_columns() const { return num_columns_ ; }

      value_type const& operator() ( size_type i, size_type j ) const {
        assert( i>=0 && i<num_rows_ ) ;
        assert( j>=0 && j<num_columns_ ) ;
        return ptr_[ i*stride_rows_ + j*stride_columns_ ] ;
      }

      value_type& operator() ( size_type i, size_type j ) {
        assert( i>=0 && i<num_rows_ ) ;
        assert( j>=0 && j<num_columns_ ) ;
        return ptr_[ i*stride_rows_ + j*stride_columns_ ] ;
      }

      double_strided_matrix& operator=( double_strided_matrix const& that ) {
        assert( that.num_rows()==num_rows() ) ;
        assert( that.num_columns()==num_columns() ) ;
        return assign( *this, that ) ;
      }

    public:
      template <typename E>
      double_strided_matrix operator=( E const& that ) {
        assert( that.num_rows()==num_rows() ) ;
        assert( that.num_columns()==num_columns() ) ;
        return assign( *this, that ) ;
      }

      template <typename E>
      double_strided_matrix& operator+=( E const& that ) {
        plus_assign( current_backend(), *this, that ) ;
        return *this ;
      }

      template <typename E>
      double_strided_matrix& operator-=( E const& that ) {
        minus_assign( current_backend(), *this, that ) ;
        return *this ;
      }

      template <typename E>
      double_strided_matrix& operator*=( E const& that ) {
        multiplies_assign( current_backend(), *this, that ) ;
        return *this ;
      }

      template <typename E>
      double_strided_matrix& operator/=( E const& that ) {
        divides_assign( current_backend(), *this, that ) ;
        return *this ;
      }

    public:
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
      T* ptr_ ;
      S  stride_rows_ ;
      S  stride_columns_ ;
      S  num_rows_ ;
      S  num_columns_ ;
  } ;

  template <typename T, typename S, typename O>
  struct glas_concept< double_strided_matrix<T,S,O> >
  : DoubleStridedDenseMatrix
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
      return strided_vector<T,S>( m.ptr()+m.stride_columns()*col, m.num_rows(), m.stride_rows() )( r ) ;
    }
  } ;


  template <typename T, typename S, typename O, typename I, typename R>
  struct matrix_selection< double_strided_matrix<T,S,O>, I, R,
        typename std::enable_if< is<DenseVector,R>::value && std::is_integral<I>::value >::type > {
    typedef typename vector_selection< strided_vector<T,S>, R>::result_type result_type ;

    static result_type apply( double_strided_matrix<T,S,O> m, I row, R const& r ) {
      assert( row<m.num_rows() && row>=0 ) ;
      return strided_vector<T,S>( m.ptr()+m.stride_rows()*row, m.num_columns(), m.stride_columns() )( r ) ;
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
