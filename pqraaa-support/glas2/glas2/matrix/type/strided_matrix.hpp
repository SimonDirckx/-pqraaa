//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_matrix_type_strided_matrix_hpp
#define glas2_matrix_type_strided_matrix_hpp

#include <glas2/matrix/algorithm/assign.hpp>
#include <glas2/matrix/algorithm/select.hpp>
#include <glas2/vector/type/range.hpp>
#include <glas2/vector/type/range_from_end.hpp>
#include <glas2/vector/type/slice.hpp>
#include <glas2/vector/type/all.hpp>
#include <glas2/vector/type/contiguous_vector.hpp>
#include <glas2/vector/type/strided_vector.hpp>
#include <glas2/matrix/type/indirect_matrix.hpp>
#include <glas2/matrix/view/matrix_selection.hpp>
#include <glas2/matrix/concept/strided_dense_matrix.hpp>
#include <glas2/concept/concept.hpp>
#include <glas2/concept/is.hpp>
#include <type_traits>
#include <cassert>

namespace glas2 {

  namespace detail_strided_matrix {
    template <typename S>
    S address( row_major, S i, S j, S s ) { return i*s+j ; }

    template <typename S>
    S address( column_major, S i, S j, S s ) { return j*s+i ; }

  } // namespace detail_strided_matrix

  template <typename T, typename S, typename O>
  class strided_matrix {
    public:
      typedef T              value_type ;
      typedef S              size_type ;
      typedef O              orientation ;
      typedef general_matrix structure ;

    public:
      strided_matrix()
      : ptr_(0)
      , stride_(0)
      , num_rows_(0)
      , num_columns_(0)
      {}

      strided_matrix( T* ptr, S s, S m, S n )
      : ptr_(ptr)
      , stride_(s)
      , num_rows_( m )
      , num_columns_( n )
      {}

      // Copy reference !!
      strided_matrix( strided_matrix const& that )
      : ptr_(that.ptr_)
      , stride_(that.stride_)
      , num_rows_( that.num_rows_ )
      , num_columns_( that.num_columns_ )
      {}

    public:
      void reset( T* ptr, S s, S m, S n ) {
        ptr_ = ptr ;
        stride_ = s ;
        num_rows_ = m ;
        num_columns_ = n ;
      }

    public:
      value_type* ptr() { return ptr_ ; }
      value_type const* const ptr() const { return ptr_ ; }
      size_type stride() const{ return stride_ ;}

      // Matrix
      size_type num_rows() const { return num_rows_ ; }
      size_type num_columns() const { return num_columns_ ; }

      value_type const* address ( size_type i, size_type j ) const {
        assert( i>=0 && i<=num_rows_ ) ;
        assert( j>=0 && j<=num_columns_ ) ;
        return ptr_+detail_strided_matrix::address(orientation(), i, j, stride_) ;
      }

      value_type* address ( size_type i, size_type j ) {
        assert( i>=0 && i<=num_rows_ ) ;
        assert( j>=0 && j<=num_columns_ ) ;
        return ptr_+detail_strided_matrix::address(orientation(), i, j, stride_) ;
      }

      value_type const& operator() ( size_type i, size_type j ) const {
        assert( i>=0 && i<num_rows_ ) ;
        assert( j>=0 && j<num_columns_ ) ;
        return ptr_[ detail_strided_matrix::address( orientation(), i, j, stride_ ) ] ;
      }

      value_type& operator() ( size_type i, size_type j ) {
        assert( i>=0 && i<num_rows_ ) ;
        assert( j>=0 && j<num_columns_ ) ;
        return ptr_[ detail_strided_matrix::address( orientation(), i, j, stride_ ) ] ;
      }

      strided_matrix operator=( strided_matrix const& that ) {
        assert( that.num_rows()==num_rows() ) ;
        assert( that.num_columns()==num_columns() ) ;
        return assign( *this, that ) ;
      }

    public:
      template <typename E>
      strided_matrix operator=( E const& that ) {
        assert( that.num_rows()==num_rows() ) ;
        assert( that.num_columns()==num_columns() ) ;
        return assign( *this, that ) ;
      }

/*      template <typename E>
      strided_matrix& operator+=( E const& that ) {
        plus_assign( current_backend(), *this, that ) ;
        return *this ;
      }

      template <typename E>
      strided_matrix& operator-=( E const& that ) {
        minus_assign( current_backend(), *this, that ) ;
        return *this ;
      }*/

      template <typename E>
      strided_matrix& operator*=( E const& that ) {
        multiplies_assign( current_backend(), *this, that ) ;
        return *this ;
      }

      template <typename E>
      strided_matrix& operator/=( E const& that ) {
        divides_assign( current_backend(), *this, that ) ;
        return *this ;
      }

    public:
      template <typename I1, typename I2>
      typename std::enable_if< !(std::is_integral<I1>::value && std::is_integral<I2>::value)
                             , typename matrix_selection< strided_matrix, I1, I2 >::result_type
                             >::type operator()( I1 const& s1, I2 const& s2 ) {
        return matrix_selection< strided_matrix, I1, I2 >::apply( *this, s1, s2 ) ;
      }

      template <typename I1, typename I2>
      typename std::enable_if< !(std::is_integral<I1>::value && std::is_integral<I2>::value)
                             , typename matrix_selection< strided_matrix, I1, I2 >::result_type
                             >::type operator()( I1 const& s1, I2 const& s2 ) const {
        return matrix_selection< strided_matrix, I1, I2 >::apply( *this, s1, s2 ) ;
      }

    private:
      T* ptr_ ;
      S  stride_ ;
      S  num_rows_ ;
      S  num_columns_ ;
  } ;

  template <typename T, typename S, typename O>
  struct glas_concept< strided_matrix<T,S,O> >
  : StridedDenseMatrix
  {};


  template <typename T, typename S, typename O>
  struct matrix_selection< strided_matrix<T,S,O>, all, all > {
    typedef strided_matrix<T,S,O> result_type ;

    static result_type apply( strided_matrix<T,S,O> m, all, all ) {
      return m ;
    }
  } ;


  template <typename T, typename S, typename R, typename I>
  struct matrix_selection< strided_matrix<T,S,column_major>, R, I,
        typename std::enable_if< is<DenseVector,R>::value && std::is_integral<I>::value >::type > {
    typedef typename vector_selection< contiguous_vector<T,S>, R>::result_type result_type ;

    static result_type apply( strided_matrix<T,S,column_major> m, R const& r, I col ) {
      assert( col<m.num_columns() && col>=0 ) ;
      return contiguous_vector<T,S>( m.address(0,col), m.num_rows() )( r ) ;
    }
  } ;


  template <typename T, typename S, typename I, typename R>
  struct matrix_selection< strided_matrix<T,S,column_major>, I, R,
        typename std::enable_if< is<DenseVector,R>::value && std::is_integral<I>::value >::type > {
    typedef typename vector_selection< strided_vector<T,S>, R>::result_type result_type ;

    static result_type apply( strided_matrix<T,S,column_major> m, I row, R const& r ) {
      assert( row<m.num_rows() && row>=0 ) ;
      return strided_vector<T,S>( m.address(row,0), m.num_columns(), m.stride() )( r ) ;
    }
  } ;


  template <typename T, typename S, typename R, typename I>
  struct matrix_selection< strided_matrix<T,S,row_major>, R, I,
        typename std::enable_if< is<DenseVector,R>::value && std::is_integral<I>::value >::type > {
    typedef typename vector_selection< strided_vector<T,S>, R>::result_type result_type ;

    static result_type apply( strided_matrix<T,S,row_major> m, R const& r, I col ) {
      assert( col<m.num_columns() && col>=0 ) ;
      return strided_vector<T,S>( m.address(0,col), m.num_rows(), m.stride() )( r ) ;
    }
  } ;


  template <typename T, typename S, typename I, typename R>
  struct matrix_selection< strided_matrix<T,S,row_major>, I, R,
        typename std::enable_if< is<DenseVector,R>::value && std::is_integral<I>::value >::type > {
    typedef typename vector_selection< contiguous_vector<T,S>, R>::result_type result_type ;

    static result_type apply( strided_matrix<T,S,row_major> m, I row, R const& r ) {
      assert( row<m.num_rows() && row>=0 ) ;
      return contiguous_vector<T,S>( m.address(row,0), m.num_columns() )( r ) ;
    }
  } ;


  template <typename T, typename S, typename R2>
  struct matrix_selection< strided_matrix<T,S,column_major>, range, R2,
        typename std::enable_if< is<DenseVector,R2>::value >::type > {
    typedef strided_matrix<T,S,column_major>                             temp_type ;
    typedef typename matrix_selection< temp_type, all, R2 >::result_type result_type ;

    static result_type apply( strided_matrix<T,S,column_major> m, range const& r, R2 const& r2 ) {
      assert( r.begin()+r.size()<=m.num_rows() ) ;
      return temp_type( &m(r.begin(),0), m.stride(), r.size(), m.num_columns() )( all(), r2 ) ;
    }
  } ;

  template <typename T, typename S, typename R2>
  struct matrix_selection< strided_matrix<T,S,column_major>, range_from_end, R2
                         , typename std::enable_if< is<DenseVector,R2>::value
                         >::type > {
    typedef strided_matrix<T,S,column_major>                            temp_type ;
    typedef typename matrix_selection< temp_type, all, R2>::result_type result_type ;

    static result_type apply( strided_matrix<T,S,column_major> m, range_from_end const& r, R2 const& r2 ) {
      assert( r.begin()<=m.num_rows()-r.from_end() ) ;
      return temp_type( &m(r.begin(),0), m.stride(), m.num_rows()-r.begin()-r.from_end(), m.num_columns() )( all(), r2 ) ;
    }
  } ;



  template <typename T, typename S, typename R1>
  struct matrix_selection< strided_matrix<T,S,row_major>, R1, range
                         , typename std::enable_if< is<DenseVector,R1>::value >::type
                         > {
    typedef strided_matrix<T,S,row_major>                                temp_type ;
    typedef typename matrix_selection< temp_type, R1, all >::result_type result_type ;

    static result_type apply( strided_matrix<T,S,row_major> m, R1 const& r1, range const& r ) {
      assert( r.begin()+r.size()<=m.num_columns() ) ;
      return temp_type( m.ptr()+r.begin()*m.stride(), m.stride(), m.num_rows(), r.size() )( r1, all() ) ;
    }
  } ;

  template <typename T, typename S, typename R1>
  struct matrix_selection< strided_matrix<T,S,row_major>, R1, range_from_end
                         , typename std::enable_if< is<DenseVector,R1>::value
                         >::type > {
    typedef strided_matrix<T,S,row_major>                            temp_type ;
    typedef typename matrix_selection< temp_type, R1, all>::result_type result_type ;

    static result_type apply( strided_matrix<T,S,row_major> m, R1 const& r1, range_from_end const& r ) {
      assert( r.begin()<=m.num_columns()-r.from_end() ) ;
      return temp_type( m.ptr()+r.begin()*m.stride(), m.stride(), m.num_rows(), m.num_columns()-r.begin()-r.from_end() )( r1, all() ) ;
    }
  } ;



  template <typename T, typename S>
  struct matrix_selection< strided_matrix<T,S,column_major>, all, range > {
    typedef strided_matrix<T,S,column_major> result_type ;

    static result_type apply( strided_matrix<T,S,column_major> m, all, range const& r ) {
      assert( r.begin()+r.size()<=m.num_columns() ) ;
      return result_type( m.ptr()+r.begin()*m.stride(), m.stride(), m.num_rows(), r.size() ) ;
    }
  } ;

  template <typename T, typename S>
  struct matrix_selection< strided_matrix<T,S,column_major>, all, range_from_end > {
    typedef strided_matrix<T,S,column_major> result_type ;

    static result_type apply( strided_matrix<T,S,column_major> m, all, range_from_end const& r ) {
      assert( r.begin()<=m.num_columns()-r.from_end() ) ;
      return result_type( m.ptr()+r.begin()*m.stride(), m.stride(), m.num_rows(), m.num_columns()-r.begin()-r.from_end() ) ;
    }
  } ;


  template <typename T, typename S>
  struct matrix_selection< strided_matrix<T,S,row_major>, range, all > {
    typedef strided_matrix<T,S,row_major> result_type ;

    static result_type apply( strided_matrix<T,S,row_major> m, range const& r, all ) {
      assert( r.begin()+r.size()<=m.num_rows() ) ;
      return result_type( m.ptr()+r.begin()*m.stride(), m.stride(), r.size(), m.num_columns() ) ;
    }
  } ;

  template <typename T, typename S>
  struct matrix_selection< strided_matrix<T,S,row_major>, range_from_end, all > {
    typedef strided_matrix<T,S,row_major> result_type ;

    static result_type apply( strided_matrix<T,S,row_major> m, range_from_end const& r, all ) {
      assert( r.begin()<=m.num_rows()-r.from_end() ) ;
      return result_type( m.ptr()+r.begin()*m.stride(), m.stride(), m.num_rows()-r.begin()-r.from_end(), m.num_columns() ) ;
    }
  } ;


  template <typename T, typename S, typename O, typename V>
  struct matrix_selection< strided_matrix<T,S,O>, all, V,
        typename std::enable_if< is<DenseVector,V>::value && !( std::is_same<V,range>::value||std::is_same<V,slice>::value||std::is_same<V,range_from_end>::value||std::is_same<V,all>::value) >::type > {
    typedef indirect_matrix<strided_matrix<T,S,O>, range,V> result_type ;

    static result_type apply( strided_matrix<T,S,O> m, all, V const& v ) {
      return result_type( m, range(0,m.num_rows()), v ) ;
    }
  } ;


  template <typename T, typename S, typename O, typename V>
  struct matrix_selection< strided_matrix<T,S,O>, V, all,
        typename std::enable_if< is<DenseVector,V>::value && !( std::is_same<V,range>::value||std::is_same<V,slice>::value||std::is_same<V,range_from_end>::value||std::is_same<V,all>::value) >::type > {
    typedef indirect_matrix<strided_matrix<T,S,O>, V, range> result_type ;

    static result_type apply( strided_matrix<T,S,O> m, V const& v, all ) {
      return result_type( m, v, range(0,m.num_columns()) ) ;
    }
  } ;



/*
  template <typename T, typename S>
  struct matrix_transform< select_tag, strided_matrix<T,S,column_major> > {
    typedef strided_matrix<T,S,row_major> result_type ;
    static result_type apply( strided_matrix<T,S,column_major> m ) {
      return result_type( m.ptr(), m.num_columns(), m.num_rows() ) ;
    }
  } ;
*/
} // namespace glas2


#endif
