//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_matrix_type_contiguous_matrix_hpp
#define glas2_matrix_type_contiguous_matrix_hpp

#include <glas2/matrix/algorithm/assign.hpp>
#include <glas2/matrix/algorithm/ops_assign.hpp>
#include <glas2/vector/type/range.hpp>
#include <glas2/vector/type/range_from_end.hpp>
#include <glas2/vector/type/slice.hpp>
#include <glas2/vector/type/all.hpp>
#include <glas2/vector/type/contiguous_vector.hpp>
#include <glas2/vector/type/strided_vector.hpp>
#include <glas2/matrix/type/strided_matrix.hpp>
#include <glas2/matrix/type/indirect_matrix.hpp>
#include <glas2/matrix/view/matrix_selection.hpp>
#include <glas2/matrix/concept/contiguous_dense_matrix.hpp>
#include <glas2/vector/concept/is_no_range_or_slice.hpp>
#include <glas2/concept/concept.hpp>
#include <glas2/concept/is.hpp>
#include <type_traits>
#include <cassert>

namespace glas2 {

  namespace detail_contiguous_matrix {
    template <typename S>
    S address( row_major, S i, S j, S nr, S nc ) { return i*nc+j ; }

    template <typename S>
    S address( column_major, S i, S j, S nr, S nc ) { return j*nr+i ; }

  } // namespace detail_contiguous_matrix

  template <typename T, typename S, typename O>
  class contiguous_matrix {
    public:
      typedef T              value_type ;
      typedef S              size_type ;
      typedef O              orientation ;
      typedef general_matrix structure ;

    public:
      contiguous_matrix()
      : ptr_(0)
      , num_rows_(0)
      , num_columns_(0)
      {}

      contiguous_matrix( T* ptr, S m, S n )
      : ptr_(ptr)
      , num_rows_( m )
      , num_columns_( n )
      {}

      // Copy reference !!
      contiguous_matrix( contiguous_matrix const& that )
      : ptr_(that.ptr_)
      , num_rows_( that.num_rows_ )
      , num_columns_( that.num_columns_ )
      {}

    public:
      void reset( T* ptr, S m, S n ) {
        ptr_ = ptr ;
        num_rows_ = m ;
        num_columns_ = n ;
      }

      friend void swap( contiguous_matrix& m1, contiguous_matrix& m2 ) {
        contiguous_matrix temp( m1.ptr_, m1.num_rows_, m1.num_columns_ ) ;
        m1.reset( m2.ptr_, m2.num_rows_, m2.num_columns_ ) ;
        m2.reset( temp.ptr_, temp.num_rows_, temp.num_columns_ ) ;
      }

    public:
      value_type* __restrict__ ptr() { return ptr_ ; }
      value_type* __restrict__ const ptr() const { return ptr_ ; }

      // Matrix
      size_type num_rows() const { return num_rows_ ; }
      size_type num_columns() const { return num_columns_ ; }

      value_type const& operator() ( size_type i, size_type j ) const {
        assert( i>=0 && i<num_rows_ ) ;
        assert( j>=0 && j<num_columns_ ) ;
        return ptr_[detail_contiguous_matrix::address(orientation(), i,j,num_rows_,num_columns_)] ;
      }

      value_type& operator() ( size_type i, size_type j ) {
        assert( i>=0 && i<num_rows_ ) ;
        assert( j>=0 && j<num_columns_ ) ;
        return ptr_[detail_contiguous_matrix::address(orientation(), i,j,num_rows_,num_columns_)] ;
      }

      value_type* address ( size_type i, size_type j ) {
        assert( i>=0 && i<=num_rows_ ) ;
        assert( j>=0 && j<=num_columns_ ) ;
        return ptr_+detail_contiguous_matrix::address(orientation(), i,j,num_rows_,num_columns_) ;
      }

      value_type const* address ( size_type i, size_type j ) const {
        assert( i>=0 && i<=num_rows_ ) ;
        assert( j>=0 && j<=num_columns_ ) ;
        return ptr_+detail_contiguous_matrix::address(orientation(), i,j,num_rows_,num_columns_) ;
      }

      contiguous_matrix& operator=( contiguous_matrix const& that ) {
        assert( that.num_rows()==num_rows() ) ;
        assert( that.num_columns()==num_columns() ) ;
        std::copy( that.ptr_, that.ptr_+num_rows_*num_columns_, ptr_ ) ;
        return *this ;
      }

    public:
      template <typename E>
      contiguous_matrix& operator=( E const& that ) {
        assign( *this, that ) ;
        return *this ;
      }
/*
      template <typename E>
      contiguous_matrix& operator+=( E const& that ) {
        plus_assign( current_backend(), *this, that ) ;
        return *this ;
      }

      template <typename E>
      contiguous_matrix& operator-=( E const& that ) {
        minus_assign( current_backend(), *this, that ) ;
        return *this ;
      }*/

      template <typename E>
      contiguous_matrix& operator*=( E const& that ) {
        multiplies_assign( current_backend(), *this, that ) ;
        return *this ;
      }

      template <typename E>
      contiguous_matrix& operator/=( E const& that ) {
        divides_assign( current_backend(), *this, that ) ;
        return *this ;
      }

    public:
      template <typename I1, typename I2>
      typename std::enable_if< !(std::is_integral<I1>::value && std::is_integral<I2>::value)
                             , typename matrix_selection< contiguous_matrix, I1, I2 >::result_type
                             >::type operator()( I1 const& s1, I2 const& s2 ) {
        return matrix_selection< contiguous_matrix, I1, I2 >::apply( *this, s1, s2 ) ;
      }

      template <typename I1, typename I2>
      typename std::enable_if< !(std::is_integral<I1>::value && std::is_integral<I2>::value)
                             , typename matrix_selection< contiguous_matrix, I1, I2 >::result_type
                             >::type operator()( I1 const& s1, I2 const& s2 ) const {
        return matrix_selection< contiguous_matrix, I1, I2 >::apply( *this, s1, s2 ) ;
      }

    private:
      T* __restrict__  ptr_ ;
      S  num_rows_ ;
      S  num_columns_ ;
  } ;


  template <typename T, typename S, typename R, typename I>
  struct matrix_selection< contiguous_matrix<T,S,column_major>, R, I,
        typename std::enable_if< is<DenseVector,R>::value && std::is_integral<I>::value >::type > {
    typedef typename vector_selection< contiguous_vector<T,S>, R>::result_type result_type ;

    static result_type apply( contiguous_matrix<T,S,column_major> m, R const& r, I col ) {
      assert( col<m.num_columns() && col>=0 ) ;
      return contiguous_vector<T,S>( m.address(0,col), m.num_rows() )( r ) ;
    }
  } ;


  template <typename T, typename S, typename I, typename R>
  struct matrix_selection< contiguous_matrix<T,S,column_major>, I, R,
        typename std::enable_if< is<DenseVector,R>::value && std::is_integral<I>::value >::type > {
    typedef typename vector_selection< strided_vector<T,S>, R>::result_type result_type ;

    static result_type apply( contiguous_matrix<T,S,column_major> m, I row, R const& r ) {
      assert( row<m.num_rows() && row>=0 ) ;
      return strided_vector<T,S>( m.address(row,0), m.num_columns(), m.num_rows() )( r ) ;
    }
  } ;


  template <typename T, typename S, typename R, typename I>
  struct matrix_selection< contiguous_matrix<T,S,row_major>, R, I,
        typename std::enable_if< is<DenseVector,R>::value && std::is_integral<I>::value >::type > {
    typedef typename vector_selection< strided_vector<T,S>, R>::result_type result_type ;

    static result_type apply( contiguous_matrix<T,S,row_major> m, R const& r, I col ) {
      assert( col<m.num_columns() && col>=0 ) ;
      return strided_vector<T,S>( m.address(0,col), m.num_rows(), m.num_columns() )( r ) ;
    }
  } ;


  template <typename T, typename S, typename I, typename R>
  struct matrix_selection< contiguous_matrix<T,S,row_major>, I, R,
        typename std::enable_if< is<DenseVector,R>::value && std::is_integral<I>::value >::type > {
    typedef typename vector_selection< contiguous_vector<T,S>, R>::result_type result_type ;

    static result_type apply( contiguous_matrix<T,S,row_major> m, I row, R const& r ) {
      assert( row<m.num_rows() && row>=0 ) ;
      return contiguous_vector<T,S>( m.address(row,0), m.num_columns() )( r ) ;
    }
  } ;


  template <typename T, typename S, typename O, typename V>
  struct matrix_selection< contiguous_matrix<T,S,O>, all, V,
        typename std::enable_if< is<DenseVector,V>::value && is_no_range_or_slice<V>::value >::type > {
    typedef indirect_matrix<contiguous_matrix<T,S,O>, range,V> result_type ;

    static result_type apply( contiguous_matrix<T,S,O> m, all, V const& v ) {
      return result_type( m, range(0,m.num_rows()), v ) ;
    }
  } ;




  template <typename T, typename S, typename O, typename V>
  struct matrix_selection< contiguous_matrix<T,S,O>, V, all,
        typename std::enable_if< is<DenseVector,V>::value && is_no_range_or_slice<V>::value >::type > {
    typedef indirect_matrix<contiguous_matrix<T,S,O>, V, range> result_type ;

    static result_type apply( contiguous_matrix<T,S,O> m, V const& v, all ) {
      return result_type( m, v, range(0,m.num_columns()) ) ;
    }
  } ;


  template <typename T, typename S, typename V>
  struct matrix_selection< contiguous_matrix<T,S,column_major>, V, range,
        typename std::enable_if< is<DenseVector,V>::value && is_no_range_or_slice<V>::value >::type > {
    typedef indirect_matrix<contiguous_matrix<T,S,column_major>, V, range> result_type ;

    static result_type apply( contiguous_matrix<T,S,column_major> m, V const& v, range r ) {
      return result_type( m, v, r ) ;
    }
  } ;



  template <typename T, typename S>
  struct matrix_selection< contiguous_matrix<T,S,column_major>, all, range > {
    typedef contiguous_matrix<T,S,column_major> result_type ;

    static result_type apply( contiguous_matrix<T,S,column_major> m, all, range const& r ) {
      assert( r.begin()+r.size()<=m.num_columns() ) ;
      return result_type( m.address(0,r.begin()), m.num_rows(), r.size() ) ;
    }
  } ;


  template <typename T, typename S>
  struct matrix_selection< contiguous_matrix<T,S,row_major>, range, all > {
    typedef contiguous_matrix<T,S,row_major> result_type ;

    static result_type apply( contiguous_matrix<T,S,row_major> m, range const& r, all ) {
      assert( r.begin()+r.size()<=m.num_rows() ) ;
      return result_type( m.address(r.begin(),0), r.size(), m.num_columns() ) ;
    }
  } ;


  template <typename T, typename S>
  struct matrix_selection< contiguous_matrix<T,S,column_major>, all, range_from_end > {
    typedef contiguous_matrix<T,S,column_major> result_type ;

    static result_type apply( contiguous_matrix<T,S,column_major> m, all, range_from_end const& r ) {
      assert( r.begin()<=m.num_columns()-r.from_end() ) ;
      return result_type( m.address(0,r.begin()), m.num_rows(), m.num_columns()-r.begin()-r.from_end() ) ;
    }
  } ;


  template <typename T, typename S>
  struct matrix_selection< contiguous_matrix<T,S,row_major>, range_from_end, all > {
    typedef contiguous_matrix<T,S,row_major> result_type ;

    static result_type apply( contiguous_matrix<T,S,row_major> m, range_from_end const& r, all ) {
      assert( r.begin()<=m.num_rows()-r.from_end() ) ;
      return result_type( m.address(r.begin(),0), m.num_rows()-r.begin()-r.from_end(), m.num_columns() ) ;
    }
  } ;


/*
  template <typename T, typename S>
  struct matrix_selection< contiguous_matrix<T,S,column_major>, range, all > {
    typedef strided_matrix<T,S,column_major> result_type ;

    static result_type apply( contiguous_matrix<T,S,column_major> m, range const& r, all ) {
      assert( r.begin()+r.size()<=m.num_rows() ) ;
      return result_type( &m(r.begin(),0), m.num_rows, r.size(), m.num_columns() ) ;
    }
  } ;
*/

  template <typename T, typename S, typename R2>
  struct matrix_selection< contiguous_matrix<T,S,column_major>, range, R2
                         , typename std::enable_if< is<DenseVector,R2>::value
                         >::type > {
    typedef strided_matrix<T,S,column_major>                            temp_type ;
    typedef typename matrix_selection< temp_type, all, R2>::result_type result_type ;

    static result_type apply( contiguous_matrix<T,S,column_major> m, range const& r, R2 const& r2 ) {
      assert( r.begin()+r.size()<=m.num_rows() ) ;
      return temp_type( m.address(r.begin(),0), m.num_rows(), r.size(), m.num_columns() )( all(), r2 ) ;
    }
  } ;


  template <typename T, typename S, typename R2>
  struct matrix_selection< contiguous_matrix<T,S,column_major>, range_from_end, R2
                         , typename std::enable_if< is<DenseVector,R2>::value
                         >::type > {
    typedef strided_matrix<T,S,column_major>                            temp_type ;
    typedef typename matrix_selection< temp_type, all, R2>::result_type result_type ;

    static result_type apply( contiguous_matrix<T,S,column_major> m, range_from_end const& r, R2 const& r2 ) {
      assert( r.begin()<=m.num_rows()-r.from_end() ) ;
      return temp_type( m.address(r.begin(),0), m.num_rows(), m.num_rows()-r.begin()-r.from_end(), m.num_columns() )( all(), r2 ) ;
    }
  } ;


  template <typename T, typename S, typename R1>
  struct matrix_selection< contiguous_matrix<T,S,row_major>, R1, range
                         , typename std::enable_if< is<DenseVector,R1>::value
                         >::type > {
    typedef strided_matrix<T,S,row_major>                            temp_type ;
    typedef typename matrix_selection< temp_type, R1, all>::result_type result_type ;

    static result_type apply( contiguous_matrix<T,S,row_major> m, R1 const& r1, range const& r ) {
      assert( r.begin()+r.size()<=m.num_columns() ) ;
      return temp_type( m.address(0,r.begin()), m.num_columns(), m.num_rows(), r.size() )( r1, all() ) ;
    }
  } ;

  template <typename T, typename S, typename R1>
  struct matrix_selection< contiguous_matrix<T,S,row_major>, R1, range_from_end
                         , typename std::enable_if< is<DenseVector,R1>::value
                         >::type > {
    typedef strided_matrix<T,S,row_major>                            temp_type ;
    typedef typename matrix_selection< temp_type, R1, all>::result_type result_type ;

    static result_type apply( contiguous_matrix<T,S,row_major> m, R1 const& r1, range_from_end const& r ) {
      assert( r.begin()<=m.num_columns()-r.from_end() ) ;
      return temp_type( m.address(0,r.begin()), m.num_columns(), m.num_rows(), m.num_columns()-r.begin()-r.from_end() )( r1, all() ) ;
    }
  } ;


  template <typename T, typename S, typename O>
  struct glas_concept< contiguous_matrix<T,S,O> >
  : ContiguousDenseMatrix
  {};

/*
  template <typename T, typename S>
  struct matrix_transform< select_tag, contiguous_matrix<T,S,column_major> > {
    typedef contiguous_matrix<T,S,row_major> result_type ;
    static result_type apply( contiguous_matrix<T,S,column_major> m ) {
      return result_type( m.ptr(), m.num_columns(), m.num_rows() ) ;
    }
  } ;
*/
} // namespace glas2


#endif
