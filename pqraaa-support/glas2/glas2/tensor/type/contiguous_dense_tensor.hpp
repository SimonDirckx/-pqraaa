//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS2 Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_tensor_type_contiguous_dense_tensor_hpp
#define glas2_tensor_type_contiguous_dense_tensor_hpp

#include <glas2/tensor/type/dense_unfolding.hpp>
#include <glas2/tensor/type/strided_dense_tensor.hpp>
#include <glas2/tensor/algorithm/assign.hpp>
#include <glas2/tensor/algorithm/fiber.hpp>
#include <glas2/tensor/algorithm/tensor_slice.hpp>
#include <glas2/tensor/algorithm/unfolding.hpp>
#include <glas2/tensor/view/tensor_selection.hpp>
#include <glas2/tensor/concept/contiguous_dense_tensor.hpp>
#include <glas2/vector/algorithm/is_equal.hpp>
#include <glas2/vector/algorithm/fill.hpp>
#include <glas2/matrix/type/double_strided_matrix.hpp>
#include <glas2/vector/type/strided_vector.hpp>
#include <glas2/vector/concept/dense_vector.hpp>
#include <glas2/vector/container/shared_vector.hpp>
#include <glas2/concept/is.hpp>
#include <type_traits>
#include <cassert>

namespace glas2 {

  template <class T, class Shape>
  class contiguous_dense_tensor
  {
  public:
    typedef T                                value_type ;

    typedef Shape                            shape_type ;
    typedef typename shape_type::value_type  size_type ;

  public:
     explicit contiguous_dense_tensor( value_type* data, shape_type const& shape )
     : data_( data )
     , shape_( shape )
     {}

     // Copy constructor
     contiguous_dense_tensor( contiguous_dense_tensor const& that )
     : data_( that.data_ )
     , shape_( that.shape_ )
     {}

  public:
     contiguous_dense_tensor& operator=( contiguous_dense_tensor const& that ) {
       assert( is_equal( that.shape(), shape_ ) ) ;
       std::copy( that.data_, that.data_ + prod(shape_), data_ ) ;
       return *this ;
     }

     template <typename E>
     contiguous_dense_tensor& operator=( E const& that ) {
       assert( is_equal( that.shape(), shape_ ) ) ;
       assign( *this, that ) ;
       return *this ;
     }

  public:
     int order() const {
       return shape_.size() ;
     }

     shape_type const& shape() const {
       return shape_ ;
     }

     value_type* ptr() const {
       return data_ ;
     }

  private:
     template <typename I>
     size_type index( I const& ind ) const {
       assert( ind.size()==order() ) ;
       size_type pos = 0 ;
       for (int i=ind.size()-1 ;i>=0; --i) {
         pos = ind(i) + pos * shape_(i) ;
       }
       return pos ;
     }

  public:
     value_type& operator[]( size_type i ) { return data_[i] ; }
     value_type const& operator[]( size_type i ) const { return data_[i] ; }

     template <typename I>
     typename std::enable_if< is<DenseVector,I>::value, value_type& >::type operator() ( I const& ind ) {
       assert( this->index(ind)<glas2::prod(shape_) ) ;
       return data_[ index(ind) ] ;
     }

     template <typename I>
     typename std::enable_if< is<DenseVector,I>::value, value_type const& >::type operator() ( I const& ind ) const {
       assert( this->index(ind)<glas2::prod(shape_) ) ;
       return data_[ index(ind) ] ;
     }

  public:
     template <typename S, typename Mode>
     typename tensor_selection< contiguous_dense_tensor, S, Mode
                              >::result_type operator()( S const& s, Mode const& m ) const {
       return tensor_selection< contiguous_dense_tensor, S, Mode >::apply( *this, s, m ) ;
     }

  private:
     value_type* data_ ;
     shape_type  shape_ ;
  } ; // contiguous_dense_tensor


  template <class T, class Shape>
  struct concept< contiguous_dense_tensor<T,Shape> >
  : ContiguousDenseTensor
  {} ;

  template <class I, class S, class T, class Shape>
  struct tensor_transform< tensor_slice_tag<I,S>, contiguous_dense_tensor<T,Shape> > {
    typedef typename contiguous_dense_tensor<T,Shape>::size_type size_type ;
    typedef double_strided_matrix<T, size_type, column_major>                  result_type ;

    static result_type apply( contiguous_dense_tensor<T,Shape>& x, I const& index, S coordinate1, S coordinate2 ) {
      assert( coordinate1>=0 && coordinate1<x.order() ) ;
      assert( coordinate2>=0 && coordinate2<x.order() ) ;
      assert( coordinate1<=coordinate2 ) ;
      assert( index(coordinate1)==0 ) ;
      assert( index(coordinate2)==0 ) ;

      size_type stride1 = 1 ;
      for (size_type i=0; i<coordinate1; ++i) stride1 *= x.shape()(i) ;

      size_type stride2 = 1 ;
      for (size_type i=0; i<coordinate2; ++i) stride2 *= x.shape()(i) ;

      size_type shift = 0 ;
      for (size_type i=x.order()-1; i>=0; --i) shift = index(i) + x.shape()(i)*shift ;

      return result_type( x.ptr()+shift, stride1, stride2, x.shape()(coordinate1), x.shape()(coordinate2) ) ;
    }
  } ;

  template <class I, class S, class T, class Shape>
  struct tensor_transform< fiber_tag<I,S>, contiguous_dense_tensor<T,Shape> > {
    typedef typename contiguous_dense_tensor<T,Shape>::size_type size_type ;
    typedef strided_vector<T, size_type>                         result_type ;

    static result_type apply( contiguous_dense_tensor<T,Shape>& x, I const& index, S coordinate ) {
      assert( coordinate>=0 && coordinate<x.order() ) ;
      assert( index(coordinate)==0 ) ;

      size_type stride = 1 ;
      for (size_type i=0; i<coordinate; ++i) stride *= x.shape()(i) ;

      size_type shift = 0 ;
      for (size_type i=x.order()-1; i>=0; --i) shift = index(i) + x.shape()(i)*shift ;

      return result_type( x.ptr()+shift, x.shape()(coordinate), stride ) ;
    }
  } ;

  template <class S, class T, class Shape>
  struct tensor_transform< unfolding_tag<S>, contiguous_dense_tensor<T,Shape> > {
    typedef typename contiguous_dense_tensor<T,Shape>::size_type size_type ;
    typedef dense_unfolding<T, size_type>                        result_type ;

    static result_type apply( contiguous_dense_tensor<T,Shape>& x, S coordinate ) {
/*      size_type row_stride = 1 ;
      for (size_type i=0; i<coordinate; ++i) row_stride *= x.shape()(i) ;
      size_type nnz = glas2::prod(x.shape()) ;
      size_type num_columns = nnz / x.shape()(coordinate) ;
      shared_vector<size_type> columns( num_columns ) ;
      shared_vector<size_type> index( num_columns ) ;
      fill( index, 0 ) ;

      size_type col = 0 ;
      for (size_type i=0; i<nnz; ++i) {
        if (index(coordinate)==0) {
          columns(col) = i ; ++col ;
        }
        index(0)++ ;
        int k=0; while ( k<x.shape().size() && index(k)==x.shape()(k) ) { index(k) = 0; ++k; index(k)++ ; }
      }

      return result_type( x.ptr(), x.shape()(coordinate), row_stride, columns ) ;
      */
      shared_vector<size_type> stride( x.order() ) ;
      stride(0) = 1 ;
      for ( size_type i=1; i<stride.size(); ++i ) stride(i) = stride(i-1) * x.shape()(i-1) ;
      return result_type( x.ptr(), x.shape(), stride, coordinate ) ;
    }
  } ;

  template <class T, class Shape, class Mode>
  struct tensor_selection< contiguous_dense_tensor<T,Shape>, range, Mode > {
    typedef strided_dense_tensor<T,Shape> result_type ;

    static result_type apply( contiguous_dense_tensor<T,Shape> x, range r, Mode m ) {
      typedef Shape shape_type ;

      assert( m < x.order() && m>=0 ) ;
      shape_type shape( x.order() ) ; shape = x.shape() ; shape(m) = r.size() ;
      shape_type stride( x.order() ) ;
      stride(0) = 1 ;
      for (int i=1; i<x.order(); ++i) stride(i) = stride(i-1) * x.shape()(i-1) ;
      return result_type( x.ptr()+stride(m)*r.begin(), shape, stride ) ;
    }
  } ;

} // Namespace glas
  

#endif
