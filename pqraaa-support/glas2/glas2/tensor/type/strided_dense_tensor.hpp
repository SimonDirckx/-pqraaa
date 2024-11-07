//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS2 Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_tensor_type_strided_dense_tensor_hpp
#define glas2_tensor_type_strided_dense_tensor_hpp

#include <glas2/tensor/type/dense_unfolding.hpp>
#include <glas2/tensor/algorithm/assign.hpp>
#include <glas2/tensor/algorithm/fiber.hpp>
#include <glas2/tensor/algorithm/tensor_slice.hpp>
#include <glas2/tensor/algorithm/unfolding.hpp>
#include <glas2/tensor/view/tensor_selection.hpp>
#include <glas2/tensor/concept/strided_dense_tensor.hpp>
#include <glas2/vector/algorithm/is_equal.hpp>
#include <glas2/matrix/type/double_strided_matrix.hpp>
#include <glas2/vector/container/shared_vector.hpp>
#include <glas2/vector/type/strided_vector.hpp>
#include <glas2/vector/concept/dense_vector.hpp>
#include <glas2/concept/is.hpp>
#include <type_traits>
#include <cassert>

namespace glas2 {

  template <class T, class Shape>
  class strided_dense_tensor
  {
  public:
    typedef T                                value_type ;

    typedef Shape                            shape_type ;
    typedef typename shape_type::value_type  size_type ;

  public:
     explicit strided_dense_tensor( value_type* data, shape_type const& shape, shape_type const& stride )
     : data_( data )
     , shape_( shape )
     , stride_( stride )
     {}

     // Copy constructor
     strided_dense_tensor( strided_dense_tensor const& that )
     : data_( that.data_ )
     , shape_( that.shape_ )
     , stride_( that.stride_ )
     {}

  public:
     strided_dense_tensor& operator=( strided_dense_tensor const& that ) {
       assert( is_equal( that.shape(), shape_ ) ) ;
       assign( *this, that ) ;
       return *this ;
     }

     template <typename E>
     strided_dense_tensor& operator=( E const& that ) {
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

     shape_type const& stride() const {
       return stride_ ;
     }

     value_type* ptr() const {
       return data_ ;
     }

  private:
     template <typename I>
     size_type index( I const& ind ) const {
       assert( ind.size()==order() ) ;
       size_type pos = 0 ;
       for (int i=0; i<ind.size(); ++i) {
         pos += ind(i) * stride_(i) ;
       }
       return pos ;
     }

  public:
     void reset( value_type* data, shape_type const& shape, shape_type const& stride ) {
       assert( order()==shape.size() ) ;
       assert( order()==stride.size() ) ;
       data_ = data ;
       shape_ = shape ;
       stride_ = stride ;
     }

  public:
     value_type& operator[]( size_type i ) { return data_[i] ; }
     value_type const& operator[]( size_type i ) const { return data_[i] ; }

     template <typename I>
     typename std::enable_if< is<DenseVector,I>::value, value_type& >::type operator() ( I const& ind ) {
       return data_[ index(ind) ] ;
     }

     template <typename I>
     typename std::enable_if< is<DenseVector,I>::value, value_type const& >::type operator() ( I const& ind ) const {
       return data_[ index(ind) ] ;
     }

  public:
     template <typename S, typename Mode>
     typename tensor_selection< strided_dense_tensor, S, Mode
                              >::result_type operator()( S const& s, Mode const& m ) const {
       return tensor_selection< strided_dense_tensor, S, Mode >::apply( *this, s, m ) ;
     }


  private:
     value_type* data_ ;
     shape_type  shape_ ;
     shape_type  stride_ ;
  } ; // strided_dense_tensor


  template <class T, class Shape>
  struct concept< strided_dense_tensor<T,Shape> >
  : StridedDenseTensor
  {} ;

  template <class I, class S, class T, class Shape>
  struct tensor_transform< tensor_slice_tag<I,S>, strided_dense_tensor<T,Shape> > {
    typedef typename strided_dense_tensor<T,Shape>::size_type size_type ;
    typedef double_strided_matrix<T, size_type, column_major>                  result_type ;

    static result_type apply( strided_dense_tensor<T,Shape>& x, I const& index, S coordinate1, S coordinate2 ) {
      assert( coordinate1>=0 && coordinate1<x.order() ) ;
      assert( coordinate2>=0 && coordinate2<x.order() ) ;
      assert( coordinate1<=coordinate2 ) ;
      assert( index(coordinate1)==0 ) ;
      assert( index(coordinate2)==0 ) ;

      size_type stride1 = x.stride()(coordinate1) ;
      size_type stride2 = x.stride()(coordinate2) ;

      size_type shift = 0 ;
      for (size_type i=0; i<x.order(); ++i) shift += index(i) * x.stride()(i) ;

      return result_type( x.ptr()+shift, stride1, stride2, x.shape()(coordinate1), x.shape()(coordinate2) ) ;
    }
  } ;

  template <class I, class S, class T, class Shape>
  struct tensor_transform< fiber_tag<I,S>, strided_dense_tensor<T,Shape> > {
    typedef typename strided_dense_tensor<T,Shape>::size_type size_type ;
    typedef strided_vector<T, size_type>                         result_type ;

    static result_type apply( strided_dense_tensor<T,Shape>& x, I const& index, S coordinate ) {
      assert( coordinate>=0 && coordinate<x.order() ) ;
      assert( index(coordinate)==0 ) ;

      size_type stride = x.stride()(coordinate) ;

      size_type shift = 0 ;
      for (size_type i=0; i<x.order(); ++i) shift += index(i) * x.stride()(i) ;

      return result_type( x.ptr()+shift, x.shape()(coordinate), stride ) ;
    }
  } ;

  template <class S, class T, class Shape>
  struct tensor_transform< unfolding_tag<S>, strided_dense_tensor<T,Shape> > {
    typedef typename strided_dense_tensor<T,Shape>::size_type size_type ;
    typedef dense_unfolding<T, size_type>                        result_type ;

    static result_type apply( strided_dense_tensor<T,Shape>& x, S coordinate ) {
      return result_type( x.ptr(), x.shape(), x.stride(), coordinate ) ;
    }
  } ;

  template <class T, class Shape, class Mode>
  struct tensor_selection< strided_dense_tensor<T,Shape>, range, Mode > {
    typedef shared_vector<typename Shape::value_type> new_shape_type ;
    typedef strided_dense_tensor<T, new_shape_type>   result_type ;

    static result_type apply( strided_dense_tensor<T,Shape> x, range r, Mode m ) {
      assert( m < x.order() && m>=0 ) ;
      new_shape_type shape( x.order() ) ; shape = x.shape() ; shape(m) = r.size() ;
      new_shape_type stride( x.order() ) ; stride = x.stride() ;
      return result_type( x.ptr()+stride(m)*r.begin(), shape, stride ) ;
    }
  } ;

} // Namespace glas
  

#endif
