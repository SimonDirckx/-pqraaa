//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas_toolbox_tensor_container_tucker_tensor_hpp
#define glas_toolbox_tensor_container_tucker_tensor_hpp

#include <glas/config.hpp>
#define GLAS_MY_INCLUDE_LEVEL GLAS_INCLUDE_CONTAINER
#include <glas/concept/check_include_level.hpp>

#include <glas/concept/round_bracketed.hpp>
#include <glas/toolbox/tensor/container/dense_tensor.hpp>
#include <glas/toolbox/tensor/concept/sizes.hpp>
#include <glas/toolbox/tensor/concept/order.hpp>
#include <glas/container/dense_matrix.hpp>
#include <glas/view/dense_tensor_unfold_view.hpp>
#include <glas/algorithm/product.hpp>
#include <glas/concept/square_bracketed.hpp>
#include <cassert>
#ifdef GLAS_WARN_CONTAINER_COPY
#ifndef NDEBUG
#include <iostream>
#include <string>
#include <limits>
#endif
#endif

#undef GLAS_INCLUDE_LEVEL
#define GLAS_INCLUDE_LEVEL GLAS_INCLUDE_CONTAINER


namespace glas {

  template <int Order, class T>
  class tucker_tensor
  {
  public:
    typedef T                                     value_type ;
    typedef T const&                              const_reference ;
    typedef T&                                    reference ;
    typedef std::ptrdiff_t                        num_rows_type ;
    typedef std::ptrdiff_t                        num_columns_type ;
    typedef value_type*                           pointer ;
    typedef value_type const*                     const_pointer ;

    typedef fixed_size_dense_vector<int>          sizes_type ;

    typedef dense_tensor< Order, T >              this_type ;

  public:
     explicit dense_tensor()
     : data_( 0 )
     {
       sizes_ = 0 ;
     }

     explicit dense_tensor( sizes_type const& sizes )
     : sizes_( sizes )
     , data_( new value_type[ product(sizes_) ] )
     {
       assert( data_ != 0 ) ;
#ifndef NDEBUG
       if (std::numeric_limits<T>::has_quiet_NaN)
         std::fill( this->storage_ptr(), this->storage_ptr()+storage_size(*this), std::numeric_limits<T>::quiet_NaN() ) ;
#endif
     }

     ~dense_tensor() {
        delete[] data_ ;
     }

     template <typename E>
     explicit dense_tensor( E const& that )
     : sizes_( sizes( that ) )
     , data_( new value_type[ storage_size(that) ] )
     {
       select_backend::assign( *this, that ) ;
     }

     // Copy constructor
     dense_tensor( dense_tensor const& that )
     : num_rows_( that.num_rows_ )
     , num_columns_( that.num_columns_ )
     , data_( new value_type[ num_rows_ * num_columns_ ] )
     {
#ifdef GLAS_WARN_CONTAINER_COPY
#ifndef NDEBUG
       std::cerr << std::string("Copy constructor of dense_tensor was called") << std::endl ;
#endif
#endif
       select_backend::assign( *this, that ) ;
     }

  public:
     dense_tensor& operator=( dense_tensor const& that ) {
       select_backend::assign( *this, that ) ;
       return *this ;
     }

     template <typename E>
     dense_tensor& operator=( E const& that ) {
       select_backend::assign( *this, that ) ;
       return *this ;
     }

  public:
     int const& order() const {
       return Order ;
     }

     sizes_type const& sizes() const {
       return sizes_ ;
     }

     pointer storage_ptr() {
       return data_ ;
     }

     const_pointer storage_ptr() const {
       return data_ ;
     }

  public:
     void resize( sizes_type const& sizes ) {
       delete [] data_ ;
       sizes_ = sizes ;
       data_ = new value_type[ storage_size(*this) ] ;
     }

  public:
     template <typename R, typename C>
     dense_tensor_unfold_view< this_type, R, C > operator() ( R const& rows, C const columns ) const {
       return dense_tensor_unfold_view< this_type const, R, C >( *this, rows, columns ) ;
     }

     template <typename R, typename C>
     dense_tensor_unfold_view< this_type, R, C > operator() ( R const& rows, C const columns ) {
       return dense_tensor_unfold_view< this_type, R, C >( *this, rows, columns ) ;
     }

  private:
     num_rows_type num_rows_ ;
     num_columns_type num_columns_ ;
     pointer     data_ ;
  } ; // dense_tensor


  template <int Order, class T>
  struct SquareBracketed< dense_tensor<Order,T> >
  : boost::mpl::true_
  {} ;

  template <int Order, class T>
  struct DenseTensorCollection< dense_tensor<Order,T> >
  : boost::mpl::true_
  {} ;

  template <int Order, class T>
  struct const_pointer< dense_tensor<Order,T> > {
    typedef T const* type ;
  } ;

} // Namespace glas
  

#endif
