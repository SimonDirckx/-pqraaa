//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas_toolbox_tensor_container_dense_tensor_hpp
#define glas_toolbox_tensor_container_dense_tensor_hpp

#include <glas/config.hpp>
#define GLAS_MY_INCLUDE_LEVEL GLAS_INCLUDE_CONTAINER
#include <glas/concept/check_include_level.hpp>

#include <glas/concept/round_bracketed.hpp>
#include <glas/toolbox/tensor/concept/shape.hpp>
#include <glas/toolbox/tensor/concept/storage_size.hpp>
#include <glas/toolbox/tensor/concept/order.hpp>
#include <glas/toolbox/tensor/view/unfolding_view.hpp>
#include <glas/algorithm/product.hpp>
#include <glas/concept/square_bracketed.hpp>
#include <glas/concept/pointer.hpp>
#include <glas/concept/const_pointer.hpp>
#include <glas/concept/storage_ptr.hpp>
#include <cassert>
#ifdef GLAS_WARN_CONTAINER_COPY
#ifndef NDEBUG
#include <iostream>
#include <vector>
#include <string>
#include <limits>
#endif
#endif

#undef GLAS_INCLUDE_LEVEL
#define GLAS_INCLUDE_LEVEL GLAS_INCLUDE_CONTAINER


namespace glas {

  template <class T, class Shape>
  class dense_tensor
  {
  public:
    typedef T                                        value_type ;
    typedef T const&                                 const_reference ;
    typedef T&                                       reference ;
    typedef value_type*                              pointer ;
    typedef value_type const*                        const_pointer ;

    typedef typename glas::value_type< Shape >::type size_type ;
    typedef Shape                                    shape_type ;

    typedef dense_tensor< T, Shape > this_type ;

  public:
     explicit dense_tensor()
     : data_( 0 )
     {}

     explicit dense_tensor( shape_type const& shape )
     : shape_( shape )
     , data_( new value_type[ product(shape_) ] )
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
     : shape_( shape( that ) )
     , data_( new value_type[ storage_size(that) ] )
     {
       select_backend::assign( *this, that ) ;
     }

     // Copy constructor
     dense_tensor( dense_tensor const& that )
     : shape_( that.shape_ )
     , data_( new value_type[ glas::product(that.shape_) ] )
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
       return shape_.size() ;
     }

     shape_type const& shape() const {
       return shape_ ;
     }

     pointer storage_ptr() {
       return data_ ;
     }

     const_pointer storage_ptr() const {
       return data_ ;
     }

  public:
     void resize( shape_type const& shape ) {
       delete [] data_ ;
       shape_ = shape ;
       data_ = new value_type[ storage_size(*this) ] ;
     }

  public:
     template <typename R, typename C>
     unfolding_view< this_type, R, C > operator() ( R const& rows, C const columns ) const {
       return unfolding_view< this_type const, R, C >( *this, rows, columns ) ;
     }

     template <typename R, typename C>
     unfolding_view< this_type, R, C > operator() ( R const& rows, C const columns ) {
       return unfolding_view< this_type, R, C >( *this, rows, columns ) ;
     }

     friend std::ostream& operator<<( std::ostream& s, dense_tensor const& that ) {
       s << std::string("[") << that.shape() << std::string("]") ;
       return s ;
     }

  private:
     shape_type shape_ ;
     pointer    data_ ;
  } ; // dense_tensor


  template <class T, class Shape>
  struct RoundBracketed< dense_tensor<T,Shape> >
  : boost::mpl::true_
  {} ;

  template <class T, class Shape>
  struct DenseTensorCollection< dense_tensor<T,Shape> >
  : boost::mpl::true_
  {} ;

  template <class T, class Shape>
  struct const_pointer< dense_tensor<T,Shape> > {
    typedef T const* type ;
  } ;

} // Namespace glas
  

#endif
