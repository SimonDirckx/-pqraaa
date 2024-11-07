//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS2 Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_tensor_container_shared_tensor_hpp
#define glas2_tensor_container_shared_tensor_hpp

#include <glas2/vector/algorithm/prod.hpp>
#include <glas2/tensor/type/contiguous_dense_tensor.hpp>
#ifndef NDEBUG
#include <algorithm>
#include <limits>
#endif

namespace glas2 {

  template <class T, class Shape>
  class shared_tensor
  : public contiguous_dense_tensor<T,Shape>
  {
  public:
    typedef T                               value_type ;
    typedef value_type*                     pointer ;
    typedef value_type const*               const_pointer ;

    typedef Shape                           shape_type ;
    typedef typename shape_type::value_type size_type ;

    typedef shared_tensor< T, Shape >        this_type ;
    typedef contiguous_dense_tensor<T,Shape> base_type ;

  public:
     explicit shared_tensor( shape_type const& shape )
     : base_type( new value_type[ glas2::prod(shape) ], shape )
     , data_( this->ptr() )
     {
       assert( data_ != 0 ) ;
#ifndef NDEBUG
       if (std::numeric_limits<T>::has_quiet_NaN)
         std::fill( this->ptr(), this->ptr()+glas2::prod(shape), std::numeric_limits<T>::quiet_NaN() ) ;
#endif
     }

     ~shared_tensor() {
        if (data_) delete[] data_ ;
     }

     // Copy constructor
     shared_tensor( shared_tensor const& that )
     : base_type( that.ptr(), that.shape() )
     , data_( 0 )
     {
#ifdef GLAS_WARN_CONTAINER_COPY
#ifndef NDEBUG
       std::cerr << std::string("Copy constructor of shared_tensor was called") << std::endl ;
#endif
#endif
       *this = that ;
     }

  public:
     base_type& operator=( shared_tensor const& that ) {
       base_type(*this) = base_type(that) ;
       return *this ;
     }

     template <typename E>
     base_type& operator=( E const& that ) {
       base_type(*this) = that ;
       return *this ;
     }

  public:
/*     void resize( shape_type const& shape ) {
       delete [] data_ ;
       data_ = new value_type[ glas2::prod(shape) ] ;
       this->resize( data_, shape ) ;
     }*/

  private:
     pointer    data_ ;
  } ; // shared_tensor


  template <class T, class Shape>
  struct concept< shared_tensor<T,Shape> >
  : concept< contiguous_dense_tensor<T,Shape> >
  {} ;

  template <typename Tag, typename T, typename Shape>
  struct tensor_transform< Tag, shared_tensor<T,Shape> >
  : tensor_transform< Tag, contiguous_dense_tensor<T,Shape> >
  {} ;


} // Namespace glas
  

#endif
