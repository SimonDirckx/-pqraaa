//  (C) Copyright Karl Meerbergen 2011.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas_toolbox_tensor_container_dynamic_dense_tensor_hpp
#define glas_toolbox_tensor_container_dynamic_dense_tensor_hpp

#include <glas/config.hpp>
#define GLAS_MY_INCLUDE_LEVEL GLAS_INCLUDE_CONTAINER
#include <glas/concept/check_include_level.hpp>

#include <glas/toolbox/tensor/concept/round_bracketed.hpp>
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

  template <typename Next>
  class dynamic_nesting {
    public:
      typedef std::ptrdiff_t size_type ;

    private:

    public:

    private:
      size_type size_ ;
      Next      next_ ;
  } ; // dynamic_nesting


  template <int Order, class T>
  class dynamic_dense_tensor
  {
  public:
    typedef T                                     value_type ;
    typedef T const&                              const_reference ;
    typedef T&                                    reference ;
    typedef std::ptrdiff_t                        size_type ;
    typedef value_type*                           pointer ;
    typedef value_type const*                     const_pointer ;

    typedef 

    typedef dynamic_dense_tensor< Order, T >   this_type ;

    typedef tensor_size<Order>

  public:
     explicit dynamic_dense_tensor()
     : num_rows_( 0 )
     , num_columns_( 0 )
     , data_( 0 )
     {}

     explicit dense_tensor( num_rows_type const& rows, num_columns_type const& columns )
     : num_rows_( rows )
     , num_columns_( columns )
     , data_( new value_type[ num_rows_ * num_columns_ ] )
     {
       assert( data_ != 0 ) ;
#ifndef NDEBUG
       if (std::numeric_limits<T>::has_quiet_NaN)
         std::fill( this->storage_ptr(), this->storage_ptr()+num_rows_ * num_columns_, std::numeric_limits<T>::quiet_NaN() ) ;
#endif
     }

     ~dense_tensor() {
        delete[] data_ ;
     }

     template <typename E>
     explicit dense_tensor( E const& that )
     : num_rows_( glas::num_rows( that ) )
     , num_columns_( glas::num_columns( that ) )
     , data_( new value_type[ num_rows_ * num_columns_ ] )
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

#ifdef GLAS_MOVE_SEMANTICS
     dense_tensor& operator=( dense_tensor&& that ) {
       if (data_) delete [] data_ ;
       data_ = that.data_ ; num_rows_ = that.num_rows_ ; num_columns_ = that.num_columns_ ;
       that.data_ = 0 ; that.num_rows_ = 0 ; that.num_columns_ = 0 ;
       return *this ;
     }
#endif

     template <typename E>
     dense_tensor& operator=( E const& that ) {
       select_backend::assign( *this, that ) ;
       return *this ;
     }

  public:
     num_rows_type num_rows() const {
       return num_rows_ ;
     }

     num_columns_type num_columns() const {
       return num_columns_ ;
     }

     typename leading_dimension_type< dense_tensor >::type leading_dimension() const {
       return detail::leading_dimension< O, num_rows_type, num_columns_type >() ( num_rows_, num_columns_ ) ;
     }

     pointer storage_ptr() {
       return data_ ;
     }

     const_pointer storage_ptr() const {
       return data_ ;
     }

  public:
     void resize( num_rows_type num_rows, num_columns_type num_columns ) {
       delete [] data_ ;
       data_ = new value_type[ num_rows * num_columns ] ;
       num_rows_ = num_rows ;
       num_columns_ = num_columns ;
     }

  public:
     const_reference operator()( num_rows_type row, num_columns_type column ) const {
       assert( row<num_rows_ ) ;
       assert( column<num_columns_ ) ;
       return detail::matrix_data<O,const_reference,const_pointer,num_rows_type, num_columns_type>() ( data_, num_rows_, num_columns_, row, column ) ;
     }

     reference operator()( num_rows_type row, num_columns_type column ) {
       assert( row<num_rows_ ) ;
       assert( column<num_columns_ ) ;
       return detail::matrix_data<O,reference,pointer,num_rows_type, num_columns_type>() ( data_, num_rows_, num_columns_, row, column ) ;
     }

  public:
     matrix_row_selector< this_type, num_rows_type > operator[]( num_rows_type r ) { return matrix_row_selector< this_type, num_rows_type >( *this, r ) ; }
     matrix_row_selector< this_type const, num_rows_type > operator[]( num_rows_type const& r ) const { return matrix_row_selector< this_type const, num_rows_type >( *this, r ) ; }

     template <typename S>
     matrix_row_selector< this_type, S > operator[]( S const& r ) { return matrix_row_selector< this_type, S >( *this, r ) ; }

     template <typename S>
     matrix_row_selector< this_type const, S > operator[]( S const& r ) const { return matrix_row_selector< this_type const, S >( *this, r ) ; }

     matrix_column_selector< this_type > operator[]( glas::all ) { return matrix_column_selector< this_type >( *this ) ; }
     matrix_column_selector< this_type const > operator[]( glas::all ) const { return matrix_column_selector< this_type const >( *this ) ; }

  public:
     dense_tensor& swap( dense_tensor& that ) {
       num_rows_type nr( num_rows_ ) ; num_rows_ = that.num_rows_ ; that.num_rows_ = nr ;
       num_columns_type nc( num_columns_ ) ; num_columns_ = that.num_columns_ ; that.num_columns_ = nc ;
       value_type* ptr( data_ ) ; data_ = that.data_ ; that.data_ = ptr ;
       return *this ;
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
