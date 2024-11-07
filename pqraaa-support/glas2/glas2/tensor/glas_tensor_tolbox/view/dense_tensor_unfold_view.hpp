//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas_toolbox_dense_tensor_view_tensor_view_hpp
#define glas_toolbox_dense_tensor_view_tensor_view_hpp

#include <glas/config.hpp>
#define GLAS_MY_INCLUDE_LEVEL GLAS_INCLUDE_VIEW
#include <glas/concept/check_include_level.hpp>

#include <glas/toolbox/tensor/concept/dense_tensor_collection.hpp>
#include <glas/toolbox/tensor/concept/sizes.hpp>
#include <glas/concept/dense_matrix_collection.hpp>
#include <glas/concept/round_bracketed.hpp>
#include <glas/toolbox/tensor/concept/sizes.hpp>
#include <cassert>

#undef GLAS_INCLUDE_LEVEL
#define GLAS_INCLUDE_LEVEL GLAS_INCLUDE_CONTAINER


namespace glas {

  template <typename Tensor, typename Rows, typename Columns>
  class dense_tensor_unfold_view
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

    typedef dense_tensor_unfold_view<T,R,C>       this_type ;

    BOOST_STATIC_ASSERT( (DenseTensorCollection<T>::value) ) ;

  public:
     explicit dense_tensor()
     : data_( 0 )
     {
       sizes_ = 0 ;
     }

     explicit dense_tensorunfold_view( T& tensor, R const& r, C const& c )
     : tensor_( tensor )
     , rows_(r)
     , columns_(c)
     {
       assert( size(r) + size(c) == order(tensor) ) ;
     }

  public:
     typedef typename value_type< typename sizes_type<T>::type >::type num_rows_type ;
     num_rows_type num_rows() const {
       return product( sizes(tensor_)[ rows_ ] ) ;
     }

     typedef typename value_type< typename sizes_type<T>::type >::type num_columns_type ;
     num_columns_type num_columns() const {
       return product( sizes(tensor_)[ columns_ ] ) ;
     }

  public:
     typedef typename reference<T>::type       reference ;
     typedef typename const_reference<T>::type const_reference ;

     reference operator() ( num_rows_type i, num_columns_type j ) {
     }

     const_reference operator() ( num_rows_type i, num_columns_type j ) const {
     }

  private:
     T&  tensor_ ;
     R const& rows_ ;
     C const& columns_ ;
  } ; // dense_tensor


  template <typename T, typename R, typename C>
  struct DenseMatrixCollection< dense_tensor_unfold_view<T,R,C> >
  : boost::mpl::true_
  {} ;

  template <typename T, typename R, typename C>
  struct RoundBrackted< dense_tensor_unfold_view<T,R,C> >
  : boost::mpl::true_
  {} ;

} // Namespace glas
  

#endif
