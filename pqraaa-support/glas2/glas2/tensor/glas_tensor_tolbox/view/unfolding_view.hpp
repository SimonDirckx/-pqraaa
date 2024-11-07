//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas_toolbox_dense_tensor_view_unfolding_view_hpp
#define glas_toolbox_dense_tensor_view_unfolding_view_hpp

#include <glas/config.hpp>
#define GLAS_MY_INCLUDE_LEVEL GLAS_INCLUDE_VIEW
#include <glas/concept/check_include_level.hpp>

#include <glas/concept/dense_matrix_collection.hpp>
#include <glas/toolbox/tensor/concept/tensor_collection.hpp>
#include <glas/toolbox/tensor/concept/shape.hpp>
#include <glas/toolbox/tensor/concept/order.hpp>
#include <glas/concept/dense_matrix_collection.hpp>
#include <glas/concept/round_bracketed.hpp>
#include <glas/concept/value_type.hpp>
#include <glas/concept/size.hpp>
#include <glas/concept/structure.hpp>
#include <glas/toolbox/tensor/concept/shape.hpp>
#include <cassert>

#undef GLAS_INCLUDE_LEVEL
#define GLAS_INCLUDE_LEVEL GLAS_INCLUDE_CONTAINER


namespace glas {

  template <typename Tensor, typename Rows, typename Columns>
  class unfolding_view
  {
  public:
    typedef Tensor                               tensor_type ;
    typedef typename Tensor::value_type          value_type ;
    typedef general_tag                          structure_type ;

    BOOST_STATIC_ASSERT( (TensorExpression<Tensor>::value) ) ;

  public:
     explicit unfolding_view( Tensor& tensor, Rows const& rows, Columns const& columns )
     : tensor_( tensor )
     , rows_(rows)
     , columns_(columns)
     {
       assert( glas::size(rows)+glas::size(columns)==glas::order(tensor) ) ;
#ifndef NDEBUG
       for (typename size_type<Rows>::type i=1; i<size(rows); ++i) {
         assert( rows[i]>rows[i-1] ) ;
       }
       for (typename size_type<Columns>::type i=1; i<size(columns); ++i) {
         assert( columns[i]>columns[i-1] ) ;
       }
#endif
     }

  public:
     
     typedef typename glas::value_type< typename shape_type<Tensor>::type >::type num_rows_type ;
     num_rows_type num_rows() const {
       return product( shape(tensor_)[ rows_ ] ) ;
     }

     typedef typename glas::value_type< typename shape_type<Tensor>::type >::type num_columns_type ;
     num_columns_type num_columns() const {
       return product( shape(tensor_)[ columns_ ] ) ;
     }

  public:
    unfolding_view& operator=( unfolding_view const& that ) {
      return select_backend::assign( *this, that ) ;
    }

    template <typename E>
    unfolding_view& operator=( E const& that ) {
      return select_backend::assign( *this, that ) ;
    }

  private:
     Tensor&        tensor_ ;
     Rows const&    rows_ ;
     Columns const& columns_ ;
  } ; // unfolding_view


  template <typename T, typename R, typename C>
  struct MatrixExpression< unfolding_view<T,R,C> >
  : TensorExpression< T >
  {} ;

  template <typename T, typename R, typename C>
  struct DenseMatrixExpression< unfolding_view<T,R,C> >
  : DenseTensorExpression< T >
  {} ;

  template <typename T, typename R, typename C>
  struct MatrixCollection< unfolding_view<T,R,C> >
  : TensorCollection< T >
  {} ;

  template <typename T, typename R, typename C>
  struct DenseMatrixCollection< unfolding_view<T,R,C> >
  : DenseTensorCollection< T >
  {} ;

  template <typename T, typename R, typename C>
  struct RoundBracketed< unfolding_view<T,R,C> >
  : boost::mpl::true_
  {} ;

  template <typename T, typename R, typename C>
  struct specialize_const_closure< unfolding_view<T,R,C> >
  : boost::mpl::true_
  {} ;

  template <typename T, typename R, typename C>
  struct const_closure_type< unfolding_view<T,R,C> >
  {
    typedef unfolding_view<typename boost::add_const<T>::type, R, C> > type ;
  } ;

  template <typename M>
  struct orientation< matrix_row_range_view<M> >
  : orientation< typename boost::remove_const<M>::type >
  {} ;

  // DenseMatrixCollection

  // StridedDenseMatrixCollection

  template <typename M>
  struct leading_dimension_type< matrix_row_range_view<M> >
  : leading_dimension_type< M >
  {} ;

  template <typename M>
  struct const_pointer< matrix_row_range_view<M> >
  : const_pointer< typename boost::remove_const<M>::type >
  {} ;

  template <typename M>
  struct pointer< matrix_row_range_view<M> >
  : boost::mpl::eval_if< is_mutable< M >
                       , pointer< typename boost::remove_const<M>::type >
                       , const_pointer< typename boost::remove_const<M>::type >
                       >
  {} ;

  template <typename M>
  struct storage_ptr_functor< matrix_row_range_view<M> const > {
    typedef typename const_pointer< typename boost::remove_const<M>::type >::type result_type ;
    result_type operator() ( matrix_row_range_view<M> const& m ) const {
      return storage_ptr( m.first_argument() ) + stride_1( m.first_argument() ) * m.second_argument().begin() ;
    }
  } ;

  template <typename M>
  struct storage_ptr_functor< matrix_row_range_view<M> > {
    typedef typename pointer< typename boost::remove_const<M>::type >::type result_type ;
    result_type operator() ( matrix_row_range_view<M>& m ) const {
      return storage_ptr( m.first_argument() ) + stride_1( m.first_argument() ) * m.second_argument().begin() ;
    }
  } ;
  template <typename M>
  struct storage_ptr_functor< matrix_row_range_view<M const> >
  : storage_ptr_functor< matrix_row_range_view<M const> const >
  {} ;

  template <typename M>
  struct leading_dimension_functor< matrix_row_range_view<M> const > {
    typename leading_dimension_type< typename boost::remove_const<M>::type >::type operator()( matrix_row_range_view<M> const& m ) const {
      return glas::leading_dimension( m.first_argument() ) ;
    }
  } ;
  template <typename M>
  struct leading_dimension_functor< matrix_row_range_view<M> >
  : leading_dimension_functor< matrix_row_range_view<M> const >
  {} ;

} // Namespace glas
  

#endif
