//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas_toolbox_dense_tensor_view_fiber_view_hpp
#define glas_toolbox_dense_tensor_view_fiber_view_hpp

#include <glas/config.hpp>
#define GLAS_MY_INCLUDE_LEVEL GLAS_INCLUDE_VIEW
#include <glas/concept/check_include_level.hpp>

#include <glas/concept/dense_matrix_collection.hpp>
#include <glas/toolbox/tensor/concept/tensor_collection.hpp>
#include <glas/toolbox/tensor/concept/shape.hpp>
#include <glas/toolbox/tensor/concept/order.hpp>
#include <glas/concept/strided_dense_vector_collection.hpp>
#include <glas/concept/round_bracketed.hpp>
#include <glas/concept/square_bracketed.hpp>
#include <glas/concept/value_type.hpp>
#include <glas/concept/size.hpp>
#include <glas/toolbox/tensor/concept/shape.hpp>
#include <boost/type_traits/add_const.hpp>
#include <boost/type_traits/remove_const.hpp>
#include <cassert>

#undef GLAS_INCLUDE_LEVEL
#define GLAS_INCLUDE_LEVEL GLAS_INCLUDE_CONTAINER


namespace glas {

  template <typename Tensor, typename Coordinates>
  class fiber_view
  {
  public:
    typedef typename boost::remove_const<Tensor>::type    tensor_type ;
    typedef Tensor&                                       tensor_reference ;
    typedef typename boost::add_const<tensor_type>::type& const_tensor_reference ;
    typedef Coordinates                                   coordinates_type ;
    typedef typename glas::size_type<tensor_type>::type   size_type ;
    typedef typename value_type<Tensor>::type             value_type ;

    BOOST_STATIC_ASSERT( (TensorExpression<Tensor>::value) ) ;

  public:
     explicit fiber_view( tensor_reference tensor, size_type index, Coordinates const& coordinates )
     : tensor_( tensor )
     , index_(index)
     , coordinates_(coordinates)
     {}

  public:
     
     size_type size() const {
       return glas::shape(tensor_)[ index_ ] ;
     }

  public:
     tensor_reference tensor() { return tensor_ ; }
     const_tensor_reference tensor() const { return tensor_ ; }

     size_type index() const {return index_;}
     Coordinates const& coordinates() const {return coordinates_;}

  public:
    fiber_view& operator=( fiber_view const& that ) {
      return select_backend::assign( *this, that ) ;
    }

    template <typename E>
    fiber_view& operator=( E const& that ) {
      return select_backend::assign( *this, that ) ;
    }

  private:
     tensor_reference   tensor_ ;
     size_type          index_ ;
     Coordinates const& coordinates_ ;
  } ; // fiber_view


  template <typename Tensor, typename Coordinates>
  struct VectorExpression< fiber_view<Tensor,Coordinates> >
  : TensorExpression< Tensor >
  {} ;

  template <typename Tensor, typename Coordinates>
  struct VectorCollection< fiber_view<Tensor,Coordinates> >
  : TensorCollection< Tensor >
  {} ;

  template <typename Tensor, typename Coordinates>
  struct StridedDenseVectorCollection< fiber_view<Tensor,Coordinates> >
  : DenseTensorCollection< Tensor >
  {} ;

/*  template <typename Tensor, typename Coordinates>
  struct RoundBracketed< fiber_view<Tensor,Coordinates> >
  : boost::mpl::true_
  {} ;

  template <typename Tensor, typename Coordinates>
  struct SquareBracketed< fiber_view<Tensor,Coordinates> >
  : boost::mpl::true_
  {} ;
*/

  // StridedDenseVectorCollection

  template <typename Tensor, typename Coordinates>
  struct const_pointer< fiber_view<Tensor,Coordinates> >
  : const_pointer< typename boost::remove_const<Tensor>::type >
  {} ;

  template <typename Tensor, typename Coordinates>
  struct stride_functor< fiber_view<Tensor,Coordinates> const > {
    typedef typename size_type< fiber_view<Tensor,Coordinates> >::type result_type ;
    result_type operator() ( fiber_view<Tensor,Coordinates> const& f ) const {
      result_type str = 1 ;
      for (typename size_type<Tensor>::type i=0; i<f.index(); ++i) str *= glas::shape(f.tensor())[i] ;
      return str ;
    }
  } ;

  template <typename Tensor, typename Coordinates>
  struct storage_ptr_functor< fiber_view<Tensor,Coordinates> const > {
    typedef typename const_pointer< Tensor >::type result_type ;
    result_type operator() ( fiber_view<Tensor,Coordinates> const& f ) const {
      result_type ptr = storage_ptr( f.tensor() ) ;
      int stride = 1 ;
      for (int j=0; j<f.index(); ++j) {
        ptr += stride*f.coordinates()[j] ;
        stride *= shape(f.tensor())[j] ;
      }
      return ptr ;
    }
  } ;

  // StridedDenseVectorCollection

  template <typename Tensor, typename Coordinates>
  struct pointer< fiber_view<Tensor,Coordinates> >
  : boost::mpl::eval_if< is_mutable< Tensor >
                       , pointer< Tensor >
                       , const_pointer< Tensor >
                       >
  {} ;

  template <typename Tensor, typename Coordinates>
  struct storage_ptr_functor< fiber_view<Tensor,Coordinates> > {
    typedef typename pointer< Tensor >::type result_type ;
    result_type operator() ( fiber_view<Tensor,Coordinates>& f ) const {
      result_type ptr = storage_ptr( f.tensor() ) ;
      int stride = 1 ;
      for (int j=0; j<f.index(); ++j) {
        ptr += stride*f.coordinates()[j] ;
        stride *= shape(f.tensor())[j] ;
      }
      return ptr ;
    }
  } ;

  template <typename Tensor, typename Coordinates>
  struct auto_bindings< fiber_view<Tensor,Coordinates> >
  : boost::mpl::true_
  {} ;

} // Namespace glas
  

#endif
