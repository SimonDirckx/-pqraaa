//  (C) Copyright Karl Meerbergen 2015.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_tensor_expression_tensor_vector_contraction_hpp
#define glas2_tensor_expression_tensor_vector_contraction_hpp

#include <glas2/vector/algorithm/inner_prod.hpp>
#include <glas2/tensor/concept/dense_tensor.hpp>
#include <glas2/vector/container/shared_vector.hpp>
#include <cassert>

namespace glas2 {

  template <typename T, typename V>
  class tensor_vector_contraction {
    public:
      typedef decltype( typename T::value_type() * typename V::value_type() ) value_type ;
      typedef decltype( typename T::size_type() + typename V::size_type() )   size_type ;
      typedef shared_vector< size_type >                                      shape_type ;

    public:
      tensor_vector_contraction( T const& t, V const& v, int mode )
      : t_( t )
      , v_( v )
      , mode_( mode )
      {
        assert( mode>=0 && mode<t.order() ) ;
        assert( t.shape()(mode)==v.size() ) ;
      }

    public:
      size_type order() const { return t_.order()-1 ; }

      shape_type shape() const {
        shape_type s(order());
        s(glas2::range(0,mode_))=t_.shape()(glas2::range(0,mode_));
        s(glas2::range(mode_,s.size()))=t_.shape()(glas2::range(mode_+1,s.size()+1));
        return s ;
      }

      template <typename I>
      value_type operator() ( I index ) const {
        assert( index.size()==order() ) ;

        shape_type s(t_.order());
        s(glas2::range(0,mode_))=index(glas2::range(0,mode_));
        s(mode_) = 0 ;
        s(glas2::range(mode_+1,s.size()))=index(glas2::range(mode_,index.size()));

        return inner_prod( glas2::fiber( t_, s, mode_ ), v_ ) ;
      }

    private:
      T   t_ ;
      V   v_ ;
      int mode_ ;
  } ;

  template <typename T, typename V>
  struct concept< tensor_vector_contraction<T,V> >
  : DenseTensor
  {} ;

} // namespace glas2

#endif
