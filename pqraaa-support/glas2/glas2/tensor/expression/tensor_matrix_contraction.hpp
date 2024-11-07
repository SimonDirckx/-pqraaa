//  (C) Copyright Karl Meerbergen 2015.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_tensor_expression_tensor_matrix_contraction_hpp
#define glas2_tensor_expression_tensor_matrix_contraction_hpp

#include <glas2/matrix/algorithm/multiply.hpp>
#include <glas2/matrix/algorithm/ops.hpp>
#include <glas2/matrix/algorithm/vec.hpp>
#include <glas2/tensor/concept/dense_tensor.hpp>
#include <glas2/vector/container/shared_vector.hpp>
//#include <glas2/backend/default_backend/matrix/matrix_vector_multiply.hpp>
#include <cassert>

namespace glas2 {

  template <typename T, typename M>
  class tensor_matrix_contraction {
    public:
      typedef decltype( typename T::value_type() * typename M::value_type() ) value_type ;
      typedef decltype( typename T::size_type() + typename M::size_type() )   size_type ;
      typedef shared_vector< size_type >                                      shape_type ;

    public:
      tensor_matrix_contraction( T const& t, M const& m, int mode1, int mode2 )
      : t_( t )
      , m_( m )
      , mode1_( mode1 )
      , mode2_( mode2 )
      {
        assert( mode1>=0 && mode1<t.order() ) ;
        assert( mode2>=0 && mode2<t.order() ) ;
        assert( t.shape()(mode1)==m.num_rows() ) ;
        assert( t.shape()(mode2)==m.num_columns() ) ;

        if (mode1_<mode2_) {
          mode_min_ = mode1_ ;
          mode_max_ = mode2_ ;
        } else {
          mode_min_ = mode2_ ;
          mode_max_ = mode1_ ;
        }
      }

    public:
      size_type order() const { return t_.order()-2 ; }

      shape_type shape() const {
        shape_type s(order());
        s(glas2::range(0,mode_min_))=t_.shape()(glas2::range(0,mode_min_));
        s(glas2::range(mode_min_,mode_max_-1))=t_.shape()(glas2::range(mode_min_+1,mode_max_));
        s(glas2::range(mode_max_-1,s.size()))=t_.shape()(glas2::range(mode_max_+1,t_.order()));
        return s ;
      }

      template <typename I>
      value_type operator() ( I index ) const {
        assert( index.size()==order() ) ;

        shape_type s(t_.order());
        s(glas2::range(0,mode_min_))=index(glas2::range(0,mode_min_));
        s(mode_min_) = 0 ;
        s(glas2::range(mode_min_+1,mode_max_))=index(glas2::range(mode_min_,mode_max_-1));
        s(mode_max_) = 0 ;
        s(glas2::range(mode_max_+1,s.size()))=index(glas2::range(mode_max_-1,index.size()));

        value_type sum = 0.0 ;
        for (size_type i=0; i<m_.num_columns(); ++i) {
          s(mode2_) = i ;
          sum += inner_prod( fiber( t_, s, mode1_ ), m_(all(), i) ) ;
        }
        return sum ;
      }

    private:
      T   t_ ;
      M   m_ ;
      int mode1_ ;
      int mode2_ ;
      int mode_min_ ;
      int mode_max_ ;
  } ;

  template <typename T, typename M>
  struct concept< tensor_matrix_contraction<T,M> >
  : DenseTensor
  {} ;

} // namespace glas2

#endif
