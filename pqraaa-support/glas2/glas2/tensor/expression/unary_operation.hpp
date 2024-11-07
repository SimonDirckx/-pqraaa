//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_tensor_expression_unary_operation_hpp
#define glas2_tensor_expression_unary_operation_hpp

#include <glas2/concept/is.hpp>
#include <glas2/concept/concept.hpp>
#include <glas2/scalar/concept/scalar.hpp>
#include <glas2/expression/unary_operation.hpp>
#include <glas2/tensor/concept/dense_tensor.hpp>
#include <glas2/tensor/concept/tensor_transform.hpp>
#include <type_traits>

namespace glas2 {

  template <typename M, typename Op>
  class unary_operation< M, Op
                        , typename std::enable_if< is<DenseTensor,M>::value>::type
                        >
  {
    public:
      unary_operation( M const& m )
      : tensor_( m )
      {}

    public:
      typedef typename M::size_type                         size_type ;
      typedef typename M::shape_type                        shape_type ;
      typedef decltype( Op() ( typename M::value_type() ) ) value_type ;

      shape_type const& shape() const { return tensor_.shape() ; }

      M tensor() const {return tensor_ ; }

      value_type operator() ( size_type i, size_type j) const { return Op()( tensor_(i,j) ) ; }

    private:
      M        tensor_ ;
  } ;

  template <typename M, typename Op>
  struct concept< unary_operation< M, Op >
                        , typename std::enable_if< is<DenseTensor,M>::value >::type
                        >
  : DenseTensor
  {} ;

  template <typename Tag, typename M, typename Op>
  struct tensor_transform< Tag, unary_operation<M,Op>
                         , typename std::enable_if< is<DenseTensor,M>::value >::type
                         >{
    typedef tensor_transform<Tag,M>                                   trans_type ;
    typedef unary_operation<typename trans_type::result_type, Op> result_type ;

    static result_type apply( unary_operation<M,Op> m ) {
      return result_type( trans_type::apply( m.tensor() ) ) ;
    }
  } ;



} // namespace glas2

#endif
