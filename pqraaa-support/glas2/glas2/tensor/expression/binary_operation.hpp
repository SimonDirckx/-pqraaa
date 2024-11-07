//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_tensor_expression_binary_operation_hpp
#define glas2_tensor_expression_binary_operation_hpp

#include <glas2/concept/is.hpp>
#include <glas2/concept/concept.hpp>
#include <glas2/scalar/concept/scalar.hpp>
#include <glas2/vector/algorithm/is_equal.hpp>
#include <glas2/expression/binary_operation.hpp>
#include <glas2/tensor/concept/dense_tensor.hpp>
#include <glas2/tensor/concept/tensor_transform.hpp>
#include <type_traits>

namespace glas2 {

  template <typename S, typename X, typename Op>
  class binary_operation< S, X, Op
                        , typename std::enable_if< is<Scalar,S>::value && is<DenseTensor,X>::value>::type
                        >
  {
    public:
      binary_operation( S const& s, X const& x )
      : scalar_(s)
      , tensor_( x )
      {}

    public:
      typedef typename X::size_type                              size_type ;
      typedef typename X::shape_type                             shape_type ;
      typedef decltype( Op() ( S(), typename X::value_type() ) ) value_type ;

      size_type const& order() const { return tensor_.order() ; }
      shape_type const& shape() const { return tensor_.shape() ; }

      value_type const& scalar() const {return scalar_ ; }
      X tensor() const {return tensor_ ; }

      value_type operator() ( size_type i, size_type j) const { return Op()( scalar_, tensor_(i,j) ) ; }

    private:
      S const& scalar_ ;
      X        tensor_ ;
  } ;

  template <typename S, typename X, typename Op>
  struct concept< binary_operation< S, X, Op >
                        , typename std::enable_if< is<Scalar,S>::value && is<DenseTensor,X>::value>::type
                        >
  : DenseTensor
  {} ;

  template <typename Tag, typename S, typename X, typename Op>
  struct tensor_transform< Tag, binary_operation<S,X,Op>
                         , typename std::enable_if< is<Scalar,S>::value && is<DenseTensor,X>::value>::type
                         >{
    typedef tensor_transform<Tag,X>                                   trans_type ;
    typedef binary_operation<S, typename trans_type::result_type, Op> result_type ;

    static result_type apply( binary_operation<S,X,Op> x ) {
      return result_type( x.scalar(), trans_type::apply( x.tensor() ) ) ;
    }
  } ;



  template <typename X, typename S, typename Op>
  class binary_operation< X, S, Op
                        , typename std::enable_if< is<Scalar,S>::value && is<DenseTensor,X>::value>::type
                        >
  {
    public:
      binary_operation( X const& x, S const& s )
      : tensor_( x )
      , scalar_(s)
      {}

    public:
      typedef typename X::size_type                              size_type ;
      typedef typename X::shape_type                             shape_type ;
      typedef decltype( Op() ( typename X::value_type(), S() ) ) value_type ;

      size_type const& order() const { return tensor_.order() ; }
      shape_type const& shape() const { return tensor_.shape() ; }

      value_type operator() ( size_type i, size_type j) const { return Op()( tensor_(i,j), scalar_ ) ; }

      value_type const& scalar() const {return scalar_ ; }
      X tensor() const {return tensor_ ; }

    private:
      X        tensor_ ;
      S const& scalar_ ;
  } ;

  template <typename X, typename S, typename Op>
  struct concept< binary_operation< X, S, Op >
                , typename std::enable_if< is<Scalar,S>::value && is<DenseTensor,X>::value>::type
                >
  : DenseTensor
  {} ;

  template <typename Tag, typename X, typename S, typename Op>
  struct tensor_transform< Tag, binary_operation<X,S,Op>
                         , typename std::enable_if< is<Scalar,S>::value && is<DenseTensor,X>::value>::type
                         >
  {
    typedef tensor_transform<Tag,X>                                   trans_type ;
    typedef binary_operation<typename trans_type::result_type, S, Op> result_type ;

    static result_type apply( binary_operation<X,S,Op> x ) {
      return result_type( trans_type::apply( x.tensor() ), x.scalar() ) ;
    }
  } ;




  template <typename X1, typename X2, typename Op>
  class binary_operation< X1, X2, Op
                        , typename std::enable_if< is<DenseTensor,X1>::value && is<DenseTensor,X2>::value>::type
                        >
  {
    public:
      binary_operation( X1 const& x1, X2 const& x2 )
      : x1_( x1 )
      , x2_( x2 )
      {}

    public:
      typedef typename X1::size_type                                                    size_type ;
      typedef typename X1::shape_type                                                   shape_type ;
      typedef decltype( Op() ( typename X2::value_type(), typename X2::value_type() ) ) value_type ;

      size_type const& order() const {
        assert( x1_.order() == x2_.order() ) ;
        return x1_.order() ;
      }

      shape_type const& shape() const {
        assert( is_equal( x1_.shape(), x2_.num_shape() ) ) ;
        return x1_.shape() ;
      }

      value_type operator() ( size_type i, size_type j) const { return Op()( x1_(i,j), x2_(i,j) ) ; }

      X1 tensor1() const {return x1_ ; }
      X2 tensor2() const {return x2_ ; }

    private:
      X1 x1_ ;
      X2 x2_ ;
  } ;

  template <typename X1, typename X2, typename Op>
  struct concept< binary_operation< X1, X2, Op >
                , typename std::enable_if< is<DenseTensor,X1>::value && is<DenseTensor,X2>::value>::type
                >
  : DenseTensor
  {} ;

  template <typename Tag, typename X1, typename X2, typename Op>
  struct tensor_transform< Tag, binary_operation<X1,X2,Op>
                         , typename std::enable_if< is<DenseTensor,X1>::value && is<DenseTensor,X2>::value>::type
                         >
  {
    typedef tensor_transform<Tag,X1>                             trans_1_type ;
    typedef tensor_transform<Tag,X2>                             trans_2_type ;
    typedef binary_operation<typename trans_1_type::result_type
                            ,typename trans_2_type::result_type
                            , Op>                                result_type ;

    static result_type apply( binary_operation<X1,X2,Op> x ) {
      return result_type( trans_1_type::apply( x.tensor1() )
                        , trans_2_type::apply( x.tensor2() )
                        ) ;
    }
  } ;

} // namespace glas2

#endif
