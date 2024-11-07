//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_vector_expression_unary_operation_hpp
#define glas2_vector_expression_unary_operation_hpp

#include <glas2/type/pass_reference.hpp>
#include <glas2/concept/is.hpp>
#include <glas2/concept/concept.hpp>
#include <glas2/expression/unary_operation.hpp>
#include <glas2/scalar/concept/scalar.hpp>
#include <glas2/vector/algorithm/vector_selection.hpp>
#include <glas2/vector/concept/dense_vector.hpp>
#include <type_traits>

namespace glas2 {

  template <typename V, typename Op>
  class unary_operation< V, Op
                        , typename std::enable_if< is<DenseVector,V>::value >::type
                        >
  {
    public:
      unary_operation( V const& v )
      : vector_( v )
      {}

      unary_operation( V const& v, Op const& op )
      : vector_( v )
      , op_( op )
      {}

    public:
      typedef typename std::decay<V>::type                                       v_type ;
      typedef typename v_type::size_type                                         size_type ;
      typedef decltype( Op()(typename V::value_type()) )   value_type ;

      size_type size() const {
        return vector_.size() ;
      }

      value_type operator() ( size_type i) const { return op_( vector_(i) ) ; }

      template <typename I>
      typename std::enable_if< !std::is_integral<I>::value, typename vector_selection< unary_operation, I >::result_type>::type operator()( I const& s ) const {
        return vector_selection< unary_operation, I >::apply( *this, s ) ;
      }

    public:
      typedef typename pass_reference<v_type>::type v_ref ;
      typename std::add_lvalue_reference< typename std::add_const<v_ref>::type>::type vector() const { return vector_; }

    private:
      v_ref  vector_ ;
      Op     op_ ;
  } ;

  template <typename V, typename Op>
  struct glas_concept< unary_operation<V,Op>
                , typename std::enable_if< is<DenseVector,V>::value>::type
                >
  : DenseVector
  {} ;


  template <typename V, typename Op, typename R>
  struct vector_selection< unary_operation<V,Op>, R, typename std::enable_if< is<glas2::DenseVector,R>::value >::type > {
    typedef typename vector_selection<V,R>::result_type select_type ;
    typedef unary_operation<select_type,Op> result_type ;

     static result_type apply( unary_operation<V,Op> const& v, R const& r ) {
      return result_type( v.vector()(r) ) ;
    }
  } ;


} // namespace glas2

#endif
