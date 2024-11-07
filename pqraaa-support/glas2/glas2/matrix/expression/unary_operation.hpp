//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_matrix_expression_unary_operation_hpp
#define glas2_matrix_expression_unary_operation_hpp

#include <glas2/concept/is.hpp>
#include <glas2/concept/concept.hpp>
#include <glas2/scalar/concept/scalar.hpp>
#include <glas2/expression/unary_operation.hpp>
#include <glas2/matrix/view/matrix_selection.hpp>
#include <glas2/matrix/concept/dense_matrix.hpp>
#include <glas2/matrix/concept/matrix_transform.hpp>
#include <type_traits>

namespace glas2 {

  template <typename M, typename Op>
  class unary_operation< M, Op
                        , typename std::enable_if< is<DenseMatrix,M>::value>::type
                        >
  {
    public:
      unary_operation( M const& m, Op const& op=Op() )
      : matrix_( m )
      , op_( op )
      {}

    public:
      typedef typename M::size_type                                       size_type ;
      //typedef typename std::result_of< Op(typename M::value_type) >::type value_type ;
      typedef decltype(Op()(typename M::value_type()))                    value_type ;
      typedef typename M::orientation                                     orientation ;

      size_type num_rows() const { return matrix_.num_rows() ; }
      size_type num_columns() const { return matrix_.num_columns() ; }

      M matrix() const {return matrix_ ; }

      value_type operator() ( size_type i, size_type j) const { return op_( matrix_(i,j) ) ; }

      template <typename I1, typename I2>
      typename std::enable_if< !(std::is_integral<I1>::value && std::is_integral<I2>::value)
                             , unary_operation< typename matrix_selection< M, I1, I2 >::result_type, Op >
                             >::type operator()( I1 const& s1, I2 const& s2 ) const {
        auto select = matrix_selection<M,I1,I2>::apply( matrix(), s1, s2 ) ;
        return unary_operation< decltype(select), Op >( select, op_ ) ;
      }

    private:
      M         matrix_ ;
      Op const& op_ ;
  } ;

  template <typename M, typename Op>
  struct glas_concept< unary_operation< M, Op >
                        , typename std::enable_if< is<DenseMatrix,M>::value >::type
                        >
  : DenseMatrix
  {} ;

  template <typename Tag, typename M, typename Op>
  struct matrix_transform< Tag, unary_operation<M,Op>
                         , typename std::enable_if< is<DenseMatrix,M>::value >::type
                         >{
    typedef matrix_transform<Tag,M>                                   trans_type ;
    typedef unary_operation<typename trans_type::result_type, Op> result_type ;

    static result_type apply( unary_operation<M,Op> m ) {
      return result_type( trans_type::apply( m.matrix() ) ) ;
    }
  } ;



} // namespace glas2

#endif
