//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_matrix_expression_matrix_scalar_expression_hpp
#define glas2_matrix_expression_matrix_scalar_expression_hpp

#include <glas2/concept/is.hpp>
#include <glas2/concept/concept.hpp>
#include <glas2/scalar/concept/scalar.hpp>
#include <glas2/matrix/concept/dense_matrix.hpp>
#include <type_traits>

namespace glas2 {

  template <typename M, typename S, typename Op>
  class matrix_scalar_expression {
    public:
      matrix_scalar_expression( M const& m, S const& s )
      : matrix_( m )
      , scalar_(s)
      {}

      static_assert( is<Scalar,S>::value, "matrix_scalar_expression: S is not a Scalar" );
      static_assert( is<DenseMatrix,M>::value, "matrix_scalar_expression: M is not a DenseMatrix" );

    public:
      typedef typename M::size_type                              size_type ;
      typedef decltype( Op() ( typename M::value_type(), S() ) ) value_type ;
      typedef typename M::orientation                            orientation ;

      size_type num_rows() const { return matrix_.num_rows() ; }
      size_type num_columns() const { return matrix_.num_columns() ; }

      value_type operator() ( size_type i, size_type j) const { return Op()( matrix_(i,j), scalar_ ) ; }

      typedef typename M::iterator iterator ;
      iterator iterate() const { return matrix_.iterate() ; }
      value_type operator[] ( size_type i) const { return Op()( matrix_[i], scalar_ ) ; }

    private:
      M        matrix_ ;
      S const& scalar_ ;
  } ;

  template <typename M, typename S, typename Op>
  struct concept< matrix_scalar_expression<M,S,Op>, typename std::enable_if< is<DenseMatrix,M>::value >::type >
  : DenseMatrix
  {} ;

} // namespace glas2

#endif
