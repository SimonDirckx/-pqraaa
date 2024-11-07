//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_matrix_expression_scalar_matrix_expression_hpp
#define glas2_matrix_expression_scalar_matrix_expression_hpp

#include <glas2/concept/is.hpp>
#include <glas2/concept/concept.hpp>
#include <glas2/matrix/concept/dense_matrix.hpp>
#include <type_traits>

namespace glas2 {

  template <typename S, typename M, typename Op>
  class scalar_matrix_expression {
    public:
      scalar_matrix_expression( S const& s, M const& m )
      : scalar_(s)
      , matrix_( m )
      {}

    public:
      typedef typename M::size_type                              size_type ;
      typedef decltype( Op() ( S(), typename M::value_type() ) ) value_type ;
      typedef typename M::orientation                            orientation ;

      size_type num_rows() const { return matrix_.num_rows() ; }
      size_type num_columns() const { return matrix_.num_columns() ; }

      value_type operator() ( size_type i, size_type j) const { return Op()( scalar_, matrix_(i,j) ) ; }

    private:
      S const& scalar_ ;
      M        matrix_ ;
  } ;

  template <typename S, typename M, typename Op>
  struct glas_concept< scalar_matrix_expression<S,M,Op>, typename std::enable_if< is<DenseMatrix,M>::value >::type >
  : DenseMatrix
  {} ;

} // namespace glas2

#endif
