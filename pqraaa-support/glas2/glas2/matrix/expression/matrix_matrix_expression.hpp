//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_matrix_expression_matrix_matrix_expression_hpp
#define glas2_matrix_expression_matrix_matrix_expression_hpp

#include <glas2/concept/is.hpp>
#include <glas2/concept/concept.hpp>
#include <glas2/scalar/concept/scalar.hpp>
#include <glas2/matrix/concept/dense_matrix.hpp>
#include <type_traits>

namespace glas2 {

  template <typename V1, typename V2, typename Op>
  class matrix_matrix_expression {
    public:
      matrix_matrix_expression( V1 const& v1, V2 const& v2 )
      : v1_( v1 )
      , v2_( v2 )
      {}

      static_assert( is<DenseMatrix,V1>::value, "matrix_matrix_expression: V1 is not a Matrix" );
      static_assert( is<DenseMatrix,V2>::value, "matrix_matrix_expression: V2 is not a Matrix" );

    public:
      typedef typename V1::size_type                                                    size_type ;
      typedef decltype( Op() ( typename V2::value_type(), typename V2::value_type() ) ) value_type ;

      size_type num_rows() const {
        assert( v1_.num_rows()==v2_.num_rows() ) ;
        return v1_.num_rows() ;
      }

      size_type num_columns() const {
        assert( v1_.num_columns()==v2_.num_columns() ) ;
        return v1_.num_columns() ;
      }

      value_type operator() ( size_type i, size_type j) const { return Op()( v1_(i,j), v2_(i,j) ) ; }

    private:
      V1 v1_ ;
      V2 v2_ ;
  } ;

  template <typename V1, typename V2, typename Op>
  struct glas_concept< matrix_matrix_expression<V1,V2,Op>, typename std::enable_if< is<DenseMatrix,V1>::value && is<DenseMatrix,V2>::value >::type >
  : DenseMatrix
  {} ;

} // namespace glas2

#endif
