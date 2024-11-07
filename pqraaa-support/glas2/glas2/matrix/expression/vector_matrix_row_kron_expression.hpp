//  (C) Copyright Karl Meerbergen 2016.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_matrix_expression_vector_matrix_row_kron_expression_hpp
#define glas2_matrix_expression_vector_matrix_row_kron_expression_hpp

#include <glas2/concept/is.hpp>
#include <glas2/concept/concept.hpp>
#include <glas2/scalar/concept/scalar.hpp>
#include <glas2/matrix/concept/dense_matrix.hpp>
#include <glas2/matrix/concept/dense_vector.hpp>
#include <type_traits>

namespace glas2 {

  template <typename V1, typename V2>
  class vector_matrix_row_kron_expression {
    public:
      vector_matrix_kron_expression( V1 const& v1, V2 const& v2 )
      : v1_( v1 )
      , v2_( v2 )
      {}

      static_assert( is<DenseVector,V1>::value, "vector_matrix_row_kron_expression: V1 is not a Vector" );
      static_assert( is<DenseMatrix,V2>::value, "vector_matrix_row_kron_expression: V2 is not a Matrix" );

    public:
      typedef decltype( typename V1::size_type() * typename V2::size_type() )   size_type ;
      typedef decltype( typename V1::value_type() * typename V2::value_type() ) value_type ;

      size_type num_rows() const {
        return v2_.num_rows() ;
      }

      size_type num_columns() const {
        return v1_.size() * v2_.num_columns() ;
      }

      value_type operator() ( size_type i, size_type j) const { return v1_(j/v2_.num_columns()) * v2_(i,j%v2_.num_columns()) ; }

    private:
      V1 v1_ ;
      V2 v2_ ;
  } ;

  template <typename V1, typename V2>
  struct glas_concept< vector_matrix_row_kron_expression<V1,V2>, typename std::enable_if< is<DenseVector,V1>::value && is<DenseMatrix,V2>::value >::type >
  : DenseMatrix
  {} ;

} // namespace glas2

#endif
