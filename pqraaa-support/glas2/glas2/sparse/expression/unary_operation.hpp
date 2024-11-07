//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_sparse_expression_unary_operation_hpp
#define glas2_sparse_expression_unary_operation_hpp

#include <glas2/concept/is.hpp>
#include <glas2/concept/concept.hpp>
#include <glas2/expression/unary_operation.hpp>
#include <glas2/sparse/concept/coordinate_sparse_matrix.hpp>
#include <glas2/matrix/concept/matrix_transform.hpp>
#include <glas2/vector/algorithm/ops.hpp>
#include <type_traits>

namespace glas2 {

  template <typename M, typename Op>
  class unary_operation< M, Op
                        , typename std::enable_if< is<CoordinateSparseMatrix,M>::value>::type
                        >
  {
    public:
      unary_operation( M const& m )
      : matrix_( m )
      {}

    public:
      typedef typename M::size_type                         size_type ;
      typedef typename M::index_type                        index_type ;
      typedef decltype( Op() ( typename M::value_type() ) ) value_type ;

      static int const index_base = M::index_base ;

      size_type num_rows() const { return matrix_.num_rows() ; }
      size_type num_columns() const { return matrix_.num_columns() ; }
      size_type num_nz() const { return matrix_.num_nz() ; }

      typedef unary_operation< typename M::data_type, Op > data_type ;
      data_type data() const { return data_type( matrix_.data() ) ; }

      typedef typename M::rows_type rows_type ;
      rows_type const& row_indices() const { return matrix_.row_indices() ; }

      typedef typename M::columns_type columns_type ;
      columns_type const& column_indices() const { return matrix_.column_indices() ; }

      size_type row(index_type i) const { return matrix_.row(i) ; }
      size_type column(index_type i) const { return matrix_.column(i) ; }

    public:
      M matrix() const {return matrix_ ; }

    private:
      M        matrix_ ;
  } ;

  template <typename M, typename Op>
  struct glas_concept< unary_operation< M, Op >
                        , typename std::enable_if< is<CoordinateSparseMatrix,M>::value >::type
                        >
  : CoordinateSparseMatrix
  {} ;

  template <typename Tag, typename M, typename Op>
  struct matrix_transform< Tag, unary_operation<M,Op>
                         , typename std::enable_if< is<CoordinateSparseMatrix,M>::value >::type
                         >{
    typedef matrix_transform<Tag,M>                               trans_type ;
    typedef unary_operation<typename trans_type::result_type, Op> result_type ;

    static result_type apply( unary_operation<M,Op> m ) {
      return result_type( trans_type::apply( m.matrix() ) ) ;
    }
  } ;


} // namespace glas2

#endif
