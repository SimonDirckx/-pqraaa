//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_sparse_expression_binary_operation_hpp
#define glas2_sparse_expression_binary_operation_hpp

#include <glas2/concept/is.hpp>
#include <glas2/concept/concept.hpp>
#include <glas2/expression/binary_operation.hpp>
#include <glas2/sparse/concept/coordinate_sparse_matrix.hpp>
#include <glas2/sparse/concept/compressed_sparse_matrix.hpp>
#include <glas2/matrix/concept/matrix_transform.hpp>
#include <glas2/vector/algorithm/ops.hpp>
#include <type_traits>

namespace glas2 {

  template <typename S, typename M, typename Op>
  class binary_operation< S, M, Op
                        , typename std::enable_if< is<Scalar,S>::value && is<CoordinateSparseMatrix,M>::value>::type
                        >
  {
    public:
      binary_operation( S const& s, M const& m )
      : scalar_(s)
      , matrix_( m )
      {}

    public:
      typedef typename M::size_type                              size_type ;
      typedef typename M::index_type                             index_type ;
      typedef decltype( Op() ( S(), typename M::value_type() ) ) value_type ;

      static int const index_base = M::index_base ;

      size_type num_rows() const { return matrix_.num_rows() ; }
      size_type num_columns() const { return matrix_.num_columns() ; }
      size_type num_nz() const { return matrix_.num_nz() ; }

      typedef binary_operation< S, typename M::data_type, Op > data_type ;
      data_type data() const { return data_type( scalar_, matrix_.data() ) ; }

      typedef typename M::rows_type rows_type ;
      rows_type const& row_indices() const { return matrix_.row_indices() ; }

      typedef typename M::columns_type columns_type ;
      columns_type const& column_indices() const { return matrix_.column_indices() ; }

      size_type row(index_type i) const { return matrix_.row(i) ; }
      size_type column(index_type i) const { return matrix_.column(i) ; }

    public:
      value_type const& scalar() const {return scalar_ ; }
      M matrix() const {return matrix_ ; }

    private:
      S const& scalar_ ;
      M        matrix_ ;
  } ;

  template <typename S, typename M, typename Op>
  struct glas_concept< binary_operation< S, M, Op >
                        , typename std::enable_if< is<Scalar,S>::value && is<CoordinateSparseMatrix,M>::value>::type
                        >
  : CoordinateSparseMatrix
  {} ;

  template <typename Tag, typename S, typename M, typename Op>
  struct matrix_transform< Tag, binary_operation<S,M,Op>
                         , typename std::enable_if< is<Scalar,S>::value && is<CoordinateSparseMatrix,M>::value>::type
                         >{
    typedef matrix_transform<Tag,M>                                   trans_type ;
    typedef binary_operation<S, typename trans_type::result_type, Op> result_type ;

    static result_type apply( binary_operation<S,M,Op> m ) {
      return result_type( m.scalar(), trans_type::apply( m.matrix() ) ) ;
    }
  } ;


  //
  // CompressedSparseMatrix
  //
  template <typename S, typename M, typename Op>
  class binary_operation< S, M, Op
                        , typename std::enable_if< is<Scalar,S>::value && is<CompressedSparseMatrix,M>::value>::type
                        >
  {
    public:
      binary_operation( S const& s, M const& m )
      : scalar_(s)
      , matrix_( m )
      {}

    public:
      typedef typename M::size_type                              size_type ;
      typedef typename M::index_type                             index_type ;
      typedef typename M::orientation                            orientation ;
      typedef decltype( Op() ( S(), typename M::value_type() ) ) value_type ;

      static int const index_base = M::index_base ;

      size_type num_rows() const { return matrix_.num_rows() ; }
      size_type num_columns() const { return matrix_.num_columns() ; }
      size_type num_nz() const { return matrix_.num_nz() ; }

      typedef binary_operation< S, typename M::data_type, Op > data_type ;
      data_type data() const { return data_type( scalar_, matrix_.data() ) ; }

      typedef typename M::compressed_index_type compressed_index_type ;
      compressed_index_type compressed_index() const { return matrix_.compressed_index() ; }

      typedef typename M::direct_index_type direct_index_type ;
      direct_index_type direct_index() const { return matrix_.direct_index() ; }

    public:
      value_type const& scalar() const {return scalar_ ; }
      M matrix() const {return matrix_ ; }

    private:
      S const& scalar_ ;
      M        matrix_ ;
  } ;

  template <typename S, typename M, typename Op>
  struct glas_concept< binary_operation< S, M, Op >
                        , typename std::enable_if< is<Scalar,S>::value && is<CompressedSparseMatrix,M>::value>::type
                        >
  : CompressedSparseMatrix
  {} ;

  template <typename Tag, typename S, typename M, typename Op>
  struct matrix_transform< Tag, binary_operation<S,M,Op>
                         , typename std::enable_if< is<Scalar,S>::value && is<CompressedSparseMatrix,M>::value>::type
                         >{
    typedef matrix_transform<Tag,M>                                   trans_type ;
    typedef binary_operation<S, typename trans_type::result_type, Op> result_type ;

    static result_type apply( binary_operation<S,M,Op> m ) {
      return result_type( m.scalar(), trans_type::apply( m.matrix() ) ) ;
    }
  } ;


/*
  template <typename M, typename S, typename Op>
  class binary_operation< M, S, Op
                        , typename std::enable_if< is<Scalar,S>::value && is<CoordinateSparseMatrix,M>::value>::type
                        >
  {
    public:
      binary_operation( M const& m, S const& s )
      : matrix_( m )
      , scalar_(s)
      {}

    public:
      typedef typename M::size_type                              size_type ;
      typedef decltype( Op() ( typename M::value_type(), S() ) ) value_type ;
      typedef typename M::orientation                            orientation ;

      size_type num_rows() const { return matrix_.num_rows() ; }
      size_type num_columns() const { return matrix_.num_columns() ; }

      value_type operator() ( size_type i, size_type j) const { return Op()( matrix_(i,j), scalar_ ) ; }

      value_type const& scalar() const {return scalar_ ; }
      M matrix() const {return matrix_ ; }

    private:
      M        matrix_ ;
      S const& scalar_ ;
  } ;

  template <typename M, typename S, typename Op>
  struct glas_concept< binary_operation< M, S, Op >
                , typename std::enable_if< is<Scalar,S>::value && is<SparseMatrix,M>::value>::type
                >
  : SparseMatrix
  {} ;

  template <typename Tag, typename M, typename S, typename Op>
  struct matrix_transform< Tag, binary_operation<M,S,Op>
                         , typename std::enable_if< is<Scalar,S>::value && is<SparseMatrix,M>::value>::type
                         >
  {
    typedef matrix_transform<Tag,M>                                   trans_type ;
    typedef binary_operation<typename trans_type::result_type, S, Op> result_type ;

    static result_type apply( binary_operation<M,S,Op> m ) {
      return result_type( trans_type::apply( m.matrix() ), m.scalar() ) ;
    }
  } ;




  template <typename M1, typename M2, typename Op>
  class binary_operation< M1, M2, Op
                        , typename std::enable_if< is<SparseMatrix,M1>::value && is<SparseMatrix,M2>::value>::type
                        >
  {
    public:
      binary_operation( M1 const& m1, M2 const& m2 )
      : m1_( m1 )
      , m2_( m2 )
      {}

    public:
      typedef typename M1::size_type                                                    size_type ;
      typedef decltype( Op() ( typename M2::value_type(), typename M2::value_type() ) ) value_type ;

      size_type num_rows() const {
        assert( m1_.num_rows()==m2_.num_rows() ) ;
        return m1_.num_rows() ;
      }

      size_type num_columns() const {
        assert( m1_.num_columns()==m2_.num_columns() ) ;
        return m1_.num_columns() ;
      }

      value_type operator() ( size_type i, size_type j) const { return Op()( m1_(i,j), m2_(i,j) ) ; }

      M1 matrix1() const {return m1_ ; }
      M2 matrix2() const {return m2_ ; }

    private:
      M1 m1_ ;
      M2 m2_ ;
  } ;

  template <typename M1, typename M2, typename Op>
  struct glas_concept< binary_operation< M1, M2, Op >
                , typename std::enable_if< is<SparseMatrix,M1>::value && is<SparseMatrix,M2>::value>::type
                >
  : SparseMatrix
  {} ;

  template <typename Tag, typename M1, typename M2, typename Op>
  struct matrix_transform< Tag, binary_operation<M1,M2,Op>
                         , typename std::enable_if< is<SparseMatrix,M1>::value && is<SparseMatrix,M2>::value>::type
                         >
  {
    typedef matrix_transform<Tag,M1>                             trans_1_type ;
    typedef matrix_transform<Tag,M2>                             trans_2_type ;
    typedef binary_operation<typename trans_1_type::result_type
                            ,typename trans_2_type::result_type
                            , Op>                                result_type ;

    static result_type apply( binary_operation<M1,M2,Op> m ) {
      return result_type( trans_1_type::apply( m.matrix1() )
                        , trans_2_type::apply( m.matrix2() )
                        ) ;
    }
  } ;
*/
} // namespace glas2

#endif
