//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_backend_default_backend_sparse_sparse_matrix_multiply_hpp
#define glas2_backend_default_backend_sparse_sparse_matrix_multiply_hpp

#include <glas2/backend/default_backend/default_backend.hpp>
#include <glas2/backend/default_backend/matrix/fill.hpp>
#include <glas2/backend/default_backend/sparse/is_coordinate_sparse_matrix_multiply.hpp>
#include <glas2/backend/default_backend/sparse/is_compressed_sparse_matrix_multiply.hpp>
#include <glas2/matrix/concept/dense_matrix.hpp>
#include <glas2/sparse/concept/coordinate_sparse_matrix.hpp>
#include <glas2/sparse/concept/compressed_sparse_matrix.hpp>
#include <glas2/sparse/type/no_value_type.hpp>
#include <glas2/concept/is.hpp>
#include <type_traits>
#include <cassert>

namespace glas2 {

  namespace detail {
    template <typename T>
    struct compute_product_element_sparse {
      template <typename V>
      auto operator() ( T sparse, V const& v ) const -> decltype(sparse*v) { return sparse * v ; }
    } ;

    template <>
    struct compute_product_element_sparse<no_value_type> {
      template <typename V>
      V const& operator() ( no_value_type , V const& v ) const { return v ; }
    } ;

  }

  template <typename S, typename M>
  struct is_coordinate_sparse_matrix_multiply< sparse_matrix_multiply<S,M> >
  : is<CoordinateSparseMatrix,S>
  {} ;

  template <typename To, typename From>
  typename std::enable_if< is<DenseMatrix,To>::value && is_coordinate_sparse_matrix_multiply<From>::value
                         , To&
                         >::type plus_assign( default_backend, To& to, From const& from )
  {
    assert( from.num_rows()==to.num_rows() );
    assert( from.num_columns()==to.num_columns() );

    auto const& matrix = from.matrix() ;
    auto const& sparse = from.sparse() ;

    for (typename From::size_type i=0; i<sparse.num_nz(); ++i) {
      //to( sparse.row(i), glas2::all() ) += sparse.data()(i) * matrix( sparse.column(i), glas2::all() ) ;
      to( sparse.row(i), glas2::all() ) += detail::compute_product_element_sparse<decltype(sparse.data()(i))>()( sparse.data()(i), matrix( sparse.column(i), glas2::all() ) ) ;
    }
    return to ;
  }


  template <typename M, typename V>
  struct is_compressed_sparse_matrix_multiply< sparse_matrix_multiply<M,V>, typename std::enable_if< is<CompressedSparseMatrix,M>::value >::type >
  : std::true_type
  {
    typedef typename M::orientation orientation ;
  } ;

  template <typename To, typename From>
  typename std::enable_if< is<DenseMatrix,To>::value && is_compressed_sparse_matrix_multiply<From>::value && std::is_same< row_major, typename is_compressed_sparse_matrix_multiply<From>::orientation >::value
                         , To&
                         >::type plus_assign( default_backend, To& to, From const& from )
  {
    assert( from.num_rows()==to.num_rows() );
    assert( from.num_columns()==to.num_columns() );

    auto const& matrix = from.matrix() ;
    auto const& sparse = from.sparse() ;

    for (typename From::size_type i=0; i<sparse.num_rows(); ++i) {
      for (typename From::size_type j=sparse.compressed_index()(i); j<sparse.compressed_index()(i+1); ++j) {
        to( i, glas2::all() ) += detail::compute_product_element_sparse<decltype(sparse.data()(j))>()( sparse.data()(j), matrix( sparse.direct_index()(j), glas2::all() ) ) ;
      }
    }
    return to ;
  }

  template <typename To, typename From>
  typename std::enable_if< is<DenseMatrix,To>::value && is_compressed_sparse_matrix_multiply<From>::value && std::is_same< column_major, typename is_compressed_sparse_matrix_multiply<From>::orientation >::value
                         , To&
                         >::type plus_assign( default_backend, To& to, From const& from )
  {
    assert( from.num_rows()==to.num_rows() );
    assert( from.num_columns()==to.num_columns() );

    auto const& matrix = from.matrix() ;
    auto const& sparse = from.sparse() ;

    for (typename From::size_type i=0; i<sparse.num_columns(); ++i) {
      for (typename From::size_type j=sparse.compressed_index()(i); j<sparse.compressed_index()(i+1); ++j) {
        to( sparse.direct_index()(j), glas2::all() ) += detail::compute_product_element_sparse<decltype(sparse.data()(j))>()( sparse.data()(j), matrix( i, glas2::all() ) ) ;
      }
    }
    return to ;
  }


  template <typename To, typename From>
  typename std::enable_if< is<DenseMatrix,To>::value && is_coordinate_sparse_matrix_multiply<From>::value
                         , To&
                         >::type minus_assign( default_backend, To& to, From const& from ) {
    assert( from.num_rows()==to.num_rows() );
    assert( from.num_columns()==to.num_columns() );

    auto const& matrix = from.matrix() ;
    auto const& sparse = from.sparse() ;

    for (typename From::size_type i=0; i<sparse.num_nz(); ++i) {
      to( sparse.row(i), glas2::all() ) -= //sparse.data()(i) * matrix( sparse.column(i), glas2::all() ) ;
       detail::compute_product_element_sparse<decltype(sparse.data()(0))>()( sparse.data()(i), matrix( sparse.column(i), glas2::all() ) ) ;
    }
    return to ;
  }


  template <typename To, typename From>
  typename std::enable_if< is<DenseMatrix,To>::value && is_compressed_sparse_matrix_multiply<From>::value && std::is_same< row_major, typename is_compressed_sparse_matrix_multiply<From>::orientation >::value
                         , To&
                         >::type minus_assign( default_backend, To& to, From const& from ) {
    assert( from.num_rows()==to.num_rows() );
    assert( from.num_columns()==to.num_columns() );

    auto const& matrix = from.matrix() ;
    auto const& sparse = from.sparse() ;

    for (typename From::size_type i=0; i<sparse.num_rows(); ++i) {
      for (typename From::size_type j=sparse.compressed_index()(i); j<sparse.compressed_index()(i+1); ++j) {
        to( i, glas2::all() ) -= detail::compute_product_element_sparse<decltype(sparse.data()(j))>()( sparse.data()(j), matrix( sparse.direct_index()(j), glas2::all() ) ) ;
      }
    }
    return to ;
  }


  template <typename To, typename From>
  typename std::enable_if< is<DenseMatrix,To>::value && is_compressed_sparse_matrix_multiply<From>::value && std::is_same< column_major, typename is_compressed_sparse_matrix_multiply<From>::orientation >::value
                         , To&
                         >::type minus_assign( default_backend, To& to, From const& from ) {
    assert( from.num_rows()==to.num_rows() );
    assert( from.num_columns()==to.num_columns() );

    auto const& matrix = from.matrix() ;
    auto const& sparse = from.sparse() ;

    for (typename From::size_type i=0; i<sparse.num_columns(); ++i) {
      for (typename From::size_type j=sparse.compressed_index()(i); j<sparse.compressed_index()(i+1); ++j) {
        to( sparse.direct_index()(j), glas2::all() ) -= detail::compute_product_element_sparse<decltype(sparse.data()(j))>()( sparse.data()(j), matrix( i, glas2::all() ) ) ;
      }
    }
    return to ;
  }


  template <typename To, typename From>
  typename std::enable_if< is<DenseMatrix,To>::value && (is_coordinate_sparse_matrix_multiply<From>::value||is_compressed_sparse_matrix_multiply<From>::value)
                         , To&
                         >::type assign( default_backend, To& to, From const& from )
  {
    assert( from.num_rows()==to.num_rows() );
    assert( from.num_columns()==to.num_columns() );

    fill( default_backend(), to, 0.0 ) ;

    return plus_assign( default_backend(), to, from ) ;
  }

} // namespace glas2

#endif
