//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_backend_default_backend_sparse_sparse_matrix_vector_multiply_hpp
#define glas2_backend_default_backend_sparse_sparse_matrix_vector_multiply_hpp

#include <glas2/backend/default_backend/default_backend.hpp>
#include <glas2/backend/default_backend/vector/fill.hpp>
#include <glas2/backend/default_backend/sparse/is_coordinate_sparse_matrix_multiply.hpp>
#include <glas2/backend/default_backend/sparse/is_compressed_sparse_matrix_multiply.hpp>
#include <glas2/matrix/concept/dense_matrix.hpp>
#include <glas2/sparse/concept/coordinate_sparse_matrix.hpp>
#include <glas2/concept/is.hpp>
#include <type_traits>
#include <cassert>
#include <iostream>
#ifdef GLAS_OPENMP
#include <omp.h>
#endif
namespace glas2 {

  template <typename M, typename V>
  struct is_coordinate_sparse_matrix_multiply< sparse_matrix_vector_multiply<M,V> >
  : is<CoordinateSparseMatrix,M>
  {} ;

  template <typename M, typename V>
  struct is_compressed_sparse_matrix_multiply< sparse_matrix_vector_multiply<M,V>, typename std::enable_if< is<CompressedSparseMatrix,M>::value >::type >
  : std::is_same< typename M::orientation, column_major >
  {
    typedef typename M::orientation orientation ;
  } ;

  template <typename To, typename From>
  typename std::enable_if< is<DenseVector,To>::value && is_coordinate_sparse_matrix_multiply<From>::value
                         , To&
                         >::type plus_assign( default_backend, To& to, From const& from )
  {
    assert( from.size()==to.size() );

    auto const& matrix = from.sparse() ;
    auto const& vector = from.vector() ;

    for (typename From::size_type i=0; i<matrix.num_nz(); ++i) {
      to( matrix.row(i) ) += matrix.data()(i) * vector(matrix.column(i)) ;
    }
    return to ;
  }

/*
  template <typename To, typename From>
  typename std::enable_if< is<DenseVector,To>::value && is_compressed_sparse_matrix_multiply<From>::value && std::is_same< row_major, typename is_compressed_sparse_matrix_multiply<From>::orientation >::value
                         , To&
                         >::type plus_assign( default_backend, To& to, From const& from )
  {
    assert( from.size()==to.size() );

    auto const& matrix = from.sparse() ;
    auto const& vector = from.vector() ;
#ifdef GLAS_OPENMP
    #pragma omp parallel for
#endif 
    for (typename From::size_type i=0; i<matrix.num_rows(); ++i) {
      for (typename From::size_type j=matrix.compressed_index()(i); j<matrix.compressed_index()(i+1); ++j) {
        to( i ) += matrix.data()(j) * vector( matrix.direct_index()(j) ) ;
      }
    }
    return to ;
  }
*/

  template <typename To, typename From>
  typename std::enable_if< is<DenseVector,To>::value && is_compressed_sparse_matrix_multiply<From>::value && std::is_same< column_major, typename is_compressed_sparse_matrix_multiply<From>::orientation >::value
                         , To&
                         >::type plus_assign( default_backend, To& to, From const& from )
  {
    assert( from.size()==to.size() );

    auto const& matrix = from.sparse() ;
    auto const& vector = from.vector() ;
    auto nc=matrix.num_columns();
    for (typename From::size_type i=0; i<nc; ++i) {
        for (typename From::size_type j=matrix.compressed_index()(i); j<matrix.compressed_index()(i+1); ++j) {
          to( matrix.direct_index()(j) ) += matrix.data()(j) * vector( i ) ;
        }
    }
    return to ;
  }


  template <typename To, typename From>
  typename std::enable_if< is<DenseVector,To>::value && is_coordinate_sparse_matrix_multiply<From>::value
                         , To&
                         >::type minus_assign( default_backend, To& to, From const& from ) {
    assert( from.size()==to.size() );

    auto const& matrix = from.sparse() ;
    auto const& vector = from.vector() ;
    for (typename From::size_type i=0; i<matrix.num_nz(); ++i) {
      to( matrix.row(i) ) -= matrix.data()(i) * vector(matrix.column(i)) ;
    }
    return to ;
  }

/*
  template <typename To, typename From>
  typename std::enable_if< is<DenseVector,To>::value && is_compressed_sparse_matrix_multiply<From>::value && std::is_same< row_major, typename is_compressed_sparse_matrix_multiply<From>::orientation >::value
                         , To&
                         >::type minus_assign( default_backend, To& to, From const& from ) {
    assert( from.size()==to.size() );

    auto const& matrix = from.sparse() ;
    auto const& vector = from.vector() ;
    for (typename From::size_type i=0; i<matrix.num_rows(); ++i) {
      for (typename From::size_type j=matrix.compressed_index()(i); j<matrix.compressed_index()(i+1); ++j) {
        to( i ) -= matrix.data()(j) * vector( matrix.direct_index()(j) ) ;
      }
    }
    return to ;
  }
*/
  template <typename To, typename From>
  typename std::enable_if< is<DenseVector,To>::value && is_compressed_sparse_matrix_multiply<From>::value && std::is_same< column_major, typename is_compressed_sparse_matrix_multiply<From>::orientation >::value
                         , To&
                         >::type minus_assign( default_backend, To& to, From const& from )
  {
    assert( from.size()==to.size() );

    auto const& matrix = from.sparse() ;
    auto const& vector = from.vector() ;
    for (typename From::size_type i=0; i<matrix.num_columns(); ++i) {
      for (typename From::size_type j=matrix.compressed_index()(i); j<matrix.compressed_index()(i+1); ++j) {
        to( matrix.direct_index()(j) ) -= matrix.data()(j) * vector( i ) ;
      }
    }
    return to ;
  }



  template <typename To, typename From>
  typename std::enable_if< is<DenseVector,To>::value && (is_coordinate_sparse_matrix_multiply<From>::value || is_compressed_sparse_matrix_multiply<From>::value)
                         , To&
                         >::type assign( default_backend, To& to, From const& from )
  {
    assert( from.size()==to.size() );

    fill( default_backend(), to, 0.0 ) ;

    return plus_assign( default_backend(), to, from ) ;
  }

} // namespace glas2

#endif
