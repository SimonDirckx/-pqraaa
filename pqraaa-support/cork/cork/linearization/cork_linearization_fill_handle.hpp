//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_linearization_cork_linearization_fill_handle_hpp
#define cork_linearization_cork_linearization_fill_handle_hpp

#include <cork/utility/timer.hpp>
#include <cork/linearization/info.hpp>
#include <cork/coefficient_matrices/coefficient_matrices4cork.hpp>
#include <glas2/vector.hpp>
#include <glas2/matrix.hpp>
#include <glas2/sparse.hpp>
#include <type_traits>
#include <exception>

namespace CORK { namespace linearization {

  template <typename Basis4CORK, typename MatrixIterator, typename CoefficientMatrices>
  class cork_linearization_fill_handle
  {
    public:
      typedef typename std::decay<Basis4CORK>::type                basis4cork_type ;
      typedef typename std::decay<CoefficientMatrices>::type       coefficient_matrices_type ;
      typedef typename std::decay<MatrixIterator>::type            matrix_iterator_type ;

      typedef typename basis4cork_type::size_type                  size_type ;

    public:
      cork_linearization_fill_handle( Basis4CORK poly, MatrixIterator matrix_iterator, CoefficientMatrices coefficient_matrices )
      : basis_( poly )
      , matrix_iterator_( matrix_iterator )
      , coefficient_matrices_( coefficient_matrices )
      {}

    public:

      //
      // Block system
      //
      // [    A - s B  ] y = [      B     ] b
      // [ (M-sN) x I  ]     [ I\otimes N ]
      //

      //
      // Assumed that A is a dense or sparse matrix
      //
      template <typename AMatrix>
      typename std::enable_if<glas2::is<glas2::DenseMatrix, AMatrix>::value>::type A( AMatrix A ) const {
        typedef typename AMatrix::value_type value_type ;
        fill(A,0.0) ;

        int n = coefficient_matrices_.num_rows() ;

        //glas2::matrix< value_type > accumulator( basis_.size(), basis_.num_terms() ) ;
        coefficient_matrices::coefficient_matrices4CORK< value_type, coefficient_matrices_type > accumulator( coefficient_matrices_ ) ;
        glas2::vector<value_type> temp( basis_.size() ) ;
        accumulator.initialize_schedule( reshape(temp, 1, basis_.size(), glas2::column_major() ) ) ;
        matrix_iterator_.schedule_a_0( accumulator, [&]( auto i ) { fill( temp, 0.0) ; temp(i) = 1.0 ; return temp ; } ) ;

        for (int i=0; i<basis_.size(); ++i) {
          auto Ar = A( glas2::range(0,n), glas2::range(i*n,i*n+n) ) ;
          coefficient_matrices_.accumulate( accumulator.accumulator()(i, glas2::all()), Ar ) ;
        }

       /*
        for (int i=0; i<basis_.size(); ++i) {
          for (int j=0; j<basis_.num_terms(); ++j) {
            if (accumulator.accumulator()(i,j)!=0.0) A( glas2::range(0,n), glas2::range(i*n,i*n+n) ) += coefficient_matrices_(j) * accumulator.accumulator()(i,j) ;
          }
        }*/

        glas2::matrix<value_type> M( basis_.size()-1, basis_.size() ) ;
        basis_.fill_M( M ) ;
        A( glas2::range_from_end(n,0), glas2::all() ) = kron( M, glas2::identity_matrix<value_type>(n,n) ) ;
      } // A()

      template <typename AMatrix>
      typename std::enable_if< glas2::is< glas2::CoordinateSparseMatrix, AMatrix >::value >::type A( AMatrix& A ) const {
        typedef typename AMatrix::value_type value_type ;
        A.reset( A.num_rows(), A.num_columns() ) ;

        int n = coefficient_matrices_.num_rows() ;

        coefficient_matrices::coefficient_matrices4CORK< value_type, coefficient_matrices_type > accumulator( coefficient_matrices_ ) ;
        glas2::vector<value_type> temp( basis_.size() ) ;
        accumulator.initialize_schedule( reshape(temp, 1, basis_.size(), glas2::column_major() ) ) ;
        matrix_iterator_.schedule_a_0( accumulator, [&]( auto i ) { fill( temp, 0.0) ; temp(i) = 1.0 ; return temp ; } ) ;

        for (int i=0; i<basis_.size(); ++i) {
          coefficient_matrices_.accumulate( accumulator.accumulator()(i, glas2::all()), A( glas2::range(0,n), glas2::range(i*n,i*n+n) ) ) ;
        }
        /*
        for (int i=0; i<basis_.size(); ++i) {
          for (int j=0; j<basis_.num_terms(); ++j) {
            if (accumulator(i,j)!=0.0) push_back( A, accumulator(i,j) * coefficient_matrices_(j), glas2::range(0,n), glas2::range(i*n,i*n+n) ) ;
          }
        }*/

        glas2::matrix<value_type> M( basis_.size()-1, basis_.size() ) ;
        basis_.fill_M( M ) ;
        for (int i=0; i<M.num_rows(); ++i) {
          for (int j=0; j<M.num_columns(); ++j) {
            if (M(i,j)!=0.0) push_back( A, M(i,j)*glas2::speye(n,n), glas2::range(i*n+n,i*n+2*n), glas2::range(j*n,j*n+n) ) ;
          }
        }
      } // A()

      template <typename BMatrix>
      typename std::enable_if< glas2::is<glas2::DenseMatrix, BMatrix>::value>::type B( BMatrix B ) const {
        typedef typename BMatrix::value_type value_type ;
        fill(B,0.0) ;

        int n = coefficient_matrices_.num_rows() ;

        coefficient_matrices::coefficient_matrices4CORK< value_type, coefficient_matrices_type > accumulator( coefficient_matrices_ ) ;
        glas2::vector<value_type> temp( basis_.size() ) ;
        accumulator.initialize_schedule( reshape(temp, 1, basis_.size(), glas2::column_major() ) ) ;
        matrix_iterator_.schedule_b_0( accumulator, [&]( auto i ) { fill( temp, 0.0) ; temp(i) = 1.0 ; return temp ; } ) ;
        /*for (int i=0; i<basis_.size(); ++i) {
          for (int j=0; j<basis_.num_terms(); ++j) {
            if (accumulator.accumulator()(i,j)!=0.0) B( glas2::range(0,n), glas2::range(i*n,i*n+n) ) += coefficient_matrices_(j) * accumulator.accumulator()(i,j) ;
          }
        }*/
        for (int i=0; i<basis_.size(); ++i) {
          auto Br = B( glas2::range(0,n), glas2::range(i*n,i*n+n) ) ;
          coefficient_matrices_.accumulate( accumulator.accumulator()(i, glas2::all()), Br ) ;
        }

        glas2::matrix<value_type> N( basis_.size()-1, basis_.size() ) ;
        basis_.fill_N( N ) ;
        B( glas2::range_from_end(n,0), glas2::all() ) = kron( N, glas2::identity_matrix<value_type>(n,n) ) ;
      } // B()

      template <typename BMatrix>
      typename std::enable_if< glas2::is<glas2::CoordinateSparseMatrix, BMatrix>::value>::type B( BMatrix& B ) const {
        typedef typename BMatrix::value_type value_type ;
        B.reset( B.num_rows(), B.num_columns() ) ;

        int n = coefficient_matrices_.num_rows() ;

        glas2::matrix< value_type > accumulator( basis_.size(), basis_.num_terms() ) ;
        glas2::vector<value_type> temp( basis_.size() ) ;
        fill( accumulator, 0 ) ;
        matrix_iterator_.schedule_b_0( accumulator, [&]( auto i ) { fill( temp, 0.0) ; temp(i) = 1.0 ; return temp ; } ) ;
        for (int i=0; i<basis_.size(); ++i) {
          for (int j=0; j<basis_.num_terms(); ++j) {
            if (accumulator(i,j)!=0.0) push_back( B, accumulator(i,j) * coefficient_matrices_(j), glas2::range(0,n), glas2::range(i*n,i*n+n) ) ;
          }
        }

        glas2::matrix<value_type> N( basis_.size()-1, basis_.size() ) ;
        basis_.fill_N( N ) ;
        for (int i=0; i<N.num_rows(); ++i) {
          for (int j=0; j<N.num_columns(); ++j) {
            if (N(i,j)!=0.0) push_back( B, N(i,j)*glas2::speye(n,n), glas2::range(i*n+n,i*n+2*n), glas2::range(j*n,j*n+n) ) ;
          }
        }
      } // B()

    private:
      Basis4CORK                       basis_ ;
      MatrixIterator                   matrix_iterator_ ;
      CoefficientMatrices              coefficient_matrices_ ;
  } ; // cork_linearization_fill_handle

} } // namespace CORK::linearization

#endif
