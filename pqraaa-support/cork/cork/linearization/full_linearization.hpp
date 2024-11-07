//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_cork_linearization_hpp
#define cork_cork_linearization_hpp

#include <glas2/vector.hpp>
#include <glas2/matrix.hpp>
#include <type_traits>
#include <limits>

namespace CORK { namespace linearization {

  template <typename Basis4CORK, typename MatrixIterator, typename CoefficientMatrices, typename LinearSolver>
  class cork_linearization
  {
    public:
      typedef typename std::decay<Basis4CORK>::type                basis4cork_type ;
      typedef typename std::decay<MatrixIterator>::type            matrix_iterator_type ;
      typedef typename std::decay<CoefficientMatrices>::type       coefficient_matrices_type ;
      typedef typename std::decay<LinearSolver>::type              linear_solver_type ;

      typedef typename basis4cork_type::grade_type                 grade_type ;
      typedef typename coefficient_matrices_type::size_type        size_type ;
      typedef typename matrix_iterator_type::value_type            value_type ;

    private:
      typedef decltype(std::abs(value_type()))                     real_value_type ;

    public:
      cork_linearization( Basis4CORK poly, MatrixIterator matrix_iterator, CoefficientMatrices coefficient_matrices, LinearSolver linear_solver )
      : basis_( poly )
      , matrix_iterator_( matrix_iterator )
      , coefficient_matrices_( coefficient_matrices )
      , linear_solver_( linear_solver )
      , drop_tol_( std::numeric_limits<real_value_type>::epsilon() )
      {
        assert( coefficient_matrices_.num_matrices()==basis_.grade()+1 ) ;
      }

    public:
      grade_type grade() const { return basis_.grade() ; }
      size_type size() const { return coefficient_matrices_.size() ; }

      basis4cork_type const& basis() const { return basis_ ; }
      coefficient_matrices_type const& coefficient_matrices() const { return coefficient_matrices_ ; }

    public:
      void shift( value_type s ) {
        basis_.shift( s ) ;
        linear_solver_.prepare_solve( s ) ;
      }

      // (A+sB)^{-1} (I x Q) vec(U)
/*      template <typename QQ, typename UUIn, typename UUOut>
      void starting_vector( QQ const& Q, int& rank, UUIn const& U_in, UUOut U_out ) {
        glas2::all all ;
        assert( Q.num_columns()>rank ) ;
        assert( U_in.num_columns()>=rank ) ;
        assert( U_out.num_columns()>=rank ) ;
        assert( U_in.num_rows()==grade() ) ;
        assert( U_out.num_rows()==grade() ) ;

        auto w = Q( all, rank ) ; // w points next to last column of Q
        {
          glas2::range r_rank(0,rank) ;
          auto Qr = Q( all, r_rank ) ;
          auto Ur = U_in( all, r_rank ) ;
          auto Zr = U_out( glas2::range_from_end(1,0), r_rank ) ;
          polynomial_matrices_.solve_upper( Ur, Zr ) ;

          glas2::shared_vector<value_type> g(rank) ;

          // Compute right-hand side for solve with P.
          glas2::fill( w, 0.0 ) ;

          polynomial_matrices_.solve_upper_no_N( Ur, Zr ) ;

          // Compute B-terms for the right-hand side.
          matrix_iterator_.initialize( Qr ) ;

          auto it_b = matrix_iterator_.begin_b() ;
          auto it_b_end = matrix_iterator_.end_b() ;
          for ( ; it_b!=it_b_end; ++it_b ) {
            if (it_b.index()!=0) it_b.schedule( polynomial_matrices_.shift() * Zr(it_b.index()-1, glas2::all()) ) ;
          }

          // Compute A-terms for the right-hand side.
          auto it_a = matrix_iterator_.begin_a() ;
          auto it_a_end = matrix_iterator_.end_a() ;
          for ( ; it_a!=it_a_end; ++it_a ) {
            if (it_a.index()!=0) {
              it_a.schedule( Zr(it_a.index()-1, glas2::all()) ) ;
            }
          }
          matrix_iterator_.multiply_add( Qr, w ) ;

          w *= -1. ;
          w += multiply( Q, U(0,glas2::all())) ;

          // Do the solve
          matrix_polynomial_.solve( w ) ;

          // Gram-Schmidt orthogonalization and possible expansion of Q.
          auto z_0 = U_out(0, r_rank) ; // First row of U will already be up to date.
          z_0 = multiply( transpose(conj(Qr)), w ) ;
          w -= multiply(Qr, z_0) ;
          g = multiply( transpose(conj(Qr)), w ) ;
          z_0 += g ;
          w -= multiply(Qr, g) ;
          auto norm_w = norm_2(w) ;
          if (norm_w!=0.0) {
            w /= norm_w ;
            U_out(0,rank) = norm_w ;
            ++rank ;
          }
        }
        assert( Q.num_columns()>=rank ) ;
        assert( U_in.num_columns()>=rank ) ;
        assert( U_out.num_columns()>=rank ) ;
        {
          glas2::range r_rank(0,rank) ;
          auto Zr = U_out( glas2::range_from_end(1,0), r_rank ) ;
          auto z_0 = U_out(0, r_rank) ;
          polynomial_matrices_.solve_lower( z_0, Zr ) ;
        }
      } // starting_vector()
*/

      // RKS step for pole set by shift().
      //
      // Argument should be an RKS quadruple
      // Input:
      //   Q
      //   U_in
      //   r
      // Output
      //   Q
      //   U_out
      //   r
      //
      template <typename QQ, typename UUIn, typename UUOut>
      void krylov_step( QQ const& Q, int& rank, UUIn const& U_in, UUOut U_out ) {
        glas2::all all ;
        assert( Q.num_columns()>rank ) ;
        assert( U_in.num_columns()>=rank ) ;
        assert( U_out.num_columns()>=rank ) ;
        assert( U_in.num_rows()==grade() ) ;
        assert( U_out.num_rows()==grade() ) ;

        auto w = Q( all, rank ) ; // w points next to last column of Q
        {
          glas2::range r_rank(0,rank) ;
          auto Qr = Q( all, r_rank ) ;
          auto Ur = U_in( all, r_rank ) ;
          auto Zr = U_out( glas2::range_from_end(1,0), r_rank ) ;
          basis_.multiply_N( Ur, Zr ) ;
          basis_.solve( Zr ) ;

          glas2::shared_vector<value_type> g(rank) ;

          // Compute right-hand side for solve with P.
          glas2::fill( w, 0.0 ) ;

          // Compute B-terms for the right-hand side.
          matrix_iterator_.initialize( Qr ) ;

          auto it_b = matrix_iterator_.begin_b() ;
          auto it_b_end = matrix_iterator_.end_b() ;
          for ( ; it_b!=it_b_end; ++it_b ) {
            it_b.schedule( Ur(it_b.index(), glas2::all()) ) ;
            if (it_b.index()!=0) it_b.schedule( basis_.shift() * Zr(it_b.index()-1, glas2::all()) ) ;
          }

          // Compute A-terms for the right-hand side.
          auto it_a = matrix_iterator_.begin_a() ;
          auto it_a_end = matrix_iterator_.end_a() ;
          for ( ; it_a!=it_a_end; ++it_a ) {
            if (it_a.index()!=0) {
              //temp_ = multiply( Qr, Zr(it_a.index()-1, glas2::all()) ) ;
              //it_a.multiply_add( temp_, w ) ;
              it_a.schedule( Zr(it_a.index()-1, glas2::all()) ) ;
            }
          }
          matrix_iterator_.multiply_add( coefficient_matrices_, Qr, w ) ;

          w *= -1.0 ;

          // Do the solve
          linear_solver_.solve( w ) ;
          auto norm_q = norm_2(w) ;

          // Gram-Schmidt orthogonalization and possible expansion of Q.
          auto z_0 = U_out(0, r_rank) ; // First row of U will already be up to date.
          z_0 = multiply( transpose(conj(Qr)), w ) ;
          w -= multiply(Qr, z_0) ;
          g = multiply( transpose(conj(Qr)), w ) ;
          z_0 += g ;
          w -= multiply(Qr, g) ;
          auto norm_w = norm_2(w) ;
          if (norm_w>drop_tol_*norm_q) {
            w /= norm_w ;
            U_out(0,rank) = norm_w ;
            ++rank ;
          }
          std::cout << Q.num_rows() << " " << rank << std::endl ;
        }
        assert( Q.num_columns()>=rank ) ;
        assert( U_in.num_columns()>=rank ) ;
        assert( U_out.num_columns()>=rank ) ;
        {
          glas2::range r_rank(0,rank) ;
          auto Zr = U_out( glas2::range_from_end(1,0), r_rank ) ;
          auto z_0 = U_out(0, r_rank) ;

          glas2::matrix< value_type > accum( Zr.num_rows(), Zr.num_columns() ) ;
          basis_.lower_solve_right_hand_side( z_0, accum ) ;
          basis_.solve( accum ) ;
          Zr += accum ;
        }
      } // RKS_step()

      template <typename UUIn, typename UUOut>
      void full_krylov_step( UUIn const& U_in, UUOut U_out ) {
        glas2::all all ;
        assert( Q.num_columns()>rank ) ;
        assert( U_in.num_columns()>=rank ) ;
        assert( U_out.num_columns()>=rank ) ;
        assert( U_in.num_rows()==grade() ) ;
        assert( U_out.num_rows()==grade() ) ;

        auto w = Q( all, rank ) ; // w points next to last column of Q
        {
          glas2::range r_rank(0,rank) ;
          auto Qr = Q( all, r_rank ) ;
          auto Ur = U_in( all, r_rank ) ;
          auto Zr = U_out( glas2::range_from_end(1,0), r_rank ) ;
          basis_.multiply_N( Ur, Zr ) ;
          basis_.solve( Zr ) ;

          glas2::shared_vector<value_type> g(rank) ;

          // Compute right-hand side for solve with P.
          glas2::fill( w, 0.0 ) ;

          // Compute B-terms for the right-hand side.
          matrix_iterator_.initialize( Qr ) ;

          auto it_b = matrix_iterator_.begin_b() ;
          auto it_b_end = matrix_iterator_.end_b() ;
          for ( ; it_b!=it_b_end; ++it_b ) {
            it_b.schedule( Ur(it_b.index(), glas2::all()) ) ;
            if (it_b.index()!=0) it_b.schedule( basis_.shift() * Zr(it_b.index()-1, glas2::all()) ) ;
          }

          // Compute A-terms for the right-hand side.
          auto it_a = matrix_iterator_.begin_a() ;
          auto it_a_end = matrix_iterator_.end_a() ;
          for ( ; it_a!=it_a_end; ++it_a ) {
            if (it_a.index()!=0) {
              //temp_ = multiply( Qr, Zr(it_a.index()-1, glas2::all()) ) ;
              //it_a.multiply_add( temp_, w ) ;
              it_a.schedule( Zr(it_a.index()-1, glas2::all()) ) ;
            }
          }
          matrix_iterator_.multiply_add( coefficient_matrices_, Qr, w ) ;

          w *= -1.0 ;

          // Do the solve
          linear_solver_.solve( w ) ;
          auto norm_q = norm_2(w) ;

          // Gram-Schmidt orthogonalization and possible expansion of Q.
          auto z_0 = U_out(0, r_rank) ; // First row of U will already be up to date.
          z_0 = multiply( transpose(conj(Qr)), w ) ;
          w -= multiply(Qr, z_0) ;
          g = multiply( transpose(conj(Qr)), w ) ;
          z_0 += g ;
          w -= multiply(Qr, g) ;
          auto norm_w = norm_2(w) ;
          if (norm_w>drop_tol_*norm_q) {
            w /= norm_w ;
            U_out(0,rank) = norm_w ;
            ++rank ;
          }
          std::cout << Q.num_rows() << " " << rank << std::endl ;
        }
        assert( Q.num_columns()>=rank ) ;
        assert( U_in.num_columns()>=rank ) ;
        assert( U_out.num_columns()>=rank ) ;
        {
          glas2::range r_rank(0,rank) ;
          auto Zr = U_out( glas2::range_from_end(1,0), r_rank ) ;
          auto z_0 = U_out(0, r_rank) ;

          glas2::matrix< value_type > accum( Zr.num_rows(), Zr.num_columns() ) ;
          basis_.lower_solve_right_hand_side( z_0, accum ) ;
          basis_.solve( accum ) ;
          Zr += accum ;
        }
      } // RKS_step()

    private:
      Basis4CORK                       basis_ ;
      MatrixIterator                   matrix_iterator_ ;
      CoefficientMatrices              coefficient_matrices_ ;
      LinearSolver                     linear_solver_ ;
      real_value_type                  drop_tol_ ;
  } ; // default_linearization


} } // namespace CORK::linearization

#endif
