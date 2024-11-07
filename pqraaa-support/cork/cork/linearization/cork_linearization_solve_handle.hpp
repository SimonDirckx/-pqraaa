//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_linearization_cork_linearization_solve_handle_hpp
#define cork_linearization_cork_linearization_solve_handle_hpp

#include <cork/utility/timer.hpp>
#include <cork/linearization/info.hpp>
#include <cork/coefficient_matrices/coefficient_matrices4cork.hpp>
#include <cork/exception/linear_solver_failure.hpp>
#include <glas2/vector.hpp>
#include <glas2/matrix.hpp>
#include <type_traits>
#include <exception>
#include <cmath>

namespace CORK { namespace linearization {

  template <typename LinearSolver, typename Basis4CORK, typename Basis4CORKHandle, typename MatrixIterator, typename CoefficientMatrices>
  class cork_linearization_solve_handle
  {
    public:
      typedef typename std::decay<Basis4CORK>::type                basis4cork_type ;
      typedef typename std::decay<Basis4CORKHandle>::type          basis_handle_type ;
      typedef typename std::decay<MatrixIterator>::type            matrix_iterator_type ;
      typedef typename std::decay<CoefficientMatrices>::type       coefficient_matrices_type ;
      typedef typename std::decay<LinearSolver>::type              linear_solver_type ;

      typedef typename basis4cork_type::size_type                  size_type ;
      typedef typename linear_solver_type::value_type              value_type ;

    private:
      typedef decltype(std::abs(value_type()))                                                         real_value_type ;
      typedef coefficient_matrices::coefficient_matrices4CORK< value_type, coefficient_matrices_type > coefficient_matrices_4cork_type ;

    public:
      cork_linearization_solve_handle( Basis4CORK poly, Basis4CORKHandle poly_handle, MatrixIterator matrix_iterator, CoefficientMatrices coefficient_matrices, LinearSolver linear_solver, info& information )
      : basis_( poly )
      , basis_handle_( poly_handle )
      , matrix_iterator_( matrix_iterator )
      , coefficient_matrices_( coefficient_matrices )
      , coefficient_matrices_4cork_( coefficient_matrices_ )
      , linear_solver_( linear_solver )
      , information_( information )
      {
        information_.linearization = std::string("CORK Linearization") ;
      }

    public:
      basis_handle_type const& basis_handle() const { return basis_handle_ ; }
      linear_solver_type const& linear_solver() const { return linear_solver_ ; }

    public:
      void shift( value_type s ) {
        basis_handle_.shift( s ) ;
        try {
          // Do this more accurately using M&N
          //linear_solver_.prepare_solve( basis_handle_ ) ;
          CORK::timer timer ;
          timer.tic() ;
          linear_solver_.prepare_solve( s ) ;
          information_.time_of_factorizations += timer.toc() ;
          ++information_.number_of_factorizations ;
        } catch (std::exception& err) {
          information_.error = err.what() ;
          throw exception::linear_solver_failure() ;
        }
      }

      value_type shift() const { return basis_handle_.shift() ; }

    public:

      //
      // Block system
      //
      // [    A - s B  ] y = [      B     ] b
      // [ (M-sN) x I  ]     [ I\otimes N ]
      //
      // Solved in two steps:
      // - solver_upper
      //   [      A - s B     ] z = [      B     ] b
      //   [  0 | ( M -s N)xI ]     [ I\otimes N ]
      // 
      // - solver_lower
      //   [  A(s) | 0  ] y = z
      //   [   *   | I  ]

      //
      // Solve Upper triangular part of linear system, in compact format
      //
      // The input vector is
      //  b = (I \otimes Q) U_in
      //
      // The output vector is
      //      (w_out )
      //  y = ( 0    ) + (I\otimes Q) U_out
      //      ( |    )
      //      ( 0    )
      // The first row of U_out will be zero.
      //
      template <typename QQ, typename UUIn, typename WOut, typename UUOut, typename Backend=glas2::default_backend>
      void solve_upper( QQ const& Q, UUIn const& U_in, WOut w_out, UUOut U_out, Backend const& backend=glas2::default_backend() ) {
        glas2::all all ;
        assert( U_in.num_columns()==Q.num_columns() ) ;
        assert( U_out.num_columns()==Q.num_columns() ) ;
        assert( U_in.num_rows()==basis_.size() ) ;
        assert( U_out.num_rows()==basis_.size() ) ;

        fill( U_out, 0.0 ) ;

        auto Zr = U_out( glas2::range_from_end(1,0), all ) ;
        //std::cout << "U_in before MN " << U_in << std::endl ;
        basis_.multiply_N( U_in, Zr ) ;
        basis_handle_.solve( Zr ) ;
      //std::cout << "Zr solve MN " << Zr << std::endl ;

        // Compute right-hand side for solve with P.
        glas2::fill( w_out, 0.0 ) ;

        coefficient_matrices_4cork_.initialize_schedule( Q, glas2::range(0,coefficient_matrices_.num_matrices() ) ) ;

        // Compute B-terms for the right-hand side.
        matrix_iterator_.schedule_b_0( coefficient_matrices_4cork_/*.accumulator()*/, [&U_in] (auto i) { return U_in(i, glas2::all()) ; } ) ;
        value_type s = shift() ;
        matrix_iterator_.schedule_b_1( coefficient_matrices_4cork_/*.accumulator()*/, [&Zr,&s] (auto i) { return s * Zr(i-1, glas2::all()) ; } ) ;

        // Compute A-terms for the right-hand side.
        matrix_iterator_.schedule_a_1( coefficient_matrices_4cork_/*.accumulator()*/, [&Zr] (auto i) { return -Zr(i-1, glas2::all()) ; } ) ;
        CORK::timer timer ;
        timer.tic() ;
        coefficient_matrices_4cork_.apply_scheduled( coefficient_matrices_, Q, w_out ) ;
        information_.time_of_matvecs += timer.toc() ;

        // Do the solve
        //std::cout << "W " << w_out << std::endl ;
        //std::cout << "Norm rhs " << norm_2(w_out) << std::endl ;
        timer.tic() ;
        linear_solver_.solve( w_out ) ;
        information_.time_of_solves += timer.toc() ;
        //std::cout << " -- norm sol " << norm_2(w_out) << std::endl ;
        ++information_.number_of_solves ;

        w_out *= basis_handle_.phi_0() ;
        fill( U_out( 0, glas2::all() ), 0.0 ) ;
      }  // solver_upper()


      template <typename UUOut>
      void solve_lower( UUOut U_out ) const {
        // Krylov step when shift and Q have the same value_type;
        // it is assumed U has the same value_type as Q
        //
        auto Zr = U_out( glas2::range_from_end(1,0), glas2::all() ) ;
        auto z_0 = U_out(0, glas2::all()) ;
        //std::cout << "Zr " << Zr << std::endl ;
        //std::cout << "z_0 " << z_0 << std::endl ;

        glas2::matrix< value_type > accum( Zr.num_rows(), Zr.num_columns() ) ;
        basis_handle_.lower_solve_right_hand_side( z_0, accum ) ;
        basis_handle_.solve( accum ) ;
        Zr -= accum ;
      } // solve_lower()


      //
      // Block system
      //
      // [    A - s B       ] u_out =  u_in
      // [ I \otimes (M-sN) ]         
      //
      // U must be column major
      //
      template <typename UU>
      void solve_without_B( UU U ) {
        glas2::all all ;
        assert( U.num_rows()==basis_.size() ) ;

        auto Zr = U( glas2::range_from_end(1,0), all ) ;
//        std::cout << "U_in before MN " << Ur << std::endl ;
        basis_handle_.solve( Zr ) ;
//      std::cout << "Zr solve MN " << Zr << std::endl ;

        auto w_out = U( 0, all ) ;

        // Compute right-hand side for solve with P.
        coefficient_matrices_4cork_.initialize_schedule( U, glas2::range(0,coefficient_matrices_.num_matrices() ) ) ;

        // Compute B-terms for the right-hand side.
        value_type s = shift() ;
        matrix_iterator_.schedule_b_1( coefficient_matrices_4cork_/*.accumulator()*/, [&Zr,&s] (auto i) { return s * Zr(i-1, glas2::all()) ; } ) ;

        // Compute A-terms for the right-hand side.
        matrix_iterator_.schedule_a_1( coefficient_matrices_4cork_/*.accumulator()*/, [&Zr] (auto i) { return -Zr(i-1, glas2::all()) ; } ) ;

        CORK::timer timer ;
        timer.tic() ;
        coefficient_matrices_4cork_.apply_scheduled( coefficient_matrices_, w_out ) ;
        information_.time_of_matvecs += timer.toc() ;

        // Do the solve
        //std::cout << "Norm w on input solver " << norm_2(w_out) << std::endl;
        // Temporary is needed because w_out does not have consecutive data locations.
        glas2::vector< value_type > w_temp( copy(w_out) ) ;
        timer.tic() ;
        linear_solver_.solve( w_temp ) ;
        information_.time_of_solves += timer.toc() ;
        w_out = basis_handle_.phi_0() * w_temp ;
        ++information_.number_of_solves ;

        glas2::matrix< value_type > accum( Zr.num_rows(), Zr.num_columns() ) ;
        basis_handle_.lower_solve_right_hand_side( w_out, accum ) ;
        basis_handle_.solve( accum ) ;
        Zr -= accum ;
      } // solve_without_B()

    private:
      Basis4CORK                       basis_ ;
      Basis4CORKHandle                 basis_handle_ ;
      MatrixIterator                   matrix_iterator_ ;
      CoefficientMatrices              coefficient_matrices_ ;
      coefficient_matrices_4cork_type  coefficient_matrices_4cork_ ;
      LinearSolver                     linear_solver_ ;
      info&                            information_ ;
  } ; // cork_linearization_solve_handle

} } // namespace CORK::linearization

#endif
