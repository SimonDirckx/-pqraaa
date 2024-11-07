//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_linearization_cork_linearization_krylov_step_handle_hpp
#define cork_linearization_cork_linearization_krylov_step_handle_hpp

#include <cork/linearization/info.hpp>
#include <cork/coefficient_matrices/coefficient_matrices4cork.hpp>
#include <glas2/vector.hpp>
#include <glas2/matrix.hpp>
#include <type_traits>
#include <exception>

namespace CORK { namespace linearization {

  template <typename LinearSolver, typename KrylovType, typename Basis4CORK, typename Basis4CORKHandle, typename MatrixIterator, typename CoefficientMatrices, typename EnableIf=void>
  class cork_linearization_krylov_step_handle
  {} ;


  template <typename LinearSolver, typename KrylovType, typename Basis4CORK, typename Basis4CORKHandle, typename MatrixIterator, typename CoefficientMatrices>
  class cork_linearization_krylov_step_handle< LinearSolver, KrylovType, Basis4CORK, Basis4CORKHandle, MatrixIterator, CoefficientMatrices
                                             , typename std::enable_if< std::is_same< typename std::decay<LinearSolver>::type::value_type, KrylovType >::value, void >::type
                                             >
  {
    public:
      typedef typename std::decay<Basis4CORK>::type                basis4cork_type ;
      typedef typename std::decay<MatrixIterator>::type            matrix_iterator_type ;
      typedef typename std::decay<CoefficientMatrices>::type       coefficient_matrices_type ;
      typedef typename std::decay<LinearSolver>::type              linear_solver_type ;

      typedef typename basis4cork_type::size_type                  size_type ;
      typedef typename linear_solver_type::value_type              value_type ;

    private:
      typedef decltype(std::abs(value_type()))                                                         real_value_type ;
      typedef coefficient_matrices::coefficient_matrices4CORK< value_type, coefficient_matrices_type > coefficient_matrices_4cork_type ;

    public:
      cork_linearization_krylov_step_handle( Basis4CORK poly, Basis4CORKHandle poly_handle, MatrixIterator matrix_iterator, CoefficientMatrices coefficient_matrices, LinearSolver linear_solver, info& information )
      : basis_( poly )
      , basis_handle_( poly_handle )
      , matrix_iterator_( matrix_iterator )
      , coefficient_matrices_( coefficient_matrices )
      , coefficient_matrices_4cork_( coefficient_matrices_ )
      , linear_solver_( linear_solver )
      , information_( information )
      {
        information_.linearization = "CORK Linearization" ;
      }

    public:
      void shift( value_type s ) {
        basis_handle_.shift( s ) ;
        try {
          // Do this more accurately using M&N
          //linear_solver_.prepare_solve( basis_handle_ ) ;
          linear_solver_.prepare_solve( s ) ;
          ++information_.number_of_factorizations ;
        } catch (std::exception& err) {
          information_.error = err.what() ;
          throw err ;
        }
      }

      value_type shift() const { return basis_handle_.shift() ; }

    public:

      // (A+sB)^{-1} (I x Q) vec(U)
/*      template <typename QQout, typename QQIn, typename UUOut, typename UUIn, typename R>
      void state_vector( QQIn const& Q_in, QQOut& Q_out, int rank_in, int& rank_out, UUIn coonst& U_in, UU U_out, R const& drop_tol ) {
        glas2::all all ;
        assert( Q_in.num_columns()>rank ) ;
        assert( Q_out.num_columns()>rank ) ;
        assert( U_in.num_columns()>=rank ) ;
        assert( U_out.num_columns()>=rank ) ;
        assert( U_in.num_rows()==basis_.size() ) ;
        assert( U_out.num_rows()==basis_.size() ) ;

        auto w = Q_out( all, rank ) ; // w points next to last column of Q
        {
          glas2::range r_rank_in(0,rank_in) ;
          auto Qr_in = Q_in( all, r_rank_in ) ;
          auto Ur_in = U_in( all, r_rank_in ) ;
          auto Zr = U_in( glas2::range_from_end(1,0), r_rank ) ;
//        std::cout << "Ur before MN " << Ur << std::endl ;
          basis_handle_.solve( Zr ) ;
//        std::cout << "Zr solve MN " << Zr << std::endl ;

          glas2::shared_vector<value_type> g(rank) ;

          // Compute right-hand side for solve with P.
          w = multiply( Qr_in, Ur_in ) ;

          coefficient_matrices_4cork_.initialize_schedule( Qr_in ) ;

          // Compute B-terms for the right-hand side.
          value_type s = shift() ;
          matrix_iterator_.schedule_b_1( coefficient_matrices_4cork_.accumulator(), [&Zr,&s] (auto i) { return s * Zr(i-1, glas2::all()) ; } ) ;

          // Compute A-terms for the right-hand side.
          matrix_iterator_.schedule_a( coefficient_matrices_4cork_.accumulator(), [&Zr] (auto i) { return -Zr(i-1, glas2::all()) ; } ) ;

          coefficient_matrices_4cork_.apply_scheduled( coefficient_matrices_, Qr, w ) ;

          // Do the solve
          //std::cout << "Norm w on input solver " << norm_2(w) << std::endl;
          linear_solver_.solve( w ) ;
          ++information_.number_of_solves ;
          auto norm_q = norm_2(w) ;
          //std::cout << "Norm w on output solver " << norm_q << std::endl;

          // Gram-Schmidt orthogonalization and possible expansion of Q.
          auto z_0 = U(0, r_rank) ; // First row of U will already be up to date.
          z_0 = multiply( transpose(conj(Qr)), w ) ;
          w -= multiply(Qr, z_0) ;
          g = multiply( transpose(conj(Qr)), w ) ;
          z_0 += g ;
          w -= multiply(Qr, g) ;
          auto norm_w = norm_2(w) ;
          //std::cout <<"norm_w " << norm_w << std::endl ;
          if (norm_w>drop_tol*norm_q) {
            w /= norm_w ;
            U_out(0,rank_out) = norm_w ;
            ++rank ;
          }
          //std::cout << Q.num_rows() << " " << rank << std::endl ;
        }
        assert( Q_out.num_columns()>=rank ) ;
        assert( U_out.num_columns()>=rank ) ;
        {
          glas2::range r_rank_out(0,rank_out) ;
          auto Zr_out = U_out( glas2::range_from_end(1,0), r_rank_out ) ;
          auto z_0_out = U_out(0, r_rank_out) ;

          glas2::matrix< value_type > accum( Zr_out.num_rows(), Zr_out.num_columns() ) ;
          basis_handle_.lower_solve_right_hand_side( z_0_out, accum ) ;
          basis_handle_.solve( accum ) ;
          Zr_out -= accum ;
        }
      } // state_vector()
      */

      //
      // Block system
      //
      // [    A - s B       ] y = [      B     ] b
      // [ I \otimes (M-sN) ]     [ I\otimes N ]
      //
      // Solved in two steps:
      // - solver_upper
      //   [  A - s B ] z = [      B     ] b
      //   [   I      ]     [ I\otimes N ]
      // 
      // - solver_lower
      //   [  A(s)  0 ... 0   ] y = z
      //   [ I \otimes (M-sN) ]

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
      template <typename QQ, typename UUIn, typename WOut, typename UUOut>
      void solve_upper( QQ const& Q, UUIn const& U_in, WOut w_out, UUOut U_out ) {
        glas2::all all ;
        assert( U_in.num_columns()==Q.num_columns() ) ;
        assert( U_out.num_columns()==Q.num_columns() ) ;
        assert( U_in.num_rows()==basis_.size() ) ;
        assert( U_out.num_rows()==basis_.size() ) ;

        fill( U_out, 0.0 ) ;

        glas2::range r_rank(0,rank) ;
        auto Qr = Q( all, r_rank ) ;
        auto Ur = U_in( all, r_rank ) ;
        auto Zr = U_out( glas2::range_from_end(1,0), r_rank ) ;
//        std::cout << "Ur before MN " << Ur << std::endl ;
        basis_.multiply_N( Ur, Zr ) ;
        basis_handle_.solve( Zr ) ;
//      std::cout << "Zr solve MN " << Zr << std::endl ;

        // Compute right-hand side for solve with P.
        glas2::fill( w, 0.0 ) ;

        coefficient_matrices_4cork_.initialize_schedule( Qr, glas2::range(1,coefficient_matrices_.num_matrices() ) ) ;

        // Compute B-terms for the right-hand side.
        matrix_iterator_.schedule_b_0( coefficient_matrices_4cork_.accumulator(), [&Ur] (auto i) { return Ur(i, glas2::all()) ; } ) ;
        value_type s = shift() ;
        matrix_iterator_.schedule_b_1( coefficient_matrices_4cork_.accumulator(), [&Zr,&s] (auto i) { return s * Zr(i-1, glas2::all()) ; } ) ;

        // Compute A-terms for the right-hand side.
        matrix_iterator_.schedule_a_1( coefficient_matrices_4cork_.accumulator(), [&Zr] (auto i) { return -Zr(i-1, glas2::all()) ; } ) ;

        coefficient_matrices_4cork_.apply_scheduled( coefficient_matrices_, Qr, w ) ;
//        std::cout << "w " << w << std::endl ;

//          w *= -1.0 ;

        // Do the solve
        //std::cout << "Norm w on input solver " << norm_2(w) << std::endl;
        linear_solver_.solve( w ) ;
        ++information_.number_of_solves ;

        fill( U_out( 0, glas2::all() ), 0.0 ) ;
      }  // solver_upper()


      template <typename UUOut>
      void solve_lower( UUOut U_out ) {
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
      template <typename QQ, typename UUIn, typename UUOut,typename R>
      void krylov_step( QQ& Q, int& rank, UUIn const& U_in, UUOut U_out, R const& drop_tol ) {
        // Krylov step when shift and Q have the same value_type;
        // it is assumed U has the same value_type as Q
        //
        glas2::all all ;
        assert( Q.num_columns()>rank ) ;
        assert( U_in.num_columns()>=rank ) ;
        assert( U_out.num_columns()>=rank ) ;
        assert( U_in.num_rows()==basis_.size() ) ;
        assert( U_out.num_rows()==basis_.size() ) ;

        fill( U_out, 0.0 ) ;

        auto w = Q( all, rank ) ; // w points next to last column of Q
        {
          glas2::range r_rank(0,rank) ;
          auto Qr = Q( all, r_rank ) ;
          auto Ur = U_in( all, r_rank ) ;
          auto Zr = U_out( glas2::range_from_end(1,0), r_rank ) ;
//        std::cout << "Ur before MN " << Ur << std::endl ;
          basis_.multiply_N( Ur, Zr ) ;
          basis_handle_.solve( Zr ) ;
//        std::cout << "Zr solve MN " << Zr << std::endl ;

          glas2::shared_vector<value_type> g(rank) ;

          // Compute right-hand side for solve with P.
          glas2::fill( w, 0.0 ) ;

          coefficient_matrices_4cork_.initialize_schedule( Qr, glas2::range(1,coefficient_matrices_.num_matrices() ) ) ;

          // Compute B-terms for the right-hand side.
          matrix_iterator_.schedule_b_0( coefficient_matrices_4cork_.accumulator(), [&Ur] (auto i) { return Ur(i, glas2::all()) ; } ) ;
          value_type s = shift() ;
          matrix_iterator_.schedule_b_1( coefficient_matrices_4cork_.accumulator(), [&Zr,&s] (auto i) { return s * Zr(i-1, glas2::all()) ; } ) ;

          // Compute A-terms for the right-hand side.
          matrix_iterator_.schedule_a_1( coefficient_matrices_4cork_.accumulator(), [&Zr] (auto i) { return -Zr(i-1, glas2::all()) ; } ) ;

          coefficient_matrices_4cork_.apply_scheduled( coefficient_matrices_, Qr, w ) ;
//          std::cout << "w " << w << std::endl ;

//          w *= -1.0 ;

          // Do the solve
          //std::cout << "Norm w on input solver " << norm_2(w) << std::endl;
          linear_solver_.solve( w ) ;
          ++information_.number_of_solves ;
          auto norm_q = norm_2(w) ;
          //std::cout << "Norm w on output solver " << norm_q << std::endl;

          // Gram-Schmidt orthogonalization and possible expansion of Q.
          auto z_0 = U_out(0, r_rank) ; // First row of U will already be up to date.
          z_0 = multiply( transpose(conj(Qr)), w ) ;
          w -= multiply(Qr, z_0) ;
          g = multiply( transpose(conj(Qr)), w ) ;
          z_0 += g ;
          w -= multiply(Qr, g) ;
          auto norm_w = norm_2(w) ;
          //std::cout <<"norm_w " << norm_w << std::endl ;
          if (norm_w>drop_tol*norm_q) {
            w /= norm_w ;
            U_out(0,rank) = norm_w ;
            ++rank ;
          }
          //std::cout << Q.num_rows() << " " << rank << std::endl ;
        }
        assert( Q.num_columns()>=rank ) ;
        assert( U_out.num_columns()>=rank ) ;
        {
          glas2::range r_rank(0,rank) ;
          auto Zr = U_out( glas2::range_from_end(1,0), r_rank ) ;
          auto z_0 = U_out(0, r_rank) ;
        //std::cout << "Zr " << Zr << std::endl ;
        //std::cout << "z_0 " << z_0 << std::endl ;

          glas2::matrix< value_type > accum( Zr.num_rows(), Zr.num_columns() ) ;
          basis_handle_.lower_solve_right_hand_side( z_0, accum ) ;
          basis_handle_.solve( accum ) ;
          Zr -= accum ;
        //std::cout << "Zr accum " << Zr << std::endl ;
        }
      } // krylov_step()

      /*template <typename QQ, typename UUIn, typename UUOut,typename R>
      typename std::enable_if< std::is_same< std::complex< typename QQ::value_type >
                                                  , value_type
                                                  >::value >::type krylov_step( QQ& Q, int& rank, UUIn const& U_in, UUOut U_out, R const& drop_tol ) {
        // Krylov step when the shift is complex and Q is real.
        // It is assumed U has the value_type "value_type".
        //
        typedef typename QQ::value_type q_value_type ;

        glas2::all all ;
        assert( Q.num_columns()>rank ) ;
        assert( U_in.num_columns()>=rank ) ;
        assert( U_out.num_columns()>=rank ) ;
        assert( U_in.num_rows()==basis_.size() ) ;
        assert( U_out.num_rows()==basis_.size() ) ;

        glas2::vector< value_type > ws( Q.num_rows() ) ;

        fill( U_out, 0.0 ) ;

        //auto w = Q( all, rank ) ; // w points next to last column of Q
        {
          glas2::range r_rank(0,rank) ;
          auto Qr = Q( all, r_rank ) ;
          auto Ur = U_in( all, r_rank ) ;
          auto Zr = U_out( glas2::range_from_end(1,0), r_rank ) ;
//        std::cout << "Ur before MN " << Ur << std::endl ;
          basis_.multiply_N( Ur, Zr ) ;
          basis_handle_.solve( Zr ) ;
//        std::cout << "Zr solve MN " << Zr << std::endl ;

          // Compute right-hand side for solve with P.
          glas2::fill( ws, 0.0 ) ;

          coefficient_matrices_4cork_.initialize_schedule( Qr, glas2::range(1,coefficient_matrices_.num_matrices() ) ) ;

          // Compute B-terms for the right-hand side.
          matrix_iterator_.schedule_b_0( coefficient_matrices_4cork_.accumulator(), [&Ur] (auto i) { return Ur(i, glas2::all()) ; } ) ;
          value_type s = shift() ;
          matrix_iterator_.schedule_b_1( coefficient_matrices_4cork_.accumulator(), [&Zr,&s] (auto i) { return s * Zr(i-1, glas2::all()) ; } ) ;

          // Compute A-terms for the right-hand side.
          matrix_iterator_.schedule_a_1( coefficient_matrices_4cork_.accumulator(), [&Zr] (auto i) { return -Zr(i-1, glas2::all()) ; } ) ;

          coefficient_matrices_4cork_.apply_scheduled( coefficient_matrices_, Qr, ws ) ;
//          std::cout << "ws " << ws << std::endl ;

//          ws *= -1.0 ;

          // Do the solve
          //std::cout << "Norm ws on input solver " << norm_2(ws) << std::endl;
          linear_solver_.solve( ws ) ;
          ++information_.number_of_solves ;
          //std::cout << "Norm w on output solver " << norm_q << std::endl;

          // Gram-Schmidt orthogonalization and possible expansion of Q.
          auto norm_q = norm_2( real(ws) ) ;

          // Orthonormalize real part
          auto w_r = Q(all,rank) ;
          glas2::vector<q_value_type> g(rank) ;
          glas2::matrix<q_value_type> H(rank+2,2) ;
          fill(H,0.0) ;

          auto h_r = H(glas2::range(0,rank),0) ;
          w_r = real( ws ) ;
          h_r = multiply( transpose(Qr), w_r ) ;
          w_r -= multiply(Qr, h_r ) ;
          g = multiply( transpose(conj(Qr)), w_r ) ;
          h_r += g ;
          w_r -= multiply(Qr, g ) ;
          auto norm_w_r = norm_2(w_r) ;
          //std::cout <<"norm_w " << norm_w << std::endl ;
          if (norm_w_r>drop_tol*norm_q) {
            w_r /= norm_w_r ;
            H(rank,0) = norm_w_r ;
            ++rank ;
          }
          auto Qr2 = Q( all, glas2::range(0,rank) ) ;

          // Orthonormalize imaginary part
          auto w_i = Q(all,rank) ;
          glas2::vector<q_value_type> gi(rank) ;

          auto h_i = H(glas2::range(0,rank),1) ;
          w_i = imag( ws ) ;
          h_i = multiply( transpose(Qr2), w_i ) ;
          w_i -= multiply(Qr2, h_i ) ;
          gi = multiply( transpose(Qr2), w_i ) ;
          h_i += gi ;
          w_i -= multiply(Qr2, gi ) ;
          auto norm_w_i = norm_2(w_i) ;
          //std::cout <<"norm_w " << norm_w << std::endl ;
          if (norm_w_i>drop_tol*norm_q) {
            w_i /= norm_w_i ;
            H(rank,1) = norm_w_i ;
            ++rank ;
          }

          // Update U_out_c with the new Q
          U_out(0,glas2::range(0,rank)) = H(glas2::range(0,rank),0) + value_type(0.,1.0) * H(glas2::range(0,rank),1) ;
        }
        assert( Q.num_columns()>=rank ) ;
        assert( U_out.num_columns()>=rank ) ;
        {
          glas2::range r_rank(0,rank) ;
          auto Zr = U_out( glas2::range_from_end(1,0), r_rank ) ;
          auto z_0 = U_out(0, r_rank) ;
        //std::cout << "Zr " << Zr << std::endl ;
        //std::cout << "z_0 " << z_0 << std::endl ;

          glas2::matrix< value_type > accum( Zr.num_rows(), Zr.num_columns() ) ;
          basis_handle_.lower_solve_right_hand_side( z_0, accum ) ;
          basis_handle_.solve( accum ) ;
          Zr -= accum ;
        //std::cout << "Zr accum " << Zr << std::endl ;
        }
      } // krylov_step()*/

    public:
      template <typename UUIn, typename UUOut>
      void full_krylov_step( UUIn const& U_in, UUOut U_out ) {
        glas2::all all ;
        assert( U_in.num_columns()==coefficient_matrices_.size() ) ;
        assert( U_out.num_columns()==coefficient_matrices_.size() ) ;
        assert( U_in.num_rows()==basis_.size() ) ;
        assert( U_out.num_rows()==basis_.size() ) ;

        auto w = U_out( 0, all ) ;
        {
          auto Zr = U_out( glas2::range_from_end(1,0), all ) ;
          basis_.multiply_N( U_in, Zr ) ;
          basis_handle_.solve( Zr ) ;

          // Compute right-hand side for solve with P.
          glas2::fill( w, 0.0 ) ;

          coefficient_matrices_4cork_.initialize_schedule( transpose(U_in), glas2::range(1,coefficient_matrices_.num_matrices() ) ) ;

          // Compute B-terms for the right-hand side.
          matrix_iterator_.schedule_b_0( coefficient_matrices_4cork_.accumulator(), [&U_in] (auto i) { return U_in(i, glas2::all()) ; } ) ;
          value_type s = shift() ;
          matrix_iterator_.schedule_b_1( coefficient_matrices_4cork_.accumulator(), [&Zr,&s] (auto i) { return s * Zr(i-1, glas2::all()) ; } ) ;

          // Compute A-terms for the right-hand side.
          matrix_iterator_.schedule_a_1( coefficient_matrices_4cork_.accumulator(), [&Zr] (auto i) { return -Zr(i-1, glas2::all()) ; } ) ;

          coefficient_matrices_4cork_.apply_schedule( coefficient_matrices_, glas2::eye(coefficient_matrices_.size(),coefficient_matrices_.size()), w ) ;

//          w *= -1.0 ;

          // Do the solve
          linear_solver_.solve( w ) ;
          ++information_.number_of_solves ;
          auto norm_q = norm_2(w) ;
        }
        {
          auto Zr = U_out( glas2::range_from_end(1,0), all ) ;
          auto z_0 = U_out(0, all) ;

          glas2::matrix< value_type > accum( Zr.num_rows(), Zr.num_columns() ) ;
          basis_handle_.lower_solve_right_hand_side( z_0, accum ) ;
          basis_handle_.solve( accum ) ;
          Zr -= accum ;
        }
      } // full_krylov_step()

    private:
      Basis4CORK                       basis_ ;
      Basis4CORKHandle                 basis_handle_ ;
      MatrixIterator                   matrix_iterator_ ;
      CoefficientMatrices              coefficient_matrices_ ;
      coefficient_matrices_4cork_type  coefficient_matrices_4cork_ ;
      LinearSolver                     linear_solver_ ;
      info&                            information_ ;
  } ; // cork_linearization_krylov_step_handle


  template <typename LinearSolver, typename KrylovType, typename Basis4CORK, typename Basis4CORKHandle, typename MatrixIterator, typename CoefficientMatrices>
  class cork_linearization_krylov_step_handle< LinearSolver, KrylovType, Basis4CORK, Basis4CORKHandle, MatrixIterator, CoefficientMatrices
                                             , typename std::enable_if< std::is_same< typename std::decay<LinearSolver>::type::value_type, std::complex<KrylovType> >::value, void >::type
                                             >
  {
    public:
      typedef typename std::decay<Basis4CORK>::type                basis4cork_type ;
      typedef typename std::decay<MatrixIterator>::type            matrix_iterator_type ;
      typedef typename std::decay<CoefficientMatrices>::type       coefficient_matrices_type ;
      typedef typename std::decay<LinearSolver>::type              linear_solver_type ;

      typedef typename basis4cork_type::size_type                  size_type ;
      typedef typename linear_solver_type::value_type              value_type ;

    private:
      typedef decltype(std::abs(value_type()))                                                         real_value_type ;
      typedef coefficient_matrices::coefficient_matrices4CORK< value_type, coefficient_matrices_type > coefficient_matrices_4cork_type ;

    public:
      cork_linearization_krylov_step_handle( Basis4CORK poly, Basis4CORKHandle poly_handle, MatrixIterator matrix_iterator, CoefficientMatrices coefficient_matrices, LinearSolver linear_solver, info& information )
      : basis_( poly )
      , basis_handle_( poly_handle )
      , matrix_iterator_( matrix_iterator )
      , coefficient_matrices_( coefficient_matrices )
      , coefficient_matrices_4cork_( coefficient_matrices_ )
      , linear_solver_( linear_solver )
      , information_( information )
      {
        information_.linearization = "CORK Linearization" ;
      }

    public:
      void shift( value_type s ) {
        basis_handle_.shift( s ) ;
        try {
          // Do this more accurately using M&N
          //linear_solver_.prepare_solve( basis_handle_ ) ;
          linear_solver_.prepare_solve( s ) ;
          ++information_.number_of_factorizations ;
        } catch (std::exception& err) {
          information_.error = err.what() ;
          throw err ;
        }
      }

      value_type shift() const { return basis_handle_.shift() ; }

    public:
      template <typename QQ, typename UUIn, typename UUOut,typename R>
      void krylov_step( QQ& Q, int& rank, UUIn const& U_in, UUOut U_out, R const& drop_tol ) {
        static_assert( std::is_same< std::complex<typename QQ::value_type>, value_type>::value, "CORK::linearization::cork_linearization_krylov_step_handle: value_type of Q does not match") ;

        // Krylov step when the shift is complex and Q is real.
        // It is assumed U has the value_type "value_type".
        //
        typedef typename QQ::value_type q_value_type ;

        glas2::all all ;
        assert( Q.num_columns()>rank ) ;
        assert( U_in.num_columns()>=rank ) ;
        assert( U_out.num_columns()>=rank ) ;
        assert( U_in.num_rows()==basis_.size() ) ;
        assert( U_out.num_rows()==basis_.size() ) ;

        glas2::vector< value_type > ws( Q.num_rows() ) ;

        fill( U_out, 0.0 ) ;

        //auto w = Q( all, rank ) ; // w points next to last column of Q
        {
          glas2::range r_rank(0,rank) ;
          auto Qr = Q( all, r_rank ) ;
          auto Ur = U_in( all, r_rank ) ;
          auto Zr = U_out( glas2::range_from_end(1,0), r_rank ) ;
//        std::cout << "Ur before MN " << Ur << std::endl ;
          basis_.multiply_N( Ur, Zr ) ;
          basis_.solve( Zr ) ;
//        std::cout << "Zr solve MN " << Zr << std::endl ;

          // Compute right-hand side for solve with P.
          glas2::fill( ws, 0.0 ) ;

          coefficient_matrices_4cork_.initialize_schedule( Qr, glas2::range(1,coefficient_matrices_.num_matrices() ) ) ;

          // Compute B-terms for the right-hand side.
          matrix_iterator_.schedule_b_0( coefficient_matrices_4cork_.accumulator(), [&Ur] (auto i) { return Ur(i, glas2::all()) ; } ) ;
          value_type s = shift() ;
          matrix_iterator_.schedule_b_1( coefficient_matrices_4cork_.accumulator(), [&Zr,&s] (auto i) { return s * Zr(i-1, glas2::all()) ; } ) ;

          // Compute A-terms for the right-hand side.
          matrix_iterator_.schedule_a_1( coefficient_matrices_4cork_.accumulator(), [&Zr] (auto i) { return -Zr(i-1, glas2::all()) ; } ) ;

          coefficient_matrices_4cork_.apply_scheduled( coefficient_matrices_, Qr, ws ) ;
//          std::cout << "ws " << ws << std::endl ;

//          ws *= -1.0 ;

          // Do the solve
          //std::cout << "Norm ws on input solver " << norm_2(ws) << std::endl;
          linear_solver_.solve( ws ) ;
          ++information_.number_of_solves ;
          //std::cout << "Norm w on output solver " << norm_q << std::endl;

          // Gram-Schmidt orthogonalization and possible expansion of Q.
          auto norm_q = norm_2( real(ws) ) ;

          // Orthonormalize real part
          auto w_r = Q(all,rank) ;
          glas2::vector<q_value_type> g(rank) ;
          glas2::matrix<q_value_type> H(rank+2,2) ;
          fill(H,0.0) ;

          auto h_r = H(glas2::range(0,rank),0) ;
          w_r = real( ws ) ;
          h_r = multiply( transpose(Qr), w_r ) ;
          w_r -= multiply(Qr, h_r ) ;
          g = multiply( transpose(conj(Qr)), w_r ) ;
          h_r += g ;
          w_r -= multiply(Qr, g ) ;
          auto norm_w_r = norm_2(w_r) ;
          //std::cout <<"norm_w " << norm_w << std::endl ;
          if (norm_w_r>drop_tol*norm_q) {
            w_r /= norm_w_r ;
            H(rank,0) = norm_w_r ;
            ++rank ;
          }
          auto Qr2 = Q( all, glas2::range(0,rank) ) ;

          // Orthonormalize imaginary part
          auto w_i = Q(all,rank) ;
          glas2::vector<q_value_type> gi(rank) ;

          auto h_i = H(glas2::range(0,rank),1) ;
          w_i = imag( ws ) ;
          h_i = multiply( transpose(Qr2), w_i ) ;
          w_i -= multiply(Qr2, h_i ) ;
          gi = multiply( transpose(Qr2), w_i ) ;
          h_i += gi ;
          w_i -= multiply(Qr2, gi ) ;
          auto norm_w_i = norm_2(w_i) ;
          //std::cout <<"norm_w " << norm_w << std::endl ;
          if (norm_w_i>drop_tol*norm_q) {
            w_i /= norm_w_i ;
            H(rank,1) = norm_w_i ;
            ++rank ;
          }

          // Update U_out_c with the new Q
          U_out(0,glas2::range(0,rank)) = H(glas2::range(0,rank),0) + value_type(0.,1.0) * H(glas2::range(0,rank),1) ;
        }
        assert( Q.num_columns()>=rank ) ;
        assert( U_out.num_columns()>=rank ) ;
        {
          glas2::range r_rank(0,rank) ;
          auto Zr = U_out( glas2::range_from_end(1,0), r_rank ) ;
          auto z_0 = U_out(0, r_rank) ;
        //std::cout << "Zr " << Zr << std::endl ;
        //std::cout << "z_0 " << z_0 << std::endl ;

          glas2::matrix< value_type > accum( Zr.num_rows(), Zr.num_columns() ) ;
          basis_.lower_solve_right_hand_side( z_0, accum ) ;
          basis_.solve( accum ) ;
          Zr -= accum ;
        //std::cout << "Zr accum " << Zr << std::endl ;
        }
      } // krylov_step_real()

    private:
      Basis4CORK                       basis_ ;
      Basis4CORKHandle                 basis_handle_ ;
      MatrixIterator                   matrix_iterator_ ;
      CoefficientMatrices              coefficient_matrices_ ;
      coefficient_matrices_4cork_type  coefficient_matrices_4cork_ ;
      LinearSolver                     linear_solver_ ;
      info&                            information_ ;
  } ; // cork_linearization_krylov_step_handle

} } // namespace CORK::linearization

#endif
