//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2015.
//
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_krylov_cork_process_hpp
#define cork_krylov_cork_process_hpp

#include <cork/options/absolute_shift_modifier.hpp>
#include <cork/options/relative_shift_modifier.hpp>
#include <cork/options/max_solver_growth_factor.hpp>
#include <cork/options/krylov_breakdown_tolerance.hpp>
#include <cork/options/debug_level.hpp>
#include <cork/options/backend.hpp>
#include <cork/krylov/info.hpp>
#include <cork/krylov/sequence_of_solve_handle.hpp>
#include <cork/krylov/cork_quadruple.hpp>
#include <cork/exception/linear_solver_failure.hpp>
#include <cork/exception/linear_solver_failure_new_shift.hpp>
#include <cmath>
#include <type_traits>

namespace CORK { namespace krylov {

  struct general_tag {} ;

  template <typename Tag, typename ShiftValueType, typename Linearization, typename ShiftGenerator, typename Quadruple, typename Options>
  class cork_process
  {} ;


  template <typename ShiftValueType, typename Linearization, typename ShiftGenerator, typename Quadruple, typename Options>
  class cork_process< general_tag, ShiftValueType, Linearization, ShiftGenerator, Quadruple, Options>
  {
    public:
      typedef typename ShiftGenerator::value_type                            shift_value_type ;
      typedef typename Quadruple::value_type                                 value_type ;
      typedef decltype( std::abs(value_type()) )                             real_type ;
      typedef typename Linearization::template solve_handle_type<value_type> solve_handle_type ;

      static_assert( std::is_convertible< shift_value_type, value_type >::value, "CORK::krylov::cork_process: ShiftValueType is not convertible to the value_type of Quadruple" ) ;

    public:
      cork_process( Linearization& linearization, ShiftGenerator& shift_generator, Quadruple& quadruple, Options const& options )
      : quadruple( quadruple )
      , linearization_( linearization )
      , shift_generator_( shift_generator )
      , solve_handles_( linearization_ )
      , options_( options )
      {
        static_cast<info&>(linearization_.information()).krylov_process = "CORK standard" ;
      }

    public:
      void initialize( shift_value_type shift ) {
        std::cout <<" WHY DO WE DO THIS?" << std::endl ;
        try {
          if (options::value_of<options::debug_level>(options_)>2)
            std::cout << "Shift in rational Krylov = " << shift << std::endl ;
          solve_handles_.shift( shift ) ;
        } catch (exception::linear_solver_failure& e) {
          shift += options::value_of< options::absolute_shift_modifier<shift_value_type> >(options_) + shift * options::value_of< options::relative_shift_modifier<shift_value_type>>(options_) ;
          linearization_.information().warnings.push_back( "A shift has changes because of option max_solver_growth_factor" ) ;
          if (options::value_of<options::debug_level>(options_)>2)
            std::cout << "Shift = " << shift << std::endl ;
          solve_handles_.shift( shift ) ;
        }
      } // initialize()

      template <typename InitialVector>
      void initial_vector( InitialVector const& initial ) {
        initial.apply( quadruple, solve_handles_, options::value_of<options::backend_key>(options_) ) ;
      } // initialize()

    public:
      void expand( int max_k ) {
        norm_h_ = norm_fro_squared( quadruple.H( glas2::range(0,quadruple.k()+1), glas2::range(0,quadruple.k()) ) ) ;

        // Make sure we do at least one iteration, otherwise there is not advance.
        assert( quadruple.k()<max_k ) ;

        for (; quadruple.k()<max_k;) {
          shift_value_type shift = shift_generator_.next_shift() ;//quadruple.k()) ;

          int old_rank = quadruple.rank ;

          solve_handles_.select( shift_generator_.solve_handle() ) ;
          real_type growth = options::value_of<options::max_solver_growth_factor<real_type>>(options_) ;
          for (int solve_trials=0; solve_trials<2 && growth>=options::value_of<options::max_solver_growth_factor<real_type>>(options_); ++solve_trials) {
            quadruple.rank = old_rank ;
            try {
            try {
              if (options::value_of<options::debug_level>(options_)>2)
                std::cout << "Shift = " << shift << std::endl ;
              solve_handles_.shift( shift ) ;
              if (options::value_of<options::debug_level>(options_)>1) std::cout << "Iteration " << quadruple.k() << "  shift = " << shift << std::endl ;

              //shift_generator.continuation_combination( quadruple.k(), quadruple.continuation_combination( glas2::range(0,quadruple.k()+1), quadruple.k() ) ) ;

              // k() is the iteration count starting from 0
              // Solve takes place in three steps:
              // 1) Solve the upper part of the U-L factorization
              solve_handles_.solve_upper( quadruple.Q_k(), quadruple.u_vector_k( quadruple.k() ), quadruple.Q(glas2::all(), quadruple.rank), quadruple.u_vector_k( quadruple.k()+1 ), options::value_of<options::backend_key>(options_) ) ;
              // 2) Update the Q and U matrices: orthonormalize column rank of Q. Modify U accordingly.
              fill( quadruple.u_vector( quadruple.k()+1 )(0, glas2::all()), 0.0) ;
              quadruple.u_vector( quadruple.k()+1 )(0,quadruple.rank) = 1.0 ;
              // Add column Q(:,rank) and update U(:,k+1) ;
              quadruple.add_column_of_Q( quadruple.k()+1 ) ;
              // 3) Solve for the lower part.
              solve_handles_.solve_lower( quadruple.u_vector_k( quadruple.k()+1 ) ) ;

              growth = norm_fro( quadruple.u_vector_k( quadruple.k()+1 ) ) ;
              if (growth>=options::value_of<options::max_solver_growth_factor<real_type>>(options_)) {
                shift += options::value_of< options::absolute_shift_modifier<shift_value_type> >(options_) + shift * options::value_of< options::relative_shift_modifier<shift_value_type>>(options_) ;
                shift_generator_.modify_shift( shift ) ;
                if (options::value_of<options::debug_level>(options_)>0) std::cout << "CORK: A shift has changed because of option max_solver_growth_factor" << std::endl ;
                linearization_.information().warnings.push_back( "A shift has changed because of option max_solver_growth_factor" ) ;
              }
            } catch (exception::linear_solver_failure& e) {
              shift += options::value_of< options::absolute_shift_modifier<shift_value_type> >(options_) + shift * options::value_of< options::relative_shift_modifier<shift_value_type>>(options_) ;
              shift_generator_.modify_shift( shift ) ;
              if (options::value_of<options::debug_level>(options_)>0) std::cout << "CORK: A shift has changed because of CORK::exception::linear_solver_failure" << std::endl ;
              linearization_.information().warnings.push_back( "A shift has changed because of CORK::exception::linear_solver_failure" ) ;
              if (solve_trials==1) throw e ;
            }
            } catch (exception::linear_solver_failure_new_shift<shift_value_type>& e) {
              shift = e.shift() ;
              shift_generator_.modify_shift( shift ) ;
              if (options::value_of<options::debug_level>(options_)>0) std::cout << "CORK: A shift has changed because of CORK::exception::linear_solver_failure_new_shift" << std::endl ;
              linearization_.information().warnings.push_back( "A shift has changed because of CORK::exception::linear_solver_failure_new_shift" ) ;
              if (solve_trials==1) throw e ;
            }
          } // For solve_trials
          if (growth>=options::value_of<options::max_solver_growth_factor<real_type>>(options_)) {
            linearization_.information().warnings.push_back( "Growth of iteration vectors is too high: results may not be accurate" ) ;
          }


          /*{
            auto multiply_handle = linearization_.template multiply_handle<value_type>() ;
            glas2::range r_rank(0, new_rank) ;
            auto Q = quadruple.Q( glas2::all(), r_rank ) ;
            glas2::vector<value_type> w_rhs(Q.num_rows()) ;
            glas2::vector<value_type> w_AB(Q.num_rows()) ;
            glas2::matrix<value_type> U_rhs(quadruple.degree(), r_rank.size()) ;
            glas2::matrix<value_type> U_AB(quadruple.degree(), r_rank.size()) ;
            auto u_k( quadruple.u_vector( quadruple.k() ) ( glas2::all(), r_rank ) ) ;
            auto u_k1( quadruple.u_vector( quadruple.k()+1 ) ( glas2::all(), r_rank ) ) ;

            multiply_handle.multiply_B( Q, u_k, w_rhs, U_rhs ) ;
            multiply_handle.multiply_AB( 1.0, -shift, Q, u_k1, w_AB, U_AB ) ;
            std::cout << "Residual linear system, first component : " << norm_2(w_AB-w_rhs) << " / " << norm_2(w_rhs) << std::endl ;
            std::cout << "Residual linear system, other components: " << norm_fro(U_AB-U_rhs) << " / " << norm_fro(U_rhs) << std::endl ;
          }*/

            // Orthonormalize vector in Krylov space
          int new_rank = quadruple.rank ;
          quadruple.next_step( new_rank, shift ) ;
          // k() has been increased by one: is now the iteration count starting from 1.

          // Check for lucky breakdown
          norm_h_ += norm_2_squared( quadruple.H( glas2::range(0,quadruple.k()), quadruple.k()-1 ) ) ;
          if (glas2::abs_squared(quadruple.H(quadruple.k(),quadruple.k()-1))<std::pow( options::value_of<options::krylov_breakdown_tolerance<real_type>>(options_),2)*norm_h_) {
            std::cout << "Breakdown bug in cork_expand() at iteration " << quadruple.k() << std::endl ;
            quadruple.back_track() ;
            break;
          }
          if (options::value_of<options::debug_level>(options_)>2) std::cout << "Iteration " << quadruple.k() << "  rank = " << old_rank << " --> " << quadruple.rank << std::endl ;
        } // end for
      } // expand()

/*      void expand_double( int max_k ) {
        typedef typename Linearization::value_type L_value_type ;

        norm_h_ = norm_fro_squared( quadruple.H( glas2::range(0,quadruple.k()+1), glas2::range(0,quadruple.k()) ) ) ;

        // Make sure we do at least one iteration, otherwise there is not advance.
        assert( quadruple.k()<max_k ) ;

        for (; quadruple.k()<max_k;) {
          int new_rank = quadruple.rank ;

          auto shift = shift_generator_.next_shift() ;//quadruple.k()) ;
          solve_handles_.select( shift_generator_.solve_handle() ) ;
          //linearization.select( shift_generator.solve_handle(quadruple.k() ) ) ;
          if (options::value_of<options::debug_level>(options_)>2)
            std::cout << "Shift = " << shift << std::endl ;
          solve_handles_.shift( shift ) ;
          quadruple.poles(quadruple.k()) = shift ;
          if (options_.debug_level>1) std::cout << "Iteration " << quadruple.k() << "  shift = " << shift << std::endl ;

          //shift_generator.continuation_combination( quadruple.k(), quadruple.continuation_combination( glas2::range(0,quadruple.k()+1), quadruple.k() ) ) ;

          // k() is the iteration count starting from 0
          auto u_k1 = quadruple.u_vector( quadruple.k()+1 ) ;
          glas2::matrix< L_value_type, glas2::row_major > u_k1_c( u_k1.num_rows(), u_k1.num_columns() ) ;
          u_k1_c = u_k1 ;
          solve_handles_.krylov_step( quadruple.Q, new_rank, quadruple.u_vector( quadruple.k() ), u_k1_c, options_.Q_drop_tol ) ;
          u_k1 = real( u_k1_c ) ;
          if ( shift.imag()!=0.0 ) {
            quadruple.u_vector( quadruple.k()+2 ) = imag( u_k1_c ) ;
            quadruple.next_step_double( new_rank, shift ) ;
          } else {
            quadruple.next_step( new_rank, real(shift) ) ;
          }

          // k() has been increased by one: is now the iteration count starting from 1.

          // Check for lucky breakdown
          norm_h_ += norm_2_squared( quadruple.H( glas2::range(0,quadruple.k()), quadruple.k()-1 ) ) ;
          if (glas2::abs_squared(quadruple.H(quadruple.k(),quadruple.k()-1))<options_.relative_tolerance*options_.relative_tolerance*norm_h_) {
            std::cout << "Breakdown bug in cork_expand() at iteration " << quadruple.k() << std::endl ;
            break;
          }
          if (options_.debug_level>2) std::cout << "Iteration " << quadruple.k() << "  rank = " << new_rank << " " << quadruple.rank << std::endl ;
        } // end for
      } // expand_double()*/

/*      real_type verify_recurrence_with_solve() {
        glas2::shared_vector< value_type > expand_vi( linearization_.size() * linearization_.size_of_basis() ) ;
        glas2::shared_vector< value_type > expand_vj( linearization_.size() * linearization_.size_of_basis() ) ;
        // Check orthogonality of the iteration vectors
        glas2::range range_rank( 0, quadruple.rank ) ;

        real_type sum;
        for (int i=0; i<quadruple.k()+1; ++i) {
          auto Ui = quadruple.u_vector(i) ;
          reshape( expand_vi.pass_ref(), linearization_.size(), linearization_.size_of_basis(), glas2::column_major() ) = 
            multiply( quadruple.Q( glas2::all(), range_rank ), transpose( Ui( glas2::all(), range_rank ) ) ) ;
          sum = std::abs( norm_2( expand_vi )-1.0 ) ;
          for (int j=0; j<quadruple.k()+1; ++j) {
            if (i!=j) {
              auto Uj = quadruple.u_vector(j) ;
              reshape( expand_vj.pass_ref(), linearization_.size(), linearization_.size_of_basis(), glas2::column_major() ) = 
                multiply( quadruple.Q( glas2::all(), range_rank ), transpose( Uj( glas2::all(), range_rank ) ) ) ;
              sum += std::abs( inner_prod( conj(expand_vj), expand_vi ) ) ;
            }
          }
        }
        std::cout << "Orthogonality of " << quadruple.k()+1 << " vectors: " << sum << std::endl ;

        if (options_.debug_level>1) std::cout << "Verify recurrence: " << std::endl ;
        glas2::shared_matrix<value_type> u_vector( linearization_.size_of_basis(), quadruple.rank+1 ) ; 
        auto const& Q( quadruple.Q ) ; 
        glas2::shared_vector< value_type > w( quadruple.Q.num_rows() ) ;

        auto shift = solve_handles_.shift() ;

        // UKH holds U * (K-\sigma H) with sigma the last shift
        // Continuation combination is I, so K - sigma * H = H * (Poles-sigma*I) + I
        glas2::shared_matrix<value_type> UKH( quadruple.U.num_rows(), quadruple.k() ) ; 
        glas2::range range_k( 0, quadruple.k() ) ;
        glas2::range range_k1( 0, quadruple.k()+1 ) ;
        UKH = multiply( quadruple.U( glas2::all(), range_k1 ), quadruple.K(range_k1, range_k) ) ;
        UKH -= shift * multiply( quadruple.U( glas2::all(), range_k1 ), quadruple.H(range_k1, range_k) ) ;
        std::cout << quadruple.K(quadruple.k(), range_k) - shift*quadruple.H(quadruple.k(), range_k) << std::endl;
        assert( norm_2( quadruple.K(quadruple.k(), range_k) - shift*quadruple.H(quadruple.k(), range_k)) < std::numeric_limits<real_type>::epsilon()*1.e1*norm_2( quadruple.K(quadruple.k(), range_k) ) ) ;

        real_type sum_rec = 0.0 ;
        for (int i=0; i<quadruple.k(); ++i) {
          int new_rank = quadruple.rank ;
          Q( glas2::all(), range_rank ) = quadruple.Q( glas2::all(), range_rank ) ;
          auto ukh_vector = reshape( UKH(glas2::all(), i), linearization_.size_of_basis(), quadruple.U.num_rows()/linearization_.size_of_basis(), glas2::row_major() ) ;

          // Solve takes place in three steps:
          // 1) Solve the upper part of the U-L factorization
          solve_handles_.solve_upper( Q, ukh_vector, w, u_vector ) ;
          // 2) Update the Q and U matrices: orthonormalize column rank of Q. Modify U accordingly.
          u_vector(0,glas2::all()) = multiply( transpose(conj(Q)), w ) ;
          // 3) Solve for the lower part.
          solve_handles_.solve_lower( u_vector ) ;

          reshape( expand_vi.pass_ref(), linearization_.size(), linearization_.size_of_basis(), glas2::column_major() ) = 
             multiply( Q, transpose( u_vector ) ) ;
          for (int j=0; j<=quadruple.k(); ++j) {
            auto Uj = quadruple.u_vector(j) ;
            reshape( expand_vj.pass_ref(), linearization_.size(), linearization_.size_of_basis(), glas2::column_major() ) = 
              multiply( quadruple.Q( glas2::all(), range_rank ), transpose( Uj( glas2::all(), range_rank ) ) ) ;
            expand_vi -= expand_vj * quadruple.H(j,i) ;
          }
          std::cout << norm_2( expand_vi(glas2::range(0,Q.num_rows()))) << " " << norm_2(expand_vi(glas2::range_from_end(Q.num_rows(),0))) << std::endl ;
          auto term = norm_2(expand_vi.pass_ref()) / norm_2( quadruple.H(glas2::range(0,quadruple.k()+1),i)) ;
          if (options_.debug_level>1) std::cout << " recurrence " << i << " : " << term << std::endl ;
          sum_rec += term ;
        }
        std::cout << "Residual of the recurrence relation: " << sum_rec << std::endl ;

        return sum + sum_rec ;
      } // verify_recurrence()*/

      auto verify_recurrence() const {
        glas2::shared_vector< value_type > expand_vi( linearization_.size() * linearization_.size_of_basis() ) ;
        glas2::shared_vector< value_type > expand_vj( linearization_.size() * linearization_.size_of_basis() ) ;
        // Check orthogonality of the iteration vectors
        glas2::range range_rank( 0, quadruple.rank ) ;

        std::cout << "Orthogonality of Q-vectors " << norm_fro( multiply( conj(transpose(quadruple.Q( glas2::all(), range_rank ))), quadruple.Q( glas2::all(), range_rank ) ) - glas2::identity_matrix<value_type>(quadruple.rank,quadruple.rank) ) << std::endl ;

        real_type sum = 0.0 ;
        for (int i=0; i<quadruple.k()+1; ++i) {
          auto Ui = quadruple.u_vector(i) ;
          reshape( expand_vi.pass_ref(), linearization_.size(), linearization_.size_of_basis(), glas2::column_major() ) = 
            multiply( quadruple.Q( glas2::all(), range_rank ), transpose( Ui( glas2::all(), range_rank ) ) ) ;
          real_type local_sum = std::abs( norm_2( expand_vi )-1.0 ) ;
          for (int j=0; j<i; ++j) {
              auto Uj = quadruple.u_vector(j) ;
              reshape( expand_vj.pass_ref(), linearization_.size(), linearization_.size_of_basis(), glas2::column_major() ) = 
                multiply( quadruple.Q( glas2::all(), range_rank ), transpose( Uj( glas2::all(), range_rank ) ) ) ;
              local_sum += std::abs( inner_prod( conj(expand_vj), expand_vi ) ) ;
          }
          if (options::value_of<options::debug_level>(options_)>1) std::cout << " orthogonality " << i << " : " << local_sum << std::endl ;
          sum += local_sum ;
        }
        std::cout << "Orthogonality of " << quadruple.k()+1 << " vectors: " << sum << std::endl ;

        if (options::value_of<options::debug_level>(options_)>1) std::cout << "Verify recurrence: " << std::endl ;
        auto Q = quadruple.Q( glas2::all(), glas2::range( 0, quadruple.rank ) ); 

        auto multiply_handle = linearization_.template multiply_handle<value_type>() ;

        // UKH holds U * (K-\sigma H) with sigma the last shift
        // Continuation combination is I, so K - sigma * H = H * (Poles-sigma*I) + I
        glas2::matrix<value_type> UH( quadruple.U.num_rows(), quadruple.k() ) ; 
        glas2::matrix<value_type> UK( quadruple.U.num_rows(), quadruple.k() ) ; 
        glas2::range range_k( 0, quadruple.k() ) ;
        glas2::range range_k1( 0, quadruple.k()+1 ) ;
        UH = multiply( quadruple.U( glas2::all(), range_k1 ), quadruple.H(range_k1, range_k) ) ;
        UK = multiply( quadruple.U( glas2::all(), range_k1 ), quadruple.K(range_k1, range_k) ) ;

        glas2::matrix<value_type> AUH( linearization_.size_of_basis(), Q.num_columns() ) ; 
        glas2::matrix<value_type> BUK( AUH.num_rows(), AUH.num_columns() ) ;
        glas2::vector<value_type> w_A( Q.num_rows() ) ;
        glas2::vector<value_type> w_B( Q.num_rows() ) ;

        real_type sum_rec = 0.0 ;
        real_type sum_rec_norm = 0.0 ;
        for (int i=0; i<quadruple.k(); ++i) {
          auto u_h_i = reshape( UH(glas2::all(),i), linearization_.size_of_basis(), UH.num_rows()/linearization_.size_of_basis(), glas2::row_major() )(glas2::all(), range_rank ) ;
          auto u_k_i = reshape( UK(glas2::all(),i), linearization_.size_of_basis(), UK.num_rows()/linearization_.size_of_basis(), glas2::row_major() )(glas2::all(), range_rank ) ;
          multiply_handle.multiply_A( Q, u_h_i, w_A, AUH ) ;
          multiply_handle.multiply_B( Q, u_k_i, w_B, BUK ) ;

          auto term = std::sqrt( norm_2_squared(w_A-w_B) + norm_fro_squared(AUH-BUK) ) ;
          auto term_norm = std::sqrt( norm_2_squared(w_A) + norm_2_squared(w_B) + norm_fro_squared(AUH) + norm_fro_squared(BUK) ) ;
          //if (options_.debug_level>1) std::cout << " recurrence " << i << " : " << term << std::endl ;
          if (options::value_of<options::debug_level>(options_)>1) std::cout << " recurrence " << i << " : " << norm_2(w_A-w_B) << " + " << norm_fro(AUH-BUK) << " / " << term_norm << std::endl ;
          sum_rec += term ;
          sum_rec_norm += term_norm ;
        }
        std::cout << "Residual of the recurrence relation: " << sum_rec / sum_rec_norm << std::endl ;

        return std::tuple( sum, sum_rec ) ;
      } // verify_recurrence()

    public:
      real_type norm_h() const { return std::sqrt(norm_h_) ; }

    public:
      Quadruple& quadruple ;

      decltype(auto) solve_handles() const { return solve_handles_ ; }

    private:
      real_type                                          norm_h_ ;
      Linearization&                                     linearization_ ;
      ShiftGenerator&                                    shift_generator_ ;
      sequence_of_solve_handle<ShiftValueType, Linearization> solve_handles_ ;
      Options const&                                     options_ ;
  } ;


  template <typename Tag, typename ShiftValueType, typename Linearization, typename ShiftGenerator, typename Quadruple, typename Options>
  decltype(auto) make_cork_process( Linearization& linearization, ShiftGenerator& shift_generator, Quadruple& quadruple, Options const& options ) {
    return cork_process<Tag,ShiftValueType,Linearization,ShiftGenerator,Quadruple,Options>( linearization, shift_generator, quadruple, options ) ;
  }


} } // namespace CORK::krylov


#endif
