//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2015.
//
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_krylov_cork_process_real_hpp
#define cork_krylov_cork_process_real_hpp

#include <cork/krylov/cork_process.hpp>
#include <cork/options/value_of.hpp>
#include <cork/options/debug_level.hpp>
#include <cork/options/backend.hpp>
#include <cork/utility/value_type_for.hpp>
#include <cmath>

namespace CORK { namespace krylov {

  struct real_matrix_complex_shift_tag {} ;

  template <typename ShiftValueType, typename Linearization, typename ShiftGenerator, typename Quadruple, typename Options>
  class cork_process< real_matrix_complex_shift_tag, ShiftValueType, Linearization, ShiftGenerator, Quadruple, Options>
  {
    public:
      typedef typename Quadruple::value_type                                     q_value_type ;
      typedef typename Linearization::template solve_handle_type<ShiftValueType> solve_handle_type ;
      typedef ShiftValueType shift_value_type ;

      typedef typename CORK::value_type_for< ShiftValueType, Linearization>::type  value_type ;
      typedef decltype( std::abs(value_type()) )                                   real_type ;

      static_assert( std::is_convertible< ShiftValueType, std::complex<q_value_type> >::value
                   && !std::is_convertible< ShiftValueType, q_value_type >::value
                   , "CORK::krylov::cork_process: ShiftValueType is not convertible to the complex value_type of Quadruple"
                   ) ;

    public:
      cork_process( Linearization& linearization, ShiftGenerator& shift_generator, Quadruple& quadruple, Options const& options )
      : quadruple( quadruple )
      , linearization_( linearization )
      , shift_generator_( shift_generator )
      , solve_handles_( linearization_ )
      , options_( options )
      {
        static_cast<info&>(linearization_.information()).krylov_process = "CORK for real pencils" ;
      }

    public:
      void initialize( shift_value_type shift ) {
        // WHY DO WE DO THIS?
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
        initial.apply( quadruple, solve_handles_, options::value_of< options::backend_key >(options_) ) ;
      } // initialize()

    public:
      // Add two vectors: Re(w) and Im(w) where w are the result of a Krylov step with the linearization
      void expand( int max_k ) {
        norm_h_ = norm_fro_squared( quadruple.H( glas2::range(0,quadruple.k()+1), glas2::range(0,quadruple.k()) ) ) ;

        // Make sure we do at least one iteration, otherwise there is not advance.
        assert( quadruple.k()<max_k ) ;

        glas2::vector< value_type > w_complex( quadruple.Q.num_rows() ) ;
        glas2::matrix< value_type > U_complex( linearization_.size_of_basis(), max_k+linearization_.size_of_basis() ) ;

        for (; quadruple.k()<max_k;) {
          auto shift = shift_generator_.next_shift() ;//quadruple.k()) ;
          if (shift.imag()!=0 && quadruple.k()>=quadruple.k_max()-1) break ;

          solve_handles_.select( shift_generator_.solve_handle() ) ;
          solve_handles_.shift( shift ) ;
          if (options::value_of<options::debug_level>(options_)>1) std::cout << "Iteration " << quadruple.k() << "  shift = " << shift << std::endl ;

          //shift_generator.continuation_combination( quadruple.k(), quadruple.continuation_combination( glas2::range(0,quadruple.k()+1), quadruple.k() ) ) ;

          // k() is the iteration count starting from 0
          // Solve takes place in three steps:
          // 1) Solve the upper part of the U-L factorization
          solve_handles_.solve_upper( quadruple.Q_k(), quadruple.u_vector_k( quadruple.k() ), w_complex, U_complex(glas2::all(),glas2::range(0,quadruple.rank)), options::value_of< options::backend_key >(options_) ) ;
          // 2) Update the Q and U matrices: orthonormalize column rank of Q. Modify U accordingly.
          // Add column Q(:,rank) and update U(:,k+1) ;

          quadruple.Q(glas2::all(), quadruple.rank) = glas2::real( w_complex ) ;
          quadruple.u_vector_k( quadruple.k()+1 ) = glas2::real(U_complex(glas2::all(),glas2::range(0,quadruple.rank))) ;
          quadruple.u_vector( quadruple.k()+1 )(0,quadruple.rank) = 1.0 ;
          quadruple.add_column_of_Q( quadruple.k()+1 ) ;

          if (shift.imag()!=0) {
            quadruple.Q(glas2::all(), quadruple.rank) = glas2::imag( w_complex ) ;
            quadruple.u_vector_k( quadruple.k()+2 ) = glas2::imag(U_complex(glas2::all(),glas2::range(0,quadruple.rank))) ;
            quadruple.u_vector( quadruple.k()+2 )(0,quadruple.rank) = 1.0 ;
            quadruple.add_column_of_Q( quadruple.k()+2 ) ;

            U_complex(glas2::all(),glas2::range(0,quadruple.rank)) = quadruple.u_vector_k( quadruple.k()+1 ) + value_type(0.,1.)*quadruple.u_vector_k( quadruple.k()+2 ) ;
          } else {
            U_complex(glas2::all(),glas2::range(0,quadruple.rank)) = quadruple.u_vector_k( quadruple.k()+1 ) ;
          }
          // 3) Solve for the lower part.
          solve_handles_.solve_lower( U_complex(glas2::all(),glas2::range(0,quadruple.rank)) ) ;
          quadruple.u_vector_k( quadruple.k()+1 ) = glas2::real(U_complex(glas2::all(),glas2::range(0,quadruple.rank))) ;
          if (shift.imag()!=0)
            quadruple.u_vector_k( quadruple.k()+2 ) = glas2::imag(U_complex(glas2::all(),glas2::range(0,quadruple.rank))) ;

          // Orthonormalize vector in Krylov space
          int new_rank = quadruple.rank ;
          if (shift.imag()!=0)
            quadruple.next_step_double( new_rank, shift ) ;
          else
            quadruple.next_step( new_rank, shift.real() ) ;
          // k() has been increased by one: is now the iteration count starting from 1.

          // Check for lucky breakdown
          norm_h_ += norm_2_squared( quadruple.H( glas2::range(0,quadruple.k()), quadruple.k()-1 ) ) ;

          // k() has been increased by one: is now the iteration count starting from 1.

          // Check for lucky breakdown
          norm_h_ += norm_2_squared( quadruple.H( glas2::range(0,quadruple.k()), quadruple.k()-1 ) ) ;
          if (glas2::abs_squared(quadruple.H(quadruple.k(),quadruple.k()-1))<std::pow(options::value_of<options::krylov_breakdown_tolerance<real_type>>(options_),2)*norm_h_) {
            std::cout << "Breakdown bug in cork_expand() at iteration " << quadruple.k() << std::endl ;
            break;
          }
          if (options::value_of<options::debug_level>(options_)>2) std::cout << "Iteration " << quadruple.k() << "  rank = " << new_rank << " " << quadruple.rank << std::endl ;
        } // end for
      } // expand()


      auto verify_recurrence() const {
        glas2::shared_vector< value_type > expand_vi( linearization_.size() * linearization_.size_of_basis() ) ;
        glas2::shared_vector< value_type > expand_vj( linearization_.size() * linearization_.size_of_basis() ) ;
        // Check orthogonality of the iteration vectors
        glas2::range range_rank( 0, quadruple.rank ) ;

        real_type sum=0.0;
        for (int i=0; i<quadruple.k()+1; ++i) {
          auto Ui = quadruple.u_vector(i) ;
          reshape( expand_vi.pass_ref(), linearization_.size(), linearization_.size_of_basis(), glas2::column_major() ) = 
            multiply( quadruple.Q( glas2::all(), range_rank ), transpose( Ui( glas2::all(), range_rank ) ) ) ;
          real_type local_sum = std::abs( norm_2( expand_vi )-1.0 ) ;
          for (int j=0; j<quadruple.k()+1; ++j) {
            if (i!=j) {
              auto Uj = quadruple.u_vector(j) ;
              reshape( expand_vj.pass_ref(), linearization_.size(), linearization_.size_of_basis(), glas2::column_major() ) = 
                multiply( quadruple.Q( glas2::all(), range_rank ), transpose( Uj( glas2::all(), range_rank ) ) ) ;
              local_sum += std::abs( inner_prod( conj(expand_vj), expand_vi ) ) ;
            }
          }
          if (options::value_of<options::debug_level>(options_)>1) std::cout << " orthogonality " << i << " : " << local_sum << std::endl ;
          sum+= local_sum ;
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
        for (int i=0; i<quadruple.k(); ++i) {
          auto u_h_i = reshape( UH(glas2::all(),i), linearization_.size_of_basis(), UH.num_rows()/linearization_.size_of_basis(), glas2::row_major() )(glas2::all(), range_rank ) ;
          auto u_k_i = reshape( UK(glas2::all(),i), linearization_.size_of_basis(), UK.num_rows()/linearization_.size_of_basis(), glas2::row_major() )(glas2::all(), range_rank ) ;
          multiply_handle.multiply_A( Q, u_h_i, w_A, AUH ) ;
          multiply_handle.multiply_B( Q, u_k_i, w_B, BUK ) ;

          auto term_norm = std::sqrt( norm_2_squared(w_A) + norm_fro_squared(AUH) ) ;
          auto term = std::sqrt( norm_2_squared(w_A-w_B) + norm_fro_squared(AUH-BUK) ) / std::sqrt( norm_2_squared(w_A) + norm_fro_squared(AUH) ) ;
          if (options::value_of<options::debug_level>(options_)>1) std::cout << " recurrence " << i << " : " << norm_2(w_A-w_B) << " + " << norm_fro(AUH-BUK) << " / " << term_norm << std::endl ;
          sum_rec += term ;
        }
        std::cout << "Residual of the recurrence relation: " << sum_rec << std::endl ;

        //return sum + sum_rec ;
        return std::tuple( sum, sum_rec ) ;
      } // verify_recurrence()

    public:
      real_type norm_h() const { return std::sqrt(norm_h_) ; }

    public:
      Quadruple& quadruple ;

      decltype(auto) solve_handles() const { return (solve_handles_) ; }

    private:
      real_type                                          norm_h_ ;
      Linearization&                                     linearization_ ;
      ShiftGenerator&                                          shift_generator_ ;
      sequence_of_solve_handle<ShiftValueType, Linearization> solve_handles_ ;
      Options const&                                     options_ ;
  } ;

} } // namespace CORK::krylov


#endif
