//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2015.
//
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_krylov_toar_process_hpp
#define cork_krylov_toar_process_hpp

#include <cork/linearization/cork_linearization_solve_handle.hpp>
#include <cork/krylov/info.hpp>
#include <cork/krylov/toar_triple.hpp>
#include <cork/options/backend.hpp>
#include <cork/options/krylov_breakdown_tolerance.hpp>
#include <cork/options/debug_level.hpp>
#include <cork/options/value_of.hpp>
#include <cmath>

namespace CORK { namespace krylov {

  template <typename ShiftValueType, typename Linearization, typename Triple, typename Options>
  class toar_process {
    public:
      typedef typename Triple::value_type                                    value_type ;
      typedef decltype( std::abs(value_type()) )                             real_type ;
      typedef typename Linearization::template solve_handle_type<ShiftValueType> solve_handle_type ;

      static_assert( std::is_convertible< ShiftValueType, value_type >::value, "CORK::krylov::toar_process: ShiftValueType is not convertible to the value_type of Triple" ) ;

    public:
      toar_process( Linearization& linearization, Triple& triple, Options const& options )
      : triple( triple )
      , linearization_( linearization )
      , solve_handle_( linearization_.template solve_handle<ShiftValueType>() )
      , options_( options )
      {
        auto& infor = static_cast<krylov::info&>(linearization_.information()) ;
        infor.krylov_process = "TOAR" ;
      }

    public:
      void initialize( ShiftValueType const& shift ) {
        solve_handle_.shift( shift ) ;
      } // initialize()

      template <typename InitialVector>
      void initial_vector( InitialVector const& initial ) {
        initial.apply( triple, solve_handle_, options::value_of<options::backend_key>(options_) ) ;
      } // initialize()

    public:
      void expand( int max_k ) {
        norm_h_ = norm_fro_squared( triple.H( glas2::range(0,triple.k()+1), glas2::range(0,triple.k()) ) ) ;

        for (; triple.k()<max_k;) {
          // k() is the iteration count starting from 0
          // 1) Solve the upper part of the U-L factorization
          solve_handle_.solve_upper( triple.Q_k(), triple.u_vector_k( triple.k() ), triple.Q(glas2::all(), triple.rank), triple.u_vector_k( triple.k()+1 ), options::value_of<options::backend_key>(options_) ) ;
          // 2) Update the Q and U matrices: orthonormalize column rank of Q. Modify U accordingly.
          triple.u_vector( triple.k()+1 )(0,triple.rank) = 1.0 ;
          // Add column Q(:,rank) and update U(:,k+1) ;
          triple.add_column_of_Q( triple.k()+1 ) ;
          // 3) Solve for the lower part.
          solve_handle_.solve_lower( triple.u_vector_k( triple.k()+1 ) ) ;

          int new_rank = triple.rank ;
          triple.next_step( new_rank ) ;
          // k() has been increased by one: is now the iteration count starting from 1.

          // Check for lucky breakdown
          norm_h_ += norm_2_squared( triple.H( glas2::range(0,triple.k()), triple.k()-1 ) ) ;
          if (glas2::abs_squared(triple.H(triple.k(),triple.k()-1))<std::pow(options::value_of<options::krylov_breakdown_tolerance<real_type>>(options_),2)*norm_h_) {
            std::cout << "Breakdown bug in toar_expand() at iteration " << triple.k() << std::endl ;
            break ;
          }
          if (options::value_of<options::debug_level>(options_)>2) std::cout << "Iteration " << triple.k() << "  rank = " << new_rank << " " << triple.rank << std::endl ;
        } // end for
      } // expand()

      real_type verify_recurrence() {
        glas2::vector< value_type > expand_vi( linearization_.size() * linearization_.size_of_basis() ) ;
        glas2::vector< value_type > expand_vj( linearization_.size() * linearization_.size_of_basis() ) ;
        // Check orthogonality of the iteration vectors
        glas2::range range_rank( 0, triple.rank ) ;

        real_type sum;
        for (int i=0; i<triple.k()+1; ++i) {
          auto Ui = triple.u_vector(i) ;
          reshape( expand_vi, linearization_.size(), linearization_.size_of_basis(), glas2::column_major() ) = 
            multiply( triple.Q( glas2::all(), range_rank ), transpose( Ui( glas2::all(), range_rank ) ) ) ;
          sum = std::abs( norm_2( expand_vi )-1.0 ) ;
          for (int j=0; j<triple.k()+1; ++j) {
            if (i!=j) {
              auto Uj = triple.u_vector(j) ;
              reshape( expand_vj, linearization_.size(), linearization_.size_of_basis(), glas2::column_major() ) = 
              multiply( triple.Q( glas2::all(), range_rank ), transpose( Uj( glas2::all(), range_rank ) ) ) ;
              sum += std::abs( inner_prod( conj(expand_vj), expand_vi ) ) ;
            }
          }
        }
        std::cout << "Orthogonality of " << triple.k()+1 << " vectors: " << sum << std::endl ;

        if (options::value_of<options::debug_level>(options_)>2) std::cout << "Verify recurrence: " << std::endl ;
        /*real_type sum_rec = 0.0 ;
        glas2::shared_matrix<value_type> u_vector( linearization_.size_of_basis(), triple.rank+1 ) ; 
        glas2::shared_matrix<value_type> Q( triple.Q.num_rows(), triple.rank+1 ) ; 
        for (int i=0; i<triple.k(); ++i) {
          int new_rank = triple.rank ;
          Q( glas2::all(), range_rank ) = triple.Q( glas2::all(), range_rank ) ;
          //linearization_.krylov_step( Q, new_rank, triple.u_vector( i ), u_vector, options_.Q_drop_tol ) ;
          solve_handle_.solve_upper( Q( glas2::all(), range_rank ), triple.u_vector_k( i ), Q(glas2::all(), triple.rank), u_vector(glas2::all(),range_rank) ) ;
          // 2) Update the Q and U matrices: orthonormalize column rank of Q. Modify U accordingly.
          u_vector(0,triple.rank) = 1.0 ;
          // Add column Q(:,rank) and update U(:,k+1) ;
          //   triple.add_column_of_Q( triple.k()+1 ) ; NOT NEEDED IN PRINCIPLE
          // 3) Solve for the lower part.
          solve_handle_.solve_lower( u_vector ) ;
          //std::cout << u_vector( glas2::all(), glas2::range(0,new_rank)) << std::endl ;

          reshape( expand_vi, linearization_.size(), linearization_.size_of_basis(), glas2::column_major() ) = 
             multiply( Q( glas2::all(), glas2::range(0,new_rank)), transpose( u_vector( glas2::all(), glas2::range(0,new_rank)) ) ) ;
          //std::cout << norm_fro( Q( glas2::all(), glas2::range(0,new_rank)) ) << " " << norm_fro( u_vector( glas2::all(), glas2::range(0,new_rank)) ) << " " << norm_2( expand_vi ) << std::endl ;
          for (int j=0; j<=triple.k(); ++j) {
            auto Uj = triple.u_vector(j) ;
            reshape( expand_vj, linearization_.size(), linearization_.size_of_basis(), glas2::column_major() ) = 
              multiply( triple.Q( glas2::all(), range_rank ), transpose( Uj( glas2::all(), range_rank ) ) ) ;
            expand_vi -= expand_vj * triple.H(j,i) ;
          }
          auto term = norm_2(expand_vi) / norm_2( triple.H(glas2::range(0,triple.k()+1),i)) ;
          if (options_.debug_level>2) std::cout << " recurrence " << i << ": " << term << std::endl ;
          sum_rec += term ;
        }*/

        auto Q = triple.Q( glas2::all(), glas2::range( 0, triple.rank ) ); 

        auto multiply_handle = linearization_.template multiply_handle<value_type>() ;

        // UKH holds U * (K-\sigma H)
        // Continuation combination is I, so K - sigma * H = H * (Poles-sigma*I) + I
        glas2::matrix<value_type> UH( triple.U.num_rows(), triple.k() ) ; 
        glas2::matrix<value_type> UK( triple.U.num_rows(), triple.k() ) ; 
        glas2::range range_k( 0, triple.k() ) ;
        glas2::range range_k1( 0, triple.k()+1 ) ;
        UH = multiply( triple.U( glas2::all(), range_k1 ), triple.H(range_k1, range_k) ) ;
        UK = multiply( triple.U( glas2::all(), range_k1 ), solve_handle_.shift() * triple.H(range_k1, range_k) ) ;
        UK += triple.U( glas2::all(), range_k ) ;

        glas2::matrix<value_type> AUH( linearization_.size_of_basis(), Q.num_columns() ) ; 
        glas2::matrix<value_type> BUK( AUH.num_rows(), AUH.num_columns() ) ;
        glas2::vector<value_type> w_A( Q.num_rows() ) ;
        glas2::vector<value_type> w_B( Q.num_rows() ) ;

        real_type sum_rec = 0.0 ;
        for (int i=0; i<triple.k(); ++i) {
          auto u_h_i = reshape( UH(glas2::all(),i), linearization_.size_of_basis(), UH.num_rows()/linearization_.size_of_basis(), glas2::row_major() )(glas2::all(), range_rank ) ;
          auto u_k_i = reshape( UK(glas2::all(),i), linearization_.size_of_basis(), UK.num_rows()/linearization_.size_of_basis(), glas2::row_major() )(glas2::all(), range_rank ) ;
          multiply_handle.multiply_A( Q, u_h_i, w_A, AUH ) ;
          multiply_handle.multiply_B( Q, u_k_i, w_B, BUK ) ;

          auto term = std::sqrt( norm_2_squared(w_A-w_B) + norm_fro_squared(AUH-BUK) ) / std::sqrt( norm_2_squared(w_A) + norm_fro_squared(AUH) ) ;
          //if (options_.debug_level>1) std::cout << " recurrence " << i << " : " << term << std::endl ;
          if (options::value_of<options::debug_level>(options_)>1) std::cout << " recurrence " << i << " : " << norm_2(w_A-w_B) << " + " << norm_fro(AUH-BUK) << std::endl ;
          sum_rec += term ;
        }
        std::cout << "Residual of the recurrence relation: " << sum_rec << std::endl ;

        return sum + sum_rec ;
      } // verify_recurrence()

    public:
      real_type norm_h() const { return std::sqrt(norm_h_) ; }

    public:
      Triple& triple ;

    private:
      real_type         norm_h_ ;
      Linearization&    linearization_ ;
      solve_handle_type solve_handle_ ;
      Options const&    options_ ;
  } ;

   
} } // namespace CORK::krylov


#endif
