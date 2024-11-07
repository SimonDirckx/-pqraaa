//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2015.
//
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_krylov_state_vector_initial_vector_hpp
#define cork_krylov_state_vector_initial_vector_hpp

#include <cork/krylov/toar_triple.hpp>
#include <glas2/matrix.hpp>
#include <glas2/vector.hpp>

namespace CORK { namespace krylov {

  template <typename RHS>
  class state_vector_initial_vector {
    public:
      inline state_vector_initial_vector( RHS const& rhs )
      : rhs_( rhs )
      , reduced_rhs_( 2, 2 )
      {
        static_assert( glas2::is< glas2::DenseVector, RHS >::value, "state_vector_initial_vector: RHS must be glas2::DenseVector" ) ;
      }

    public:
      template <typename Triple, typename SolveHandle, typename Backend>
      typename std::enable_if< std::is_same< typename Triple::value_type, typename SolveHandle::value_type >::value
                             >::type apply( Triple& triple, SolveHandle const& solve_handle, Backend const& backend ) const
      {
        fill(reduced_rhs_, 0.) ;

        triple.Q( glas2::all(), 0 ) = rhs_ ;
        solve_handle.linear_solver().solve( triple.Q( glas2::all(), 0 ) ) ;
        reduced_rhs_(0,0) = norm_2( backend, triple.Q( glas2::all(), 0 ) ) ;
        triple.Q( glas2::all(), 0 ) /= reduced_rhs_(0,0) ;
        glas2::shared_vector< typename SolveHandle::value_type > one( 1 ) ; fill(one,1.0) ;
        glas2::shared_matrix< typename SolveHandle::value_type > z0( solve_handle.basis_handle().size(), 1 ) ;

        // Evaluate functions using the MN solve
        z0(0,0) = 1.0 ;
        auto z0r = z0(glas2::range_from_end(1,0), glas2::all()) ;
        solve_handle.basis_handle().lower_solve_right_hand_side(one,z0r) ;
        solve_handle.basis_handle().solve(z0r) ;
        z0r *= -1. ;

        auto U = triple.u_vector(0) ;
        fill( U, 0.0 ) ;
        U(glas2::all(), 0) = z0(glas2::all(),0) ;
        reduced_rhs_(0,0) *= norm_2(U(glas2::all(),0)) ;

#ifndef NDEBUG
        bool ok =
#endif
        triple.initial_vector( 1 ) ;
        assert( ok ) ; // If not OK, there is a problem with Basis4CORK
        rank_ = triple.rank ;
      } // apply()


      // Only supported for Rational Krylov
      template <typename Quadruple, typename SolveHandle, typename Backend>
      typename std::enable_if< std::is_convertible< std::complex<typename Quadruple::value_type>, typename SolveHandle::value_type >::value
                               && std::is_same< typename RHS::value_type, typename Quadruple::value_type >::value
                             >::type apply( Quadruple& quadruple, SolveHandle const& solve_handle, Backend const& backend ) const
      {
        typedef std::complex<typename Quadruple::value_type> value_type ;
        fill(reduced_rhs_, 0.) ;

        glas2::vector< value_type > rhs_c( rhs_.size() ) ;
        rhs_c = rhs_ ;
        solve_handle.linear_solver().solve( rhs_c ) ;

        // Compute the first two columns of Q
        glas2::vector< value_type > first_vector(2) ;
        reduced_rhs_(0,0) = norm_2(real(rhs_c)) ;
        first_vector(0) = 1.0 ;
        first_vector(1) = 0.0 ;
        quadruple.Q( glas2::all(), 0 ) = real(rhs_c) / reduced_rhs_(0,0) ;
        int rank = 1 ;

        auto norm_i = norm_2(imag(rhs_c))/reduced_rhs_(0,0) ;
        if (norm_i>0.0) {
          quadruple.Q( glas2::all(), 1 ) = imag(rhs_c)/reduced_rhs_(0,0) ;
          auto g = inner_prod( quadruple.Q(glas2::all(),0), quadruple.Q(glas2::all(),1) ) ;
          quadruple.Q( glas2::all(), 1 ) -= quadruple.Q(glas2::all(),0) * g ;
          first_vector(0) += value_type(0.,g) ;
          g = inner_prod( quadruple.Q(glas2::all(),0), quadruple.Q(glas2::all(),1) ) ;
          quadruple.Q( glas2::all(), 1 ) -= quadruple.Q(glas2::all(),0) * g ;
          first_vector(0) += value_type(0.,g) ;
          auto norm_ii = norm_2( quadruple.Q( glas2::all(), 1 ) ) ;
          first_vector(1) = value_type(0.0, norm_ii) ;
          if (norm_ii>0.0) {
            quadruple.Q( glas2::all(), 1 ) /= norm_ii ;
            rank = 2 ;
          }
        }

        //first_vector /= first_vector(0).real() ;

        // Compute the first two columns of U
        glas2::shared_vector< typename SolveHandle::value_type > one( rank ) ; fill(one,-1.0) ;
        glas2::shared_matrix< typename SolveHandle::value_type > z0( solve_handle.basis_handle().size(), rank ) ;

        // Evaluate functions using the MN solve
        z0(0, glas2::all() ) = first_vector(glas2::range(0,rank)) ;
        auto z0r = z0(glas2::range_from_end(1,0), glas2::all()) ;
        solve_handle.basis_handle().lower_solve_right_hand_side(-z0(0,glas2::all()),z0r) ;
        solve_handle.basis_handle().solve(z0r) ;

        // Split U in real and imaginary parts.
        // Real part
        auto U0 = quadruple.u_vector(0) ;
        fill( U0, 0.0 ) ;
        U0(glas2::all(),  glas2::range(0,rank)) = real( z0(glas2::all(), glas2::range(0,rank)) ) ;
        auto norm_u0 = norm_fro( U0(glas2::all(),  glas2::range(0,rank)) ) ;
        reduced_rhs_(0,0) *= norm_u0 ;
        z0 /= norm_u0 ;

#ifndef NDEBUG
        bool ok =
#endif
        quadruple.initial_vector( rank ) ;
        assert( ok ) ; // If not OK, there is a problem with Basis4CORK

        // Imaginary part
        if (norm_2( imag( z0(glas2::all(),0) ) )!=0.0) {
          auto U1 = quadruple.u_vector(1) ;
          fill( U1, 0.0 ) ;
          U1(glas2::all(), glas2::range(0,rank)) = imag( z0(glas2::all(),glas2::range(0,rank)) ) ;
          static_cast<typename Quadruple::triple_type&>(quadruple).next_step( rank ) ;
          quadruple.K(0,0) = solve_handle.shift().imag() ;
          quadruple.K(1,0) = 0.0 ;
          quadruple.K(glas2::range(0,2),0) += solve_handle.shift().real() * quadruple.H(glas2::range(0,2),0) ;
          fill( quadruple.K(glas2::range_from_end(2,0),0), 0.0 ) ;

          reduced_rhs_(glas2::all(),1) = reduced_rhs_(0,0)*quadruple.H(glas2::range(0,2),0) ;
        }
        rank_ = quadruple.rank ;
      } // apply()

      auto const& reduced_rhs_real() const { return reduced_rhs_ ; }

      template <typename ValueType>
      typename std::enable_if< std::is_arithmetic<ValueType>::value, glas2::shared_vector<ValueType> >::type reduced_rhs() const {
        glas2::shared_vector<ValueType> reduced_rhs_complex( rank_ ) ;
        reduced_rhs_complex( glas2::range(0,rank_) ) = reduced_rhs_( glas2::range(0,rank_), 0 ) ;
        assert( norm_2( reduced_rhs_( glas2::range(0,rank_), 1) )==0.0 ) ;
        return reduced_rhs_complex ;
      }

      template <typename ValueType>
      typename std::enable_if< std::is_arithmetic<typename ValueType::value_type>::value, glas2::shared_vector<ValueType> >::type reduced_rhs() const {
        glas2::shared_vector<ValueType> reduced_rhs_complex( rank_ ) ;
        reduced_rhs_complex = reduced_rhs_( glas2::range(0,rank_), 0 ) + ValueType(0.,1.) * reduced_rhs_( glas2::range(0,rank_), 1 ) ;
        return reduced_rhs_complex ;
      }

    private:
      RHS const&                                        rhs_ ; 
      mutable glas2::matrix< typename RHS::value_type > reduced_rhs_ ;
      mutable int                                       rank_ ;
  } ;
   
} } // namespace CORK::krylov


#endif
