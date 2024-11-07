//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2015.
//
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_krylov_structured_initial_vector_hpp
#define cork_krylov_structured_initial_vector_hpp

#include <cork/krylov/toar_triple.hpp>
#include <glas2/matrix.hpp>
#include <glas2/vector.hpp>

namespace CORK { namespace krylov {

  template <typename X>
  class structured_initial_vector {
    public:
      inline structured_initial_vector( X const& x )
      : x_( x )
      , norm_x_( norm_2(x) )
      {
        static_assert( glas2::is< glas2::DenseVector, X >::value, "structured_initial_vector: X must be glas2::DenseVector" ) ;
      }

    public:
      template <typename Triple, typename SolveHandle, typename Backend>
      decltype(std::abs(typename std::decay<X>::type::value_type())) apply( Triple& triple, SolveHandle const& solve_handle, Backend const& backend ) const
      {
        static_assert( std::is_convertible< typename SolveHandle::value_type, typename Triple::value_type >::value, "CORK::krylov::structured_initial_vector: value types of Triple and SolveHandle are not compatible" ) ;
        triple.Q( glas2::all(), 0 ) = x_ / norm_x_ ;
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
        auto norm_u = norm_2( U(glas2::all(), 0) ) ;

#ifndef NDEBUG
        bool ok =
#endif
        triple.initial_vector( 1 ) ;
        assert( ok ) ; // If not OK, there is a problem with Basis4CORK
        return norm_x_ * norm_u ;
      } // apply()

    private:
      X const&                                     x_ ; 
      decltype(std::abs(typename std::decay<X>::type::value_type())) norm_x_ ;
  } ;
   
} } // namespace CORK::krylov


#endif
