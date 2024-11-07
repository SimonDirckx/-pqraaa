//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2015.
//
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_krylov_user_defined_initial_vector_hpp
#define cork_krylov_user_defined_initial_vector_hpp

#include <cork/vector.hpp>
#include <cork/matrix.hpp>
#include <cork/exception/initial_vector.hpp>
#include <cork/krylov/toar_triple.hpp>
#include <glas2/matrix.hpp>
#include <glas2/vector.hpp>
#include <functional>

namespace CORK { namespace krylov {

  template <typename ValueType>
  class user_defined_initial_vector {
    public:
      typedef ValueType                                         value_type ;
      typedef std::function<void(CORK::matrix<value_type>)>     function_type ;

    public:
      inline user_defined_initial_vector( bool structured, value_type const& lambda, function_type const& init_function )
      : structured_( structured )
      , lambda_( lambda )
      , init_function_( init_function )
      {}

    public:
      // For solvers with the same real or complex value type as the triple.
      template <typename Triple, typename SolveHandle, typename Backend>
      typename std::enable_if< !std::is_same<std::complex<value_type>,typename SolveHandle::value_type>::value>::type
      apply( Triple& triple, SolveHandle const& solve_handle, Backend const& backend ) const {
        auto initial_vec =  triple.Q(glas2::all(), glas2::range(0,triple.degree()));
        init_function_(initial_vec) ;
        auto U = triple.u_vector(0) ;

        int rank ;
        if (structured_) {
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
          rank = 1 ;
        } else {
          rank = initial_vec.num_columns() ;
          U(glas2::all(), glas2::range(0, rank) ) = glas2::identity_matrix<typename Triple::value_type>(U.num_rows(), rank);
        }
        glas2::fill( U(glas2::all(), glas2::range_from_end(rank,0)), 0.0 ) ;

#ifndef NDEBUG
        bool ok =
#endif
        triple.initial_vector( rank, false ) ;
        if (triple.rank <= 0) { throw CORK::exception::initial_vector() ; }
        assert( ok ) ; // If not OK, there is a problem with Basis4CORK
      }

    public:
      // For complex valued solver with real value type for the triple.
      template <typename Triple, typename SolveHandle, typename Backend>
      typename std::enable_if< std::is_same<std::complex<value_type>,typename SolveHandle::value_type>::value>::type
      apply( Triple& triple, SolveHandle const& solve_handle, Backend const& backend ) const {
        auto initial_vec =  triple.Q(glas2::all(), glas2::range(0,triple.degree()));
        init_function_(initial_vec) ;
        auto U = triple.u_vector(0) ;

        int rank ;
        if (structured_) {
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
          U(glas2::all(), 1) = glas2::real( z0(glas2::all(),0) ) ;
          if (norm_2(glas2::imag( z0(glas2::all(),0) ))!=0.0)
            std::cout << "CORK: WARNING: only the real part of the structured vector is taken as initial vector" << std::endl ;
          rank = 1 ;
        } else {
          rank = triple.degree() ;
          U(glas2::all(), glas2::range(0, triple.degree()) ) = glas2::identity_matrix<typename Triple::value_type>(U.num_rows(),  triple.degree());
        }

#ifndef NDEBUG
        bool ok =
#endif
        triple.initial_vector( rank, false ) ;
        if (triple.rank <= 0) { throw CORK::exception::initial_vector() ; }
        assert( ok ) ; // If not OK, there is a problem with Basis4CORK
      }

    private:
      bool                 structured_ ;
      value_type           lambda_ ;
      function_type        init_function_;
 //     function_type const& init_function_;
  } ;
   
} } // namespace CORK::krylov


#endif
