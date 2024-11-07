//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2015.
//
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_krylov_random_initial_vector_hpp
#define cork_krylov_random_initial_vector_hpp

#include <cork/krylov/toar_triple.hpp>
#include <glas2/matrix.hpp>
#include <glas2/vector.hpp>

namespace CORK { namespace krylov {

  class random_initial_vector {
    public:
      inline random_initial_vector( int rank )
      : rank_( rank )
      {}

    public:
      template <typename Triple, typename SolveHandle, typename Backend>
      void apply( Triple& triple, SolveHandle const&, Backend const& backend ) const {
        triple.rank = rank_ ;

        glas2::randomize( triple.Q(glas2::all(), glas2::range(0,rank_)) ) ;
        auto U = triple.u_vector(0) ;
        glas2::fill( U(glas2::all(), glas2::range_from_end(rank_,0)), 0.0 ) ;
        glas2::randomize( U(glas2::all(), glas2::range(0, rank_) ) ) ;

#ifndef NDEBUG
        bool ok =
#endif
        triple.initial_vector( rank_, false ) ;
        assert( ok ) ; // If not OK, there is a problem with Basis4CORK
      }

    private:
      int rank_ ; 
  } ;
   
} } // namespace CORK::krylov


#endif
