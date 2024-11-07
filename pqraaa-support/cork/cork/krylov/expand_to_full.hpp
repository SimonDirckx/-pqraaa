//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_krylov_expand_to_full_hpp
#define cork_krylov_expand_to_full_hpp

#include <cork/krylov/toar_triple.hpp>
#include <glas2/matrix.hpp>
#include <glas2/vector.hpp>

namespace CORK { namespace krylov {

  // Expand Krylov vectors in compact format to full format
  template <typename Triple, typename VV>
  void expand_to_full( Triple const& triple, VV V, int k ) {
    auto n = triple.Q.num_rows() ;
    auto r = triple.rank ;
    auto d = triple.degree() ;
    glas2::all all ;
    glas2::range range_r(0,r) ;
    glas2::range range_k1(0,k+1) ;

    assert( V.num_rows()==n*d ) ;
    assert( V.num_columns()==k+1 ) ;

    for (int i=0; i<d; ++i) {
      V(glas2::range(i*n, i*n+n), all ) = multiply( triple.Q(all,range_r), triple.u_block(i)(range_r,range_k1) ) ;
    }
  } // expand_to_full()


  // Expand Krylov vectors in compact format to full format
  template <typename Triple, typename VV>
  void expand_to_full( Triple const& triple, VV V ) {
    expand_to_full( triple, V, triple.k() ) ;
  } // expand_to_full()
   
} } // namespace CORK::krylov


#endif
