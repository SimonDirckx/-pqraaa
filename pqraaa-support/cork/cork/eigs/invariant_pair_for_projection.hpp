//  (C) Copyright Karl Meerbergen 2020.
//
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_eigs_invariant_pair_for_projection_hpp
#define cork_eigs_invariant_pair_for_projection_hpp

#include <cork/eigs/invariant_pair.hpp>

namespace CORK { namespace eigs {

  //
  // Invariant pair for explicit projection
  //
  template <typename T>
  struct invariant_pair_for_projection_type {
  } ;

  template <typename T, typename Options, typename PT>
  struct invariant_pair_for_projection_type< krylov::cork_quadruple<T,Options,PT> > {
    typedef typename krylov::cork_quadruple<T,Options,PT>::sharedmatrix_type matrix_type ;
    typedef invariant_pair< matrix_type const, matrix_type const, matrix_type const > type ;
  } ;

  template <typename T, typename Options, typename InnerProd>
  struct invariant_pair_for_projection_type< krylov::toar_triple<T,Options,InnerProd> > {
    typedef typename krylov::toar_triple<T,Options,InnerProd>::sharedmatrix_type matrix_type ;
    typedef glas2::shared_matrix< T > k_type ;
    typedef invariant_pair< matrix_type const, k_type const, matrix_type const > type ;
  } ;

  template <typename T, typename Options, typename PT>
  auto make_invariant_pair_for_projection( krylov::cork_quadruple<T,Options,PT> const& quad, int order ) {
    return typename invariant_pair_for_projection_type< krylov::cork_quadruple<T,Options,PT> >::type( quad.Q, quad.K, quad.H, order ) ;
  }


  template <typename T, typename Options, typename InnerProd>
  auto make_invariant_pair_for_projection( krylov::toar_triple<T,Options,InnerProd> const& quad, int order ) {
    glas2::shared_matrix< T > K( quad.H.num_columns(), quad.H.num_columns() ) ;
    K = glas2::identity_matrix< T >(K.num_rows(), K.num_columns) ;
    return typename invariant_pair_for_projection_type< krylov::toar_triple<T,Options,InnerProd> >::type( quad.Q, K, quad.H, order ) ;
  }

} } // namespace CORK::eigs

#endif
