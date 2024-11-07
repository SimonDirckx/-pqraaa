//  (C) Copyright Karl Meerbergen 2021.
//
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_eigen_pairs_hpp
#define cork_eigen_pairs_hpp

#include <cork/eigs/cork_eigen_pairs.hpp>
#include <cork/eigs/do_cork.hpp>
#include <cork/eigs/do_qz.hpp>
#include <cork/eigs/qz_eigen_pairs.hpp>
#include <cork/utility/ref.hpp>
#include <cork/utility/tuple.hpp>

namespace CORK {

  // Only for CORK
  template <typename CorkResult>
  auto eigen_pairs( CorkResult result ) {
    return CORK::eigs::cork_eigen_pairs( pop_front(CORK::deref(result)) ) ;
  } // eigen_pairs()

  /* With two arguments:
     template <typename Problem, typename CorkResult>
     auto eigen_pairs( Problem const& problem, CorkResult result ) ;
  */

  // Implementations for CORK and QZ
  namespace detail {
    template <typename Problem, typename CorkResult>
    auto eigen_pairs( eigs::do_cork, Problem const& problem, CorkResult result ) {
      return CORK::eigs::cork_eigen_pairs( problem, result ) ;
    } // eigen_pairs()

    template <typename Problem, typename CorkResult>
    auto eigen_pairs( eigs::do_qz, Problem const& problem, CorkResult result ) {
      return CORK::eigs::qz_eigen_pairs( problem, result ) ;
    } // eigen_pairs()
  }

  template <typename Problem, typename CorkResult>
  auto eigen_pairs( Problem const& problem, CorkResult result ) {
    return detail::eigen_pairs( typename std::tuple_element<0,typename CORK::deref_type<CorkResult>::type>::type(), problem, pop_front(CORK::deref(result)) ) ;
  } // eigen_pairs()

} // namespace CORK


#endif
