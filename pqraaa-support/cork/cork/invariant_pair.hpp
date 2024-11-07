//  (C) Copyright Karl Meerbergen 2021.
//
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_invariant_pair_hpp
#define cork_invariant_pair_hpp

#include <cork/eigs/do_cork.hpp>
#include <cork/eigs/do_qz.hpp>
#include <cork/eigs/cork_invariant_pair.hpp>
#include <cork/eigs/qz_invariant_pair.hpp>
#include <cork/utility/ref.hpp>
#include <cork/utility/tuple.hpp>

namespace CORK {

 template <typename CorkResult>
 auto invariant_pair( CorkResult& cork_result ) {
   return CORK::eigs::cork_invariant_pair( pop_front(CORK::deref(cork_result)) ) ;
 } // invariant_pair()


  // Implementations for CORK and QZ
  namespace detail {
    template <typename Problem, typename CorkResult>
    auto invariant_pair( eigs::do_cork, Problem const& problem, CorkResult result ) {
      return CORK::eigs::cork_invariant_pair( problem, result ) ;
    } // invariant_pair()

    template <typename Problem, typename CorkResult>
    auto invariant_pair( eigs::do_qz, Problem const& problem, CorkResult result ) {
      return CORK::eigs::qz_invariant_pair( problem, result ) ;
    } // invariant_pair()
  }

  template <typename Problem, typename CorkResult>
  auto invariant_pair( Problem const& problem, CorkResult result ) {
    return detail::invariant_pair( typename std::tuple_element<0,typename CORK::deref_type<CorkResult>::type>::type(), problem, pop_front(CORK::deref(result)) ) ;
  } // invariant_pair()

} // namespace CORK


#endif
