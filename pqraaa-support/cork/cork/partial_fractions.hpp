//  (C) Copyright Karl Meerbergen 2021.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_partial_fractions_hpp
#define cork_partial_fractions_hpp

#include <cork/vector.hpp>
#include <cork/basis/union_of_bases.hpp>
#include <cork/basis4cork/union_of_bases.hpp>
#include <cork/matrix_iterator/union_of_bases.hpp>
#include <cork/basis/monomial.hpp>
#include <cork/basis4cork/monomial.hpp>
#include <cork/matrix_iterator/monomial.hpp>
#include <cork/basis/partial_fractions.hpp>
#include <cork/basis4cork/partial_fractions.hpp>
#include <cork/matrix_iterator/partial_fractions.hpp>

namespace CORK {

  template <typename P, typename W>
  auto partial_fractions( int degree, P const& poles, W const& weights ) {
    basis::monomial monomial( degree ) ;
    basis::partial_fractions< P const&, W const& > pf( poles, weights ) ;
    return basis::make_union_of_bases( monomial, pf ) ;
  }

} // namespace CORK::basis

#endif
