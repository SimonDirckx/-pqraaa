//  (C) Copyright Karl Meerbergen 2021.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_approximation_sv_aaa_2_rational_newton_hpp
#define cork_approximation_sv_aaa_2_rational_newton_hpp

#include <cork/approximation/sv_aaa_approximation.hpp>
#include <cork/approximation/sv_aaa_poles.hpp>
#include <cork/approximation/potential_theory_approximation.hpp>
#include <cork/utility/matlab.hpp>
#include <cork/lapack/eig.hpp>
#include <glas2/vector.hpp>
#include <glas2/matrix.hpp>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace CORK { namespace approximation {

  template <typename AAAApproximation, typename DiscreteDomain, typename Options>
  decltype(auto) sv_aaa_2_rational_newton( AAAApproximation const& aaa, DiscreteDomain const& domain, Options const& options ) {
    typedef typename AAAApproximation::value_type            value_type ;
    typedef decltype(std::abs(value_type()))         real_value_type ;
    typedef std::complex< real_value_type >          complex_value_type ;

    basis::rational_newton rat_newton( aaa.nodes(), aaa.poles(), scalings ) ;

    potential_theory_scaling( domain, repr.basis() ) ;
    return repr ;
  } // sv_aaa_2_rational_newton()

} } // namespace CORK::approximation

#endif
