//  (C) Copyright Karl Meerbergen 2021.
//
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_mor_aaa_linearize_hpp
#define cork_mor_aaa_linearize_hpp

#ifndef CORK_CORK2
#include <cork/eigs/cork3.hpp>
#else
#include <cork/eigs/cork2.hpp>
#endif
#include <cork/linearization/union_of_bases.hpp>
#include <cork/linearization/barycentric_rational.hpp>
#include <cork/linearization/combined.hpp>
#include <cork/linearization/cork_linearization.hpp>
#include <cork/linearization/monomial.hpp>
#include <cork/matrix_valued_function/matrix_polynomial_with_shadow.hpp>
#include <cork/matrix_valued_function/sv_aaa.hpp>
#include <cork/matrix_valued_function/norm_est.hpp>
#include <cork/options/value_of.hpp>
#include <cork/options/debug_level.hpp>
#include <cork/utility/value_type_for.hpp>
#include <tuple>
#include <typeinfo>

namespace CORK {

  template <typename NEP, typename R, typename Options>
  decltype (auto) aaa_linearize( NEP const& problem, R const& f_max, Options const& options ) {
    std::vector< std::string > warnings ;

    auto const& nep_krylov = matrix_valued_function::select_nep_from_pair< typename decltype(shift_selector)::value_type, std::false_type >( problem ) ;
    typedef typename std::decay<decltype(nep_krylov)>::type nep_krylov_type ;

    typedef typename nep_krylov_type::template value_type_for< typename decltype(shift_selector)::value_type > value_type ;
    typedef decltype(std::abs(value_type()))                                                                   real_type ;
    typedef std::complex< real_type >                                                                          eig_value_type ;

    auto const& nep_eig = matrix_valued_function::select_nep_from_pair< eig_value_type, std::true_type >( problem ) ;
    typedef typename std::decay<decltype(nep_eig)>::type nep_eig_type ;
    matrix_valued_function::norm_est norm_est(nep_eig) ;

    auto PEP = make_matrix_polynomial_with_shadow( CORK::matrix_valued_function::sv_aaa( nep_krylov, aaa_domain, norm_est, options ), nep_krylov ) ;
    if (options::value_of<options::debug_level>(options)>0)
      std::cout << "Degree of rational AAA polynomial " << PEP.num_terms() << std::endl ;

    eigs::info information ;
    auto linearization = CORK::linearization::make_cork_linearization( PEP, information ) ;

    return linearization
  } // aaa_linearize()

} // namespace CORK

#endif
