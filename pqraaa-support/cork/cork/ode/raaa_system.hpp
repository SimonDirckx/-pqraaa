//  (C) Copyright Karl Meerbergen 2021.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_matrix_ode_raaa_system_hpp
#define cork_matrix_ode_raaa_system_hpp

#include <cork/linearization/union_of_bases.hpp>
#include <cork/linearization/barycentric_rational_real_strong.hpp>
#include <cork/linearization/combined.hpp>
#include <cork/linearization/cork_linearization.hpp>
#include <cork/linearization/monomial.hpp>
#include <cork/matrix_valued_function/select_nep_from_pair.hpp>
#include <cork/matrix_valued_function/norm_est.hpp>
#include <cork/matrix_valued_function/sv_raaa.hpp>
#include <cork/matrix_valued_function/use_shadowed_nonlinear_matrix.hpp>
#include <cork/matrix_valued_function/matrix_polynomial_with_shadow.hpp>
#include <cork/eigs/info.hpp>
#include <memory>

namespace CORK { namespace ode {

  //!
  //! Basis: type of polynomial basis, can be evaluated in any value_type in principle, unless the basis does not support this
  //! CoefficientMatrices: matrix coefficients, can have any value_type
  //! T: is the type of the matrix polynomial. Usually this is determined by the type of the matrices, but sometimes, the basis may have a specific
  //!    type so that T should be different. This allows the user to set the type manually if it is different from the default.
  //!

  template <typename ValueType, typename NEP, typename Domain, typename Options>
  auto raaa_system( NEP const& problem, Domain const& domain, Options const& options ) {
    std::vector< std::string > warnings ;

    auto const& nep_time = matrix_valued_function::select_nep_from_pair< ValueType, std::false_type >( problem ) ;

    auto const& nep_aaa = matrix_valued_function::select_nep_from_pair< typename Domain::value_type, std::true_type >( problem ) ;
    matrix_valued_function::norm_est norm_est(nep_aaa) ;

    auto PEP_aaa = CORK::matrix_valued_function::sv_raaa( nep_time, nep_aaa, domain, norm_est, options ) ;

    auto PEP = matrix_valued_function::use_shadowed_nonlinear_matrix( options, PEP_aaa, nep_time ) ;
//   auto PEP = make_matrix_polynomial_with_shadow( CORK::matrix_valued_function::sv_raaa( nep_time, nep_aaa, domain, norm_est, options ), nep_time ) ;
//    auto PEP = CORK::matrix_valued_function::sv_raaa( nep_time, nep_aaa, domain, norm_est, options ) ;
    if (options::value_of<options::debug_level>(options)>0)
      std::cout << "Degree of rational AAA polynomial " << PEP.num_terms() << std::endl ;

    std::shared_ptr<linearization::info> information( new linearization::info ) ;
    CORK::linearization::cork_linearization<typename std::decay<decltype(PEP)>::type> linearization( PEP, *information ) ;

    return std::tuple( linearization, information ) ;
  } // raaa_system()

} } // namespace CORK::ode

#endif
