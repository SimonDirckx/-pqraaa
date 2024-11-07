//  (C) Copyright Karl Meerbergen 2022.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_approximation_svd_aaa_ls_hpp
#define cork_approximation_svd_aaa_ls_hpp

#include <cork/basis/monomial.hpp>
#include <cork/basis/barycentric_rational_strong.hpp>
#include <cork/approximation/svd_aaa.hpp>
#include <cork/approximation/least_squares.hpp>

namespace CORK { namespace approximation {


  template <typename TestSet, typename TestValues, typename FuncValues, typename FindWorstIndex, typename Options>
  auto SVD_AAA_LS( TestSet const& test_set, TestValues const& test_values, FuncValues const& fun_values, FindWorstIndex const& find_worst_index, int degree, Options const& options ) {
    auto aaa_repr = SVD_AAA( test_set, test_values, find_worst_index, options ) ;

    glas2::shared_vector< typename TestSet::value_type> nodes( copy( aaa_repr.nodes() ) ) ;
    glas2::shared_vector< typename decltype(aaa_repr.weights())::value_type > weights( copy( aaa_repr.weights() ) ) ;
    basis::barycentric_rational_strong< decltype(weights), decltype(nodes) > basis( weights, nodes ) ;

    // Perform least squares step
    return least_squares( test_set, fun_values, basis::monomial(degree), basis, options, true ) ;
  } // SVD_AAA_LS()


} } // namespace CORK::approximation

#endif
