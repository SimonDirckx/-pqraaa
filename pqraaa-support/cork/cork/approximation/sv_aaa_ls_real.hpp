//  (C) Copyright Karl Meerbergen 2021.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_approximation_sv_aaa_ls_real_hpp
#define cork_approximation_sv_aaa_ls_real_hpp

#include <cork/approximation/sv_aaa_real.hpp>
#include <cork/approximation/real_least_squares.hpp>
#include <cork/basis/barycentric_rational_real_strong.hpp>

namespace CORK { namespace approximation {


  template <typename TestSet, typename TestValues, typename FuncValues, typename FindWorstIndex, typename Options>
  auto SV_AAA_LS_real_discretized( TestSet const& test_set, TestValues const& test_values, FuncValues const& fun_values, FindWorstIndex const& find_worst_index, int degree, Options const& options ) {
    auto aaa_repr = SV_AAA_real( test_set, test_values, find_worst_index, options ) ;

    glas2::shared_vector< typename TestSet::value_type> nodes( copy( aaa_repr.nodes() ) ) ;
    glas2::shared_vector< typename decltype(aaa_repr.weights())::value_type > weights( copy( aaa_repr.weights() ) ) ;
    basis::barycentric_rational_real_strong basis( weights, nodes ) ;

    // Perform least squares step
    return real_least_squares( test_set, fun_values, basis::monomial(degree), basis, options, false ) ;
  } // SV_AAA_LS_real_discretized()

  
  template <typename AAAFunctionSequence, typename FunctionSequence, typename Domain, typename FindWorstIndex, typename Options>
  auto SV_AAA_LS_real( AAAFunctionSequence const& aaa_fun_sequence, FunctionSequence const& fun_sequence, FindWorstIndex const& find_worst_index, Domain const& domain, int degree, Options const& options ) {
    typedef typename Domain::value_type                                               argument_value_type ;
    typedef typename FunctionSequence::template value_type_for< argument_value_type > value_type ;

    // Construct the test set
    auto test_set_ini = domain.discretize( std::max(2*options::value_of<options::max_degree>(options),
                                                 options::value_of<options::number_of_sample_points>(options)) ) ;
    auto test_set1 = unique_test_set( test_set_ini, []( auto const& x, auto const& y ) { return x==y || std::conj(x)==y ; } ) ;

    // Evaluate and scale the functions
    glas2::matrix<value_type> aaa_values( test_set1.size(), aaa_fun_sequence.num_terms() ) ;
    glas2::matrix<value_type> test_values( test_set1.size(), fun_sequence.num_terms() ) ;

    for (typename glas2::vector<value_type>::size_type i=0; i<test_set1.size(); ++i) {
      aaa_fun_sequence.evaluate( test_set1(i), aaa_values(i,glas2::all()) ) ;
      if (options::value_of<options::debug_level>(options)>0) {
        if (is_infinite(glas2::norm_2(test_values(i,glas2::all()))))
          throw exception::rational_approximation( "AAA approximation: infinite function values in sample point" ) ;
        if (std::isnan(glas2::norm_2(test_values(i,glas2::all()))))
          throw exception::rational_approximation( "AAA approximation: NaN values in sample point" ) ;
      }
      fun_sequence.evaluate( test_set1(i), test_values(i,glas2::all()) ) ;
      if (options::value_of<options::debug_level>(options)>0) {
        if (is_infinite(glas2::norm_2(test_values(i,glas2::all()))))
          throw exception::rational_approximation( "AAA approximation: infinite function values in sample point" ) ;
        if (std::isnan(glas2::norm_2(test_values(i,glas2::all()))))
          throw exception::rational_approximation( "AAA approximation: NaN values in sample point" ) ;
      }
    }

    return SV_AAA_LS_real_discretized( test_set1, aaa_values, test_values, find_worst_index, degree, options ) ;
  } // SV_AAA_LS_real()


  struct SV_AAA_LS_real_selector {
    template <typename TestSet, typename AAAValues, typename TestValues, typename FindWorstIndex, typename Filter, typename Options>
    auto operator() ( TestSet const& test_set, AAAValues const& aaa_values, TestValues const& test_values, FindWorstIndex const& find_worst_index, Filter const& filter, int degree, Options const& options ) const {
      return SV_AAA_LS_real_discretized( test_set, aaa_values, test_values, find_worst_index, degree, options ) ;
    }
  } ; // struct SV_AAA_filtered_LS_real_discretized_selector

  
} } // namespace CORK::approximation

#endif
