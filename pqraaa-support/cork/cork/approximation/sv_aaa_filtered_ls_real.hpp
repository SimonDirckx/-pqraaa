//  (C) Copyright Karl Meerbergen 2022.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_approximation_sv_aaa_filtered_ls_real_hpp
#define cork_approximation_sv_aaa_filtered_ls_real_hpp

#include <cork/approximation/sv_aaa_real.hpp>
#include <cork/approximation/aaa_real_poles.hpp>
#include <cork/approximation/real_least_squares.hpp>
#ifdef CORK_USE_INVERSE_NEWTON
#include <cork/basis/accumulated_partial_fractions_real.hpp>
#else
#include <cork/basis/partial_fractions_real.hpp>
#endif
#include <cmath>

namespace CORK { namespace approximation {


  template <typename TestSet, typename AAATestValues, typename TestValues, typename FindWorstIndex, typename FilterPoles, typename Options>
  auto SV_AAA_filtered_LS_real_discretized( TestSet const& test_set, AAATestValues const& aaa_test_values, TestValues const& test_values, FindWorstIndex const& find_worst_index, FilterPoles const& filter_poles, int degree, Options const& options ) {
    auto aaa_repr = SV_AAA_real( test_set, aaa_test_values, find_worst_index, options ) ;

    // Compute poles
    auto poles = aaa_real_poles( aaa_repr ) ;

    typedef typename decltype(poles)::value_type value_type ;
    typedef decltype(std::abs(value_type())) real_value_type ;

    if (options::value_of<options::debug_level>(options)>2) std::cout << "AAA-LS: Poles before filtering " << poles << std::endl ;
    auto filtered_poles = filter_poles( poles ) ;
    if (options::value_of<options::debug_level>(options)>2) std::cout << "AAA-LS: Poles after filtering " << filtered_poles << std::endl ;
    // Sort: real poles come first
    {
      glas2::vector<int> ind( filtered_poles.size() ) ; ind =  glas2::range(0,filtered_poles.size()) ;
      glas2::vector< real_value_type > poles_i( glas2::copy(glas2::abs(imag(filtered_poles))) ) ;
      glas2::sort( poles_i, ind ) ;
      glas2::vector<value_type> copy_poles( glas2::copy(filtered_poles) ) ;
      filtered_poles = copy_poles(ind) ;
    }

    // Determine weights
    glas2::shared_vector<value_type> copy_poles( glas2::copy(filtered_poles) ) ;
    glas2::shared_vector<real_value_type> weights( copy_poles.size() ) ;
    fill( weights, 1.0 ) ;
#ifdef CORK_USE_INVERSE_NEWTON
    basis::accumulated_partial_fractions_real basis2(weights, copy_poles) ;
#else
    CORK::basis::partial_fractions_real basis2(copy_poles,weights) ;
#endif
    {
      glas2::vector<real_value_type> norms( weights.size() ) ;
      glas2::vector<value_type> eval( 1+copy_poles.size() ) ;
      fill( norms, 0.0 ) ;
      for (int i=0; i<test_set.size(); ++i) {
        basis2.evaluate( test_set(i), eval ) ;
        for (int j=0; j<norms.size(); ++j) {
          auto eval_abs = std::abs(eval(j+1)) ;
          if (eval_abs>norms(j)) norms(j) = eval_abs ;
        }
      }
      for (int i=0; i<copy_poles.size(); ++i) {
        if (std::imag(copy_poles(i))!=0.) {
          assert( std::abs(copy_poles(i)-std::conj(copy_poles(i+1))) == 0.0 ) ;
          norms(i) = norms(i+1) = std::max( norms(i), norms(i+1) ) ;
          ++i ;
        }
      }
      weights = 1./norms ;
    }

    // Perform least squares step
    return real_least_squares( test_set, test_values, basis::monomial(degree), basis2, options, true ) ;
  } // SV_AAA_filtered_LS_real_discretized()

  
  template <typename AAAFunctionSequence, typename FunctionSequence, typename Domain, typename FindWorstIndex, typename FilterPoles, typename Options>
  auto SV_AAA_filtered_LS_real( AAAFunctionSequence const& aaa_fun_sequence, FunctionSequence const& fun_sequence, FindWorstIndex const& find_worst_index, FilterPoles const& filter_poles, Domain const& domain, int degree, Options const& options ) {
    typedef typename Domain::value_type                                               argument_value_type ;
    typedef typename FunctionSequence::template value_type_for< argument_value_type > value_type ;

    // Construct the test set
    auto test_set_ini = domain.discretize( std::max(2*options::value_of<options::max_degree>(options),
                                                 options::value_of<options::number_of_sample_points>(options)) ) ;
    auto test_set1 = unique_test_set( test_set_ini, []( auto const& x, auto const& y ) { return x==y || std::conj(x)==y ; } ) ;

    // Evaluate and scale the functions
    glas2::matrix<value_type> aaa_test_values( test_set1.size(), aaa_fun_sequence.num_terms() ) ;
    glas2::matrix<value_type> test_values( test_set1.size(), fun_sequence.num_terms() ) ;

    for (typename glas2::vector<value_type>::size_type i=0; i<test_set1.size(); ++i) {
      aaa_fun_sequence.evaluate( test_set1(i), aaa_test_values(i,glas2::all()) ) ;
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

    return SV_AAA_filtered_LS_real_discretized( test_set1, aaa_test_values, test_values, find_worst_index, filter_poles, degree, options ) ;
  } // SV_AAA_filtered_LS_real()


  struct SV_AAA_filtered_LS_real_selector {
    template <typename TestSet, typename AAATestValues, typename TestValues, typename FindWorstIndex, typename Filter, typename Options>
    auto operator() ( TestSet const& test_set, AAATestValues const& aaa_test_values, TestValues const& test_values, FindWorstIndex const& find_worst_index, Filter const& filter, int degree, Options const& options ) const {
      return SV_AAA_filtered_LS_real_discretized( test_set, aaa_test_values, test_values, find_worst_index, filter, degree, options ) ;
    }
  } ; // struct SV_AAA_filtered_LS_real_discretized_selector

  
} } // namespace CORK::approximation

#endif
