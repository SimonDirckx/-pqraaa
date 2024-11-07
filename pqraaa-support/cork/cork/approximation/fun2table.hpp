//  (C) Copyright Karl Meerbergen 2022.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_approximation_fun2table_hpp
#define cork_approximation_fun2table_hpp

#include <cork/exception/rational_approximation.hpp>
#include <cork/approximation/unique_test_set.hpp>
#include <cork/options/number_of_sample_points.hpp>
#include <cork/options/value_of.hpp>
#include <glas2/matrix.hpp>

namespace CORK { namespace approximation {

  template <typename FunctionSequence, typename Domain, typename Options>
  auto fun2table( FunctionSequence const& fun_sequence, Domain const& domain, Options const& options ) {
    typedef typename Domain::value_type                                               argument_value_type ;
    typedef typename FunctionSequence::template value_type_for< argument_value_type > result_value_type ;

    // Construct the test set
    auto test_set_ini = domain.discretize( std::max(3,options::value_of<options::number_of_sample_points>(options)) ) ;
    auto test_set = unique_test_set( test_set_ini ) ;
    //test_set_ini.resize( 0 ) ;

    // Evaluate and scale the functions
    glas2::shared_matrix<result_value_type> test_values( test_set.size(), fun_sequence.num_terms() ) ;
    for (typename decltype(test_set)::size_type i=0; i<test_set.size(); ++i) {
      fun_sequence.evaluate( test_set(i), test_values(i,glas2::all()) ) ;
      auto ff = norm_inf( test_values(i,glas2::all()) ) ;
      if (std::isnan(ff) || std::isinf(ff)) throw exception::rational_approximation("Functions evaluate to Inf or NaN") ;
    }

    return std::tuple( test_set, test_values ) ;
  } // fun2table()


} } // namespace CORK::approximation

#endif
