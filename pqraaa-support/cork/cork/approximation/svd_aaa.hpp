//  (C) Copyright Karl Meerbergen 2018.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_approximation_svd_aaa_hpp
#define cork_approximation_svd_aaa_hpp

#include <cork/exception/rational_approximation.hpp>
#include <cork/approximation/sv_aaa.hpp>
#include <cork/approximation/unique_test_set.hpp>
#include <cork/options/max_degree.hpp>
#include <cork/options/debug_level.hpp>
#include <cork/options/number_of_sample_points.hpp>
#include <cork/options/value_of.hpp>
#include <cork/options/aaa_cheap_correction_tolerance.hpp>
#include <cork/options/aaa_function_drop_tol.hpp>
#include <cork/utility/matlab.hpp>
#include <cork/options/value_of.hpp>
//#include <fstream>

namespace CORK { namespace approximation {

  template <typename Points, typename FunctionValues, typename FindWorstIndex, typename Options>
  typename std::enable_if< glas2::is<glas2::DenseMatrix, FunctionValues>::value
                         , SV_AAA_approximation<typename Points::value_type, typename FunctionValues::value_type>
                         >::type SVD_AAA( Points const& test_set, FunctionValues test_values, FindWorstIndex const& find_worst_index, Options const& options, bool scale_functions=true ) {
    typedef typename Points::value_type         argument_value_type ;
    typedef typename FunctionValues::value_type result_value_type ;
    typedef decltype(std::abs(result_value_type()))    real_type ;

    assert( test_set.size()>=3 ) ;
    assert( test_values.num_rows()==test_set.size() ) ;

    glas2::shared_matrix<result_value_type> function_scaling(test_values.num_columns(), test_values.num_columns() ) ;
    function_scaling = glas2::eye(test_values.num_columns(), test_values.num_columns() );

    glas2::vector<real_type> norm_f( test_values.num_columns() ) ;

    if (scale_functions) {
    // Evaluate and scale the functions
    for (typename FunctionValues::size_type j=0; j<test_values.num_columns(); ++j) {
        norm_f(j) = norm_inf( test_values(glas2::all(),j) ) ;
        if (std::isnan(norm_f(j)) || std::isinf(norm_f(j))) throw exception::rational_approximation("Functions evaluate to Inf or NaN") ;
        if (norm_f(j)!=0.0) test_values(glas2::all(),j) /= norm_f(j) ;
      /*  if (norm_1(test_values(glas2::all(),j))<1.1) {
            std::stringstream ss ;
            ss<< "AAA: possible singularity in point for function number " << j ;
            throw std::runtime_error( ss.str() ) ;
        }*/
      }
    }

    real_type abs_tol = options::value_of<options::aaa_stop_criterion<real_type>>(options).tolerance() ;

    int n_functions = test_values.num_columns();
    // Evaluate the functions
    {
      glas2::vector<real_type>  svd_val(test_values.num_columns()) ;
      int info = boost::numeric::bindings::lapack::gesvd( 'O', 'A', test_values, svd_val, function_scaling, function_scaling ) ;
      assert(info>=0) ;
      if (info>0) throw exception::rational_approximation("GESVD failed in determining SVD of functions values: use aaa_algorithm::sv_aaa() instead") ;
      for ( ; svd_val(n_functions-1)<abs_tol; --n_functions ) ;
      //for (int i=0; i<n_functions; ++i) function_scaling(i,glas2::all()) *= svd_val(i) ;
      for (int i=0; i<n_functions; ++i) test_values(glas2::all(),i) *= svd_val(i) ;
    }
    if (options::value_of<options::debug_level>(options)>0) std::cout << "Number of functions is reduced from " << test_values.num_columns() << " to " << n_functions << std::endl ;
    /*{
    std::ofstream ff("v2.m") ;
    ff << CORK::matlab( test_values( glas2::all(), glas2::range(0,n_functions) ), "v2" ) ;
    }*/

    // Change the find_worst_index for the SVD
    auto find_worst_index_svd = find_worst_index.transform( transpose(function_scaling(glas2::range(0,n_functions),glas2::all())) ) ;
    auto repr_r = SV_AAA( test_set, test_values( glas2::all(), glas2::range(0,n_functions) ), find_worst_index_svd, options, false ) ;

    SV_AAA_approximation<argument_value_type,result_value_type> repr( test_values.num_columns(), repr_r.n() ) ;
    repr.n() = repr_r.n() ;
    repr.nodes() = repr_r.nodes() ;
    repr.weights() = repr_r.weights() ;
    repr.coefficients() = multiply( repr_r.coefficients(), function_scaling(glas2::range(0,n_functions),glas2::all()) ) ;
    if (scale_functions) {
      for (int i=0; i<test_values.num_columns(); ++i) {
        repr.coefficients()(glas2::all(),i) *= norm_f(i) ;
      }
    }

    return repr ;
  } // SVD_AAA()


  /*
  template <typename FunctionSequence, typename Domain, typename FindWorstIndex, typename Options>
  typename std::enable_if< !glas2::is< glas2::DenseVector, FunctionSequence>::value
                         , SV_AAA_approximation<typename Domain::value_type, typename FunctionSequence::template value_type_for<typename Domain::value_type> >
                         >::type SVD_AAA( FunctionSequence const& fun_sequence, FindWorstIndex const& find_worst_index, Domain const& domain,  Options const& options ) {
    typedef typename Domain::value_type                                               argument_value_type ;
    typedef typename FunctionSequence::template value_type_for< argument_value_type > result_value_type ;

    // Construct the test set
    auto test_set_ini = domain.discretize( std::max(3,options::value_of<options::number_of_sample_points>(options)) ) ;
    auto test_set = unique_test_set( test_set_ini ) ;
    //test_set_ini.resize( 0 ) ;

    // Evaluate and scale the functions
    glas2::matrix<result_value_type> test_values( test_set.size(), fun_sequence.num_terms() ) ;
    for (typename decltype(test_set)::size_type i=0; i<test_set.size(); ++i) {
      fun_sequence.evaluate( test_set(i), test_values(i,glas2::all()) ) ;
    }
    return SVD_AAA( test_set, test_values, find_worst_index, options ) ;
  } // SVD_AAA()


  template <typename FunctionSequence, typename FindWorstIndex, typename Domain>
  typename std::enable_if< !glas2::is< glas2::DenseVector, FunctionSequence>::value
                         , SV_AAA_approximation<typename Domain::value_type, typename FunctionSequence::template value_type_for<typename Domain::value_type> >
                         >::type SVD_AAA( FunctionSequence const& fun_sequence, FindWorstIndex const& find_worst_index, Domain const& domain ) {
    return SVD_AAA( fun_sequence, find_worst_index, domain, std::tuple<>() ) ;
  } // SVD_AAA
  */

  template <typename Points, typename FunctionValues, typename FindWorstIndex>
  typename std::enable_if< glas2::is< glas2::DenseVector, Points>::value
                         && glas2::is< glas2::DenseMatrix, FunctionValues>::value
                         , SV_AAA_approximation<typename Points::value_type, typename FunctionValues::value_type>
                         >::type SVD_AAA( Points const& points, FunctionValues const& function_values, FindWorstIndex const& find_worst_index ) {
    return SVD_AAA( points, function_values, find_worst_index, std::tuple<>() ) ;
  } // SVD_AAA

} } // namespace CORK::approximation

#endif
