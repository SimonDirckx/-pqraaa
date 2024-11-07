//  (C) Copyright Karl Meerbergen 2022.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_approximation_sv_rat_accumulated_real_least_squares_hpp
#define cork_approximation_sv_rat_accumulated_real_least_squares_hpp

#include <cork/approximation/sv_aaa_triple.hpp>
#include <cork/approximation/set_valued_function.hpp>
#include <cork/basis/accumulated_partial_fractions_real.hpp>
#include <cork/basis/monomial.hpp>
#include <cork/basis/union_of_bases.hpp>
#include <cork/options/value_of.hpp>
#include <cork/options/debug_level.hpp>
#include <glas2/matrix.hpp>
#include <glas2/vector.hpp>
#include <glas2/bindings.hpp>
#include <boost/numeric/bindings/lapack/driver/gels.hpp>
#include <cassert>
#include <cmath>
#include <limits>
#include <iomanip>

namespace CORK { namespace approximation {
  

  template <typename SamplePoints, typename FunctionValues, typename Poles, typename Options>
  auto SV_rat_accumulated_real_least_squares( SamplePoints const& test_set, FunctionValues const& test_values, Poles const& poles, int degree, Options const& options ) {
    typedef typename FunctionValues::value_type  value_type ;
    typedef decltype(std::abs(value_type()))     real_value_type ;

    // Count the number of complex points
    int n_points = test_set.size() ;
    for (int i=0;i<test_set.size();++i) {
      n_points += test_set(i).imag()!=0.0 ;
    }
    assert( n_points>poles.size()+1 ) ;

    glas2::shared_vector< real_value_type > weights( poles.size() ) ;
    glas2::shared_vector< typename Poles::value_type> poles_c( poles.size() ) ;
    int i_left = 0 ; int i_right=poles.size() ;
    for (int i=0; i<poles.size(); ++i) {
      if (poles(i).imag()==0.0) {
        poles_c(i_left) = poles(i) ; ++i_left ;
      } else {
        --i_right ; poles_c(i_right) = poles(i+1) ;
        --i_right ; poles_c(i_right) = poles(i) ;
        ++i ;
      }
    }
    assert(i_left==i_right) ;

    std::cout << "poles " << poles_c << std::endl ;
    fill( weights, 1.0 ) ;

    basis::monomial basis1(degree) ;
    basis::accumulated_partial_fractions_real basis2(weights, poles_c) ;
    auto basis = basis::make_union_of_bases( basis1, basis2 ) ;

    // Set the matrix for least squares approximation
    glas2::matrix< real_value_type > cauchy( n_points, basis.num_terms() ) ;
    glas2::vector< value_type > coefs( basis.num_terms() ) ;
    int j_v = 0 ;
    for (int j=0; j<test_set.size(); ++j) {
      basis.evaluate( test_set(j), coefs ) ;
      if (test_set(j).imag()==0) {
        assert( norm_2( glas2::imag( coefs ) ) == 0.0 ) ;
        cauchy( j_v, glas2::all() ) = glas2::real( coefs ) ;
      } else {
        cauchy( j_v, glas2::all() ) = glas2::real( coefs ) ;
        ++j_v ;
        cauchy( j_v, glas2::all() ) = glas2::imag( coefs ) ;
      }
      ++j_v ;
    }
    assert( j_v==cauchy.num_rows() ) ;

    // Determine optimal weights
    // First column are all ones.
    for (int i=degree+1; i<cauchy.num_columns(); ++i) {
      real_value_type nrm = norm_inf( cauchy(glas2::all(), i) ) ;
      if (poles_c(i-degree-1).imag()!=0.) nrm += norm_inf( cauchy(glas2::all(), i+1) ) ;

      if (nrm!=0) basis.basis_2().weights()(i-degree-1) = 1./nrm ;
      cauchy(glas2::all(), i) *= basis.basis_2().weights()(i-degree-1) ;
      if (poles_c(i-1).imag()!=0.) {
        ++i ;
        basis.basis_2().weights()(i-degree-1) = basis.basis_2().weights()(i-degree-2) ;
        cauchy(glas2::all(), i) *= basis.basis_2().weights()(i-degree-1) ;
      }
    }
    std::cout << "weights " << basis.basis_2().weights() << std::endl ;

    // Set right-hand sides
    glas2::matrix< real_value_type > rhs( n_points, test_values.num_columns() ) ;
    j_v = 0 ;
    for (int j=0; j<test_set.size(); ++j) {
      rhs(j_v, glas2::all()) = glas2::real( test_values(j, glas2::all()) ) ;
      if (test_set(j).imag()!=0.0) {
        ++j_v ;
        rhs(j_v, glas2::all()) = glas2::imag( test_values(j, glas2::all()) ) ;
      }
      ++j_v ;
    }

    // Solve least squares problem
    //std::cout << "cauchy " << cauchy << std::endl ;
    //std::cout << "rhs " << rhs << std::endl ;
#ifndef NDEBUG
    int info =
#endif
      boost::numeric::bindings::lapack::gels( cauchy, rhs ) ;
    real_value_type res_norm = glas2::norm_fro(rhs(glas2::range_from_end(cauchy.num_columns(),0), glas2::all())) ;
    if (options::value_of<options::debug_level>(options)>1)
      std::cout << "AAA-LS residual norm " << res_norm << std::endl ;
    //std::cout << "sol " << rhs << std::endl ;
    assert( info==0 ) ;

    glas2::matrix< real_value_type > coefficients( poles.size()+1, rhs.num_columns() ) ;

    return set_valued_function< value_type, decltype(basis) >( basis, rhs( glas2::range(0,poles_c.size()+1), glas2::all() ) ) ;
  } // SV_rat_accumulated_real_least_squares()


  struct SV_rat_accumulated_real_least_squares_selector {
    template <typename SamplePoints, typename FunctionValues, typename Poles, typename Options>
    auto operator() ( SamplePoints const& test_set, FunctionValues const& test_values, Poles const& poles, int degree, Options const& options ) const {
      return SV_rat_accumulated_real_least_squares( test_set, test_values, poles, degree, options ) ;
    }
  } ; // struct SV_rat_accumulated_real_least_squares_selector

} } // namespace CORK::approximation

#endif
