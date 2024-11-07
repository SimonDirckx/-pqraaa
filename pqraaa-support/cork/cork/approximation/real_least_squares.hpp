//  (C) Copyright Karl Meerbergen 2021.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_approximation_real_least_squares_hpp
#define cork_approximation_real_least_squares_hpp

#include <cork/approximation/sv_aaa_triple.hpp>
#include <cork/approximation/set_valued_function.hpp>
#include <cork/basis/partial_fractions_real.hpp>
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
  

  template <typename SamplePoints, typename FunctionValues, typename Poly, typename Rat, typename Options>
  auto real_least_squares( SamplePoints const& test_set, FunctionValues const& test_values, Poly const& polynomial, Rat const& rat, Options const& options, bool with_constant_term ) {
    typedef typename FunctionValues::value_type  value_type ;
    typedef decltype(std::abs(value_type()))     real_value_type ;

    // Count the number of complex points
    int n_points = test_set.size() ;
    for (int i=0;i<test_set.size();++i) {
      n_points += test_set(i).imag()!=0.0 ;
    }

    auto basis = basis::make_union_of_bases( polynomial, rat ) ;

    // Set the matrix for least squares approximation
    glas2::matrix< real_value_type > cauchy( n_points, basis.num_terms()-1+with_constant_term ) ;
    glas2::vector< value_type > coefs( basis.num_terms() ) ;
    glas2::range_from_end r(1-with_constant_term, 0) ;
    int j_v = 0 ;
    for (int j=0; j<test_set.size(); ++j) {
      basis.evaluate( test_set(j), coefs ) ;
      if (test_set(j).imag()==0) {
        assert( norm_2( glas2::imag( coefs ) ) == 0.0 ) ;
        cauchy( j_v, glas2::all() ) = glas2::real( coefs(r) ) ;
      } else {
        cauchy( j_v, glas2::all() ) = glas2::real( coefs(r) ) ;
        ++j_v ;
        cauchy( j_v, glas2::all() ) = glas2::imag( coefs(r) ) ;
      }
      ++j_v ;
    }
    assert( j_v==cauchy.num_rows() ) ;

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

    if (!with_constant_term) {
      glas2::matrix< real_value_type > temp( basis.num_terms()-1, rhs.num_columns() ) ;
      temp = rhs( glas2::range(0, basis.num_terms()-1), glas2::all() ) ;
      rhs( glas2::range(1, basis.num_terms()), glas2::all() ) = temp ;
      fill( rhs( 0, glas2::all() ), 0. ) ;
    }
    return set_valued_function< real_value_type, decltype(basis) >( basis, rhs( glas2::range(0,basis.num_terms()), glas2::all() ) ) ;
  } // real_least_squares()

} } // namespace CORK::approximation

#endif
