//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_matrix_valued_function_difference_hpp
#define cork_matrix_valued_function_difference_hpp

#include <cassert>
#include <cork/matrix_valued_function/matrix_polynomial.hpp>
#include <cork/matrix_valued_function/matrix_polynomial.hpp>
#include <cork/basis/union_of_bases.hpp>
#include <cork/basis/barycentric_rational.hpp>
#include <cork/options/debug_level.hpp>
#include <cork/options/value_of.hpp>
#include <cork/options/number_of_sample_points.hpp>
#include <cork/coefficient_matrices/combined.hpp>
#include <cork/utility/norm_est.hpp>
#include <cmath>

namespace CORK { namespace matrix_valued_function {

  template <typename NonlinearMatrixOne, typename NonlinearMatrixTwo, typename Domain, typename Options>
  decltype (auto) difference ( NonlinearMatrixOne const& one, NonlinearMatrixTwo const& two, Domain const& domain, Options const& options, int n_points=0 ) {
    typedef typename Domain::value_type        argument_type ;
    typedef typename std::common_type< typename NonlinearMatrixOne::template value_type_for<argument_type>
                                     , typename NonlinearMatrixTwo::template value_type_for<argument_type>
                                     >::type   value_type ;
    typedef decltype( std::abs(value_type()) ) real_type ;

    if (!n_points) n_points = options::value_of<options::number_of_sample_points>(options) ;

    assert( one.size()==two.size() ) ;
    auto points = domain.discretize(  n_points ) ;
    real_type nrm_one ;
    real_type max_diff = 0.0 ;
      /*glas2::vector<value_type> w( one.size() ) ;
      glas2::vector<value_type> v( one.size() ) ;
      fill(v,1.0) ;*/
      //std::cout << one.coefficient_matrices().combinations() << std::endl;

    for (int i=0; i<points.size(); ++i) {
/*      fill(w,0.0) ;
      one.multiply_add( points(i), v, w ) ;
      std::cout << "one: " << w << std::endl ;
      fill(w,0.0) ;
      two.multiply_add( points(i), v, w ) ;
      std::cout << "two: " << w << std::endl ;
*/
      auto nrm = CORK::norm_est<value_type>( one.size(), [&](auto const& v, auto& w ) {
                fill( w, 0.0 ) ;
                one.multiply_add( points(i), v, w ) ;
          } ) ;
      nrm_one = nrm.first ;
      nrm = CORK::norm_est<value_type>( one.size(), [&](auto const& v, auto& w ) {
                fill( w, 0.0 ) ;
                one.multiply_add( points(i), v, w ) ;
                w *= -1.0 ;
                two.multiply_add( points(i), v, w ) ;
          } ) ;
      if (options::value_of<options::debug_level>(options)>1) std::cout << "Error on rational approximation for " << points(i) << " --> " << nrm.first << " / " << nrm_one << std::endl ;
      if (nrm.first/nrm_one>max_diff) max_diff = nrm.first/nrm_one ;
    }
    return max_diff ;
  } // difference()

} } // namespace CORK::matrix_valued_function

#endif
