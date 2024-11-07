//  (C) Copyright Karl Meerbergen 2018.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_approximation_symm_aaa_hpp
#define cork_approximation_symm_aaa_hpp

#include <cork/approximation/aaa.hpp>
#include <cork/approximation/aaa_approximation.hpp>
#include <cork/domain/upper_half_plane.hpp>
//#include <glas2/matrix.hpp>
//#include <glas2/vector.hpp>
//#include <glas2/bindings.hpp>
//#include <boost/numeric/bindings/lapack/driver/gesvd.hpp>
//#include <cassert>

namespace CORK { namespace approximation {

  template <typename T, typename F, typename Domain>
  aaa_approximation<T> symm_AAA( F const& fun, Domain const& domain, aaa_options<T> const& options=aaa_options<T>() ) {
    using namespace glas2::ops ;

    typedef T                                value_type ;
//    typedef decltype(std::abs(value_type())) real_type ;

    // Construct the test set for the upper part of the domain and make the AAA approximation
    glas2::vector<value_type> test_set( std::max(3,options.n_points) ) ;
    glas2::vector<value_type> test_values( test_set.size() ) ;
    aaa_test_set( fun, make_upper_half_plane(domain), test_set, test_values ) ;

    glas2::vector<int> ind(options.max_degree) ;
    aaa_approximation<T> repr_upper = aaa_approximate<T>( test_set, test_values, ind, options ) ;

    // We copy the nodes to the lower part and recompute the weights.
/*    int n_real = sum( imag(test_set)==glas2::zeros<real_type>(test_set.size()) ) ;
    int n_complex = test_set.size() - n_real ;
    int d_real = sum( imag(repr_upper.nodes())==0. ) ;
    int d_complex = repr_upper.n() - d_real ;

    glas2::vector<value_type> full_test_set( n_real-d_real + 2*(n_complex-d_complex) ) ;
    glas2::vector<value_type> full_test_values( full_test_set.size() ) ;
    glas2::vector<bool> is_not_node( test_set.size() ) ;
    fill( is_not_node, true ) ;
    fill( is_not_node(ind(glas2::range(0,repr_upper.n()))), false ) ;

    typename decltype(full_test_set)::size_type i_full = 0 ;
    for (typename decltype(test_set)::size_type i=0; i<test_set.size(); ++i) {
      if (is_not_node(i)) {
        if (imag(test_set(i))==0.0) {
          full_test_set(i_full) = test_set(i) ;
          full_test_values(i_full) = test_values(i) ;
          ++i_full ;
        } else {
          full_test_set(i_full) = test_set(i) ;
          full_test_values(i_full) = test_values(i) ;
          ++i_full ;
          full_test_set(i_full) = conj(test_set(i)) ;
          full_test_values(i_full) = conj(test_values(i)) ;
          ++i_full ;
        }
      }
    }
*/
    aaa_approximation<T> repr( repr_upper.n() ) ;
    for (typename decltype(test_set)::size_type i=0; i<repr_upper.n(); ++i) {
      if (imag(repr_upper.nodes()(i))==0.0) {
        repr.add_node( repr_upper.weights()(i), repr_upper.nodes()(i), repr_upper.coefficients()(i) ) ;
      } else {
        repr.add_node( repr_upper.weights()(i), repr_upper.nodes()(i), repr_upper.coefficients()(i) ) ;
        repr.add_node( conj(repr_upper.weights()(i)), conj(repr_upper.nodes()(i)), conj(repr_upper.coefficients()(i)) ) ;
      }
    }
/*
    // Compute Loewner matrix:
    glas2::matrix< value_type > L( full_test_values.size(), repr.n() ) ;
    auto nodes( repr.nodes() ) ;
    auto values( repr.coefficients() ) ;
    for (int i=0; i<repr.n(); ++i) {
      L( glas2::all(), i ) = (full_test_values - values(i)) / (full_test_set-nodes(i)) ;
    }

    // Compute SVD of L

    glas2::vector<real_type>  svd_val(repr.n()) ;
    glas2::matrix<value_type> svd_u(1,1) ;
    int info = boost::numeric::bindings::lapack::gesvd( 'N', 'O', L, svd_val, svd_u, svd_u ) ;
    assert(info==0) ;
    std::cout << svd_val << std::endl;
    repr.weights() = conj(L(repr.n()-1,glas2::all())) ;

    for (int i=1; i<repr.n(); ++i) {
      if (svd_val(i)<options.tolerance*svd_val(0)/sqrt(real_type(L.num_rows()))) {
        repr.n() = i ;
        break ;
      }
    }
    std::cout << L(glas2::range(2,6), glas2::all() ) << std::endl ;
      std::cout << "d = " << repr.n() << std::endl ;

    // Remove small weights
    real_type norm_w = norm_2( repr.weights() ) ;
    for (int i=0; i<repr.weights().size(); ++i) {
      if (norm_2( repr.weights()( glas2::range_from_end(i,0) ) ) < options.tolerance * norm_w ) {
        repr.n() = i ;
        break ;
      } 
    }
*/
    return std::move(repr) ;
  } // symm_AAA()

} } // namespace CORK::approximation

#endif
