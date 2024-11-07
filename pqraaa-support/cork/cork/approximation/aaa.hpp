//  (C) Copyright Karl Meerbergen 2018.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_approximation_aaa_hpp
#define cork_approximation_aaa_hpp

#include <cork/approximation/aaa_options.hpp>
#include <cork/approximation/aaa_approximation.hpp>
#include <cork/options/max_degree.hpp>
#include <cork/options/number_of_sample_points.hpp>
#include <cork/options/aaa_stop_criterion.hpp>
#include <cork/options/debug_level.hpp>
#include <cork/options/value_of.hpp>
#include <glas2/matrix.hpp>
#include <glas2/vector.hpp>
#include <glas2/bindings.hpp>
#include <boost/numeric/bindings/lapack/driver/gesvd.hpp>
#include <boost/numeric/bindings/lapack/computational/potrf.hpp>
#include <cassert>

namespace CORK { namespace approximation {

  template <typename F, typename Domain, typename V>
  void aaa_test_set( F const& fun, Domain const& domain, V& test_set, V& test_values ) {
    auto ts = domain.discretize(test_set.size()) ;
    test_set = ts ;
    for (typename V::size_type i=0; i<test_set.size(); ++i) {
      test_values(i) = fun( test_set(i) ) ;
    }
  } // aaa_test_set()

  template <typename T, typename V, typename I, typename Options>
  aaa_approximation<T> aaa_approximate( V& test_set, V& test_values, I& ind, Options const& options ) {
    typedef T                                value_type ;
    typedef decltype(std::abs(value_type())) real_type ;
    assert( ind.size() == options::value_of<options::max_degree>(options) ) ;
    
    glas2::matrix< value_type > C(test_set.size(), options::value_of<options::max_degree>(options));
    glas2::matrix< value_type > Q(test_set.size(), options::value_of<options::max_degree>(options));
    glas2::matrix< value_type > S(options::value_of<options::max_degree>(options), options::value_of<options::max_degree>(options)); fill(S,0.0) ;
    glas2::matrix< value_type > H(options::value_of<options::max_degree>(options), options::value_of<options::max_degree>(options)); fill(H,0.0) ;
    glas2::vector<value_type> N(test_set.size()) ;
    glas2::vector<value_type> D(test_set.size()) ;
    glas2::vector<value_type> R(test_set.size()) ;
    glas2::vector<real_type>  svd_val(options::value_of<options::max_degree>(options)) ;
    
    aaa_approximation<T> repr( options::value_of<options::max_degree>(options) ) ;
    repr.reset( 0 ) ;

    // Compute the AAA approximation
    real_type norm_f = norm_inf( test_values ) ;
    fill( R, sum( test_values ) / real_type(test_values.size()) ) ;
    repr.error() = norm_inf(R) / norm_f ;
    
    if (options::value_of<options::debug_level>(options)>0) std::cout << "Degree 0 has approximation error " << repr.error() << std::endl ;

    for (int d=0; d<options::value_of<options::max_degree>(options) && repr.error()>options::value_of<options::aaa_stop_criterion<real_type>>(options).tolerance(); ++d) {

      ind(d) = max_ind( abs( test_values - R ) ) ;
      //std::cout << "Max for index "<< ind(d) << " at point " << test_set(ind(d)) << " with value " << test_values(ind(d)) << " and residual "
      //          << test_values(ind(d)) - R(ind(d)) << std::endl ;
      real_type r_max = std::abs( test_values(ind(d)) - R(ind(d)) ) ;
      repr.error() = r_max/norm_f ;
      if (options::value_of<options::debug_level>(options)>0) std::cout << "Degree " << d << " has approximation error " << repr.error() << std::endl ;
      if (r_max<options::value_of<options::aaa_stop_criterion<real_type>>(options).tolerance()*norm_f) break ;

      // Add interpolation point
      repr.add_node( 0.0, test_set(ind(d)), test_values(ind(d)) ) ;

      // Set column of Cauchy matrix
      C(glas2::all(), d) = glas2::ones<real_type>(test_set.size()) / (test_set-repr.nodes()(d)) ;
      fill( C(ind(glas2::range(0,d+1)),d), 0.0 ) ;

      // Column d of Loewner matrix
      auto L = Q(glas2::all(), d) ;
      auto Qi = Q(glas2::all(), glas2::range(0,d)) ;
      L = C(glas2::all(),d) * test_values - C(glas2::all(),d)*test_values(ind(d)) ;

      // Update H and S to compensate for the removal of the rows
      glas2::vector< value_type > temp( d ) ;
      glas2::range d_r(0,d) ;
      if (d>0) {
        temp(d_r) = multiply( transpose(S(d_r,d_r)), Qi(ind(d),d_r) ) ;
        glas2::matrix<value_type> ee( d, d ) ;
        ee = glas2::identity_matrix<value_type>(d,d) - outer_prod(conj(temp(d_r)),temp(d_r)) ;
        auto Si( glas2::upper(ee) ) ;
        int info = boost::numeric::bindings::lapack::potrf( Si ) ;
        assert(info==0) ;
        // H(d_r,d_r) = Si * H(d_r,d_r) ;
        inplace_multiply( Si, H(d_r,d_r) ) ;
        // S(d_r,d_r) = S(d_r,d_r) * inv(Si) ;
        inplace_solve( S(d_r,d_r), Si ) ;
      }
      S(d,d) = 1.0 ;
      fill( Qi(ind(glas2::range(0,d+1)),glas2::all()), 0.0 ) ;
    
      real_type nv = norm_2(L);
      temp = multiply( transpose(conj(Qi)), L ) ;
      H(d_r,d) = multiply( transpose(conj(S(d_r,d_r))), temp ) ;
      temp = multiply( S(d_r,d_r), H(d_r,d) ) ;
      L -= multiply(Qi, temp ) ;
      H(d,d) = norm_2(L) ;
      // Reorthogonalization is necessary for higher precision

      glas2::vector<value_type> h_new( d ) ;

      if (d>0) {
      for (int ot = 0; (ot<3) && (std::abs(H(d,d)) < 0.707 * nv); ++ot) {
          temp = multiply( transpose(conj(Qi)), L ) ;
          h_new = multiply( transpose(conj(S(d_r,d_r))), temp ) ;
          temp = multiply( S(d_r,d_r), h_new ) ;
          L -= multiply(Qi, temp ) ;
          H(d_r,d) += h_new ;

          nv = std::abs(H(d,d)) ;
          H(d,d) = norm_2(L) ;
        }
      }
      L /= H(d,d) ;

      glas2::matrix<value_type> svd_u(1,1) ;
      glas2::matrix<value_type> svd_vh(d+1,d+1) ;
      svd_vh = H(glas2::range(0,d+1),glas2::range(0,d+1)) ;
      auto svd_val_d = svd_val(glas2::range(0,d+1)) ;
      int info = boost::numeric::bindings::lapack::gesvd( 'N', 'O', svd_vh, svd_val_d, svd_u, svd_u ) ;
      assert(info==0) ;
      // truncate small singular values.
      repr.weights() = conj(svd_vh(d,glas2::all())) ;

      // Get the rational approximation
      glas2::range d_r1(0,d+1) ;
      N = multiply(C(glas2::all(),d_r1), repr.weights()*repr.coefficients()) ; // Numerator
      D = multiply(C(glas2::all(),d_r1), repr.weights()) ;                     // Denominator
      R = N/D;
      R(ind(d_r1)) = test_values(ind(d_r1));
//      std::cout << "Residual norm " << norm_2(test_values-R) << std::endl ;
    }
/*
    for (int i=1; i<repr.n(); ++i) {
      if (svd_val(i)<options.tolerance*svd_val(0)/sqrt(real_type(test_values.size()))) {
        repr.n() = i ;
        break ;
      }
    }
      std::cout << "d = " << repr.n() << std::endl ;*/

    // Remove small weights
    real_type norm_w = norm_2( repr.weights() ) ;
    for (int i=0; i<repr.weights().size(); ++i) {
      if (norm_2( repr.weights()( glas2::range_from_end(i,0) ) ) < options::value_of<options::aaa_stop_criterion<real_type>>(options).tolerance() * norm_w ) {
        repr.n() = i ;
        break ;
      } 
    }

    return std::move(repr);
  } // aaa_approximate()

  template <typename T, typename F, typename Domain>
  aaa_approximation<T> AAA( F const& fun, Domain const& domain ) {
    return AAA<T>( fun, domain, std::tuple<>() ) ;
  }


  template <typename T, typename F, typename Domain, typename Options>
  aaa_approximation<T> AAA( F const& fun, Domain const& domain, Options const& options ) {
    typedef T                                value_type ;

    // Construct the test set
    glas2::vector<value_type> test_set( std::max(3,options::value_of<options::number_of_sample_points>(options)) ) ;
    glas2::vector<value_type> test_values( test_set.size() ) ;
    aaa_test_set( fun, domain, test_set, test_values ) ;

    // Compute AAA
    glas2::vector<int> ind(options::value_of<options::max_degree>(options)) ;
    return aaa_approximate<T>( test_set, test_values, ind, options ) ;
/*
    // Construct the test set
    glas2::vector<value_type> test_set( std::max(3,options.n_points) ) ;
    domain.discretize( test_set ) ;

    glas2::vector<int> ind(options.max_degree) ;
    glas2::matrix< value_type > C(test_set.size(), options.max_degree);
    glas2::matrix< value_type > Q(test_set.size(), options.max_degree);
    glas2::matrix< value_type > S(options.max_degree, options.max_degree); fill(S,0.0) ;
    glas2::matrix< value_type > H(options.max_degree, options.max_degree); fill(H,0.0) ;
    glas2::vector<value_type> N(test_set.size()) ;
    glas2::vector<value_type> D(test_set.size()) ;
    glas2::vector<value_type> R(test_set.size()) ;

    aaa_approximation<T> repr( options.max_degree ) ;
    repr.reset( 0 ) ;

    // Evaluate the functions
    glas2::vector<value_type> test_values( test_set.size() ) ;
    for (typename glas2::vector<value_type>::size_type i=0; i<test_set.size(); ++i) {
      test_values(i) = fun( test_set(i) ) ;
    }

    // Compute the AAA approximation
    real_type norm_f = norm_inf( test_values ) ;
    fill( R, sum( test_values ) / real_type(test_values.size()) ) ;
    repr.error() = norm_inf(R) / norm_f ;

    for (int d=0; d<options.max_degree && repr.error()>options.tolerance; ++d) {
      ind(d) = max_ind( abs( test_values - R ) ) ;
      real_type r_max = std::abs( test_values(ind(d)) - R(ind(d)) ) ;
      if (r_max<options.tolerance*norm_f) break ;

      // Add interpolation point
      repr.add_node( 0.0, test_set(ind(d)), test_values(ind(d)) ) ;

      // Set column of Cauchy matrix
      C(glas2::all(), d) = glas2::ones<real_type>(test_set.size()) / (test_set-repr.nodes()(d)) ;
      fill( C(ind(glas2::range(0,d+1)),d), 0.0 ) ;

      // Column d of Loewner matrix
      auto L = Q(glas2::all(), d) ;
      auto Qi = Q(glas2::all(), glas2::range(0,d)) ;
      L = C(glas2::all(),d) * test_values - C(glas2::all(),d)*test_values(ind(d)) ;

      // Update H and S to compensate for the removal of the rows
      glas2::vector< value_type > temp( d ) ;
      glas2::range d_r(0,d) ;
      if (d>0) {
        temp(d_r) = multiply( transpose(S(d_r,d_r)), Qi(ind(d),d_r) ) ;
        glas2::matrix<value_type> ee( d, d ) ;
        ee = glas2::identity_matrix<value_type>(d,d) - outer_prod(conj(temp(d_r)),temp(d_r)) ;
        auto Si( glas2::upper(ee) ) ;
        int info = boost::numeric::bindings::lapack::potrf( Si ) ;
        assert(info==0) ;
        // H(d_r,d_r) = Si * H(d_r,d_r) ;
        inplace_multiply( Si, H(d_r,d_r) ) ;
        // S(d_r,d_r) = S(d_r,d_r) * inv(Si) ;
        inplace_solve( S(d_r,d_r), Si ) ;
      }
      S(d,d) = 1.0 ;
      fill( Qi(ind(glas2::range(0,d+1)),glas2::all()), 0.0 ) ;
    
      real_type nv = norm_2(L);
      temp = multiply( transpose(conj(Qi)), L ) ;
      H(d_r,d) = multiply( transpose(conj(S(d_r,d_r))), temp ) ;
      temp = multiply( S(d_r,d_r), H(d_r,d) ) ;
      L -= multiply(Qi, temp ) ;
      H(d,d) = norm_2(L) ;
      // Reorthogonalization is necessary for higher precision

      glas2::vector<value_type> h_new( d ) ;

      if (d>0) {
      for (int ot = 0; (ot<3) && (std::abs(H(d,d)) < 0.707 * nv); ++ot) {
          temp = multiply( transpose(conj(Qi)), L ) ;
          h_new = multiply( transpose(conj(S(d_r,d_r))), temp ) ;
          L -= multiply(Qi, h_new) ;
          H(d_r,d) += h_new ;

          nv = std::abs(H(d,d)) ;
          H(d,d) = norm_2(L) ;
        }
      }
      L /= H(d,d) ;

      glas2::vector<real_type>  svd_val(d+1) ;
      glas2::matrix<value_type> svd_u(1,1) ;
      glas2::matrix<value_type> svd_vh(d+1,d+1) ;
      svd_vh = H(glas2::range(0,d+1),glas2::range(0,d+1)) ;
      int info = boost::numeric::bindings::lapack::gesvd( 'N', 'O', svd_vh, svd_val, svd_u, svd_u ) ;
      assert(info==0) ;
      repr.weights() = conj(svd_vh(d,glas2::all())) ;

      // Get the rational approximation
      glas2::range d_r1(0,d+1) ;
      N = multiply(C(glas2::all(),d_r1), repr.weights()*repr.coefficients()) ; // Numerator
      D = multiply(C(glas2::all(),d_r1), repr.weights()) ;                     // Denominator
      R = N/D;
      R(ind(d_r1)) = test_values(ind(d_r1));
//      std::cout << "Residual norm " << norm_2(test_values-R) << std::endl ;
    }

    // Remove small weights
    real_type norm_w = norm_2( repr.weights() ) ;
    for (int i=0; i<repr.weights().size(); ++i) {
      if (norm_2( repr.weights()( glas2::range_from_end(i,0) ) ) < options.tolerance * norm_w ) {
        repr.n() = i ;
        break ;
      } 
    }

    return std::move(repr) ;*/
  } // AAA()

} } // namespace CORK::approximation

#endif
