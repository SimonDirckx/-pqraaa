//  (C) Copyright Karl Meerbergen 2018.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_approximation_symm_aaa_hpp
#define cork_approximation_symm_aaa_hpp

#include <cork/approximation/aaa_approximation.hpp>
#include <cork/options/max_degree.hpp>
#include <cork/options/debug_level.hpp>
#include <cork/options/number_of_sample_points.hpp>
#include <cork/options/value_of.hpp>
#include <cork/options/aaa_stop_criterion.hpp>
#include <cork/options/aaa_cheap_correction_tolerance.hpp>
#include <glas2/matrix.hpp>
#include <glas2/vector.hpp>
#include <glas2/bindings.hpp>
#include <boost/numeric/bindings/lapack/driver/gesvd.hpp>
#include <boost/numeric/bindings/lapack/computational/potrf.hpp>
#include <cassert>

namespace CORK { namespace approximation {

  template <typename T, typename F, typename Domain, typename Options>
  aaa_approximation<T> symm_AAA( F const& fun, Domain const& domain, Options const& options ) {
    typedef T                                value_type ;
    typedef decltype(std::abs(value_type())) real_type ;

    // Construct the test set
    int half_set_size = std::max(4,options::value_of<options::number_of_sample_points>(options)+1)/ 2 ;
    glas2::vector<value_type> test_set( half_set_size*2 ) ;
    auto half = domain.discretize( half_set_size ) ;
    test_set(glas2::range(0,half_set_size)) = half ;
    test_set(glas2::range_from_end(half_set_size,0)) = conj(half) ;

    int max_degree = options::value_of<options::max_degree>(options) ;
    glas2::vector<int> ind(max_degree) ;
    glas2::matrix< value_type > C(test_set.size(), max_degree);
    glas2::matrix< value_type > Q(test_set.size(), max_degree);
    glas2::matrix< value_type > S(max_degree, max_degree); fill(S,0.0) ;
    glas2::matrix< value_type > H(max_degree, max_degree); fill(H,0.0) ;
    glas2::vector<value_type> N(test_set.size()) ;
    glas2::vector<value_type> D(test_set.size()) ;
    glas2::vector<value_type> R(test_set.size()) ;

    real_type tolerance = options::value_of<options::aaa_stop_criterion<real_type>>(options).tolerance() ;

    aaa_approximation<T> repr( max_degree ) ;
    repr.reset( 0 ) ;

    // Evaluate the functions
    glas2::vector<value_type> test_values( test_set.size() ) ;
    for (typename glas2::vector<value_type>::size_type i=0; i<half_set_size; ++i) {
      test_values(i) = fun( test_set(i) ) ;
      test_values(half_set_size+i) = conj( test_values(i) ) ;
    }

    // Compute the AAA approximation
    real_type norm_f = norm_inf( test_values ) ;
    fill( R, sum( test_values ) / real_type(test_values.size()) ) ;
    repr.error() = norm_inf(R) / norm_f ;

    for (int d=0; d<max_degree && repr.error()>tolerance; /*++d*/) {
      ind(d) = max_ind( abs( test_values(glas2::range(0,half_set_size)) - R(glas2::range(0,half_set_size)) ) ) ;
      real_type r_max = std::abs( test_values(ind(d)) - R(ind(d)) ) ;
      if (r_max<tolerance*norm_f) break ;

      int todo = 0;
      if (imag(test_set(ind(d)))==0) {
        // Add one interpolation point
        assert( imag(test_values(ind(d)))==0.0 ) ;
        repr.add_node( 0.0, test_set(ind(d)), test_values(ind(d)) ) ;

        // Set column of Cauchy matrix
        C(glas2::all(), d) = glas2::ones<real_type>(test_set.size()) / (test_set-repr.nodes()(d)) ;
        fill( C(ind(glas2::range(0,d+1)),d), 0.0 ) ;
        todo = 1 ;
      } else {
        // Add two complex conjugate interpolation points
        repr.add_node( 0.0, test_set(ind(d)), test_values(ind(d)) ) ;
        ind(d+1) = ind(d) + half_set_size ;
        repr.add_node( 0.0, test_set(ind(d+1)), test_values(ind(d+1)) )  ;
        assert( test_set(ind(d+1))==conj(test_set(ind(d))) ) ;
        assert( test_values(ind(d+1))==conj(test_values(ind(d))) ) ;

        // Set column of Cauchy matrix
        C(glas2::all(), d) = glas2::ones<real_type>(test_set.size()) / (test_set-repr.nodes()(d)) ;
        C(glas2::all(), d+1) = glas2::ones<real_type>(test_set.size()) / (test_set-repr.nodes()(d+1)) ;
        fill( C(ind(glas2::range(0,d+2)),glas2::range(d,d+2)), 0.0 ) ;
        todo = 2 ;
      }
      if (d+todo>max_degree) {
        repr.n() = d ;
        break ;
      }

      // Update H and S to compensate for the removal of rows
      {
      auto Qd = Q(glas2::all(), glas2::range(0,d)) ;
      glas2::matrix< value_type > temp( todo, d ) ;
      glas2::range d_r(0,d) ;
      if (d>0) {
        glas2::matrix<value_type> ee( d, d ) ;
        temp = multiply( Qd(ind(glas2::range(d,d+todo)), d_r), S(d_r,d_r) ) ;
        ee = glas2::identity_matrix<value_type>(d,d) ;
        ee -= multiply( conj(transpose(temp)), temp ) ;
        auto Si( glas2::upper(ee) ) ;
        int info = boost::numeric::bindings::lapack::potrf( Si ) ;
        assert(info==0) ;
        // H(d_r,d_r) = Si * H(d_r,d_r) ;
        inplace_multiply( Si, H(d_r,d_r) ) ;
        // S(d_r,d_r) = S(d_r,d_r) * inv(Si) ;
        inplace_solve( S(d_r,d_r), Si ) ;
      }
      S(d,d) = 1.0 ;
      fill( Qd(ind(glas2::range(0,d+todo)),glas2::all()), 0.0 ) ;
      }

      for (int i_todo=0; i_todo<todo; ++i_todo, ++d) {
        // Column d_i_todo of Loewner matrix
        glas2::range d_r(0,d) ;
        auto Qi = Q(glas2::all(), glas2::range(0,d)) ;
        auto L = Q(glas2::all(), d) ;
        L = C(glas2::all(),d) * test_values - C(glas2::all(),d)*test_values(ind(d)) ;

        glas2::vector< value_type > temp( d ) ;
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
      }
      // d has been increased here.

      glas2::vector<real_type>  svd_val(d) ;
      glas2::matrix<value_type> svd_u(1,1) ;
      glas2::matrix<value_type> svd_vh(d,d) ;
      svd_vh = H(glas2::range(0,d),glas2::range(0,d)) ;
      int info = boost::numeric::bindings::lapack::gesvd( 'N', 'O', svd_vh, svd_val, svd_u, svd_u ) ;
      std::cout << svd_val << std::endl ;
      assert(info==0) ;
      repr.weights() = conj(svd_vh(d-1,glas2::all())) ;

      // Get the rational approximation
      glas2::range d_r1(0,d) ;
      N = multiply(C(glas2::all(),d_r1), repr.weights()*repr.coefficients()) ; // Numerator
      D = multiply(C(glas2::all(),d_r1), repr.weights()) ;                     // Denominator
      R = N/D;
      R(ind(d_r1)) = test_values(ind(d_r1));
//      std::cout << "Residual norm " << norm_2(test_values-R) << std::endl ;
    }

    // Remove small weights
    real_type norm_w = norm_2( repr.weights() ) ;
    for (int i=0; i<repr.weights().size(); ++i) {
      if (norm_2( repr.weights()( glas2::range_from_end(i,0) ) ) < tolerance * norm_w ) {
        repr.n() = i ;
        break ;
      }
      if (i<repr.weights().size()-1 && repr.nodes()(i)==conj(repr.nodes()(i+1))) ++i ;
    }

    return std::move(repr) ;
  } // AAA()

} } // namespace CORK::approximation

#endif
