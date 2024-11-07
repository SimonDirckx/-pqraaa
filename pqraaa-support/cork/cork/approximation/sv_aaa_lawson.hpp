//  (C) Copyright Karl Meerbergen 2018.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_approximation_sv_aaa_hpp
#define cork_approximation_sv_aaa_hpp

#include <cork/approximation/aaa_options.hpp>
#include <cork/approximation/sv_aaa_approximation.hpp>
#include <cork/basis/barycentric_rational.hpp>
#include <glas2/matrix.hpp>
#include <glas2/vector.hpp>
#include <glas2/bindings.hpp>
#include <boost/numeric/bindings/lapack/driver/gesvd.hpp>
#include <boost/numeric/bindings/lapack/computational/pocon.hpp>
#include <boost/numeric/bindings/lapack/computational/potrf.hpp>
#include <cassert>
#include <iomanip>

namespace CORK { namespace approximation {

  template <typename T, typename FunctionSequence, typename Domain>
  typename std::enable_if< !glas2::is< glas2::DenseVector, FunctionSequence>::value, SV_AAA_approximation<T> >::type SV_AAA( FunctionSequence const& fun_sequence, Domain const& domain ) {
    return SV_AAA<T>( fun_sequence, domain, aaa_options< decltype(std::abs(T())) >() ) ;
  } // SV_AAA

  template <typename Points, typename FunctionValues>
  typename std::enable_if< glas2::is< glas2::DenseVector, Points>::value
                         && glas2::is< glas2::DenseMatrix, FunctionValues>::value
                         , SV_AAA_approximation<typename FunctionValues::value_type>
                         >::type SV_AAA( Points const& points, FunctionValues const& function_values ) {
    return SV_AAA( points, function_values, aaa_options< decltype(std::abs(typename FunctionValues::value_type())) >() ) ;
  } // SV_AAA

  template <typename T, typename FunctionSequence, typename Domain, typename Options>
  typename std::enable_if< !glas2::is< glas2::DenseVector, FunctionSequence>::value, SV_AAA_approximation<T> >::type SV_AAA( FunctionSequence const& fun_sequence, Domain const& domain, Options const& options ) {
    typedef T                                value_type ;

    // Construct the test set
    glas2::vector<value_type> test_set = domain.discretize( std::max(3,options.n_points) ) ;

    // Evaluate and scale the functions
    glas2::matrix<value_type> test_values( test_set.size(), fun_sequence.num_terms() ) ;
    for (typename glas2::vector<value_type>::size_type i=0; i<test_set.size(); ++i) {
      fun_sequence.evaluate( test_set(i), test_values(i,glas2::all()) ) ;
    }
    return SV_AAA( test_set, test_values, options ) ;
  } // SV_AAA()

  template <typename Points, typename FunctionValues, typename Options>
  typename std::enable_if< glas2::is<glas2::DenseMatrix, FunctionValues>::value
                         , SV_AAA_approximation<typename FunctionValues::value_type>
                         >::type SV_AAA( Points const& test_set, FunctionValues test_values, Options const& options, bool scale_functions=true ) {
    typedef typename FunctionValues::value_type value_type ;
    typedef decltype(std::abs(value_type()))    real_type ;
    typedef int                                 size_type ;

    // Check if points are close to each other

    assert( test_set.size()>=3 ) ;
    assert( test_values.num_rows()==test_set.size() ) ;

    glas2::vector<int> ind_C(options.max_degree) ;
    glas2::matrix< value_type > C(test_set.size(), options.max_degree);
    glas2::matrix< value_type > Q(test_set.size()*test_values.num_columns(), options.max_degree);
    glas2::matrix< value_type > S(options.max_degree, options.max_degree); fill(S,0.0) ;
    glas2::matrix< value_type > H(options.max_degree, options.max_degree); fill(H,0.0) ;
    glas2::vector<value_type> D(test_set.size()) ;
    glas2::matrix<value_type> R(test_set.size(), test_values.num_columns() ) ;

    SV_AAA_approximation<value_type> repr( test_values.num_columns(), options.max_degree ) ;
    repr.reset( 0 ) ;

    glas2::vector<real_type> norm_f( test_values.num_columns() ) ;

    if (scale_functions) {
      // Evaluate and scale the functions
      for (typename FunctionValues::size_type j=0; j<test_values.num_columns(); ++j) {
        norm_f(j) = norm_inf( test_values(glas2::all(),j) ) ;
        if (norm_f(j)!=0.0) test_values(glas2::all(),j) /= norm_f(j) ;
        if (norm_1(test_values(glas2::all(),j))<1.1) {
            std::stringstream ss ;
            ss<< "AAA: possible singularity in point for function number " << j ;
            throw std::runtime_error( ss.str() ) ;
        }
      }
    }

    // Compute the SV_AAA approximation
    repr.error() =  abs( vec(test_values)( max_ind( abs( vec(test_values) ) ) ) ) ;
    fill( R, 0.0 ) ;

    for (int d=0; d<options.max_degree && repr.error()>options.tolerance; ++d) {
      // Find index with largest residual
      /*real_type r_max =0.0;
      for (int i=0; i<test_values.num_rows(); ++i) {
        real_type temp_nrm = norm_2( test_values(i,glas2::all()) - R(i,glas2::all()) ) ;
        if (temp_nrm>r_max) {
          r_max = temp_nrm ;
          ind_C(d) = i ;
        }
      }*/
      int ind_Q = max_ind( abs( vec(test_values - R) ) ) ;
      //if (d==0) ind_Q = 79156;
      ind_C(d) = ind_Q % test_set.size() ;
      real_type r_max = std::abs( vec(test_values)(ind_Q) - vec(R)(ind_Q) ) ;
      repr.error() = r_max ;
      if (options.debug_level>0) std::cout << "Degree " << d << " has approximation error " << repr.error() << std::endl ;
      if (r_max<options.tolerance) break ;

      // Add interpolation point
      repr.add_node( 0.0, test_set(ind_C(d)), test_values(ind_C(d),glas2::all()) ) ;
      assert( repr.n() == d+1 ) ;

      // Set column of Cauchy matrix
      C(glas2::all(), d) = glas2::ones<real_type>(test_set.size()) / (test_set-repr.nodes()(d)) ;
      fill( C(ind_C(d), glas2::all()), 0.0 ) ;
      fill( C(ind_C(glas2::range(0,d+1)),d), 0.0 ) ;

      // Column d of Loewner matrix
      auto L = Q(glas2::all(), d) ;
      auto Qi = Q(glas2::all(), glas2::range(0,d)) ;
      for (size_type i=0; i<test_values.num_columns(); ++i) {
        L( glas2::range(i*test_set.size(),(i+1)*test_set.size())) = C(glas2::all(),d) * test_values(glas2::all(),i) - C(glas2::all(),d) * test_values(ind_C(d),i) ;
      }

      bool cheaply_correct = options.fast ;

      // Update H and S to compensate for the removal of the rows
      glas2::vector< value_type > temp( d ) ;
      glas2::range d_r(0,d) ;
      if (d>0) {
        if (cheaply_correct) {
          glas2::matrix<value_type> ee( d, d ) ;
          ee = glas2::identity_matrix<value_type>(d,d) ;
          // Do this for all functions
          for (int i=0; i<test_values.num_columns(); ++i) {
            temp(d_r) = multiply( transpose(S(d_r,d_r)), Qi(i*test_set.size()+ind_C(d),d_r) ) ;
            ee -= outer_prod(conj(temp(d_r)),temp(d_r)) ;
          }
          real_type Si_norm = norm_1(ee) ;
          auto Si( glas2::upper(ee) ) ;
          int info = boost::numeric::bindings::lapack::potrf( Si ) ;
          if (info==0) {
            // We should perform the following operation here, but this may correct H when
            // cheaply_correct is set to false further.
            // H(d_r,d_r) = Si * H(d_r,d_r) ;
            // S(d_r,d_r) = S(d_r,d_r) * inv(Si) ;
            inplace_solve( S(d_r,d_r), Si ) ;

            // Check condition number
            real_type Si_cond ;
            info = boost::numeric::bindings::lapack::pocon( Si, Si_norm, Si_cond ) ;
            assert( info==0 ) ;
            if (Si_cond<options.cheap_correction_tolerance) cheaply_correct = false ;
          } else {
            cheaply_correct = false ;
          }
          if (cheaply_correct) {
            // H(d_r,d_r) = Si * H(d_r,d_r) ;
            inplace_multiply( Si, H(d_r,d_r) ) ;
          }
        }

        // Remove rows
        for (int i=0; i<test_values.num_columns(); ++i) {
          fill( Qi(i*test_set.size()+ind_C(glas2::range(0,d+1)),glas2::all()), 0.0 ) ;
        }

        if (!cheaply_correct) {
          glas2::matrix<value_type> Si( d, d ) ;
          fill(Si,0.0) ;
          for (int i=0; i<d; ++i) {
            auto temp_i = temp(glas2::range(0,i)) ;

            temp_i = multiply( transpose(conj(Q(glas2::all(),glas2::range(0,i)))), Q(glas2::all(),i) ) ;
            Si(glas2::range(0,i),i) = temp_i ;
            Q(glas2::all(),i) -= multiply( Q(glas2::all(),glas2::range(0,i)), temp_i ) ;

            temp_i = multiply( transpose(conj(Q(glas2::all(),glas2::range(0,i)))), Q(glas2::all(),i) ) ;
            Si(glas2::range(0,i),i) += temp_i ;
            Q(glas2::all(),i) -= multiply( Q(glas2::all(),glas2::range(0,i)), temp_i ) ;

            Si(i,i) = norm_2( Q(glas2::all(),i) ) ;
            Q(glas2::all(),i) /= Si(i,i) ;
          }

          inplace_multiply( glas2::upper(Si), H(d_r,d_r) ) ;
          S(d_r,d_r) = glas2::identity_matrix<value_type>(d,d) ;
        }
      }
      S(d,d) = 1.0 ;
    
      //std::cout << "Qj " << L << std::endl ;
      // Orthogonalize L against Qi
      real_type nv = norm_2(L);
      temp = multiply( transpose(conj(Qi)), L ) ;
      H(d_r,d) = multiply( transpose(conj(S(d_r,d_r))), temp ) ;
      temp = multiply( S(d_r,d_r), H(d_r,d) ) ;
      L -= multiply(Qi, temp ) ;
      H(d,d) = norm_2(L) ;
      // Reorthogonalization is necessary for higher precision

      glas2::vector<value_type> h_new( d ) ;

      if (d>0) {
        for (int ot = 0; (ot<2) && (glas2::real(H(d,d)) < 0.707 * nv); ++ot) {
          temp = multiply( transpose(conj(Qi)), L ) ;
          h_new = multiply( transpose(conj(S(d_r,d_r))), temp ) ;
          temp = multiply( S(d_r,d_r), h_new ) ;
          L -= multiply(Qi, temp ) ;
          H(d_r,d) += h_new ;

          nv = glas2::real(H(d,d)) ;
          assert( glas2::imag(H(d,d))==0.0 ) ;
          H(d,d) = norm_2(L) ;
        }
      }
      L /= H(d,d) ;



      glas2::vector<real_type>  svd_val(d+1) ;
      glas2::matrix<value_type> svd_u(1,1) ;
      glas2::matrix<value_type> svd_vh(d+1,d+1) ;
      svd_vh = H(glas2::range(0,d+1),glas2::range(0,d+1)) ;
#ifndef NDEBUG
      int info =
#endif
        boost::numeric::bindings::lapack::gesvd( 'N', 'O', svd_vh, svd_val, svd_u, svd_u ) ;
      //std::cout << "svd : " << svd_val << std::endl ;
      assert(info==0) ;
      repr.weights() = conj(svd_vh(d,glas2::all())) ;
      //std::cout << "nodes " << repr.nodes() << std::endl ;
      //std::cout << "weights " << repr.weights() << std::endl ;
      //std::cout << "coefficients " << repr.coefficients() << std::endl ;

      // Get the rational approximation
      glas2::range d_r1(0,d+1) ;
      D = multiply(C(glas2::all(),d_r1), repr.weights()) ;                     // Denominator
      for (int i=0; i<test_values.num_columns(); ++i) {
        R(glas2::all(),i) = multiply(C(glas2::all(),d_r1), repr.weights()*repr.coefficients()(glas2::all(),i)) ; // Numerator
        R(glas2::all(),i) /= D;
      }
      R(ind_C(d_r1), glas2::all()) = test_values(ind_C(d_r1), glas2::all());
    } // AAA loop

    // Improve with Lawson
    if (options.lawson>0) {
      glas2::vector< real_type > lawson_weights( test.set.size(), test_values.num_columns() ) ;
      fill( lawson_weights, 1.0 ) ;

      glas2::vector<real_type>  svd_val(repr.n()) ;
      glas2::matrix<value_type> svd_u(1,1) ;
      auto svd_vt = L( glas2::all(), glas2::range(0,repr.n()) ) ;

      for (int ii=0; ii<L.num_rows(); ++ii) {
        svd_vt(ii, glas2::all() ) *= lawson_weights( ii ) ;
      }

      for (int i=0; options.lawson; ++i) {
        boost::numeric::bindings::lapack::gesvd( 'N', 'O', L, svd_val, svd_u, svd_u ) ;

        D = multiply(C(glas2::all(),d_r1), repr.weights()) ;                     // Denominator
        for (int j=0; j<test_values.num_columns(); ++j) {
          R(glas2::all(),j) = multiply(C(glas2::all(),d_r1), repr.weights()*repr.coefficients()(glas2::all(),j)) ; // Numerator
          R(glas2::all(),j) /= D;
        }
        R -= test_values ;

        for (int j=0; i<test_set.size(); ++j) {
        lawson_weights = lawson_weights * abs(R) ;
      }
      repr.weights() = conj(svd_vh(d,glas2::all())) ;
    } // Lawson loop

    if (scale_functions) {
      for (int i=0; i<test_values.num_columns(); ++i) {
        repr.coefficients()(glas2::all(),i) *= norm_f(i) ;
      }
    }

    // Sort following decreasing weights
    {
      glas2::vector< value_type > 
    }

    // Remove small weights
    real_type norm_w = norm_2( repr.weights() ) ;
    for (int i=0; i<repr.weights().size(); ++i) {
      if (norm_2( repr.weights()( glas2::range_from_end(i,0) ) ) < options.tolerance * norm_w ) {
        repr.n() = i ;
        break ;
      } 
    }

    return repr ;
  } // SV_AAA()

} } // namespace CORK::approximation

#endif
