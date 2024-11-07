//  (C) Copyright Karl Meerbergen 2018.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_approximation_qr_aaa_hpp
#define cork_approximation_qr_aaa_hpp

#include <cork/approximation/aaa_options.hpp>
#include <cork/approximation/sv_aaa_approximation.hpp>
#include <cork/basis/barycentric_rational.hpp>
#include <cork/lapack/qr.hpp>
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
  typename std::enable_if< !glas2::is< glas2::DenseVector, FunctionSequence>::value, SV_AAA_approximation<T> >::type QR_AAA( FunctionSequence const& fun_sequence, Domain const& domain ) {
    return QR_AAA<T>( fun_sequence, domain, aaa_options< decltype(std::abs(T())) >() ) ;
  } // QR_AAA

  template <typename Points, typename FunctionValues>
  typename std::enable_if< glas2::is< glas2::DenseVector, Points>::value
                         && glas2::is< glas2::DenseMatrix, FunctionValues>::value
                         , SV_AAA_approximation<typename FunctionValues::value_type>
                         >::type QR_AAA( Points const& points, FunctionValues const& function_values ) {
    return QR_AAA( points, function_values, aaa_options< decltype(std::abs(typename FunctionValues::value_type())) >() ) ;
  } // QR_AAA

  template <typename T, typename FunctionSequence, typename Domain, typename Options>
  typename std::enable_if< !glas2::is< glas2::DenseVector, FunctionSequence>::value, SV_AAA_approximation<T> >::type QR_AAA( FunctionSequence const& fun_sequence, Domain const& domain, Options const& options ) {
    typedef T                                value_type ;
    typedef decltype(std::abs(value_type())) real_type ;
    typedef int                              size_type ;

    // Construct the test set
    glas2::vector<value_type> test_set = domain.discretize( std::max(3,options.n_points) ) ;

    // Evaluate and scale the functions
    glas2::matrix<value_type> test_values( test_set.size(), fun_sequence.num_terms() ) ;
    for (typename glas2::vector<value_type>::size_type i=0; i<test_set.size(); ++i) {
      fun_sequence.evaluate( test_set(i), test_values(i,glas2::all()) ) ;
    }
    return QR_AAA( test_set, test_values, options ) ;
  } // QR_AAA()

  template <typename Points, typename FunctionValues, typename Options>
  typename std::enable_if< glas2::is<glas2::DenseMatrix, FunctionValues>::value
                         , SV_AAA_approximation<typename FunctionValues::value_type>
                         >::type QR_AAA( Points const& test_set, FunctionValues test_values, Options const& options ) {
    typedef typename FunctionValues::value_type value_type ;
    typedef decltype(std::abs(value_type()))    real_type ;
    typedef int                                 size_type ;

    assert( test_set.size()>=3 ) ;
    assert( test_values.num_rows()==test_set.size() ) ;

    glas2::vector<int> ind_Q(options.max_degree) ;
    glas2::vector<int> ind_C(options.max_degree) ;
    glas2::matrix< value_type > C(test_set.size(), options.max_degree);
    glas2::matrix< value_type > S(options.max_degree, options.max_degree); fill(S,0.0) ;
    glas2::matrix< value_type > H(options.max_degree, options.max_degree); fill(H,0.0) ;
    glas2::vector<value_type> D(test_set.size()) ;
    glas2::matrix<value_type> function_scaling(test_values.num_columns(), test_values.num_columns() ) ;

    SV_AAA_approximation<value_type> repr( test_values.num_columns(), options.max_degree ) ;
    repr.reset( 0 ) ;

    glas2::vector<real_type> norm_f( test_values.num_columns() ) ;

    int n_functions = test_values.num_columns();
    // Evaluate the functions
    {
      glas2::vector<real_type>  svd_val(test_values.num_columns()) ;
      glas2::matrix<value_type> temp(1,1) ;
#ifndef NDEBUG
      int info =
#endif
      boost::numeric::bindings::lapack::gesvd( 'O', 'A', test_values, svd_val, function_scaling, function_scaling ) ;
      for ( ; svd_val(n_functions-1)<options.function_drop_tol*svd_val(0); --n_functions ) ;
      for (int i=0; i<n_functions; ++i) function_scaling(i,glas2::all()) *= svd_val(i) ;
    }
    if (options.debug_level>0) std::cout << "Number of functions is reduced from " << test_values.num_columns() << " to " << n_functions << std::endl ;
    glas2::matrix<value_type> R(test_set.size(), n_functions ) ;
    glas2::matrix< value_type > Q(test_set.size()*n_functions, options.max_degree);
    auto test_values_r( test_values(glas2::all(), glas2::range(0,n_functions)) ) ;

    // Compute the SV_AAA approximation
    repr.error() =  abs( vec(test_values_r)( max_ind( abs( vec(test_values_r) ) ) ) ) ;
    fill( R, 0.0 ) ;
    if (options.debug_level>0) std::cout << "Degree 0 has approximation error " << repr.error() << std::endl ;

    real_type last_error = 1.0 ;

    for (int d=0; d<options.max_degree && repr.error()>options.tolerance; ++d) {
      // Find index with largest residual
      ind_Q(d) = max_ind( abs( vec(test_values_r - R) ) ) ;
      //if (d==0) ind_Q(d) = 79156;
      ind_C(d) = ind_Q(d) % test_set.size() ;
      real_type r_max = std::abs( vec(test_values_r)(ind_Q(d)) - vec(R)(ind_Q(d)) ) ;
      repr.error() = r_max ;
      if (options.debug_level>0) std::cout << "Degree " << d << " has approximation error " << repr.error() << std::endl ;
      if (r_max<options.tolerance) break ;
      //if (r_max<options.tolerance || r_max>1.1*last_error) break ;
      last_error = r_max ;

      // Add interpolation point
      repr.add_node( 0.0, test_set(ind_C(d)), test_values(ind_C(d),glas2::all()) ) ; // Add test_values because of dimensions of repr.coefficients
      assert( repr.n() == d+1 ) ;

      // Set column of Cauchy matrix
      C(glas2::all(), d) = glas2::ones<real_type>(test_set.size()) / (test_set-repr.nodes()(d)) ;
      fill( C(ind_C(glas2::range(0,d+1)),d), 0.0 ) ;

      // Column d of Loewner matrix
      auto L = Q(glas2::all(), d) ;
      auto Qi = Q(glas2::all(), glas2::range(0,d)) ;
      for (size_type i=0; i<test_values_r.num_columns(); ++i) {
        L( glas2::range(i*test_set.size(),(i+1)*test_set.size())) = C(glas2::all(),d) * test_values_r(glas2::all(),i) - C(glas2::all(),d) * test_values_r(ind_C(d),i) ;
      }

      bool cheaply_correct = true ;

      // Update H and S to compensate for the removal of the rows
      glas2::vector< value_type > temp( d ) ;
      glas2::range d_r(0,d) ;
      if (d>0) {
        if (cheaply_correct) {
          glas2::matrix<value_type> ee( d, d ) ;
          ee = glas2::identity_matrix<value_type>(d,d) ;
          // Do this for all functions
          for (int i=0; i<test_values_r.num_columns(); ++i) {
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
        for (int i=0; i<test_values_r.num_columns(); ++i) {
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
      //std::cout << "H " << svd_vh << std::endl ;
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
      for (int i=0; i<test_values_r.num_columns(); ++i) {
        R(glas2::all(),i) = multiply(C(glas2::all(),d_r1), repr.weights()*repr.coefficients()(glas2::all(),i)) ; // Numerator
        R(glas2::all(),i) /= D;
      }
      R(ind_C(d_r1), glas2::all()) = test_values_r(ind_C(d_r1), glas2::all());
    } // End of main loop

    glas2::vector< value_type > temp_coefs( n_functions ) ;
    for (int i=0; i<repr.n(); ++i) {
      temp_coefs = repr.coefficients()(i, glas2::range(0,n_functions)) ;
      repr.coefficients()(i,glas2::all()) = multiply( transpose(function_scaling(glas2::range(0,n_functions),glas2::all())), temp_coefs ) ;
    }

    // Remove small weights
    real_type norm_w = norm_2( repr.weights() ) ;
    for (int i=0; i<repr.weights().size(); ++i) {
      if (norm_2( repr.weights()( glas2::range_from_end(i,0) ) ) < options.tolerance * norm_w ) {
        repr.n() = i ;
        break ;
      } 
    }

    // Scale weights to make norm one of basis functions
/*    basis::barycentric_rational< decltype(repr.weights()), decltype(repr.nodes()) > basis( repr.weights(), repr.nodes() ) ;
    glas2::vector<value_type> eval_basis( basis.num_terms() ) ;
    real_type max_val = 0.0 ;
    for (int i=0; i<test_set.size(); ++i) {
      basis.evaluate( test_set(i), eval_basis ) ;
      max_val = std::max( max_val, norm_inf(eval_basis(glas2::range_from_end(1,0))) ) ;
    }
    if (options.debug_level>0) std::cout << "AAA: scale basis functions with factor " << max_val << std::endl ;
    repr.weights() /= max_val ;
    repr.coefficients() *= max_val ;
*/
    return repr ;
  } // QR_AAA()

} } // namespace CORK::approximation

#endif
