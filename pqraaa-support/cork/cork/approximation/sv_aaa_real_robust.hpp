//  (C) Copyright Karl Meerbergen 2018.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_approximation_sv_aaa_real_robust_hpp
#define cork_approximation_sv_aaa_real_robust_hpp

#include <cork/approximation/aaa_options.hpp>
#include <cork/approximation/sv_aaa_approximation_real.hpp>
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
  decltype (auto) SV_AAA_real_robust( FunctionSequence const& fun_sequence, Domain const& domain ) {
    return SV_AAA_real_robust<T>( fun_sequence, domain, aaa_options< decltype(std::abs(T())) >() ) ;
  } // SV_AAA

  template <typename T, typename FunctionSequence, typename Domain, typename Options>
  decltype(auto) SV_AAA_real_robust( FunctionSequence const& fun_sequence, Domain const& domain, Options const& options ) {
    typedef T                                value_type ;
    typedef decltype(std::abs(value_type())) real_value_type ;

    // Construct the test set
    glas2::vector<value_type> test_set1 = domain.discretize( std::max(3,options.n_points) ) ;

    //test_set1(glas2::range_from_end(1,0)) = glas2::linspace(-4.0, 6.0, test_set1.size()-1) ;
    //test_set1(0) = value_type(0.0,1.0) ;

    //
    // Sort the test_set following real and complex points. Complex points come first.
    //

    // Count the number of complex points
    int n_c = 0.0 ;
    for (int i=0;i<test_set1.size();++i) {
      n_c += test_set1(i).imag()!=0.0 ;
    }

    // Reorder the points.
    int n_points = n_c*2 + (test_set1.size()-n_c) ;
    glas2::vector<value_type> test_set( test_set1.size() ) ;
    int ii = 0 ;
    for (int i=0;i<test_set1.size();++i) {
      if (test_set1(i).imag()!=0.0) {
        test_set(ii) = test_set1(i) ;
        ++ii ;
      }
    }
    for (int i=0;i<test_set1.size();++i) {
      if (test_set1(i).imag()==0.0) {
        test_set(ii) = test_set1(i) ;
        ++ii ;
      }
    }
    assert( ii==test_set.size() ) ;

    int complex_offset = test_set1.size() ;
    //std::cout << "test_set " << test_set << std::endl ;
    assert( complex_offset+n_c==n_points ) ;

    glas2::vector<int> ind_C(options.max_degree+1) ; // We add one position as sentinel for complex supports.
    glas2::matrix< real_value_type > Q(n_points*fun_sequence.num_terms(), options.max_degree+1); fill(Q,0.0) ; // We add one column as sentinel for complex supports.
    glas2::matrix< real_value_type > S(options.max_degree, options.max_degree); S = glas2::eye( S.num_rows(), S.num_columns() ) ;
    glas2::matrix< real_value_type > H(options.max_degree+1, options.max_degree+1); fill(H,0.0) ; // Sentinel

    SV_AAA_approximation_real<value_type> repr( fun_sequence.num_terms(), options.max_degree ) ;
    repr.reset( 0 ) ;

    glas2::vector<real_value_type> norm_f( fun_sequence.num_terms() ) ;

    // Evaluate and scale the functions
    glas2::matrix<value_type> test_values( test_set.size(), fun_sequence.num_terms() ) ;
    for (typename glas2::vector<value_type>::size_type i=0; i<test_set.size(); ++i) {
      fun_sequence.evaluate( test_set(i), test_values(i,glas2::all()) ) ;
    }
    for (typename FunctionSequence::size_type j=0; j<fun_sequence.num_terms(); ++j) {
      norm_f(j) = norm_inf( test_values(glas2::all(),j) ) ;
      if (norm_f(j)!=0.0) test_values(glas2::all(),j) /= norm_f(j) ;
    }
    //std::cout << "test_values " << test_values << std::endl ;

    // The approximation is now zero. So, the error is the function itself
    repr.error() = norm_inf( vec(test_values) ) ;
    if (options.debug_level>0) std::cout << "Degree " << 0 << " has approximation error " << repr.error() << std::endl ;

    assert( options.max_degree>0 ) ;

    for (int d=0; d<options.max_degree && repr.error()>options.tolerance; ) {
      assert( repr.n()==d ) ;
      //
      // Determine maximum values of the residual
      //
      glas2::vector<value_type> local_values( fun_sequence.num_terms() ) ;
      real_value_type r_max = 0.0 ;
      for (int i=0; i<test_set.size(); ++i) {
        // Evaluate representation
        repr.eval( test_set(i), local_values ) ;
        local_values -= test_values( i, glas2::all() ) ;
        //std::cout << "R " << local_values << std::endl;
        int local_ind = max_ind( abs( local_values ) ) ;
        if (r_max<std::abs( local_values(local_ind) )) {
          r_max = std::abs( local_values(local_ind) ) ;
          ind_C(d) = i ;
        }
      }

      repr.error() = r_max ;
      if (options.debug_level>0) std::cout << "Degree " << d << " has approximation error " << repr.error() << std::endl ;
      if (r_max<options.tolerance) break ;

      bool is_complex = test_set(ind_C(d)).imag()!=0.0 ;
      if (is_complex) {
        ind_C(d+1) = ind_C(d)+complex_offset ;
      }

      // Add interpolation point
      if (is_complex) {
        repr.add_node( 0.0, test_set(ind_C(d)), real(test_values(ind_C(d),glas2::all())) ) ;
        repr.add_node( 0.0, conj(test_set(ind_C(d))), imag(test_values(ind_C(d),glas2::all())) ) ;
      } else {
        repr.add_node( 0.0, test_set(ind_C(d)), real(test_values(ind_C(d),glas2::all())) ) ;
      }
      assert( !is_complex && (repr.n() == d+1) || is_complex && (repr.n() == d+2) ) ;

      int new_cols = 1+is_complex ;

      // Reshape columns of L as matrices.
      auto Qi = reshape( Q(glas2::all(), d ), n_points, fun_sequence.num_terms(), glas2::column_major() ) ;
      auto Qi1 = reshape( Q(glas2::all(), d+1 ), n_points, fun_sequence.num_terms(), glas2::column_major() ) ;

      // Set columns of L
      if (is_complex) {
        for (int i=0; i<test_set.size(); ++i) {
          if (i!=ind_C(d)) {
            glas2::vector<value_type> term_rr = glas2::copy( ( test_values(i, glas2::all()) - test_values(ind_C(d),glas2::all()) ) / (test_set(i) - test_set(ind_C(d)) ) ) ;
            //std::cout << "C term_rr " << term_rr << std::endl ;
            if (test_set(i).imag()==0.) {
              Qi( i, glas2::all() ) = 2.0*glas2::real( term_rr ) ;
              Qi1( i, glas2::all() ) = -2.0*glas2::imag( term_rr ) ;
            } else {
              glas2::vector<value_type> term_rc = glas2::copy( ( test_values(i, glas2::all()) - conj(test_values(ind_C(d),glas2::all())) ) / (test_set(i) - conj(test_set(ind_C(d))) ) ) ;
              glas2::vector<value_type> term_1 = copy( term_rr + term_rc ) ;
              glas2::vector<value_type> term_2 = copy( term_rr - term_rc ) ;
              Qi( i, glas2::all() ) = std::sqrt(0.5) * glas2::real( term_1 ) ;
              Qi( i+complex_offset, glas2::all() ) = std::sqrt(0.5) * glas2::imag (term_1 ) ;
              Qi1( i, glas2::all() ) = -std::sqrt(0.5) * glas2::imag( term_2 ) ;
              Qi1( i+complex_offset, glas2::all() ) = std::sqrt(0.5) * glas2::real (term_2 ) ;
            }
          }
        }
      } else {
        for (int i=0; i<test_set.size(); ++i) {
          if (i!=ind_C(d)) {
            glas2::vector<value_type> term_rr = glas2::copy( ( test_values(i, glas2::all()) - test_values(ind_C(d),glas2::all()) ) / (test_set(i) - test_set(ind_C(d)) ) ) ;
            //std::cout << "R term_rr " << term_rr << std::endl ;
            if (test_set(i).imag()==0.) {
              assert( norm_2( imag( term_rr ) )==0.0 ) ;
              Qi( i, glas2::all() ) = glas2::real( term_rr ) ;
            } else {
              Qi( i, glas2::all() ) = std::sqrt(2.0)*glas2::real( term_rr ) ;
              Qi( i+complex_offset, glas2::all() ) = -std::sqrt(2.0)*glas2::imag( term_rr ) ;
            }
          }
        }
      }
      //std::cout << "Qj " << Q(glas2::all(),d) << std::endl ;
      //std::cout << "0.5*Qj " << 0.5*Q(glas2::all(),d) << std::endl ;
      //std::cout << "Qj1 " << Q(glas2::all(),d+1) << std::endl ;
      //std::cout << "0.5*Qj1 " << 0.5*Q(glas2::all(),d+1) << std::endl ;

      //bool efficient = false ;
      if (!options.fast) {
      // Update H and S to compensate for the removal of the rows
      glas2::range d_r(0,d) ;
      if (d>0) {
        glas2::vector< real_value_type > temp( d ) ;
        glas2::matrix<real_value_type> ee( d, d ) ;
        ee = glas2::identity_matrix<real_value_type>(d,d) ;
        // Do this for all functions
        for (int i=0; i<fun_sequence.num_terms(); ++i) {
          for (int ii=0; ii<new_cols; ++ii) {
            temp(d_r) = multiply( transpose(S(d_r,d_r)), Q(i*n_points+ind_C(d+ii),d_r) ) ;
            ee -= outer_prod(temp(d_r),temp(d_r)) ;
          }
        }
        //std::cout << "ee " << std::setprecision(16) << ee << std::endl ;
        real_value_type Si_norm = norm_1(ee) ;
        auto Si( glas2::upper(ee) ) ;
        int info = boost::numeric::bindings::lapack::potrf( Si ) ;
        if (info==0) {
            // We should perform the following operation here, but this may correct H when
            // cheaply_correct is set to false further.
            // H(d_r,d_r) = Si * H(d_r,d_r) ;
            // S(d_r,d_r) = S(d_r,d_r) * inv(Si) ;
          inplace_solve( S(d_r,d_r), Si ) ;

          // Check condition number
          real_value_type Si_cond ;
          info = boost::numeric::bindings::lapack::pocon( Si, Si_norm, Si_cond ) ;
          assert( info==0 ) ;
            // H(d_r,d_r) = Si * H(d_r,d_r) ;
          inplace_multiply( Si, H(d_r,d_r) ) ;
          //std::cout << "S" << S(d_r,d_r) << std::endl ;
        }

      }

      // Remove rows in d-th and d+1-st columns of Q.
      glas2::matrix<int> selection( d+new_cols, fun_sequence.num_terms() ) ;
      selection( glas2::all(), 0 ) = ind_C( glas2::range(0,d+new_cols) ) ;
      for (int i=1; i<fun_sequence.num_terms(); ++i) {
        selection( glas2::all(), i ) = selection( glas2::all(), i-1 ) + n_points ;
      }
      fill( Q( vec(selection), glas2::range(0,d+new_cols)), 0.0 ) ;
    
      // Gram-Schmidt of columns d..d+new_cols-1 of Q.
      //
      for (int j=0; j<new_cols; ++j) {
        glas2::vector< real_value_type > temp(d+j) ;
        glas2::range r_j(0,d+j) ;
        auto Q_rj = Q( glas2::all(), r_j ) ;
        auto Qj = Q( glas2::all(), d+j ) ;
      //std::cout << "Qj " << Qj << std::endl ;
        real_value_type nv = norm_2( Qj );
        temp = multiply( transpose(conj(Q_rj)), Qj ) ;
        H(r_j,d+j) = multiply( transpose(conj(S(r_j,r_j))), temp ) ;
        temp = multiply( S(r_j,r_j), H(r_j,d+j) ) ;
        Qj -= multiply(Q_rj, temp ) ;
        //std::cout << "Qj " << Qj << std::endl ;
        H(d+j,d+j) = norm_2( Qj ) ;

        // Reorthogonalization is necessary for higher precision
        glas2::vector<real_value_type> h_new( d+j ) ;

        if (d+j>0) {
          for (int ot = 0; (ot<2) && (H(d+j,d+j) < 0.707 * nv); ++ot) {
            temp = multiply( transpose(conj(Q_rj)), Qj ) ;
            h_new = multiply( transpose(conj(S(r_j,r_j))), temp ) ;
            temp = multiply( S(r_j,r_j), h_new ) ;
            Qj -= multiply(Q_rj, temp ) ;
            H(r_j,d+j) += h_new ;
  
            nv = H(d+j,d+j) ;
            H(d+j,d+j) = norm_2( Qj ) ;
          }
        }
        Qj /= H(d+j,d+j) ;
        //std::cout << "H = " << H(glas2::range(0,d+j+1), glas2::range(0,d+j+1)) << std::endl ;
      } // Gram-Schmidt done
      //std::cout << "Q " << nice_print(Q(glas2::all(),glas2::range(0,d+new_cols) ))<< std::endl ;


      glas2::vector<real_value_type>  svd_val(d+new_cols) ;
      glas2::matrix<real_value_type> svd_u(1,1) ;
      glas2::matrix<real_value_type> svd_vh(d+new_cols,d+new_cols) ;
      svd_vh = H(glas2::range(0,d+new_cols),glas2::range(0,d+new_cols)) ;
      //std::cout << "H " << svd_vh << std::endl ;
#ifndef NDEBUG
      int info =
#endif
        boost::numeric::bindings::lapack::gesvd( 'N', 'O', svd_vh, svd_val, svd_u, svd_u ) ;
      assert(info==0) ;
      //std::cout << "svd : " << svd_val << std::endl ;

      assert( repr.n()==d+new_cols ) ;
      int ii = 0 ;
      for (int i=0; i<d+new_cols; ++i) {
        if (ind_C(i)<test_set.size()) {
          if (test_set(ind_C(i)).imag()==0) {
            repr.weights()(ii) = svd_vh(d+new_cols-1,i) ;
            ++ii ;
          }  else {
            repr.weights()(ii) = value_type( svd_vh(d+new_cols-1,i), svd_vh(d+new_cols-1,i+1) ) ;
            ++ii ;
            repr.weights()(ii) = conj(repr.weights()(ii-1)) ;
            ++ii ;
            ++i ;
          }
        }
      }
      } else {
      // Remove rows in d-th and d+1-st columns of Q.
      glas2::matrix<int> selection( d+new_cols, fun_sequence.num_terms() ) ;
      selection( glas2::all(), 0 ) = ind_C( glas2::range(0,d+new_cols) ) ;
      for (int i=1; i<fun_sequence.num_terms(); ++i) {
        selection( glas2::all(), i ) = selection( glas2::all(), i-1 ) + n_points ;
      }
      fill( Q( vec(selection), glas2::range(0,d+new_cols)), 0.0 ) ;

        glas2::vector<real_value_type>  svd_val(d+new_cols) ;
        glas2::matrix<real_value_type> svd_u(1,1) ;
        glas2::matrix<real_value_type> svd_L( Q.num_rows(), d+new_cols ) ;
        svd_L = Q( glas2::all(), glas2::range(0,d+new_cols) ) ;
        int info = boost::numeric::bindings::lapack::gesvd( 'N', 'O', svd_L, svd_val, svd_u, svd_u ) ;
      assert(info==0) ;
      //std::cout << "svd : " << svd_val << std::endl ;

      assert( repr.n()==d+new_cols ) ;
      int ii = 0 ;
      for (int i=0; i<d+new_cols; ++i) {
        if (ind_C(i)<test_set.size()) {
          if (test_set(ind_C(i)).imag()==0) {
            repr.weights()(ii) = svd_L(d+new_cols-1,i) ;
            ++ii ;
          }  else {
            repr.weights()(ii) = value_type( svd_L(d+new_cols-1,i), svd_L(d+new_cols-1,i+1) ) ;
            ++ii ;
            repr.weights()(ii) = conj(repr.weights()(ii-1)) ;
            ++ii ;
            ++i ;
          }
        }
      }
      }
      //std::cout << "nodes " << repr.nodes() << std::endl ;
      //std::cout << "weights " << repr.weights() << std::endl ;
      //std::cout << "coefficients " << repr.coefficients() << std::endl ;

      d += new_cols ;
    }

    for (int i=0; i<fun_sequence.num_terms(); ++i) {
      repr.coefficients()(glas2::all(),i) *= norm_f(i) ;
    }

    // Remove small weights
    real_value_type norm_w = norm_2( repr.weights() ) ;
    //std::cout << "weights " << repr.weights() << std::endl ;
    for (int i=0; i<repr.weights().size(); ++i) {
      if (norm_2( repr.weights()( glas2::range_from_end(i,0) ) ) < options.tolerance * norm_w ) {
        repr.n() = i ;
        break ;
      } 
    }

    // Sort nodes and weights so that real nodes come first.
    glas2::vector< value_type > nodes( repr.n() ) ;
    glas2::vector< value_type > weights( repr.n() ) ;
    glas2::matrix< real_value_type, glas2::row_major > coefficients( repr.coefficients().num_rows(), repr.coefficients().num_columns() ) ;
    nodes = repr.nodes() ;
    weights = repr.weights() ;
    coefficients = repr.coefficients() ;
    int i_real_end = 0;
    int i_complex_end = nodes.size() ;
    for (int i=0; i<nodes.size(); ++i) {
      if (nodes(i).imag()==0.0) {
        repr.nodes()( i_real_end ) = nodes(i) ;
        repr.weights()( i_real_end ) = weights(i) ;
        repr.coefficients()( i_real_end, glas2::all() ) = coefficients(i, glas2::all()) ;
        ++i_real_end ;
      } else {
        assert( i<nodes.size()-1) ;
        assert( std::conj(nodes(i))==nodes(i+1) ) ;
        assert( std::conj(weights(i))==weights(i+1) ) ;
        repr.nodes()( i_complex_end-2 ) = nodes(i) ;
        repr.nodes()( i_complex_end-1 ) = nodes(i+1) ;
        repr.weights()( i_complex_end-2 ) = weights(i) ;
        repr.weights()( i_complex_end-1 ) = weights(i+1) ;
        repr.coefficients()( i_complex_end-2, glas2::all() ) = coefficients(i, glas2::all()) ;
        repr.coefficients()( i_complex_end-1, glas2::all() ) = coefficients(i+1, glas2::all()) ;
        i_complex_end -= 2 ;
        ++i ;
      }
    }
    assert( i_real_end==i_complex_end ) ;

    return repr ;
  } // SV_AAA()

} } // namespace CORK::approximation

#endif
