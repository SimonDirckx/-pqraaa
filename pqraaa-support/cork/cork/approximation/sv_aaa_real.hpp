//  (C) Copyright Karl Meerbergen 2018.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_approximation_sv_aaa_real_hpp
#define cork_approximation_sv_aaa_real_hpp

#include <cork/approximation/sv_aaa_approximation_real.hpp>
#include <cork/approximation/sv_aaa_2_partial_fractions.hpp>
#include <cork/approximation/sv_aaa_poles_real.hpp>
#include <cork/approximation/unique_test_set.hpp>
#include <cork/exception/rational_approximation.hpp>
#include <cork/options/max_degree.hpp>
#include <cork/options/debug_level.hpp>
#include <cork/options/number_of_sample_points.hpp>
#include <cork/options/aaa_stop_criterion.hpp>
#include <cork/options/compute_poles.hpp>
#include <cork/options/value_of.hpp>
#include <cork/options/aaa_fast.hpp>
#include <cork/options/aaa_always_use_infinite_support_point.hpp>
#include <cork/options/aaa_max_stagnation_iterations.hpp>
#include <cork/options/aaa_cheap_correction_tolerance.hpp>
#include <cork/utility/is_infinite.hpp>
#include <glas2/matrix.hpp>
#include <glas2/vector.hpp>
#include <glas2/bindings.hpp>
#include <boost/numeric/bindings/lapack/driver/gesvd.hpp>
#include <boost/numeric/bindings/lapack/computational/pocon.hpp>
#include <boost/numeric/bindings/lapack/computational/potrf.hpp>
#include <cassert>
#include <cmath>
#include <limits>
#include <iomanip>

namespace CORK { namespace approximation {
  
  template <typename TestSet, typename TestValues, typename IndC, typename QQ, typename Repr>
  int add_column_L( TestSet const& test_set, TestValues const& test_values, int d, IndC& ind_C, int n_points, int n_funs, int complex_offset, QQ& Q, Repr& repr ) {
    typedef typename TestValues::value_type value_type ;
    typedef decltype(std::abs(value_type())) real_value_type ;

      bool is_complex = test_set(ind_C(d)).imag()!=0.0 ;
      if (is_complex) {
        ind_C(d+1) = ind_C(d)+complex_offset ;
      }

      // Add interpolation point
      if (is_complex) {
        repr.add_node( 0.0, test_set(ind_C(d)), glas2::real(test_values(ind_C(d),glas2::all())) ) ;
        repr.add_node( 0.0, std::conj(test_set(ind_C(d))), glas2::imag(test_values(ind_C(d),glas2::all())) ) ;
      } else {
        assert( norm_inf(imag(test_values(ind_C(d),glas2::all()))) == 0.0 ) ;
        repr.add_node( 0.0, test_set(ind_C(d)), real(test_values(ind_C(d),glas2::all())) ) ;
      }
      assert( !is_complex && (repr.n() == d+1) || is_complex && (repr.n() == d+2) ) ;

      int new_cols = 1+is_complex ;

      // Reshape columns of L as matrices.
      auto Qi = reshape( Q(glas2::all(), d ), n_points, n_funs, glas2::column_major() ) ;
      auto Qi1 = reshape( Q(glas2::all(), d+1 ), n_points, n_funs, glas2::column_major() ) ;

      // Set columns of L
      if (is_complex) {
        for (int i=0; i<test_set.size(); ++i) {
          glas2::vector<value_type> term_rr( test_values.num_columns() ) ;
          glas2::vector<value_type> term_rc( test_values.num_columns() ) ;
          if (i!=ind_C(d)) {
            if (is_infinite(std::abs(test_set(ind_C(d))))) {
              term_rr = test_values(i, glas2::all()) - test_values(ind_C(d),glas2::all()) ;
              term_rc = test_values(i, glas2::all()) - glas2::conj(test_values(ind_C(d),glas2::all())) ;
              //term_rr = test_values(i, glas2::all()) ;
            } else if (std::abs(test_set(i))==std::numeric_limits<real_value_type>::infinity()) {
              fill(term_rr,0.0); 
              fill(term_rc,0.0); 
              //term_rr = test_values(i, glas2::all()) - test_values(ind_C(d),glas2::all()) ;
              //term_rc = test_values(i, glas2::all()) - glas2::conj(test_values(ind_C(d),glas2::all())) ;
            } else {
              term_rr = ( test_values(i, glas2::all()) - test_values(ind_C(d),glas2::all()) ) / (test_set(i) - test_set(ind_C(d)) ) ;
              term_rc = ( test_values(i, glas2::all()) - conj(test_values(ind_C(d),glas2::all())) ) / (test_set(i) - conj(test_set(ind_C(d))) ) ;
            }
            //std::cout << "C term_rr " << term_rr << std::endl ;
            if (test_set(i).imag()==0.) {
              Qi( i, glas2::all() ) = 2.0*glas2::real( term_rr ) ;
              Qi1( i, glas2::all() ) = -2.0*glas2::imag( term_rr ) ;
            } else {
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
          glas2::vector<value_type> term_rr( test_values.num_columns() ) ;
          if (i!=ind_C(d)) {
            if (is_infinite(std::abs(test_set(ind_C(d))))) {
              term_rr = test_values(i, glas2::all()) - test_values(ind_C(d),glas2::all()) ;
            } else if (is_infinite(std::abs(test_set(i)))) {
              fill( term_rr, 0.0 ) ;
              //term_rr = test_values(i, glas2::all()) - test_values(ind_C(d),glas2::all()) ;
            //std::cout << "Inf term_rr " << test_set(i) << " " << test_values(i, glas2::all()) << " " << test_values(ind_C(d),glas2::all()) << " " << term_rr << std::endl ;
            } else {
              term_rr = ( test_values(i, glas2::all()) - test_values(ind_C(d),glas2::all()) ) / (test_set(i) - test_set(ind_C(d)) ) ;
            //std::cout << "Fin term_rr " << test_values(i, glas2::all()) << " " << test_values(ind_C(d),glas2::all()) << " " << term_rr << std::endl ;
            }
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

      return new_cols ;
  } // add_column_L()


  template <typename SamplePoints, typename FunctionValues, typename FindWorstIndex, typename Options>
  typename std::enable_if< glas2::is<glas2::DenseMatrix, FunctionValues>::value
                         , SV_AAA_approximation_real<typename FunctionValues::value_type>
                         >::type
              SV_AAA_real( SamplePoints const& test_set1, FunctionValues const& test_values1, FindWorstIndex const& find_worst_index, Options const& options ) {
    typedef typename FunctionValues::value_type  value_type ;
    typedef decltype(std::abs(value_type()))     real_value_type ;

    // Count the number of complex points
    int n_c = 0.0 ;
    for (int i=0;i<test_set1.size();++i) {
      n_c += test_set1(i).imag()!=0.0 ;
    }

    // Reorder the points.
    int n_points = n_c*2 + (test_set1.size()-n_c) ;
    glas2::vector<value_type> test_set( test_set1.size() ) ;
    glas2::matrix<value_type> test_values( test_values1.num_rows(), test_values1.num_columns() ) ;
    glas2::matrix<value_type> R( test_values1.num_rows(), test_values1.num_columns() ) ;
    int ii = 0 ;
    for (int i=0;i<test_set1.size();++i) {
      if (test_set1(i).imag()!=0.0) {
        test_set(ii) = test_set1(i) ;
        test_values(ii,glas2::all()) = test_values1(i,glas2::all()) ;
        ++ii ;
      }
    }
    for (int i=0;i<test_set1.size();++i) {
      if (test_set1(i).imag()==0.0) {
        test_set(ii) = test_set1(i) ;
        test_values(ii,glas2::all()) = test_values1(i,glas2::all()) ;
        ++ii ;
      }
    }
    assert( ii==test_set.size() ) ;

    int complex_offset = test_set1.size() ;
    //std::cout << "test_set " << test_set << std::endl ;
    assert( complex_offset+n_c==n_points ) ;

    int max_degree = options::value_of<options::max_degree>(options) ;
    glas2::vector<int> ind_C(max_degree+1) ; // We add one position as sentinel for complex supports.
    glas2::matrix< real_value_type > Q(n_points*test_values.num_columns(), max_degree+1); fill(Q,0.0) ; // We add one column as sentinel for complex supports.
    glas2::matrix< real_value_type > S(max_degree, max_degree); S = glas2::eye( S.num_rows(), S.num_columns() ) ;
    glas2::matrix< real_value_type > H(max_degree+1, max_degree+1); fill(H,0.0) ; // Sentinel

    SV_AAA_approximation_real<value_type> repr( test_values.num_columns(), max_degree ) ;
    repr.reset( 0 ) ;

    glas2::vector<real_value_type> norm_f( test_values.num_columns() ) ;

    // Scale the functions
    for (typename FunctionValues::size_type j=0; j<test_values.num_columns(); ++j) {
      norm_f(j) = norm_inf( test_values(glas2::all(),j) ) ;
      if (norm_f(j)!=0.0) test_values(glas2::all(),j) /= norm_f(j) ;
    }
    //std::cout << "test_values " << test_values << std::endl ;

    // The approximation is now zero. So, the error is the function itself
    repr.error() = norm_inf( vec(test_values) ) ;
    int number_stagnation_iterations = 0 ;

    assert( max_degree>0 ) ;
    real_value_type tolerance = std::max(0.,options::value_of<options::aaa_stop_criterion<real_value_type>>(options).tolerance()) ;

    for (int d=0; d<max_degree && repr.error()>tolerance; ) {
      assert( repr.n()==d ) ;
      //
      // Determine maximum values of the residual
      //
      glas2::vector<value_type> local_values( test_values.num_columns() ) ;
      real_value_type r_max = 0.0 ;
      for (int i=0; i<test_set.size(); ++i) {
        // Evaluate representation
        repr.eval( test_set(i), R( i, glas2::all() ) ) ;
        //R( i,  -= test_values( i, glas2::all() ) ;
        //std::cout << "R " << local_values << std::endl;
        /*auto local = norm_inf( local_values ) ;
        if (r_max<local) {
          r_max = local ;
           ind_C(d) = i ;
        }*/
      }
      R -= test_values ;
      auto [find_worst_i, find_worst_r] = find_worst_index( test_set, R ) ;
      r_max = find_worst_r ;
      if (d==0 && options::value_of<options::aaa_always_use_infinite_support_point>(options)) {
        auto it = std::find( test_set.begin(), test_set.end(), std::numeric_limits<real_value_type>::infinity() ) ;
        if (it!=test_set.end()) ind_C(d) = it-test_set.begin() ;
        else ind_C(d) = find_worst_i ;
        std::cout << "infinity " << it-test_set.begin() << " " << ind_C(d) <<  " " << test_set(ind_C(d)) << std::endl ;
      } else {
        ind_C(d) = find_worst_i ;
      }

      real_value_type previous_error = repr.error() ;
      repr.error() = r_max ;
      if (options::value_of<options::debug_level>(options)>0) std::cout << "Degree " << d << " has approximation error " << repr.error() << std::endl ;
      if (r_max<tolerance) break ;

      if (d>0 && repr.error()>previous_error && repr.error()<std::sqrt(std::numeric_limits<real_value_type>::epsilon())) {
        ++number_stagnation_iterations ;
        if (number_stagnation_iterations>options::value_of<options::aaa_max_stagnation_iterations>(options)) {
          repr.warnings().push_back( std::string("Stagnation detecting in AAA") ) ;
          if (options::value_of<options::debug_level>(options)>1) std::cout << "Stagnation detecting in AAA: " << number_stagnation_iterations << std::endl ;
          break ;
        }
        // Check for Froissart doublets
        auto pf = sv_aaa_2_partial_fractions( repr ) ;
        real_value_type norm_coef = norm_fro( pf.coefficients() ) ;
        std::vector< int > froissart_points ;
        for (int i=0; i<pf.n(); ++i) {
          if (norm_2(pf.coefficients()(i,glas2::all()))<std::numeric_limits<real_value_type>::epsilon()*norm_coef) {
            froissart_points.push_back( glas2::max_ind( - abs_squared( pf.nodes()(i)-repr.nodes() ) ) ) ;
          }
        }
        if (froissart_points.size()>0) {
          repr.warnings().push_back( std::string("Froissart doublets detected") ) ;
          if (options::value_of<options::debug_level>(options)>1) std::cout << froissart_points.size() << " Froissart doublets detected." << std::endl ;
          break ;
        }
      }

      int new_cols = add_column_L( test_set, test_values, d, ind_C, n_points, test_values.num_columns(), complex_offset, Q, repr ) ;
        //std::cout << "Qj " << Q(glas2::all(),d) << std::endl ;
        //std::cout << "0.5*Qj " << 0.5*Q(glas2::all(),d) << std::endl ;
        //std::cout << "Qj1 " << Q(glas2::all(),d+1) << std::endl ;
        //std::cout << "0.5*Qj1 " << 0.5*Q(glas2::all(),d+1) << std::endl ;

      //bool efficient = false ;
      if (!options::value_of<options::aaa_fast>(options)) {
      // Update H and S to compensate for the removal of the rows
      glas2::range d_r(0,d) ;
      if (d>0) {
        glas2::vector< real_value_type > temp( d ) ;
        glas2::matrix<real_value_type> ee( d, d ) ;
        ee = glas2::identity_matrix<real_value_type>(d,d) ;
        // Do this for all functions
        for (int i=0; i<test_values.num_columns(); ++i) {
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
      glas2::matrix<int> selection( d+new_cols, test_values.num_columns() ) ;
      selection( glas2::all(), 0 ) = ind_C( glas2::range(0,d+new_cols) ) ;
      for (int i=1; i<test_values.num_columns(); ++i) {
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
        if (H(d+j,d+j)!=0.0) Qj /= H(d+j,d+j) ;
        //std::cout << "H = " << H(glas2::range(0,d+j+1), glas2::range(0,d+j+1)) << std::endl ;
      } // Gram-Schmidt done
      //std::cout << "Q " << nice_print(Q(glas2::all(),glas2::range(0,d+new_cols) ))<< std::endl ;


      glas2::vector<real_value_type> svd_val(d+new_cols) ;
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
      //std::cout << "smallest singular value " << svd_val(svd_val.size()-1) << std::endl ;

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
      glas2::matrix<int> selection( d+new_cols, test_values.num_columns() ) ;
      selection( glas2::all(), 0 ) = ind_C( glas2::range(0,d+new_cols) ) ;
      for (int i=1; i<test_values.num_columns(); ++i) {
        selection( glas2::all(), i ) = selection( glas2::all(), i-1 ) + n_points ;
      }
      fill( Q( vec(selection), glas2::range(0,d+new_cols)), 0.0 ) ;

        glas2::vector<real_value_type>  svd_val(d+new_cols) ;
        glas2::matrix<real_value_type> svd_u(1,1) ;
        glas2::matrix<real_value_type> svd_L( Q.num_rows(), d+new_cols ) ;
        svd_L = Q( glas2::all(), glas2::range(0,d+new_cols) ) ;
#ifndef NDEBUG
        int info =
#endif
          boost::numeric::bindings::lapack::gesvd( 'N', 'O', svd_L, svd_val, svd_u, svd_u ) ;
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

    // Remove Froissart doublets
    auto pf = sv_aaa_2_partial_fractions( repr ) ;
    real_value_type norm_coef = norm_fro( pf.coefficients() ) ;
    std::vector< int > froissart_points ;
    for (int i=0; i<pf.n(); ++i) {
      if (norm_2(pf.coefficients()(i,glas2::all()))<std::numeric_limits<real_value_type>::epsilon()*norm_coef) {
        froissart_points.push_back( glas2::max_ind( - abs_squared( pf.nodes()(i)-repr.nodes() ) ) ) ;
      }
    }
    if (froissart_points.size()>0) {
      int d = repr.n() ;
      std::sort( froissart_points.begin(), froissart_points.end() ) ;
      froissart_points.push_back(repr.n()) ; // sentinel element.
      int ii = froissart_points[0];
      for (typename std::vector< int >::size_type i=1; i<froissart_points.size(); ++i) {
        for ( ; ii+i<froissart_points[i]; ++ii ) {
          ind_C(ii) = ind_C(ii+i) ;
        }
        --d ;
      }
      assert( repr.n()==d+froissart_points.size()-1) ; // -1 because of sentinel element
      {
        repr.n() = 0 ;
        for (int i=0; i<d; ) {
          if (ind_C(i)<test_set.size())
            i += add_column_L( test_set, test_values, i, ind_C, n_points, test_values.num_columns(), complex_offset, Q, repr ) ;
          else
            ++i ;
        }
      }
      for (int i=0; i<test_values.num_columns(); ++i) {
        fill( Q(i*test_set.size()+ind_C(glas2::range(0,d)),glas2::range(0,d)), 0.0 ) ;
      }
      auto Qd = Q(glas2::all(), glas2::range(0,d)) ;
      glas2::matrix<real_value_type> svd_u(d,d) ;
      glas2::vector<real_value_type>  svd_val(d) ;
      boost::numeric::bindings::lapack::gesvd( 'N', 'A', Qd, svd_val, svd_u, svd_u ) ;
      ii = 0 ;
      for (int i=0; i<d; ++i) {
        if (ind_C(i)<test_set.size()) {
          if (test_set(ind_C(i)).imag()==0) {
            repr.weights()(ii) = svd_u(d-1,i) ;
            ++ii ;
          }  else {
            repr.weights()(ii) = value_type( svd_u(d-1,i), svd_u(d-1,i+1) ) ;
            ++ii ;
            repr.weights()(ii) = conj(repr.weights()(ii-1)) ;
            ++ii ;
            ++i ;
          }
        }
      }
      repr.warnings().push_back( std::string("Froissart doublets detected") ) ;
      if (options::value_of<options::debug_level>(options)>1) {
        std::cout << froissart_points.size() << " Froissart doublets detected." << std::endl ;
        std::cout << "Reduction of the degree from " << d+froissart_points.size() << " to " << repr.n() << std::endl ;
      }

      // Recompute error
      for (int i=0; i<test_set.size(); ++i) {
        repr.eval( test_set(i), R( i, glas2::all() ) ) ;
      }
      R -= test_values ;
      int ind_Q = max_ind( abs( vec(R) ) ) ;
      repr.error() = std::abs( vec(R)(ind_Q) ) ;
    } // if (froissart_points.size()>0)

    for (int i=0; i<test_values.num_columns(); ++i) {
      repr.coefficients()(glas2::all(),i) *= norm_f(i) ;
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

    // Compute poles
    if (options::value_of<options::compute_poles>(options)) {
      auto poles = sv_aaa_poles_real( repr ) ;
      std::cout << "Poles " << poles << std::endl ;
    }

    return repr ;
  } // SV_AAA_real()

  template <typename FunctionSequence, typename Domain, typename FindWorstIndex, typename Options>
  typename std::enable_if< !glas2::is<glas2::DenseMatrix, FunctionSequence>::value
                         , SV_AAA_approximation_real<typename FunctionSequence::template value_type_for< typename Domain::value_type > >
                         >::type
                SV_AAA_real( FunctionSequence const& fun_sequence, FindWorstIndex const& find_worst_index, Domain const& domain, Options const& options ) {
    typedef typename Domain::value_type                                               argument_value_type ;
    typedef typename FunctionSequence::template value_type_for< argument_value_type > value_type ;

    // Construct the test set
    auto test_set_ini = domain.discretize( std::max(2*options::value_of<options::max_degree>(options),
                                                 options::value_of<options::number_of_sample_points>(options)) ) ;
    auto test_set1 = unique_test_set( test_set_ini, []( auto const& x, auto const& y ) { return x==y || std::conj(x)==y ; } ) ;

    // Evaluate and scale the functions
    glas2::matrix<value_type> test_values( test_set1.size(), fun_sequence.num_terms() ) ;
    for (typename glas2::vector<value_type>::size_type i=0; i<test_set1.size(); ++i) {
      fun_sequence.evaluate( test_set1(i), test_values(i,glas2::all()) ) ;
      if (options::value_of<options::debug_level>(options)>0) {
        if (is_infinite(glas2::norm_2(test_values(i,glas2::all()))))
          throw exception::rational_approximation( "AAA approximation: infinite function values in sample point" ) ;
        if (std::isnan(glas2::norm_2(test_values(i,glas2::all()))))
          throw exception::rational_approximation( "AAA approximation: NaN values in sample point" ) ;
      }
    }

    return SV_AAA_real( test_set1, test_values, find_worst_index, options ) ;
  } // SV_AAA_real()

  template <typename FunctionSequence, typename Domain>
  decltype (auto) SV_AAA_real( FunctionSequence const& fun_sequence, Domain const& domain ) {
    return SV_AAA_real( fun_sequence, domain, std::tuple<>() ) ;
  } // SV_AAA_real()

} } // namespace CORK::approximation

#endif
