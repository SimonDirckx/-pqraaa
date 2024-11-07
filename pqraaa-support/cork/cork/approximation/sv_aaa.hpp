//  (C) Copyright Karl Meerbergen 2018.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_approximation_sv_aaa_hpp
#define cork_approximation_sv_aaa_hpp

#include <cork/exception/rational_approximation.hpp>
#include <cork/approximation/sv_aaa_approximation.hpp>
#include <cork/approximation/unique_test_set.hpp>
#include <cork/approximation/sv_aaa_2_partial_fractions.hpp>
#include <cork/options/max_degree.hpp>
#include <cork/options/aaa_max_stagnation_iterations.hpp>
#include <cork/options/aaa_stop_criterion.hpp>
#include <cork/options/debug_level.hpp>
#include <cork/options/number_of_sample_points.hpp>
#include <cork/options/value_of.hpp>
#include <cork/options/aaa_fast.hpp>
#include <cork/options/aaa_cheap_correction_tolerance.hpp>
#include <glas2/matrix.hpp>
#include <glas2/vector.hpp>
#include <glas2/bindings.hpp>
#include <boost/numeric/bindings/lapack/driver/gesvd.hpp>
#include <boost/numeric/bindings/lapack/computational/pocon.hpp>
#include <boost/numeric/bindings/lapack/computational/potrf.hpp>
#include <cassert>
#include <tuple>
#include <limits>
#include <iomanip>

namespace CORK { namespace approximation {


  template <typename Points, typename FunctionValues, typename Options, typename FindWorstIndex>
  typename std::enable_if< glas2::is<glas2::DenseMatrix, FunctionValues>::value
                         , SV_AAA_approximation<typename Points::value_type, typename FunctionValues::value_type>
                         >::type SV_AAA( Points const& test_set, FunctionValues test_values, FindWorstIndex const& find_worst_index, Options const& options, bool scale_functions=false ) {
    typedef typename Points::value_type         argument_type ;
    typedef typename FunctionValues::value_type value_type ;
    typedef decltype(std::abs(value_type()))    real_type ;
    typedef int                                 size_type ;

    // Check if points are close to each other

    assert( test_set.size()>=3 ) ;
    assert( test_values.num_rows()==test_set.size() ) ;

    int LDQ = test_set.size() ;

    glas2::vector<int> ind_C(options::value_of<options::max_degree>(options)) ;
    glas2::matrix< value_type > C(test_set.size(), options::value_of<options::max_degree>(options));
    glas2::matrix< value_type > Q(LDQ*test_values.num_columns(), options::value_of<options::max_degree>(options));
    glas2::matrix< value_type > S(options::value_of<options::max_degree>(options), options::value_of<options::max_degree>(options)); fill(S,0.0) ;
    glas2::matrix< value_type > H(options::value_of<options::max_degree>(options), options::value_of<options::max_degree>(options)); fill(H,0.0) ;
    glas2::vector<value_type> D(test_set.size()) ;
    glas2::matrix<value_type> R(test_set.size(), test_values.num_columns() ) ;

    real_type tolerance = options::value_of<options::aaa_stop_criterion<real_type>>(options).tolerance() ;
    tolerance = std::max(tolerance, 0.0 ) ;
    int number_stagnation_iterations = 0 ;

    SV_AAA_approximation<argument_type,value_type> repr( test_values.num_columns(), options::value_of<options::max_degree>(options) ) ;
    repr.reset( 0 ) ;

    glas2::vector<real_type> norm_f( test_values.num_columns() ) ;

    if (scale_functions) {
      // Evaluate and scale the functions
      for (typename FunctionValues::size_type j=0; j<test_values.num_columns(); ++j) {
        norm_f(j) = norm_inf( test_values(glas2::all(),j) ) ;
        if (std::isnan(norm_f(j)) || std::isinf(norm_f(j))) {
          if (options::value_of<options::debug_level>(options)>1) {
            std::cout << "Function values in the discrete points\n" ;
            std::cout << "Argument ->  Function values:\n" ;
            for (int i=0; i<test_set.size(); ++i) {
              std::cout << test_set(i) << " -> " << test_values(i, glas2::all() ) << std::endl ;
            }
          }
          throw exception::rational_approximation("Functions evaluate to Inf or NaN") ;
        }
        if (norm_f(j)!=0.0) test_values(glas2::all(),j) /= norm_f(j) ;
      }
    }
    for (typename FunctionValues::size_type j=0; j<test_values.num_columns(); ++j) {
      if (norm_1(test_values(glas2::all(),j))<1.1*norm_inf(test_values(glas2::all(),j))) {
        std::stringstream ss ;
        ss<< "AAA: possible singularity in point for function number " << j ;
        repr.warnings().push_back( std::string( ss.str() ) ) ;
      }
    }

    // Compute the SV_AAA approximation
    repr.error() =  abs( vec(test_values)( max_ind( abs( vec(test_values) ) ) ) ) ;
    //fill( R, 0.0 ) ;
    for (typename FunctionValues::size_type j=0; j<test_values.num_columns(); ++j) {
      fill( R(glas2::all(),j), glas2::sum(test_values(glas2::all(),j)) / real_type(test_values.num_rows()) ) ;
    }

    for (; repr.n()<options::value_of<options::max_degree>(options) && repr.error()>tolerance; ) {
      int d = repr.n() ;

      // Find index with largest residual
      /*real_type r_max =0.0;
      for (int i=0; i<test_values.num_rows(); ++i) {
        real_type temp_nrm = norm_2( test_values(i,glas2::all()) - R(i,glas2::all()) ) ;
        if (temp_nrm>r_max) {
          r_max = temp_nrm ;
          ind_C(d) = i ;
        }
      }*/

      real_type previous_error = repr.error() ;
      R -= test_values ;
      auto [find_worst_i, find_worst_r] = find_worst_index( test_set, R ) ;
      ind_C(d) = find_worst_i ;
      repr.error() = find_worst_r ;
      //int ind_Q = max_ind( abs( vec(R) ) ) ;
      //ind_C(d) = ind_Q % test_set.size() ;
      //repr.error() = std::abs( vec(R)(ind_Q) ) ;
      if (options::value_of<options::debug_level>(options)>1) std::cout << "Degree " << d << " has approximation error " << repr.error() << std::endl ;
      if (repr.error()<tolerance) break ;

      if (d>0 && repr.error()>previous_error && repr.error()<std::sqrt(std::numeric_limits<real_type>::epsilon())) {
        ++number_stagnation_iterations ;
        if (number_stagnation_iterations>options::value_of<options::aaa_max_stagnation_iterations>(options)) {
          repr.warnings().push_back( std::string("Stagnation detecting in AAA") ) ;
          if (options::value_of<options::debug_level>(options)>1) std::cout << "Stagnation detecting in AAA: " << number_stagnation_iterations << std::endl ;
          break ;
        }
        std::cout << "d = " << d << " aaa.n() " << repr.n() << std::endl;
        // Check for Froissart doublets
        auto pf = sv_aaa_2_partial_fractions( repr ) ;
        real_type norm_coef = norm_fro( pf.coefficients() ) ;
        std::vector< int > froissart_points ;
        for (int i=0; i<pf.n(); ++i) {
          if (norm_2(pf.coefficients()(i,glas2::all()))<std::numeric_limits<real_type>::epsilon()*norm_coef) {
            froissart_points.push_back( glas2::max_ind( - abs_squared( pf.nodes()(i)-repr.nodes() ) ) ) ;
          }
        }
        if (froissart_points.size()>0) {
          repr.warnings().push_back( std::string("Froissart doublets detected") ) ;
          if (options::value_of<options::debug_level>(options)>1) std::cout << froissart_points.size() << " Froissart doublets detected." << std::endl ;
          break ;
        }
      } else {
        // Add interpolation point
        repr.add_node( 0.0, test_set(ind_C(d)), test_values(ind_C(d),glas2::all()) ) ;
        if (options::value_of<options::debug_level>(options)>2) std::cout << "New support point " << test_set(ind_C(d)) << std::endl ;
        assert( repr.n() == d+1 ) ;

        // Set column of Cauchy matrix
        // Infinity cannot be both an active sample point and a support point
        if (std::abs(repr.nodes()(d))==std::numeric_limits<real_type>::infinity())
          C(glas2::all(), d) = test_set ;
        else
          for (int i=0; i<test_set.size(); ++i) {
            if (std::abs(test_set.size())==std::numeric_limits<real_type>::infinity())
              C(i, d) = test_set(i) ;
            else
              C(i, d) = 1.0 / (test_set(i)-repr.nodes()(d)) ;
          }
        fill( C(ind_C(d), glas2::all()), 0.0 ) ;
        fill( C(ind_C(glas2::range(0,d+1)),d), 0.0 ) ;

        // Column d of Loewner matrix
        auto L = Q(glas2::all(), d) ;
        auto Qi = Q(glas2::all(), glas2::range(0,d)) ;
        for (size_type i=0; i<test_values.num_columns(); ++i) {
          L( glas2::range(i*LDQ,(i+1)*LDQ) ) = C(glas2::all(),d) * test_values(glas2::all(),i) - C(glas2::all(),d) * test_values(ind_C(d),i) ;
        }
        //std::cout << "L " << L << std::endl;

        bool cheaply_correct = options::value_of<options::aaa_fast>(options) ;

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
              if (Si_cond<options::value_of<options::aaa_cheap_correction_tolerance<real_type>>(options)) cheaply_correct = false ;
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
        //std::cout << H(glas2::range(0,d+1),glas2::range(0,d+1)) << std::endl ;

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
        assert( repr.n() == d+1 ) ;
      } // New point added

      d = repr.n()-1 ;

      glas2::vector<real_type>  svd_val(d+1) ;
      glas2::matrix<value_type> svd_u(1,1) ;
      glas2::matrix<value_type> svd_vh(d+1,d+1) ;
      svd_vh = H(glas2::range(0,d+1),glas2::range(0,d+1)) ;
      //std::cout << svd_vh << std::endl ;
      int info = boost::numeric::bindings::lapack::gesvd( 'N', 'O', svd_vh, svd_val, svd_u, svd_u ) ;
      assert(info>=0) ;
      if (info>0) throw exception::rational_approximation("GESVD failed in determining AAA weights") ;
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
    } // end for d

    // Remove Froissart doublets
    std::vector< int > froissart_points ;
    if (repr.n()>0) {
      auto pf = sv_aaa_2_partial_fractions( repr ) ;
      real_type norm_coef = norm_fro( pf.coefficients() ) ;
      for (int i=0; i<pf.n(); ++i) {
        if (norm_2(pf.coefficients()(i,glas2::all()))<std::numeric_limits<real_type>::epsilon()*norm_coef) {
          froissart_points.push_back( glas2::max_ind( - abs_squared( pf.nodes()(i)-repr.nodes() ) ) ) ;
        }
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
          //C( glas2::all(), ii ) = C( glas2::all(), ii+i ) ;
        }
        --d ;
      }
      assert( repr.n()==d+froissart_points.size()-1) ; // -1 because of sentinel element
      {
        repr.n() = 0 ;
        for (int i=0; i<d; ++i) {
          repr.add_node( 0.0, test_set(ind_C(i)), test_values(ind_C(i),glas2::all()) ) ;
        }
      }
      for (size_type j=0; j<d; ++j) {
        C(glas2::all(), j) = glas2::ones<real_type>(test_set.size()) / (test_set-repr.nodes()(j)) ;
        for (size_type i=0; i<test_values.num_columns(); ++i) {
          Q( glas2::range(i*test_set.size(),(i+1)*test_set.size()), j ) = C(glas2::all(),j) * test_values(glas2::all(),i) - C(glas2::all(),j) * test_values(ind_C(j),i) ;
        }
      }
      for (int i=0; i<test_values.num_columns(); ++i) {
        fill( Q(i*test_set.size()+ind_C(glas2::range(0,d)),glas2::range(0,d)), 0.0 ) ;
      }
      auto Qd = Q(glas2::all(), glas2::range(0,d)) ;
      glas2::matrix<value_type> svd_u(d,d) ;
      glas2::vector<real_type>  svd_val(d) ;
      boost::numeric::bindings::lapack::gesvd( 'N', 'A', Qd, svd_val, svd_u, svd_u ) ;
      repr.weights() = conj(svd_u(d-1,glas2::all())) ;
      repr.warnings().push_back( std::string("Froissart doublets detected") ) ;
      if (options::value_of<options::debug_level>(options)>1) {
        std::cout << froissart_points.size() << " Froissart doublets detected." << std::endl ;
        std::cout << "Reduction of the degree from " << d+froissart_points.size() << " to " << repr.n() << std::endl ;
      }

      // Recompute error
      glas2::range d_r1(0,d) ;
      D = multiply(C(glas2::all(),d_r1), repr.weights()) ;                     // Denominator
      for (int i=0; i<test_values.num_columns(); ++i) {
        R(glas2::all(),i) = multiply(C(glas2::all(),d_r1), repr.weights()*repr.coefficients()(glas2::all(),i)) ; // Numerator
        R(glas2::all(),i) /= D;
      }
      R(ind_C(d_r1), glas2::all()) = test_values(ind_C(d_r1), glas2::all());
      R -= test_values ;
      int ind_Q = max_ind( abs( vec(R) ) ) ;
      repr.error() = std::abs( vec(R)(ind_Q) ) ;
    } // if (froissart_points.size()>0)

    if (scale_functions) {
      for (int i=0; i<test_values.num_columns(); ++i) {
        repr.coefficients()(glas2::all(),i) *= norm_f(i) ;
      }
    }

    // Remove small weights
/*    real_type norm_w = norm_2( repr.weights() ) ;
    for (int i=0; i<repr.weights().size(); ++i) {
      if (norm_2( repr.weights()( glas2::range_from_end(i,0) ) ) < options::value_of<options::rational_approximation_tolerance<real_type>>(options) * norm_w ) {
        repr.n() = i ;
        break ;
      } 
    }
*/
    return repr ;
  } // SV_AAA()


  template <typename Points, typename FunctionValues, typename FindWorstIndex>
  typename std::enable_if< glas2::is< glas2::DenseVector, Points>::value
                         && glas2::is< glas2::DenseMatrix, FunctionValues>::value
                         , SV_AAA_approximation<typename Points::value_type, typename FunctionValues::value_type>
                         >::type SV_AAA( Points const& points, FunctionValues const& function_values, FindWorstIndex const& find_worst_index ) {
    return SV_AAA( points, function_values, find_worst_index, std::tuple<>() ) ;
  } // SV_AAA()

  /*
  template <typename FunctionSequence, typename Domain, typename Options, typename FindWorstIndex>
  typename std::enable_if< !glas2::is< glas2::DenseVector, FunctionSequence>::value, SV_AAA_approximation<typename Domain::value_type, typename FunctionSequence::template value_type_for<typename Domain::value_type> > >::type SV_AAA( FunctionSequence const& fun_sequence, FindWorstIndex const& find_worst_index, Domain const& domain, Options const& options ) {
    //typedef T                                value_type ;
    typedef typename Domain::value_type                                               argument_value_type ;
    typedef typename FunctionSequence::template value_type_for< argument_value_type > result_value_type ;

    // Construct the test set
    auto test_set_ini = domain.discretize( std::max( 2*options::value_of<options::max_degree>(options), options::value_of<options::number_of_sample_points>(options)) ) ;
    auto test_set = unique_test_set( test_set_ini ) ;
    //test_set_ini.resize( 0 ) ;

    // Evaluate and scale the functions
    glas2::shared_matrix<result_value_type> test_values( test_set.size(), fun_sequence.num_terms() ) ;
    for (typename decltype(test_set)::size_type i=0; i<test_set.size(); ++i) {
      //if (std::abs(test_set(i))!=std::numeric_limits<decltype(std::abs(test_set(i)))>::infinity())
        fun_sequence.evaluate( test_set(i), test_values(i,glas2::all()) ) ;
      //else
        //fill( test_values(i,glas2::all()), 0.0 ) ;
    }
    return SV_AAA( test_set, test_values, find_worst_index, options ) ;
  } // SV_AAA()

  template <typename FunctionSequence, typename Domain, typename FindWorstIndex>
  typename std::enable_if< !glas2::is< glas2::DenseVector, FunctionSequence>::value, SV_AAA_approximation<typename Domain::value_type, typename FunctionSequence::template value_type_for<typename Domain::value_type> > >::type SV_AAA( FunctionSequence const& fun_sequence, FindWorstIndex const& find_worst_index,  Domain const& domain ) {
    return SV_AAA( fun_sequence, find_worst_index, domain, std::tuple<>() ) ;
  } // SV_AAA()
  */

} } // namespace CORK::approximation

#endif
