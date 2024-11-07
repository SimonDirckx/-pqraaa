//  (C) Copyright Karl Meerbergen 2018.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_approximation_sv_aaa_hpp
#define cork_approximation_sv_aaa_hpp

#include <cork/approximation/sv_aaa_approximation.hpp>
#include <cork/approximation/unique_test_set.hpp>
#include <cork/approximation/sv_aaa_2_partial_fractions.hpp>
#include <cork/basis/barycentric_rational.hpp>
#include <cork/options/max_degree.hpp>
#include <cork/options/rational_approximation_tolerance.hpp>
#include <cork/options/debug_level.hpp>
#include <cork/options/domain_num_points.hpp>
#include <cork/options/value_of.hpp>
#include <cork/options/fast_aaa.hpp>
#include <cork/options/aaa_cheap_correction_tolerance.hpp>
#include <glas2/matrix.hpp>
#include <glas2/vector.hpp>
#include <glas2/bindings.hpp>
#include <boost/numeric/bindings/lapack/driver/gesvd.hpp>
#include <boost/numeric/bindings/lapack/computational/pocon.hpp>
#include <boost/numeric/bindings/lapack/computational/potrf.hpp>
#include <cassert>
#include <iomanip>

namespace CORK { namespace approximation {


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

    glas2::vector<int> ind_C(options::value_of<options::max_degree>(options)) ;
    glas2::matrix< value_type > C(test_set.size(), options::value_of<options::max_degree>(options));
    glas2::matrix< value_type > Q(test_set.size()*test_values.num_columns(), options::value_of<options::max_degree>(options));
    glas2::matrix< value_type > S(options::value_of<options::max_degree>(options), options::value_of<options::max_degree>(options)); fill(S,0.0) ;
    glas2::matrix< value_type > H(options::value_of<options::max_degree>(options), options::value_of<options::max_degree>(options)); fill(H,0.0) ;
    glas2::vector<value_type> D(test_set.size()) ;
    glas2::matrix<value_type> R(test_set.size(), test_values.num_columns() ) ;

    SV_AAA_approximation<value_type> repr( test_values.num_columns(), options::value_of<options::max_degree>(options) ) ;
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
            //throw std::runtime_error( ss.str() ) ;
            repr.warnings().push_back( std::string( ss.str() ) ) ;
        }
      }
    }

    // Compute the SV_AAA approximation
    repr.error() =  abs( vec(test_values)( max_ind( abs( vec(test_values) ) ) ) ) ;
    //fill( R, 0.0 ) ;
    for (typename FunctionValues::size_type j=0; j<test_values.num_columns(); ++j) {
      fill( R(glas2::all(),j), glas2::sum(test_values(glas2::all(),j)) / real_type(test_values.num_rows()) ) ;
    }

    for (; repr.n()<options::value_of<options::max_degree>(options) && repr.error()>options::value_of<options::rational_approximation_tolerance<real_type>>(options); ) {
      // Find index with largest residual
      /*real_type r_max =0.0;
      for (int i=0; i<test_values.num_rows(); ++i) {
        real_type temp_nrm = norm_2( test_values(i,glas2::all()) - R(i,glas2::all()) ) ;
        if (temp_nrm>r_max) {
          r_max = temp_nrm ;
          ind_C(d) = i ;
        }
      }*/
      int d = repr.n() ;
      int ind_Q = max_ind( abs( vec(test_values - R) ) ) ;
      //if (d==0) ind_Q = 79156;
      ind_C(d) = ind_Q % test_set.size() ;
      real_type r_max = std::abs( vec(test_values)(ind_Q) - vec(R)(ind_Q) ) ;
      real_type previous_error = repr.error() ;
      repr.error() = r_max ;
      if (options::value_of<options::debug_level>(options)>0) std::cout << "Degree " << d << " has approximation error " << repr.error() << std::endl ;
      if (r_max<options::value_of<options::rational_approximation_tolerance<real_type>>(options)) break ;

      // Check for Froissart doublets if error has increased
      if (d>0 && 0.*repr.error()>previous_error) {
        // Check for Froissart doublets
        auto pf = sv_aaa_2_partial_fractions( repr ) ;
        real_type norm_coef = norm_fro( pf.coefficients() ) ;
        std::vector< int > froissart_points ;
        for (int i=0; i<pf.n(); ++i) {
          if (norm_2(pf.coefficients()(i,glas2::all()))<options::value_of<options::rational_approximation_tolerance<real_type>>(options)*norm_coef) {
            froissart_points.push_back( glas2::max_ind( - abs_squared( pf.nodes()(i)-repr.nodes() ) ) ) ;
          }
        }
        if (froissart_points.size()>0) {
          std::sort( froissart_points.begin(), froissart_points.end() ) ;
          froissart_points.push_back(repr.n()) ; // sentinel element.
          int ii = froissart_points[0];
          for (typename std::vector< int >::size_type i=1; i<froissart_points.size(); ++i) {
            for ( ; ii+i<froissart_points[i]; ++ii ) {
              ind_C(ii) = ind_C(ii+i) ;
              C( glas2::all(), ii ) = C( glas2::all(), ii+i ) ;
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
          assert( d==repr.n() ) ;
          for (typename std::vector< int >::size_type i=0; i<froissart_points.size()-1; ++i) {
            C(ind_C(froissart_points[i]), glas2::range(0,d)) = glas2::ones<real_type>(repr.n()) / (test_set(froissart_points[i])-repr.nodes()) ;
          }
          for (size_type j=0; j<d; ++j) {
            //C(glas2::all(), j) = glas2::ones<real_type>(test_set.size()) / (test_set-repr.nodes()(j)) ;
            for (size_type i=0; i<test_values.num_columns(); ++i) {
              Q( glas2::range(i*test_set.size(),(i+1)*test_set.size()), j ) = C(glas2::all(),j) * test_values(glas2::all(),i) - C(glas2::all(),j) * test_values(ind_C(j),i) ;
            }
          }
          for (int i=0; i<test_values.num_columns(); ++i) {
            fill( Q(i*test_set.size()+ind_C(glas2::range(0,d)),glas2::range(0,d)), 0.0 ) ;
          }
          repr.warnings().push_back( std::string("Froissart points detected") ) ;
          if (options::value_of<options::debug_level>(options)>0) std::cout << "Froissart  points detected: " << froissart_points.size() << std::endl ;
        }
      } else {
        // Add interpolation point
        repr.add_node( 0.0, test_set(ind_C(d)), test_values(ind_C(d),glas2::all()) ) ;
        if (options::value_of<options::debug_level>(options)>1) std::cout << "New support point " << test_set(ind_C(d)) << std::endl ;
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
        //std::cout << "L " << L << std::endl;

        bool cheaply_correct = options::value_of<options::fast_aaa>(options) ;

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
      } // New point added

      glas2::vector<real_type>  svd_val(repr.n()) ;
      glas2::matrix<value_type> svd_u(1,1) ;
      glas2::matrix<value_type> svd_vh(repr.n(),repr.n()) ;
      svd_vh = H(glas2::range(0,repr.n()),glas2::range(0,repr.n())) ;
      //std::cout << svd_vh << std::endl ;
#ifndef NDEBUG
      int info =
#endif
        boost::numeric::bindings::lapack::gesvd( 'N', 'O', svd_vh, svd_val, svd_u, svd_u ) ;
      assert(info==0) ;
      repr.weights() = conj(svd_vh(repr.n()-1,glas2::all())) ;
      //std::cout << "nodes " << repr.nodes() << std::endl ;
      //std::cout << "weights " << repr.weights() << std::endl ;
      //std::cout << "coefficients " << repr.coefficients() << std::endl ;

      // Get the rational approximation
      glas2::range d_r1(0,repr.n()) ;
      D = multiply(C(glas2::all(),d_r1), repr.weights()) ;                     // Denominator
      for (int i=0; i<test_values.num_columns(); ++i) {
        R(glas2::all(),i) = multiply(C(glas2::all(),d_r1), repr.weights()*repr.coefficients()(glas2::all(),i)) ; // Numerator
        R(glas2::all(),i) /= D;
      }
      R(ind_C(d_r1), glas2::all()) = test_values(ind_C(d_r1), glas2::all());
    } // end for d

    // Remove Froissart doublets
    auto pf = sv_aaa_2_partial_fractions( repr ) ;
    real_type norm_coef = norm_fro( pf.coefficients() ) ;
    std::vector< int > froissart_points ;
    for (int i=0; i<pf.n(); ++i) {
      if (norm_2(pf.coefficients()(i,glas2::all()))<options::value_of<options::rational_approximation_tolerance<real_type>>(options)*norm_coef) {
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
          C( glas2::all(), ii ) = C( glas2::all(), ii+i ) ;
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
      assert( repr.n()==d) ; // -1 because of sentinel element
      for (typename std::vector< int >::size_type i=0; i<froissart_points.size()-1; ++i) {
        C(ind_C(froissart_points[i]), glas2::range(0,d)) = glas2::ones<real_type>(repr.n()) / (test_set(froissart_points[i])-repr.nodes()) ;
      }
      for (size_type j=0; j<d; ++j) {
        //C(glas2::all(), j) = glas2::ones<real_type>(test_set.size()) / (test_set-repr.nodes()(j)) ;
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
      repr.warnings().push_back( std::string("Froissart  points detected") ) ;
      if (options::value_of<options::debug_level>(options)>0) std::cout << "Froissart  points detected: " << froissart_points.size() << std::endl ;

      // Recompute error
      glas2::range d_r1(0,d) ;
      D = multiply(C(glas2::all(),d_r1), repr.weights()) ;                     // Denominator
      for (int i=0; i<test_values.num_columns(); ++i) {
        R(glas2::all(),i) = multiply(C(glas2::all(),d_r1), repr.weights()*repr.coefficients()(glas2::all(),i)) ; // Numerator
        R(glas2::all(),i) /= D;
      }
      R(ind_C(d_r1), glas2::all()) = test_values(ind_C(d_r1), glas2::all());
      int ind_Q = max_ind( abs( vec(test_values - R) ) ) ;
      real_type r_max = std::abs( vec(test_values)(ind_Q) - vec(R)(ind_Q) ) ;
      repr.error() = r_max ;
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
  } // SV_AAA()


  template <typename T, typename FunctionSequence, typename Domain>
  typename std::enable_if< !glas2::is< glas2::DenseVector, FunctionSequence>::value, SV_AAA_approximation<T> >::type SV_AAA( FunctionSequence const& fun_sequence, Domain const& domain ) {
    return SV_AAA<T>( fun_sequence, domain, std::tuple<>() ) ;
  } // SV_AAA()

  template <typename Points, typename FunctionValues>
  typename std::enable_if< glas2::is< glas2::DenseVector, Points>::value
                         && glas2::is< glas2::DenseMatrix, FunctionValues>::value
                         , SV_AAA_approximation<typename FunctionValues::value_type>
                         >::type SV_AAA( Points const& points, FunctionValues const& function_values ) {
    return SV_AAA( points, function_values, std::tuple<>() ) ;
  } // SV_AAA()

  template <typename T, typename FunctionSequence, typename Domain, typename Options>
  typename std::enable_if< !glas2::is< glas2::DenseVector, FunctionSequence>::value, SV_AAA_approximation<T> >::type SV_AAA( FunctionSequence const& fun_sequence, Domain const& domain, Options const& options ) {
    typedef T                                value_type ;

    // Construct the test set
    auto test_set_ini = domain.discretize( std::max(3,options::value_of<options::domain_num_points>(options)) ) ;
    auto test_set = unique_test_set( test_set_ini ) ;
    test_set_ini.resize( 0 ) ;

    // Evaluate and scale the functions
    glas2::shared_matrix<value_type> test_values( test_set.size(), fun_sequence.num_terms() ) ;
    for (typename glas2::vector<value_type>::size_type i=0; i<test_set.size(); ++i) {
      fun_sequence.evaluate( test_set(i), test_values(i,glas2::all()) ) ;
    }
    return SV_AAA( test_set, test_values, options ) ;
  } // SV_AAA()

} } // namespace CORK::approximation

#endif
