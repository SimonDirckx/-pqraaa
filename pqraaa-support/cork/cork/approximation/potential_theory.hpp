//  (C) Copyright Karl Meerbergen 2019.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_approximation_potential_theory_hpp
#define cork_approximation_potential_theory_hpp

#include <cork/utility/matlab.hpp>
#include <cork/lapack/eig.hpp>
#include <cork/basis/rational_newton.hpp>
#include <cork/approximation/potential_theory_approximation.hpp>
#include <cork/options/value_of.hpp>
#include <cork/options/rational_approximation_tolerance.hpp>
#include <cork/exception/rational_approximation.hpp>
#include <glas2/matrix.hpp>
#include <glas2/vector.hpp>
#include <glas2/bindings/vector.hpp>
#include <glas2/bindings/matrix.hpp>
#include <boost/numeric/bindings/lapack/driver/gesv.hpp>
#include <cassert>
#include <cmath>
#include <limits>
#include <exception>
#include <sstream>

namespace CORK { namespace approximation {

  // Given: border, basis.poles
  template <typename Sigma, typename Xi, typename Basis>
  void potential_theory_leja_bagby( Sigma const& sigma, Xi const& xi, Basis& basis ) {
    typedef typename Basis::template value_type_for< typename Sigma::value_type >          value_type ;
    typedef decltype(std::abs(value_type()))    real_type ;

    glas2::vector< real_type > evals_sigma( sigma.size() ) ;
    glas2::vector< real_type > evals_xi( xi.size() ) ;

    // First node
    basis.nodes()(0) = sigma(0) ;
    evals_sigma = glas2::abs( sigma - basis.nodes()(0) ) ;
    evals_xi = glas2::abs( xi - basis.nodes()(0) ) ;
    real_type nrm = std::max( norm_inf( evals_sigma ), norm_inf( evals_xi ) ) ;
    if (nrm==0) throw exception::rational_approximation("NLEIGS: domain and singularity set overlap") ;
    evals_sigma /= nrm ;
    evals_xi /= nrm ;

    for (int i=1; i<basis.nodes().size(); ++i) {
      //std::cout << " i = " << i << " : " << evals_sigma(glas2::range(0,std::min<int>(10,evals_sigma.size()))) << " , " << evals_xi(glas2::range(0,std::min<int>(10,evals_sigma.size()))) << std::endl ;
      basis.nodes()(i) = sigma( glas2::max_ind(evals_sigma) ) ;
      basis.poles()(i-1) = xi( glas2::max_ind(-evals_xi) ) ;

//      evals_sigma /= evals_sigma(glas2::max_ind(evals_sigma)) ;
//      evals_xi /= evals_xi(glas2::max_ind(evals_xi)) ;

      evals_sigma = evals_sigma * glas2::abs( sigma - basis.nodes()(i) ) ;
      assert( prod( sigma - basis.poles()(i-1) )!=0.0);
      if (prod( glas2::abs( sigma - basis.poles()(i-1) ) )==0.) throw exception::rational_approximation("NLEIGS: domain and singularity set overlap") ;
      evals_sigma /= glas2::abs( sigma - basis.poles()(i-1) ) ;

      evals_xi = evals_xi * glas2::abs( xi - basis.nodes()(i) ) ;
      for (int j=0; j<evals_xi.size(); ++j) {
        real_type t = std::abs( xi(j) - basis.poles()(i-1) ) ;
        if (t==0.) evals_xi(j) = std::numeric_limits<real_type>::infinity() ;
        else evals_xi(j) /= t ;
      }
    }
    // Change the last pole to infinity.
    //basis.poles()(basis.poles().size()-1) = std::numeric_limits<real_type>::infinity() ;
  } // potential_theory_leja_bagby()


  // Given: border, basis.poles, basis.nodes
  template <typename Border, typename Basis>
  void potential_theory_scaling( Border const& border, Basis& basis ) {
    typedef typename Basis::template value_type_for< typename Border::value_type >       value_type ;
    typedef decltype(std::abs(value_type())) real_type ;

    glas2::vector< value_type > evals( border.size() ) ;

    // Degree 0:
    evals = glas2::ones<value_type>(evals.size()) ;
    if (basis.poles()(basis.grade()-1)!=std::numeric_limits<real_type>::infinity()) {
      assert( prod( glas2::abs(border - basis.poles()(basis.grade()-1)) )!=0.0 ) ;
      evals /= border-basis.poles()(basis.grade()-1) ;
    }
    basis.scaling()(0) = glas2::max( glas2::abs( evals ) ) ;
    evals /= basis.scaling()(0) ;

    // Scale the Newton polynomials
    for (typename Basis::size_type i=1; i<basis.num_terms(); ++i) {
      evals = evals * (border-basis.nodes()(i-1)) ;
      if (basis.poles()(i-1)!=std::numeric_limits<real_type>::infinity()) {
        assert( prod( glas2::abs(border - basis.poles()(i-1)) )!=0.0 ) ;
        evals /= basis.poles()(i-1)-border ;
      }
      basis.scaling()(i) = glas2::max( glas2::abs( evals ) ) ;
      evals /= basis.scaling()(i) ;
    }
  } // potential_theory_scaling()


  template <typename Basis, typename Function, typename Coefficients>
  void potential_theory_coefficients( Basis& basis, Function const& f, Coefficients& coef ) {
    assert( coef.size()==basis.nodes().size() ) ;

    typedef typename Coefficients::value_type       value_type ;
    typedef decltype(std::abs(value_type())) real_type ;
    glas2::matrix<value_type> H( basis.nodes().size(), basis.nodes().size() ) ;
    glas2::matrix<value_type> K( basis.nodes().size(), basis.nodes().size() ) ;
   // glas2::matrix<value_type> M( basis.nodes().size(), basis.nodes().size() ) ;
    fill(H,0.0) ;
    diagonal(H,-1) = basis.scaling()(glas2::range(1,basis.nodes().size())) ;
    diagonal(H,0) = basis.nodes();

    fill(K,0.0) ;
    fill( diagonal(K,0), 1.0 ) ;
    diagonal(K,-1) = basis.scaling()(glas2::range(1,basis.nodes().size())) * basis.poles()( glas2::range(0,basis.nodes().size()-1) ) ;
  //  std::cout << "res " << multiply(transpose(H),coef)-2.*multiply(transpose(K),coef) << std::endl;

    /*fill(M,0.0) ;
    diagonal(M,0) = basis.nodes();
    diagonal(M,-1) = basis.scaling()(glas2::range_from_end(1,0)) * (basis.poles()-basis.nodes()(glas2::range_from_end(0,1))) ;
    diagonal(M,-2) = -basis.scaling()(glas2::range_from_end(1,1)) * basis.scaling()(glas2::range_from_end(2,0)) * basis.poles()(glas2::range_from_end(1,0)) ;
    std::cout << "res " << multiply(transpose(M),coef)-2.*coef << std::endl;*/

    // Compute the matrix function times vector
    glas2::vector< value_type > rhs( H.num_rows() ) ;
    fill(rhs,0.0) ; rhs(0) = basis.scaling()(0) ;
    std::cout << CORK::matlab( H, "H" ) << std::endl ;
    std::cout << CORK::matlab( K, "K" ) << std::endl ;

    // Compute matrix function from eigenvalues.
    typedef std::complex<real_type> complex_type ;
    glas2::matrix<complex_type> X(H.num_rows(),H.num_columns());
    glas2::matrix<complex_type> Y(H.num_rows(),H.num_columns());
    glas2::vector<complex_type> e(H.num_columns()); // H = Y^-* a X^{-1}, K = Y^-* b X^{-1}, K^-1 = X b^-1 Y', H K^-1 = Y^-* e Y'
    //int info = lapack::eig( M, X, e ) ; // M = X e X^{-1} , f(M) = X f(e) X^{-1}
    int info = lapack::eig( H, K, X, Y, e ) ; // M = X e X^{-1} , f(M) = X f(e) X^{-1}
    assert(info==0) ;
    std::cout << "eig " << e << std::endl ;

    coef = multiply(transpose(conj(Y)),rhs) ;
    for (int i=0; i<e.size(); ++i) {
      coef(i) *= f( e(i) ) ;
    }

    glas2::vector<int> pivots(rhs.size());
    X = transpose(conj(Y)) ;
    info = boost::numeric::bindings::lapack::gesv( X, pivots, coef ) ;
    assert(info==0) ;
  } // potential_theory_coefficients()



  template <typename T, typename FunctionSequence, typename Basis, typename Border, typename Options>
  decltype (auto) potential_theory( FunctionSequence const& f, Basis const& basis, Border const& border, Options const& options ) {
    typedef T       value_type ;
    typedef decltype(std::abs(value_type())) real_type ;

    glas2::matrix<value_type> coefs( basis.grade()+1, f.num_terms() ) ;
    glas2::vector<value_type> values( f.num_terms() ) ;

    //
    // Compute divided difference from a matrix function on H * inv(K)
    //
    // Here, we use the eigendecomposition
    //
    //  Y' H X = alpha, Y' K X = beta and, so,  Y' H inv(K) inv(Y') = e with e eigenvalues of H inv(K)
    //
    //  f(Hinv(K)) = inv(Y') f(e) Y'
    //
    glas2::matrix<value_type> H( basis.nodes().size(), basis.nodes().size() ) ;
    glas2::matrix<value_type> K( basis.nodes().size(), basis.nodes().size() ) ;
    fill(H,0.0) ;
    fill(K,0.0) ;
    fill( diagonal(K,0), 1.0 ) ;
    diagonal(K,-1) = basis.scaling()(glas2::range_from_end(1,0)) ;
    diagonal(H,0) = basis.nodes();
    diagonal(H,-1) = basis.scaling()(glas2::range_from_end(1,0)) * basis.poles() ;
  //  std::cout << "res " << multiply(transpose(H),coef)-2.*multiply(transpose(K),coef) << std::endl;

    // Compute matrix function from eigenvalues.
    typedef std::complex<real_type> complex_type ;
    glas2::matrix<complex_type> X(H.num_rows(),H.num_columns());
    glas2::matrix<complex_type> Y(H.num_rows(),H.num_columns());
    glas2::vector<complex_type> e(H.num_columns()); // H = Y^-* a X^{-1}, K = Y^-* b X^{-1}, K^-1 = X b^-1 Y', H K^-1 = Y^-* e Y'
    //int info = lapack::eig( M, X, e ) ; // M = X e X^{-1} , f(M) = X f(e) X^{-1}
    int info = lapack::eig( H, K, X, Y, e ) ; // M = X e X^{-1} , f(M) = X f(e) X^{-1}
    assert(info==0) ;

    for (int i=0; i<e.size(); ++i) {
      f.evaluate( e(i), values ) ;
      coefs(i ,glas2::all()) = values * basis.scaling()(0) * std::conj(Y(0,i));
    }

    glas2::vector<int> pivots(coefs.num_rows());
    X = transpose(conj(Y)) ;
    info = boost::numeric::bindings::lapack::gesv( X, pivots, coefs ) ;
    assert(info>=0) ;
    std::stringstream ss ;
    ss << "lapack::GESV failed in CORK::approximation::potential_theory().  Lapack info = " << info ;
    if (info>0) throw std::runtime_error(ss.str()) ;

    //
    // Truncate the polynomial following relative tolerance
    //
    assert( options::value_of<options::rational_approximation_tolerance<real_type>>(options)<1. ) ;
    real_type error_level = 1.0 ;
    // Truncate approximation
    glas2::vector<real_type> norms( coefs.num_columns() ) ;
    fill(norms, 0.0) ;
    int grade = coefs.num_rows()-1 ;
    for (int i=0; i<coefs.num_rows(); ++i) {
      error_level = glas2::norm_inf(coefs(i,glas2::all())/norms ) ;
      if (error_level<=options::value_of<options::rational_approximation_tolerance<real_type>>(options)) {
        grade = i ;
        break ;
      }
      for (int j=0; j<norms.size(); ++j) norms(j) = std::max( norms(j), std::abs(coefs(i,j))) ;
    }

    std::cout << "Potential theory: degree of approximation is " << grade << " with error level " << error_level << std::endl ;

    // Fix degree of polynomial:
    glas2::range r_grade(0,grade) ;
    glas2::range r1_grade(0,grade+1) ;

    //
    // Build a new approximation with the last pole at infinity.
    //

    // Set last pole at infty.

    // Build data structure for the approximation
    potential_theory_approximation< value_type > repr( grade, f.num_terms() ) ;
    repr.nodes() = basis.nodes()(r1_grade) ;
    repr.poles() = basis.poles()(r_grade) ;
    repr.poles()(grade-1) = std::numeric_limits<real_type>::infinity() ;
    repr.scaling() = basis.scaling()(r1_grade) ;
    potential_theory_scaling( border, repr.basis() ) ;

    auto Hr = H(r1_grade,r1_grade) ;
    auto Kr = K(r1_grade,r1_grade) ;
    auto er = e(r1_grade) ;
    auto pivotsr = pivots(r1_grade) ;
    auto Xr = X(r1_grade,r1_grade) ;
    auto Yr = Y(r1_grade,r1_grade) ;
    auto coefs_r = coefs(r1_grade,glas2::all()) ;

    fill(Hr,0.0) ;
    fill(Kr,0.0) ;
    fill( diagonal(Kr,0), 1.0 ) ;
    diagonal(Kr,-1) = repr.scaling()(glas2::range(1,grade+1)) ;
    diagonal(Hr,0) = repr.nodes();
    diagonal(Hr,-1) = repr.scaling()(glas2::range(1,grade+1)) * repr.poles() ;
    // Correction for pole at inf.
    K(grade,grade-1) = 0.0 ;
    H(grade,grade-1) = repr.scaling()(grade) ;

    info = lapack::eig( Hr, Kr, Xr, Yr, er ) ; // M = X e X^{-1} , f(M) = X f(e) X^{-1}
    assert(info==0) ;

    for (int i=0; i<er.size(); ++i) {
      f.evaluate( er(i), values ) ;
      coefs(i ,glas2::all()) = values /* repr.scaling()(0)*/ * std::conj(Yr(0,i)); // scaling(0) == 1.0
    }

    Xr = transpose(conj(Yr)) ;
    info = boost::numeric::bindings::lapack::gesv( Xr, pivotsr, coefs_r ) ;
    assert(info==0) ;

    //
    // Build data structure for the approximation
    //
    repr.coefficients() = coefs(r1_grade, glas2::all()) ;
    return repr ;
  } // potential_theory()


} } // namespace CORK::approximation

#endif
