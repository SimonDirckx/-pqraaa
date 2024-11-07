//  (C) Copyright Karl Meerbergen 2018.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_approximation_sv_aaa_2_partial_fractions_hpp
#define cork_approximation_sv_aaa_2_partial_fractions_hpp

#include <cork/approximation/sv_aaa_approximation.hpp>
#include <cork/basis4cork/barycentric_rational_strong.hpp>
#include <cork/exception/lapack_error.hpp>
#include <cork/lapack/eig.hpp>
#include <glas2/vector.hpp>
#include <glas2/matrix.hpp>
#include <glas2/bindings/matrix.hpp>
#include <boost/numeric/bindings/lapack/driver/gels.hpp>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace CORK { namespace approximation {

  struct nofilter {
    template <typename T>
    bool operator() ( T const& t ) const {
      return false ;
    }
  } ;

  template <typename AAAApproximation, typename Filter=nofilter>
  decltype(auto) sv_aaa_2_partial_fractions( AAAApproximation const& aaa, Filter const& filter=nofilter(), bool make_coefs=false ) {
    typedef typename AAAApproximation::value_type            value_type ;
    typedef decltype(std::abs(value_type()))         real_value_type ;
    typedef std::complex< real_value_type >          complex_value_type ;

    glas2::range range_n( 0, aaa.n() ) ;

    basis::barycentric_rational_strong< decltype(aaa.weights()), decltype(aaa.nodes()) > basis( aaa.weights(), aaa.nodes() ) ;
    basis4CORK::basis4CORK< decltype(basis) > basis4cork(basis) ;

    // Set up linearization
    glas2::matrix< complex_value_type > M( aaa.n()-1, aaa.n() ) ;
    glas2::matrix< complex_value_type > N( M.num_rows(), M.num_columns() ) ;
    glas2::matrix< complex_value_type > N_copy( M.num_rows(), M.num_rows() ) ;
    //glas2::matrix< value_type > M_copy( M.num_rows(), M.num_columns() ) ;
    glas2::matrix< complex_value_type > X( M.num_rows(), M.num_rows() ) ;
    glas2::matrix< complex_value_type > Y( M.num_rows(), M.num_rows() ) ;

    basis4cork.fill_M( M ) ;
    basis4cork.fill_N( N ) ;

    // Solve eigenvalue problem
    glas2::vector<complex_value_type> e( M.num_rows() ) ;
    glas2::vector<complex_value_type> e_beta( e.size() ) ;
    N_copy = N(glas2::all(), glas2::range_from_end(1,0)) ;
    //M_copy = M ;
    auto Ne = N(glas2::all(), glas2::range_from_end(1,0)) ;
    auto Me = M(glas2::all(), glas2::range_from_end(1,0)) ;
    int info = boost::numeric::bindings::lapack::ggev( 'V', 'V', Me, Ne, e, e_beta, Y, X ) ;
    assert( info>= 0) ;
    if (info>0) throw std::runtime_error("QZ method failed for SV_AAA 2 PF") ;

    // Select wanted finite eigenvalues
    glas2::vector< int > selection_all( e.size() ) ;
    int n_pf = 0 ;
    for (int i=0; i<e.size(); ++i) {
      if (!filter(e(i)/e_beta(i))) {
        selection_all( n_pf ) = i ;
        ++n_pf ;
      }
    }
    auto selection( selection_all(glas2::range(0,n_pf)) ) ;

    SV_AAA_approximation< complex_value_type, complex_value_type > pf( aaa.coefficients().num_columns(), n_pf ) ;

    pf.n() = selection.size() ;
    pf.nodes() = e(selection) / e_beta(selection) ;
    fill( pf.weights(), 1.0 ) ;

    assert( pf.nodes().end()==std::unique( pf.nodes().begin(), pf.nodes().end() ) ) ;

    if (make_coefs) {
      glas2::matrix<complex_value_type> coef_mat( aaa.n(), n_pf+1 ) ;
      glas2::matrix<complex_value_type> coef_rhs( aaa.n(), aaa.coefficients().num_columns() ) ;
      coef_rhs = aaa.coefficients() ;
      for (int i=0; i<aaa.n(); ++i) {
        coef_mat(i,0) = 1. ;
        for (int j=0; j<n_pf; ++j) {
          coef_mat(i,j+1) = 1./(aaa.nodes()(i) - pf.nodes()(j)) ;
        }
      }
      info = boost::numeric::bindings::lapack::gels( coef_mat, coef_rhs ) ;
      assert( info>=0 ) ;
      if (info!=0) {
        throw exception::lapack_error( "LAPACK GELS: matrix does not have full rank" ) ;
      }
      //std::cout << "coefs " << coef_rhs( glas2::range(0,n_pf+1), glas2::all() ) << std::endl ;
      pf.coefficients() = coef_rhs( glas2::range(1,n_pf+1), glas2::all() ) ;
    }

    return pf ;
  } // sv_aaa_2_partial_fractions()

} } // namespace CORK::approximation

#endif
