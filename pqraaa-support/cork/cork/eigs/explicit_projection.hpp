#ifndef cork_3_eigs_explicit_projection_hpp
#define cork_3_eigs_explicit_projection_hpp

#include <cork/eigs/nep_residual_norm.hpp>
#include <cork/eigs/info.hpp>
#include <cork/eigs/invariant_pair_for_projection.hpp>
#include <cork/matrix_valued_function/nonlinear_matrix.hpp>
#include <cork/coefficient_matrices/any_glas.hpp>
#include <cork/coefficient_matrices/combined.hpp>
#include <cork/linearization/cork_linearization.hpp>
#include <cork/options/value_of.hpp>
#include <cork/options/keep_quadruple.hpp>
#include <cork/options/stop_criterion.hpp>
#include <cork/exception/lapack_error.hpp>
#include <cork/lapack/eig.hpp>
#include <cork/lapack/order_schur.hpp>
#include <cork/lapack/schur.hpp>
#include <glas2/matrix.hpp>
#include <glas2/vector.hpp>
#include <cmath>
#include <limits>
#include <memory>

namespace CORK { namespace eigs {

  template <typename MatrixPolynomial, typename Quadruple, typename Geometry, typename Options>
  typename std::enable_if< !coefficient_matrices::is_combined< typename MatrixPolynomial::coefficient_matrices_type >::value
                         , typename invariant_pair_for_projection_type<Quadruple>::type
                         >::type
  explicit_projection( MatrixPolynomial const& pep, Quadruple& quad, Geometry const& eigenvalue_selector, Options const& options, eigs::info& information ) {
    typedef typename std::common_type< typename Quadruple::value_type, typename MatrixPolynomial::coefficient_matrices_type::value_type>::type value_type ;
    typedef typename std::conditional< std::is_arithmetic<value_type>::value, std::complex<value_type>, value_type>:: type                     complex_value_type ;
    typedef decltype( std::abs(value_type()) )                                                                                                 real_value_type ;

    typedef glas2::shared_matrix<value_type> proj_matrix_type ;

    auto coefficient_matrices = pep.coefficient_matrices() ;

    // Estimate norms of coefficient matrices
    // We assume that the coefficient_matrices are not combined
    glas2::vector< real_value_type > a_norms( coefficient_matrices.num_matrices() ) ;
    glas2::vector< value_type > random_x( pep.size() ) ;
    randomize( random_x ) ; random_x /= norm_2( random_x ) ;
    glas2::vector< value_type > a_x( pep.size() ) ;
    for (int i=0; i<coefficient_matrices.num_matrices(); ++i) {
      fill( a_x, 0.0 ) ;
      coefficient_matrices.multiply_add( i, random_x, a_x ) ;
      a_norms(i) = norm_2( a_x ) ;
    }

    // Make coefficient matrices of projection
    // We assume that the coefficient_matrices are not combined
    glas2::matrix< value_type > AQ( pep.size(), quad.Q_k().num_columns() ) ;
    std::vector< proj_matrix_type > coefficient_matrices_proj_sequence ;
    for (int i=0; i<coefficient_matrices.num_matrices(); ++i) {
      fill(AQ, 0. ) ;
      for (int j=0; j<quad.Q_k().num_columns(); ++j)
        coefficient_matrices.multiply_add( i, quad.Q_k()(glas2::all(),j), AQ(glas2::all(),j) ) ;
      coefficient_matrices_proj_sequence.push_back( proj_matrix_type( quad.rank, quad.rank ) ) ;
      coefficient_matrices_proj_sequence.back() = multiply( transpose(quad.Q_k()), AQ ) ;
    }

    // Make matrix polynomial and linearization
    coefficient_matrices::any_glas< decltype(coefficient_matrices_proj_sequence) > coefs_proj( coefficient_matrices_proj_sequence ) ;
    matrix_valued_function::nonlinear_matrix pep_proj( pep.basis(), coefs_proj ) ;

    linearization::cork_linearization< decltype(pep_proj) > proj_lin( pep_proj, information ) ;

    // Compute Schur form
    auto fill_handle = proj_lin.fill_handle() ;
    glas2::matrix<value_type> AA( proj_lin.num_rows(), proj_lin.num_columns() ) ;
    glas2::matrix<value_type> BB( proj_lin.num_rows(), proj_lin.num_columns() ) ;
    fill_handle.A( AA ) ;
    fill_handle.B( BB ) ;
    glas2::matrix<value_type> Xf( proj_lin.num_rows(), proj_lin.num_columns() ) ;
    glas2::matrix<value_type> Yf( proj_lin.num_rows(), proj_lin.num_columns() ) ;
    glas2::vector<complex_value_type> e( proj_lin.num_rows() ) ;
    int ierr = lapack::schur( AA, BB, Xf, Yf, e ) ;
    if (ierr!=0) {
      std::stringstream s ;
      s << "zGGES failed with info = " << ierr ;
      throw exception::lapack_error( s.str() ) ;
    }

    glas2::matrix<value_type> AAc( copy(AA) ) ;
    glas2::matrix<value_type> BBc( copy(BB) ) ;
    glas2::matrix<complex_value_type> Z( e.size(), e.size() ) ;
    ierr = CORK::lapack::eig( AAc, BBc, Z, e ) ;
    if (ierr!=0) {
      std::stringstream s ;
      s << "zGGEV failed with info = " << ierr ;
      throw exception::lapack_error( s.str() ) ;
    }
    //std::cout << "e = " << e << std::endl ;

    assert( quad.Q_k().num_columns()*proj_lin.size_of_basis()==Z.num_rows() ) ;
    glas2::vector<complex_value_type> X_Z( quad.rank ) ;
    glas2::vector<complex_value_type> Q_X_Z( quad.Q.num_rows() ) ;
    glas2::vector<complex_value_type> temp1( quad.Q.num_rows() ) ;
    glas2::vector<real_value_type> resid( e.size() ) ;
    real_value_type norm_A ;

    // Keep eigenvalues with small residual norm.
    glas2::vector<int> keep( e.size() ) ;
    glas2::vector<value_type> coefs( pep.basis().num_terms() ) ;
    int i_keep = 0 ;
    for (int i=0; i<e.size(); ++i) {
      Z(glas2::all(), i) /= norm_2( Z(glas2::all(), i) ) ;
      X_Z = multiply( Xf(glas2::range(0, quad.rank), glas2::all()), Z(glas2::all(), i) ) ;
      Q_X_Z = multiply( quad.Q_k(), X_Z(glas2::range(0, quad.rank)) ) ;
      nep_residual_norm( pep, e(i), Q_X_Z, temp1, resid(i) ) ;

      pep.basis().evaluate( e(i), coefs ) ;
      norm_A = inner_prod( a_norms, glas2::abs(coefs) ) ;
      if (resid(i) < options::value_of<options::stop_criterion<real_value_type>>(options).relative_tolerance()*norm_A) {
        keep(i_keep) = i; ++i_keep ;
      }
    }
    //std::cout << "e = " << e << std::endl ;
    //std::cout << "resid = " << resid << std::endl ;
    //std::cout << "e(keep) = " << e( keep(glas2::range(0,i_keep)) ) << std::endl ;
    i_keep = detail::sort_ritz_values_final( eigenvalue_selector, e(keep(glas2::range(0,i_keep))), Z(glas2::all(), keep(glas2::range(0,i_keep))) ) ;
    std::cout << "e(keep) = " << e( keep(glas2::range(0,i_keep)) ) << std::endl ;
    std::cout << "resid(keep) = " << resid( keep(glas2::range(0,i_keep)) ) << std::endl ;

    // Copy information to Invariant Pair
    typedef typename invariant_pair_for_projection_type<Quadruple>::type inv_pair_type ;
    inv_pair_type inv_pair( make_invariant_pair_for_projection( quad, 0 ) ) ;
    if (i_keep<=quad.k_max() && options::value_of< options::keep_quadruple >( options ) ) {
      inv_pair.order( i_keep ) ;
    } else {
      typename inv_pair_type::x_type x( quad.Q.num_rows(), i_keep ) ;
      typename inv_pair_type::s_type s( i_keep, i_keep ) ;
      typename inv_pair_type::t_type t( i_keep, i_keep ) ;
      std::cout << "New inv pair" << std::endl ;
      inv_pair_type inv_pair_new( x, s, t, i_keep ) ;
      inv_pair = std::move( inv_pair_new ) ;
    }

    if (i_keep>0) {
      int schur_size ;
      ierr = lapack::order_schur( AA, BB, Xf, Yf, e, keep(glas2::range(0,i_keep)), schur_size ) ;
      if (ierr==1 && schur_size<i_keep) {
        information.warnings.push_back( "Schur ordering failed: less eigenvalues are retained" ) ;
      } else if (ierr!=0) {
        std::stringstream s ;
        s << "zTGSEN failed with info = " << ierr ;
        throw exception::lapack_error( s.str() ) ;
      }
      if (schur_size==0) {
         information.error = "Failure in reordering of the Schur form using LAPACK: *TGSEN" ;
         information.number_converged_and_wanted = 0 ;
         information.number_converged = 0 ;
      }

      glas2::range r_keep(0,i_keep) ;
      glas2::matrix< value_type > buf(8,i_keep) ;
      for (int i=0; i<quad.Q.num_rows(); i+=buf.num_rows()) {
        int n_rows = std::min(buf.num_rows(),quad.Q.num_rows()-i) ;
        buf( glas2::range(0,n_rows), glas2::all() ) = multiply( quad.Q_k()(glas2::range(i,i+n_rows), glas2::all()), Xf(glas2::range(0,quad.rank), r_keep) ) ;
        inv_pair.vectors()(glas2::range(i,i+n_rows), glas2::all()) = buf( glas2::range(0,n_rows), glas2::all() ) ;
      }

      inv_pair.matrix_S() = AA( r_keep, r_keep ) ;
      inv_pair.matrix_T() = BB( r_keep, r_keep ) ;
    }

    return std::move( inv_pair ) ;
  } // explicit_projection()

} } // namespace CORK::eigs

#endif
