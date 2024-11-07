//  (C) Copyright Karl Meerbergen 2023.
//
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_eigs_qz_invariant_pair_hpp
#define cork_eigs_qz_invariant_pair_hpp

#include <cork/lapack/schur.hpp>
#include <cork/lapack/order_schur.hpp>
#include <cork/exception/lapack_error.hpp>

namespace CORK { namespace eigs {

 template <typename NEP, typename QZResult>
 auto qz_invariant_pair( NEP const& nep, QZResult const& qz_result ) {
   auto const& qz_pair = std::get<0>(qz_result) ;
   auto const& eigenvalue_selector = std::get<2>(qz_result) ;
   auto const& information = std::get<eigs::info>(qz_result) ;

   typedef typename std::decay<decltype(qz_pair)>::type::value_type value_type ;
   typedef glas2::shared_matrix< value_type >                       matrix_type ;
   typedef decltype(std::abs(value_type()))                         real_value_type ;
   typedef std::complex<real_value_type>                            complex_value_type ;

   int n_wanted = -1 ;
   glas2::vector< complex_value_type > eigen_values( qz_pair.size() ) ;

   glas2::matrix< int > order( 1, eigen_values.size() ) ;
   glas2::matrix< complex_value_type > eigen_vectors_lin( eigen_values.size(), eigen_values.size() ) ;
   glas2::matrix< value_type > A_copy( eigen_values.size(), eigen_values.size() ) ; A_copy = qz_pair.A_ ;
   glas2::matrix< value_type > B_copy( eigen_values.size(), eigen_values.size() ) ; B_copy = qz_pair.B_ ;

   A_copy = qz_pair.A_ ; B_copy = qz_pair.B_ ;
   glas2::matrix< value_type > schur_vec_X( eigen_values.size(), eigen_values.size() ) ;
   glas2::matrix< value_type > schur_vec_Y( eigen_values.size(), eigen_values.size() ) ;

   int ierr = lapack::schur( A_copy, B_copy, schur_vec_X, schur_vec_Y, eigen_values ) ;
   assert (ierr>=0) ;
   if (ierr>0) {
      std::ostringstream err_str ;
      err_str << "CORK algorithm Schur decomposition failed (This may happen in the complex case due to a bug in LAPACK's QZ.\nCORK: LAPACK INFO = " << ierr ;
      throw exception::lapack_error(err_str.str()) ;
   }

   // Sort eigenvalues and keep wanted ones.
   order(0,glas2::all()) = glas2::range(0,order.num_columns() ) ;
   n_wanted = detail::sort_ritz_values_final( eigenvalue_selector, eigen_values, order ) ;

   int schur_size ;
   ierr = lapack::order_schur( A_copy, B_copy, schur_vec_X, schur_vec_Y, eigen_values, order(0,glas2::all()), schur_size ) ;
   if (ierr==1 && schur_size<n_wanted) {
     //information.warnings.push_back( "Schur ordering failed: less eigenvalues are retained" ) ;
   }
   n_wanted = std::min(schur_size,n_wanted) ;
   if (schur_size==0) {
      //information.error = "Failure in reordering of the Schur form using LAPACK: *TGSEN" ;
      n_wanted = information.number_converged ;
   }

   matrix_type X( nep.size(), n_wanted ) ;
   matrix_type S( n_wanted, n_wanted ) ;
   matrix_type T( n_wanted, n_wanted ) ;

   S = A_copy( glas2::range(0,n_wanted), glas2::range(0,n_wanted) ) ;
   T = B_copy( glas2::range(0,n_wanted), glas2::range(0,n_wanted) ) ;
   X = schur_vec_X( glas2::range(0,X.num_rows()), glas2::range(0,n_wanted) ) ;
   return std::make_tuple( X, S, T) ;
 } // qz_invariant_pair()


} } // namespace CORK::eigs


#endif
