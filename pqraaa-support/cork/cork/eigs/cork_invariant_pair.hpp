//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2015.
//
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_eigs_cork_invariant_pair_hpp
#define cork_eigs_cork_invariant_pair_hpp

#include <cork/eigs/invariant_pair.hpp>

namespace CORK { namespace eigs {

 template <typename CorkResult>
 auto cork_invariant_pair( CorkResult const& cork_result ) {
   auto const& quadruple = std::get<0>(cork_result) ;
   auto const& eigenvalue_selector = std::get<2>(cork_result) ;
   auto const& information = std::get<3>(cork_result) ;

   typedef typename std::decay<decltype(quadruple)>::type::value_type value_type ;
   typedef glas2::shared_matrix< value_type > matrix_type ;
   typedef decltype(std::abs(value_type()))     real_value_type ;
   typedef std::complex<real_value_type>        complex_value_type ;

   int n_wanted = -1 ;
   glas2::vector< complex_value_type > eigen_values( information.number_converged ) ;

   glas2::matrix< int > order( 1, information.number_converged ) ;
   if (information.number_converged>0) {
      glas2::matrix< value_type > schur_matrix_copy_H( information.number_converged, information.number_converged ) ;
      glas2::matrix< value_type > schur_matrix_copy_K( information.number_converged, information.number_converged ) ;
      glas2::matrix< complex_value_type > eigen_vectors_small( information.number_converged, information.number_converged ) ;
      //glas2::vector< complex_value_type > full_eigenvector( quadruple.degree()*quadruple.rank ) ;

      schur_matrix_copy_H = quadruple.H( glas2::range(0,information.number_converged), glas2::range(0,information.number_converged) ) ;
      schur_matrix_copy_K = quadruple.K( glas2::range(0,information.number_converged), glas2::range(0,information.number_converged) ) ;
      CORK::lapack::eig( schur_matrix_copy_K, schur_matrix_copy_H, eigen_vectors_small, eigen_values ) ;

      // Sort eigenvalues and keep wanted ones.
      order(0,glas2::all()) = glas2::range(0,order.num_columns() ) ;
      n_wanted = detail::sort_ritz_values_final( eigenvalue_selector, eigen_values, order ) ;
   }
   glas2::matrix< value_type > schur_matrix_copy_H( information.number_converged, information.number_converged ) ;
   glas2::matrix< value_type > schur_matrix_copy_K( information.number_converged, information.number_converged ) ;
   schur_matrix_copy_H = quadruple.H( glas2::range(0,information.number_converged), glas2::range(0,information.number_converged) ) ;
   schur_matrix_copy_K = quadruple.K( glas2::range(0,information.number_converged), glas2::range(0,information.number_converged) ) ;
   glas2::matrix< value_type > schur_vec_X( information.number_converged, information.number_converged ) ;
   glas2::matrix< value_type > schur_vec_Y( information.number_converged, information.number_converged ) ;
   schur_vec_X = schur_vec_Y = glas2::identity_matrix<value_type>( information.number_converged, information.number_converged ) ;
   int schur_size ;
   int ierr = lapack::order_schur( schur_matrix_copy_K, schur_matrix_copy_H, schur_vec_X, schur_vec_Y, eigen_values, order(0,glas2::range(0,n_wanted)), schur_size ) ;
   if (ierr==1 && schur_size<n_wanted) {
     //information.warnings.push_back( "Schur ordering failed: less eigenvalues are retained" ) ;
   }
   n_wanted = std::min(schur_size,n_wanted) ;
   if (schur_size==0) {
      //information.error = "Failure in reordering of the Schur form using LAPACK: *TGSEN" ;
      n_wanted = information.number_converged ;
   }

   glas2::matrix<value_type> combinations( quadruple.rank, n_wanted ) ; 
   combinations = multiply( quadruple.u_block_k(0)(glas2::all(),glas2::range(0,information.number_converged)), schur_vec_Y(glas2::all(),glas2::range(0,n_wanted)) ) ;
  
   matrix_type X( quadruple.Q.num_rows(), n_wanted ) ;
   matrix_type S( n_wanted, n_wanted ) ;
   matrix_type T( n_wanted, n_wanted ) ;

   S = schur_matrix_copy_K( glas2::range(0,n_wanted), glas2::range(0,n_wanted) ) ;
   T = schur_matrix_copy_H( glas2::range(0,n_wanted), glas2::range(0,n_wanted) ) ;
   X = multiply( quadruple.Q_k(), combinations ) ;
   return std::make_tuple( X, S, T) ;
   //return invariant_pair< matrix_type, matrix_type, matrix_type >( std::move( X), std::move(S), std::move(T), information.number_converged ) ;
 } // cork_invariant_pair()


} } // namespace CORK::eigs


#endif
