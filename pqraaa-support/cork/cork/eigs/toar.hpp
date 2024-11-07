//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2015.
//
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_eigs_toar_hpp
#define cork_eigs_toar_hpp

#include <cork/eigs/detail/sort_ritz_values.hpp>
#include <cork/eigs/detail/implicit_restart.hpp>
#include <cork/eigs/info.hpp>
#include <cork/krylov/toar_process.hpp>
#include <cork/krylov/toar_triple.hpp>
#include <cork/options/value_of.hpp>
#include <cork/options/debug_level.hpp>
#include <cork/options/stop_criterion.hpp>
#include <cork/options/test_recurrence.hpp>
#include <cork/options/filter_infinite_eigenvalue.hpp>
#include <cork/options/max_solves.hpp>
#include <cork/lapack/schur.hpp>
#include <cork/lapack/order_schur.hpp>
#include <cork/lapack/eig.hpp>
#include <cork/lapack/qr.hpp>
#include <cork/utility/matlab.hpp>
#include <cork/utility/value_type_for.hpp>
#include <glas2/bindings/matrix.hpp>
#include <glas2/bindings/vector.hpp>
#include <limits>

namespace CORK { namespace eigs {

 //
 // Important: arguments are passed by reference.
 template <typename Linearization, typename Geometry, typename Shift, typename InitialVector, typename Triple, typename StopCriterion, typename Options>
 info toar( Linearization& linearization, Geometry const& eigenvalue_selector, Shift const& shift, InitialVector const& initial_vector, Triple& triple,
            StopCriterion const& stop_criterion, Options const& options ) {
   typedef typename value_type_for< typename Geometry::value_type, Linearization >::type value_type ;
   typedef decltype(std::abs(value_type())) real_value_type ;
   typedef std::complex<real_value_type>    complex_value_type ;

   int stagnate = 0 ;
   glas2::vector< complex_value_type >  theta_values( triple.k_max() ) ;                        // Ritz values of the shift-and-invert matrix,
                                                                                                // i.e., the eigenvalues of the Hessenberg matrix
   glas2::vector< complex_value_type >  lambda_values( triple.k_max() ) ;                       // Ritz values of the pencil, are a map of the theta values.
   glas2::matrix< complex_value_type >  eigen_vectors( triple.k_max(), triple.k_max() ) ;       // Eigenvectors of the Hessenberg matrix
   glas2::matrix< value_type >          schur_vectors( triple.k_max()+1, triple.k_max()+1 ) ;   // Schur vectors of the Hessenberg matrix
   glas2::vector< real_value_type >     eigen_resid( triple.k_max() ) ;                         // Norms of the residual of the Ritz values of the shift-and-invert matrix

   auto& information = static_cast<eigs::info&>( linearization.information() ) ;

   krylov::toar_process< value_type, Linearization, Triple, Options > process( linearization, triple, options ) ;

   // Set initial vector for Arnoldi
   process.initialize( shift ) ;
   process.initial_vector( initial_vector ) ;

   // Number of converged and wanted eigenvalues
   information.number_converged = 0 ;
   information.number_converged_and_wanted = 0 ;

   information.maximum_subspace_dimension = 0 ;
   information.minimum_subspace_dimension = triple.k_max() ;

   // Number of converged eigenvalues, included wanted and unwanted ones.
   // n_converged_unsorted >= information.number_converged

   bool converged = false ;
   for (; !converged && information.number_of_solves+triple.k_max()<options::value_of<options::max_solves>(options)+triple.k();) {
     // -------------------------------------
     // Expand the Krylov space using Arnoldi
     // -------------------------------------
     //
     // Expand from order triple.k() to order triple.k_max()
     //
     int max_k = std::min<int>( triple.k_max(), options::value_of<options::max_solves>(options)-information.number_of_solves+triple.k() ) ;
     process.expand( max_k ) ;
     real_value_type norm_H = process.norm_h() ;
     information.maximum_subspace_dimension = std::max<int>( information.maximum_subspace_dimension, triple.k() ) ;
     information.minimum_subspace_dimension = std::min<int>( information.maximum_subspace_dimension, triple.k() ) ;

     if (options::value_of<options::test_recurrence>(options)) process.verify_recurrence() ;

     if (options::value_of<options::filter_infinite_eigenvalue>(options)>0) {
       detail::implicit_qr( triple, options, schur_vectors ) ;
       if (options::value_of<options::test_recurrence>(options)) process.verify_recurrence() ;
     }

     // ---------------------------
     // Compute the Schur form of H
     // ---------------------------
     //
     //  beta is the right bottom corner value of H
     value_type beta = triple.H( triple.k(), triple.k()-1 ) ;

     glas2::range range_k(0,triple.k()) ;
     auto schur_matrix = triple.H( range_k, range_k ) ;
     auto schur_vec = schur_vectors( range_k, range_k ) ;
     auto theta_val = theta_values( range_k ) ;
     lapack::schur( schur_matrix, schur_vec, theta_val ) ;

     // -----------------------------------------------------------------
     // Sort the Ritz values and the Schur form
     // Let p = sort_ritz_values.keep()
     // Let c = information.number_converged
     //
     // 1) Push the p wanted Ritz values to the left upper corner of H
     // 2) Compute the residual norms of Ritz pairs and compute the
     //    number of converged Ritz pairs.
     // 3) Push the c converged ritz values to the left upper corner of H
     // 4) Lock the Schur form
     // 5) Apply all transformations to the Krylov vectors
     // -----------------------------------------------------------------
     //
     detail::sort_ritz_values< Geometry, decltype(options) > sort_ritz_values( eigenvalue_selector, options, information, stagnate, diagonal(schur_matrix,-1) ) ;

     // 1) Sort eigenvalues depending on the selection criterion, eigenvalue_selector
     auto lambda_val = lambda_values( range_k ) ;
     lambda_val = shift + glas2::ones<real_value_type>(theta_val.size()) / theta_val ;
     glas2::shared_vector<int> order( triple.k() ) ;

     sort_ritz_values.step_one( lambda_val, order, eigenvalue_selector.n_wanted_max() ) ;

     glas2::range range_keep( 0, sort_ritz_values.keep() ) ;
     auto schur_mat = triple.H( range_k, range_k ) ;
     {
       int schur_size ;
       int ierr = lapack::order_schur( schur_mat, schur_vec, theta_val, order(range_keep), schur_size ) ;
       if (ierr==1 && schur_size<sort_ritz_values.keep()) {
         information.warnings.push_back( "Schur ordering failed: less eigenvalues are retained" ) ;
         sort_ritz_values.keep( schur_size ) ;
       }
       if (schur_size==0) {
          information.error = "Failure in reordering of the Schur form using LAPACK: *TGSEN" ;
          information.number_converged_and_wanted = 0 ;
          information.number_converged = 0 ;
          return information ;
       }
     }

     auto schur_resid = triple.H( triple.k(), range_keep ) ;
     schur_resid = beta * schur_vec( triple.k()-1, range_keep ) ;

     // 2) Compute eigenvalues and eigenvectors of the Schur matrix and check convergence
     glas2::matrix< value_type > schur_mat2( sort_ritz_values.keep(), sort_ritz_values.keep() ) ;
     glas2::vector< real_value_type > resid( sort_ritz_values.keep() ) ;
     schur_mat2 = triple.H( range_keep, range_keep ) ;
     auto theta_val2 = theta_values( range_keep ) ;
     auto eigen_vec = eigen_vectors( range_keep, range_keep ) ;
     CORK::lapack::eig( schur_mat2, eigen_vec, theta_val2 ) ;
     resid = glas2::abs( multiply( transpose(eigen_vec), schur_resid ) ) ;

     // theta's are eigenvalues of shift-and-invert operator. Transform to get lambda's
     auto lambda_val2 = lambda_values( range_keep ) ;
     lambda_val2 = shift + glas2::ones<real_value_type>(theta_val2.size()) / theta_val2 ;

     auto order2 = order( range_keep ) ;
     converged = sort_ritz_values.step_two( triple, eigen_vec, schur_vec(glas2::all(), range_keep), lambda_val2, order2, resid, stop_criterion, norm_H ) ;
     //converged = sort_ritz_values.step_two( lambda_val2, order2, resid, stop_tolerance ) ;

     // 3) Reorder Schur form so that converged eigenvalues come first
     // Recompute the residual terms.
     glas2::matrix< value_type > permutation( sort_ritz_values.keep(), sort_ritz_values.keep() ) ;
     permutation = glas2::eye( sort_ritz_values.keep(), sort_ritz_values.keep() ) ;
     glas2::matrix< value_type > temp( triple.k(), sort_ritz_values.keep() ) ;

     auto schur_mat3 = schur_mat( range_keep, range_keep ) ;
     {
       int schur_size ;
       int ierr = lapack::order_schur( schur_mat3, permutation, theta_val2, order2(glas2::range(0,information.number_converged)), schur_size ) ;
       if (ierr==1) {
          information.error = "Failure in reordering of the Schur form using LAPACK: *TGSEN" ;
          information.number_converged_and_wanted = schur_size ;
          information.number_converged = schur_size ;
          temp = schur_vec( range_k, range_keep ) ;
          schur_vec( range_k, range_keep ) = multiply( temp, permutation ) ;
          triple.transform_vectors( schur_vectors( glas2::range(0,triple.k()+1), glas2::range(0,schur_size+1) ) ) ;
          return information ;
       }
     }

     // Apply ordering to the Schur vectors
     // (Use eigen_vectors for temporary storage)
     temp = schur_vec( range_k, range_keep ) ;
     schur_vec( range_k, range_keep ) = multiply( temp, permutation ) ;
     triple.H( sort_ritz_values.keep(), range_keep ) = beta * schur_vec( triple.k()-1, range_keep ) ;
     fill( schur_resid, 0.0 ) ; // Wipe out old residual because this can corrupt the recurrence relation.
     
     // 4) Lock converged vectors
     fill( triple.H( sort_ritz_values.keep(), glas2::range(0,information.number_converged) ), 0.0 ) ;
     //std::cout << triple.H( glas2::range(0,dim_keep+1), range_keep ) << std::endl;

     // 5) Reduction of the Krylov vectors
     fill( schur_vectors(triple.k(),glas2::all()), 0.0 ) ;
     fill( schur_vectors(glas2::all(),sort_ritz_values.keep()), 0.0 ) ;
     schur_vectors( triple.k(), sort_ritz_values.keep() ) = 1.0 ;
     triple.transform_vectors( schur_vectors( glas2::range(0,triple.k()+1), glas2::range(0, sort_ritz_values.keep()+1) ) ) ;

     ++information.number_of_restarts ;
     information.minimum_subspace_dimension = std::min<int>( information.minimum_subspace_dimension, triple.k() ) ;

     if (options::value_of<options::test_recurrence>(options)) process.verify_recurrence() ;
     if (options::value_of<options::debug_level>(options)>0) std::cout << "---------- New restart --------" << std::endl ;
   } // end for

   return information ;
 } // toar()

 // GLas arguments are passed by copy.
 template <typename Linearization, typename Geometry, typename Shift, typename Triple, typename EigenValues, typename EigenVectors>
 int toar_eigen( Linearization& linearization, Geometry const& eigenvalue_selector, Shift const& shift, Triple& triple, EigenValues ritz_values, EigenVectors eigen_vectors ) {
   typedef typename Triple::value_type      value_type ;
   typedef decltype(std::abs(value_type())) real_value_type ;
   typedef std::complex<real_value_type>    complex_value_type ;

   auto& information = static_cast<eigs::info&>( linearization.information() ) ;

   assert( information.number_converged>=ritz_values.size() ) ;
   assert( triple.Q.num_rows()==eigen_vectors.num_rows() ) ;
   assert( ritz_values.size()==eigen_vectors.num_columns() ) ;

   if (information.number_converged==0) return 0 ;

   glas2::vector< complex_value_type > eigen_values( information.number_converged ) ;
   glas2::matrix< value_type > schur_matrix_copy( information.number_converged, information.number_converged ) ;
   glas2::matrix< complex_value_type > eigen_vectors_small( information.number_converged, information.number_converged ) ;

   // Compute eigenvalues
   schur_matrix_copy = triple.H( glas2::range(0,information.number_converged), glas2::range(0,information.number_converged) ) ;
   CORK::lapack::eig( schur_matrix_copy, eigen_vectors_small, eigen_values ) ;

   eigen_values = shift + glas2::ones<real_value_type>(eigen_values.size()) / eigen_values ;

   // Sort eigenvalues and keep wanted ones.
   int n_wanted = detail::sort_ritz_values_final( eigenvalue_selector, eigen_values, eigen_vectors_small ) ;
   ritz_values = eigen_values( glas2::range(0,ritz_values.size() ) ) ;

   // Compute eigenvectors
   glas2::matrix< complex_value_type > combinations( triple.rank, information.number_converged ) ;
   combinations = multiply( triple.U( glas2::range(0,triple.rank), glas2::range(0,information.number_converged) ), eigen_vectors_small ) ;
   eigen_vectors = multiply( triple.Q(glas2::all(),glas2::range(0,triple.rank)), combinations(glas2::all(), glas2::range(0,eigen_vectors.num_columns())) ) ;
   for (int i=0; i<eigen_vectors.num_columns(); ++i) {
     eigen_vectors( glas2::all(), i ) /= norm_2( eigen_vectors( glas2::all(), i) ) ;
   }

   return n_wanted ;
 } // toar_eigen()
   
} } // namespace CORK::eigs


#endif
