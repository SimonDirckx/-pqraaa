//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2015.
//
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_eigs_cork_hpp
#define cork_eigs_cork_hpp

#include <cork/utility/value_type_for.hpp>
#include <cork/eigs/detail/sort_ritz_values.hpp>
#include <cork/eigs/detail/implicit_restart.hpp>
#include <cork/eigs/explicit_projection.hpp>
#include <cork/eigs/invariant_pair.hpp>
#include <cork/eigs/eigen_pairs.hpp>
//#include <cork/eigs/eigen_pairs.hpp>
#include <cork/krylov/cork_process.hpp>
#include <cork/krylov/cork_process_real.hpp>
#include <cork/krylov/cork_quadruple.hpp>
#include <cork/options/explicit_projection.hpp>
#include <cork/options/max_solves.hpp>
#include <cork/options/test_recurrence.hpp>
#include <cork/options/filter_infinite_eigenvalue.hpp>
#include <cork/options/value_of.hpp>
#include <cork/eigs/info.hpp>
#include <cork/lapack/schur.hpp>
#include <cork/lapack/order_schur.hpp>
#include <cork/lapack/eig.hpp>
#include <cork/utility/timer.hpp>
#include <glas2/bindings/matrix.hpp>
#include <glas2/bindings/vector.hpp>
#include <limits>

namespace CORK { namespace eigs {

 //
 // Important: arguments are passed by reference.
 template <typename Linearization, typename EigSelector, typename ShiftGenerator,typename InitialVector, typename Quadruple, typename StopCriterion, typename Options>
 info cork( Linearization& linearization, EigSelector const& eigenvalue_selector, ShiftGenerator & shift_generator, InitialVector const& initial_vector, Quadruple& quadruple,
            StopCriterion const& stop_criterion, Options const& options ) {

   typedef typename value_type_for< typename ShiftGenerator::value_type, Linearization >::type value_type ;
   typedef decltype(std::abs(value_type()))                                             real_value_type ;
   typedef std::complex<real_value_type>                                                complex_value_type ;

   eigs::info& information( static_cast<eigs::info&>( linearization.information() ) ) ;
   int stagnate = 0 ;
   glas2::shared_vector< complex_value_type >  lambda_values( quadruple.k_max() ) ;                            // Ritz values of the pencil
   glas2::shared_matrix< complex_value_type >  eigen_vectors( quadruple.k_max(), quadruple.k_max() ) ;         // Eigenvectors of the Hessenberg matrix
   glas2::shared_matrix< value_type >          schur_vectors_X( quadruple.k_max()+1, quadruple.k_max()+1 ) ;   // Schur vectors of the Hessenberg matrix
   glas2::shared_matrix< value_type >          schur_vectors_Y( quadruple.k_max()+1, quadruple.k_max()+1 ) ;   // Schur vectors of the Hessenberg matrix
   glas2::shared_vector< value_type >          beta_H( quadruple.k_max() ) ;                                   // Copy of last row of H
   glas2::shared_vector< value_type >          beta_K( quadruple.k_max() ) ;                                   // Copy of last row of K

   typedef typename std::conditional< std::is_same< std::complex<typename Quadruple::value_type>, value_type>::value
                                    , krylov::real_matrix_complex_shift_tag
                                    , krylov::general_tag
                                    >::type krylov_tag ;

   krylov::cork_process< krylov_tag, value_type
                       , Linearization, ShiftGenerator, Quadruple, Options
                       > process( linearization, shift_generator, quadruple, options ) ;

   CORK::timer timer ;
   timer.tic() ;

   // Set initial vector for Arnoldi
   process.initialize( shift_generator.shift(0) ) ;
   process.initial_vector( initial_vector ) ;

   // Number of converged and wanted eigenvalues
   information.number_converged = 0 ;
   information.number_converged_and_wanted = 0 ;
   information.number_of_solves = 0 ;

   information.maximum_subspace_dimension = 0 ;
   information.minimum_subspace_dimension = quadruple.k_max() ;

   // Number of converged eigenvalues, included wanted and unwanted ones.
   // information.number_converged >= information.number_converged_and_wanted
   bool restarted = false ;

   bool converged = false ;
   while ( !converged && information.number_of_solves+quadruple.k_max() < options::value_of<options::max_solves>(options)+quadruple.k() ) {
     /*if (restarted) {
     //std::cout << "new U " << quadruple.u_vector_k(0) << std::endl ;
     restarted = false ;
     }*/
     // ----------------------------------------------
     // Expand the Krylov space using rational Arnoldi
     // ----------------------------------------------
     //
     // Expand from order quadruple.k() to order quadruple.k_max()
     //
     int max_k = std::min<int>( quadruple.k_max(), options::value_of<options::max_solves>(options)-information.number_of_solves+quadruple.k() ) ;
     max_k = std::min<int>( max_k, linearization.size()*linearization.size_of_basis() ) ;
     /**/
     //if (restarted) --max_k ;
     restarted = false ;
     /**/
     process.expand( max_k ) ;
     real_value_type norm_H = process.norm_h() ;
     information.maximum_subspace_dimension = std::max<int>( information.maximum_subspace_dimension, quadruple.k() ) ;

     if (options::value_of<options::test_recurrence>(options)) process.verify_recurrence() ;

     if (options::value_of<options::filter_infinite_eigenvalue>(options)>0) {
       detail::implicit_qz( quadruple, options, schur_vectors_X ) ;
       if (options::value_of<options::test_recurrence>(options)) process.verify_recurrence() ;
     }

     // ----------------------------------------------------------
     // Compute the generalized Schur form of pencil K - \lambda*H
     // ----------------------------------------------------------

     glas2::range range_k(0,quadruple.k()) ;
     glas2::range range_l(0,information.number_converged) ;
     glas2::range range_k_l(information.number_converged,quadruple.k()) ;

     fill( schur_vectors_X, 0.0 ) ; glas2::fill(diagonal(schur_vectors_X),1.0) ;
     fill( schur_vectors_Y, 0.0 ) ; glas2::fill(diagonal(schur_vectors_Y),1.0) ;

     //
     // beta is the right bottom corner value of H
     auto beta_H_k = beta_H( range_k ) ;
     auto beta_K_k = beta_K( range_k ) ;
     beta_H_k = quadruple.H( quadruple.k(), range_k_l ) ;
     beta_K_k = quadruple.K( quadruple.k(), range_k_l ) ;

     auto lambda_val = lambda_values( range_k_l ) ;
     {
       auto schur_matrix_H = quadruple.H( range_k_l, range_k_l ) ;
       auto schur_matrix_K = quadruple.K( range_k_l, range_k_l ) ;
       auto schur_vec_X = schur_vectors_X( range_k_l, range_k_l ) ;
       auto schur_vec_Y = schur_vectors_Y( range_k_l, range_k_l ) ;
       lapack::schur( schur_matrix_K, schur_matrix_H, schur_vec_X, schur_vec_Y, lambda_val ) ;
       glas2::matrix<value_type> temp( range_l.size(), range_k_l.size() ) ;
       temp = quadruple.H( range_l, range_k_l ) ;
       quadruple.H( range_l, range_k_l ) = multiply( temp, schur_vec_X ) ;
       temp = quadruple.K( range_l, range_k_l ) ;
       quadruple.K( range_l, range_k_l ) = multiply( temp, schur_vec_X ) ;
     }

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
     detail::sort_ritz_values< EigSelector, decltype(options) >
        sort_ritz_values( eigenvalue_selector, options, information, stagnate
                        , glas2::abs_squared(diagonal(quadruple.K( range_k, range_k ),-1))+glas2::abs_squared(diagonal(quadruple.K( range_k, range_k ),-1))
                        ) ;

     glas2::shared_vector< int > order( quadruple.k() ) ;

     int eig_wanted = eigenvalue_selector.n_wanted_max();
     sort_ritz_values.step_one( lambda_val, order, eig_wanted ) ;
     int kept = sort_ritz_values.keep() ;
     int n_locked = information.number_converged ;

     // 0) Modify shifts if shift_generator supports this.
     int n_wanted = std::min( eig_wanted, sort_ritz_values.wanted_first( lambda_val, order ) ) ;
     shift_generator.update_shifts( lambda_val(order(glas2::range(0,n_wanted))) ) ;

     // 1) Sort eigenvalues depending on the selection criterion, eigenvalue_selector
     // First reordering
     // Locked and wanted Ritz values first. Locked but unwanted are pushed towards the right.
     // It may be that the ordering fails: in that case, do not purge unwanted eigenvalues.
     glas2::vector<int> order_unlocked( order.size() ) ;
     int keep_unlocked = 0 ;
     {
       glas2::vector<int> order_locked( std::min<int>(n_locked, order.size()) ) ;
       int keep_locked = 0 ;
       for (int i=0; i<sort_ritz_values.keep(); ++i) {
         if (n_locked>order(i)) {
           order_locked(keep_locked) = order(i) ;
           ++keep_locked ;
         } else {
           order_unlocked(keep_unlocked) = order(i) ;
           ++keep_unlocked ;
         }
       }
       auto schur_matrix_H_l = quadruple.H( range_l, range_l ) ;
       auto schur_matrix_K_l = quadruple.K( range_l, range_l ) ;
       auto schur_vec_X_l = schur_vectors_X( range_l, range_l ) ;
       auto schur_vec_Y_l = schur_vectors_Y( range_l, range_l ) ;
       auto lambda_val = lambda_values( range_l ) ;
       glas2::range range_keep( 0, keep_locked ) ;
       int schur_size ;
       int ierr = lapack::order_schur( schur_matrix_K_l, schur_matrix_H_l, schur_vec_X_l, schur_vec_Y_l, lambda_val, order_locked(range_keep), schur_size ) ;
       if (ierr==1 && schur_size<keep_locked) {
         information.warnings.push_back( "Schur ordering failed: all eigenvalues are retained" ) ;
       } else {
         n_locked = schur_size ;
       }
     }
     range_l = glas2::range(0,n_locked) ;
     range_k_l = glas2::range(n_locked,quadruple.k()) ;

     // Second reordering
     // Other Ritz values.
     auto schur_matrix_H_u = quadruple.H( range_k_l, range_k_l ) ;
     auto schur_matrix_K_u = quadruple.K( range_k_l, range_k_l ) ;
     auto schur_vec_X_u = schur_vectors_X( range_k_l, range_k_l ) ;
     auto schur_vec_Y_u = schur_vectors_Y( range_k_l, range_k_l ) ;
     auto lambda_val_u = lambda_values( range_k_l ) ;

     // It may be that the ordering fails: in that case, reduce the number of eigenvalues to keep
     {
       glas2::range range_keep( 0, keep_unlocked ) ;
       int schur_size ;
       int ierr = lapack::order_schur( schur_matrix_K_u, schur_matrix_H_u, schur_vec_X_u, schur_vec_Y_u, lambda_val_u, order_unlocked(range_keep), schur_size ) ;
       if (ierr==1 && schur_size<kept) {
         information.warnings.push_back( "Schur ordering failed: less eigenvalues are retained" ) ;
         sort_ritz_values.keep( n_locked+schur_size ) ;
       }
     }
     kept = sort_ritz_values.keep() ;
     glas2::range range_keep( 0, kept ) ;
     range_k_l = glas2::range(n_locked,kept) ;

     auto schur_resid = quadruple.H( quadruple.k(), range_keep ) ;
     //schur_resid = beta * schur_vec_X( quadruple.k()-1, range_keep ) ;
     quadruple.H( quadruple.k(), range_keep ) = multiply( transpose(schur_vectors_X(range_k, range_keep)), beta_H_k ) ;
     glas2::vector< real_value_type > resid( quadruple.k() ) ;

     // 2) Compute eigenvalues, eigenvectors of the unlocked Schur matrix and compute residual norms for unlocked Schur matrices
     //    that are kept in the subspace.
     {
       glas2::matrix< value_type > schur_mat_H2_l( range_k_l.size(), range_k_l.size() ) ;
       glas2::matrix< value_type > schur_mat_K2_l( range_k_l.size(), range_k_l.size() ) ;
       glas2::matrix< complex_value_type > eigen_vec_l( range_k_l.size(), range_k_l.size() ) ;
       auto lambda_val2_l( lambda_values( range_k_l ) ) ;
       CORK::lapack::eig( schur_mat_K2_l, schur_mat_H2_l, eigen_vec_l, lambda_val2_l ) ;
       resid(range_k_l) = glas2::abs( multiply( transpose(eigen_vec_l), schur_resid(range_k_l) ) ) ;
     }

     auto order2 = order( range_keep ) ;
     if (stop_criterion.eigenvector_needed()) {
       // Compute eigenvalues, eigenvectors of the entire Schur pair. This is needed for the convergence check
       // if this is based on Ritz pairs.

       glas2::matrix< value_type > schur_mat_H2( range_k.size(), range_k.size() ) ;
       glas2::matrix< value_type > schur_mat_K2( range_k.size(), range_k.size() ) ;
       glas2::matrix< complex_value_type > eigen_vec_k( range_k.size(), range_k.size() ) ;
       auto lambda_val_k( lambda_values( range_k ) ) ;
       CORK::lapack::eig( schur_mat_K2, schur_mat_H2, eigen_vec_k, lambda_val_k ) ;

       converged = sort_ritz_values.step_two( quadruple, eigen_vec_k, schur_vectors_Y(glas2::all(), range_k), lambda_val_k, order2, resid, stop_criterion, norm_H ) ;
       converged = 0 ;
     } else {
       auto order2_k_l = order( glas2::range( n_locked, kept ) ) ;
       order2( range_l ) = range_l ;
       glas2::matrix< complex_value_type > eigen_vec_k_l( range_k_l.size(), range_k_l.size() ) ;
       converged = n_locked + sort_ritz_values.step_two( quadruple, eigen_vec_k_l, schur_vectors_Y(glas2::all(), range_k_l), lambda_values( range_k_l )
                                                       , order2_k_l, resid(range_k_l), stop_criterion, norm_H
                                                       ) ;
       converged = 0 ;
       // Keep order of locked Schur.
       order2( range_l ) = range_l ;
     }

     // Reorder Schur form so that converged eigenvalues come first
     // Recompute the residual terms.
     glas2::shared_matrix< value_type > permutation_X( kept, kept ) ;
     permutation_X = glas2::eye( kept, kept ) ;

     glas2::shared_matrix< value_type > permutation_Y( kept, kept ) ;
     permutation_Y = glas2::eye( kept, kept ) ;

     // Reordering the Schur form of the unlocked part. Newly locked part is pushed to the left top corner.
     auto schur_mat_H3 = quadruple.H( range_keep, range_keep ) ;
     auto schur_mat_K3 = quadruple.K( range_keep, range_keep ) ;
     glas2::shared_matrix< value_type > temp( quadruple.k(), sort_ritz_values.keep() ) ;
     {
       auto lambda_val_k( lambda_values( range_keep ) ) ;
       int schur_size ;
       int ierr = lapack::order_schur( schur_mat_K3, schur_mat_H3, permutation_X, permutation_Y, lambda_val_k, order2(glas2::range(n_locked,converged)), schur_size ) ;
       if (ierr==1) {
          information.error = "Failure in reordering of the Schur form using LAPACK: *TGSEN" ;
          kept = schur_size ;
          information.number_converged_and_wanted = n_locked+schur_size ; // This is wrong
          information.number_converged = schur_size ;
          temp = schur_vectors_Y( range_k, range_keep ) ;
          schur_vectors_Y( glas2::range(0,quadruple.k()+1), range_keep ) = multiply( temp, permutation_Y ) ;
          fill( schur_vectors_Y( glas2::all(), kept ), 0.0 ) ; schur_vectors_Y( quadruple.k(), kept ) = 1.0 ;
          quadruple.transform_vectors( schur_vectors_Y( glas2::range(0,quadruple.k()+1), glas2::range(0,kept+1) ) ) ;
          return information ;
       }
     }

     // 3) Apply ordering to the Schur vectors
     // (Use eigen_vectors for temporary storage)
     temp = schur_vectors_X( range_k, range_keep ) ;
     schur_vectors_X( range_k, range_keep ) = multiply( temp, permutation_X ) ;

     temp = schur_vectors_Y( range_k, range_keep ) ;
     schur_vectors_Y( range_k, range_keep ) = multiply( temp, permutation_Y ) ;

     // Wipe out old residual because this can corrupt the recurrence relation.
     fill( schur_resid, 0.0 ) ;
     fill( quadruple.K( quadruple.k(), range_keep ), 0.0 ) ;

     //quadruple.H( kept, range_keep ) = beta * schur_vec_X( quadruple.k()-1, range_keep ) ;
     //quadruple.K( kept, range_keep ) = quadruple.poles(quadruple.k()-1) * quadruple.H( kept, range_keep ) ;
     quadruple.H( kept, range_keep ) = multiply( transpose(schur_vectors_X(range_k, range_keep)), beta_H_k ) ;
     quadruple.K( kept, range_keep ) = multiply( transpose(schur_vectors_X(range_k, range_keep)), beta_K_k ) ;
     
     // 4) Lock converged vectors
     int to_be_locked = 0 ;
     /*for (int i=0; i<sort_ritz_values.converged_unsorted(); ++i) {
       if (norm_2(quadruple.H( kept, glas2::range(0,i+1) ) ) < stop_tolerance) {
         to_be_locked = i+1 ;
       }
     }*/
     to_be_locked = information.number_converged ;
     fill( quadruple.H( kept, glas2::range(0,to_be_locked) ), 0.0 ) ;
     fill( quadruple.K( kept, glas2::range(0,to_be_locked) ), 0.0 ) ;

     fill( quadruple.H( glas2::range_from_end(kept+1,0), range_keep ), 0.0 ) ;
     fill( quadruple.K( glas2::range_from_end(kept+1,0), range_keep ), 0.0 ) ;

     // 5) Reduction of the Krylov vectors
     fill( schur_vectors_Y(quadruple.k(),glas2::all()), 0.0 ) ;
     fill( schur_vectors_Y(glas2::all(),kept), 0.0 ) ;
     schur_vectors_Y( quadruple.k(), kept ) = 1.0 ;
     quadruple.transform_vectors( schur_vectors_Y( glas2::range(0,quadruple.k()+1), glas2::range(0,kept+1) ) ) ;
     /**/
     //std::cout << "restart U " << quadruple.u_vector_k(0) << std::endl ;
     restarted = true ;
     /**/

     ++information.number_of_restarts ;
     information.minimum_subspace_dimension = std::min<int>( information.minimum_subspace_dimension, quadruple.k() ) ;

     if (options::value_of<options::test_recurrence>(options)) process.verify_recurrence() ;
     if (options::value_of<options::debug_level>(options)>0) std::cout << "---------- New restart --------" << std::endl ;
   } // end for


   information.total_time += timer.toc() ;

   return information ;
 } // cork()


/*
 // GLas arguments are passed by copy.
 template <typename Linearization, typename EigSelector, typename Quadruple, typename EigenValues, typename EigenVectors>
 void cork_eigen( Linearization& linearization, EigSelector const& eigenvalue_selector, Quadruple& quadruple, EigenValues ritz_values, EigenVectors eigen_vectors ) {
   typedef typename Quadruple::value_type     value_type ;
   typedef decltype(std::abs(value_type()))   real_value_type ;
   typedef std::complex<real_value_type>      complex_value_type ;

   auto& information = static_cast<eigs::info&>( linearization.information() ) ;

   assert( information.number_converged>=ritz_values.size() ) ;
   assert( quadruple.Q.num_rows()==eigen_vectors.num_rows() ) ;
   assert( ritz_values.size()==eigen_vectors.num_columns() ) ;

   if (information.number_converged==0) return ;

   glas2::vector< complex_value_type > eigen_values( information.number_converged ) ;
   glas2::shared_matrix< value_type > schur_matrix_copy_H( information.number_converged, information.number_converged ) ;
   glas2::shared_matrix< value_type > schur_matrix_copy_K( information.number_converged, information.number_converged ) ;
   glas2::shared_matrix< complex_value_type > eigen_vectors_small( information.number_converged, information.number_converged ) ;
   glas2::shared_matrix< complex_value_type > eigen_vectors_small_times_H( information.number_converged, information.number_converged ) ;
   glas2::shared_matrix< complex_value_type > combinations( quadruple.rank, information.number_converged ) ;
   //glas2::vector< complex_value_type > full_eigenvector( quadruple.degree()*quadruple.rank ) ;

   schur_matrix_copy_H = quadruple.H( glas2::range(0,information.number_converged), glas2::range(0,information.number_converged) ) ;
   schur_matrix_copy_K = quadruple.K( glas2::range(0,information.number_converged), glas2::range(0,information.number_converged) ) ;
   CORK::lapack::eig( schur_matrix_copy_K, schur_matrix_copy_H, eigen_vectors_small, eigen_values ) ;
   eigen_vectors_small_times_H = multiply( quadruple.H( glas2::range(0,information.number_converged), glas2::range(0,information.number_converged) )
                                         , eigen_vectors_small
                                         ) ;
   / *for (int i=0; i<information.number_converged; ++i) {
     full_eigenvector = multiply( quadruple.U( glas2::range(0,quadruple.degree()*quadruple.rank), glas2::range(0,information.number_converged) ), eigen_vectors_small_times_H( glas2::all(), i ) ) ;
     real_value_type nrm_u = norm_2(full_eigenvector( glas2::range(0,quadruple.rank)));
     int index = 0 ;
     for (int j=1; j<quadruple.degree()-1; ++i) {
       if (norm_2(full_eigenvector( glas2::range(j*quadruple.rank,(j+1)*quadruple.rank)))>nrm_u) {
         nrm_u = norm_2(full_eigenvector( glas2::range(j*quadruple.rank,(j+1)*quadruple.rank))) ;
         index = j ;
       }
       combinations(glas2::all(),i) = full_eigenvector(glas2::range(index*quadruple.rank,(index+1)*quadruple.rank)) ;
       std::cout << "selected " << index << std::endl ;
     }
   }* /
   combinations = multiply( quadruple.U( glas2::range(0,quadruple.rank), glas2::range(0,information.number_converged) ), eigen_vectors_small_times_H ) ;

   // Sort eigenvalues and keep wanted ones.
   detail::sort_ritz_values_final( eigenvalue_selector, eigen_values, combinations ) ;

   eigen_vectors = multiply( quadruple.Q(glas2::all(),glas2::range(0,quadruple.rank)), combinations(glas2::all(), glas2::range(0,eigen_vectors.num_columns()) ) ) ;
   for (int i=0; i<eigen_vectors.num_columns(); ++i) {
     eigen_vectors( glas2::all(), i ) /= norm_2( eigen_vectors( glas2::all(), i) ) ;
   }
   ritz_values = eigen_values( glas2::range(0,ritz_values.size() ) ) ;
 } // cork_eigen()
 */  
} } // namespace CORK::eigs


#endif
