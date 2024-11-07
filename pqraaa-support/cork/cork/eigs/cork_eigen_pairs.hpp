//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2015.
//
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_eigs_cork_eigen_pairs_hpp
#define cork_eigs_cork_eigen_pairs_hpp

#include <cork/matrix_valued_function/select_nep_from_pair.hpp>
#include <cork/eigs/eigen_residual_norms.hpp>
#include <cork/eigs/detail/sort_ritz_values.hpp>
#include <cork/lapack/eig.hpp>
#include <tuple>

namespace CORK { namespace eigs {


  template <typename CorkResult>
  decltype(auto) cork_eigen_pairs( CorkResult const& cork_result ) {
    typedef typename std::tuple_element<0,CorkResult>::type  quadruple_type ;
    typedef typename quadruple_type::value_type  value_type ;
    typedef decltype(std::abs(value_type()))     real_value_type ;
    typedef std::complex<real_value_type>        complex_value_type ;

    auto const& quadruple = std::get<0>(cork_result) ;
    auto const& eigenvalue_selector = std::get<2>(cork_result) ;
    auto const& information = std::get<3>(cork_result) ;



//    std::cout << "U1 in cork_eigen_pairs: " << quadruple.U(glas2::all(), 0) << std::endl;

    glas2::vector< complex_value_type > eigen_values( information.number_converged ) ;
    glas2::matrix< complex_value_type > combinations( quadruple.rank, information.number_converged ) ;

    int n_wanted = -1 ;

    if (information.number_converged>0) {
      glas2::matrix< value_type > schur_matrix_copy_H( information.number_converged, information.number_converged ) ;
      glas2::matrix< value_type > schur_matrix_copy_K( information.number_converged, information.number_converged ) ;
      glas2::matrix< complex_value_type > eigen_vectors_small( information.number_converged, information.number_converged ) ;
      glas2::matrix< complex_value_type > eigen_vectors_small_times_H( information.number_converged, information.number_converged ) ;
      //glas2::vector< complex_value_type > full_eigenvector( quadruple.degree()*quadruple.rank ) ;

      schur_matrix_copy_H = quadruple.H( glas2::range(0,information.number_converged), glas2::range(0,information.number_converged) ) ;
      schur_matrix_copy_K = quadruple.K( glas2::range(0,information.number_converged), glas2::range(0,information.number_converged) ) ;
      CORK::lapack::eig( schur_matrix_copy_K, schur_matrix_copy_H, eigen_vectors_small, eigen_values ) ;
 //     std::cout << "U2 in cork_eigen_pairs: " << quadruple.U(glas2::all(), 0) << std::endl;

      eigen_vectors_small_times_H = multiply( quadruple.H( glas2::range(0,information.number_converged), glas2::range(0,information.number_converged) )
                                            , eigen_vectors_small) ;
      combinations = multiply( quadruple.U( glas2::range(0,quadruple.rank), glas2::range(0,information.number_converged) )
                             , eigen_vectors_small_times_H ) ;

//      std::cout << "eigen_vectors_small_0: " << eigen_vectors_small( glas2::all(), 0) << std::endl;
//      std::cout << "eigen_vectors_small_times_H0: " << eigen_vectors_small_times_H( glas2::all(), 0) << std::endl;
//      std::cout << "U: " << quadruple.U( glas2::all(), 0) << std::endl;
//      std::cout << "combination0: " << combinations( glas2::all(), 0) << std::endl;
      // Sort eigenvalues and keep wanted ones.
      n_wanted = detail::sort_ritz_values_final( eigenvalue_selector, eigen_values, combinations ) ;
   }

   n_wanted = std::max(0,n_wanted) ;

   glas2::shared_matrix< complex_value_type > eigen_vectors( quadruple.Q.num_rows(), n_wanted ) ;
   // Hier worden de eigenvectoren berekend, dus op deze code baseren
   // eerst vermenigvuldigen met een u block (functie van cork_quadruple) (kleine matrix, kan op voorhand) en dan met Q (eigenvector per eigenvector)
   if (n_wanted>0) {
     eigen_vectors = multiply ( quadruple.Q(glas2::all(),glas2::range(0,quadruple.rank))
                            , combinations(glas2::all(), glas2::range(0,n_wanted) ) ) ;
     for (int i=0; i<eigen_vectors.num_columns(); ++i) {
       //if (i == 0) { std::cout << "eigen_vectors_small_times_H: " << eigen_vectors_small_times_H( glas2::all(), i) << std::endl; }
//       if (i == 0) { std::cout << "combination0: " << combinations( glas2::all(), i) << std::endl; }
//       if (i == 0) { std::cout << "eigvec0: " << eigen_vectors( glas2::all(), i) << std::endl; }
        eigen_vectors( glas2::all(), i ) /= norm_2( eigen_vectors( glas2::all(), i) ) ;
      }
    }

    glas2::shared_vector< complex_value_type > eigen_values_nwanted( n_wanted ) ;
    eigen_values_nwanted = eigen_values(glas2::range(0,n_wanted)) ;

    //return eigen_pairs<decltype(eigen_vectors), decltype(eigen_values_nwanted)>( eigen_vectors, eigen_values_nwanted ) ;
    return std::tuple( eigen_values_nwanted, eigen_vectors ) ;
  } // cork_eigen_pairs()


  template <typename Problem, typename CorkResult>
  decltype(auto) cork_eigen_pairs( Problem const& problem, CorkResult const& cork_result ) {
    auto [eigen_values, eigen_vectors] = cork_eigen_pairs( cork_result ) ;

    typedef typename std::tuple_element<0,CorkResult>::type  quadruple_type ;
    typedef typename quadruple_type::value_type  value_type ;
    typedef decltype(std::abs(value_type()))     real_value_type ;
    typedef std::complex<real_value_type>        complex_value_type ;

    auto const& nep_eig = matrix_valued_function::select_nep_from_pair< complex_value_type, std::false_type >( CORK::deref(problem) ) ;

    glas2::shared_vector< real_value_type > residual_norms( eigen_values.size() ) ;
    eigen_residual_norms( nep_eig, eigen_values, eigen_vectors, residual_norms ) ;

    return std::make_tuple( eigen_values, eigen_vectors, residual_norms ) ;
 } // cork_eigen_pairs()


} } // namespace CORK::eigs


#endif
