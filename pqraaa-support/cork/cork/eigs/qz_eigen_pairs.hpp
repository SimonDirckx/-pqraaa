//  (C) Copyright Karl Meerbergen 2023.
//
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_eigs_qz_eigen_pairs_hpp
#define cork_eigs_qz_eigen_pairs_hpp

#include <cork/matrix_valued_function/select_nep_from_pair.hpp>
#include <cork/eigs/eigen_residual_norms.hpp>
#include <cork/eigs/detail/sort_ritz_values.hpp>
#include <cork/lapack/eig.hpp>
#include <tuple>

namespace CORK { namespace eigs {

  template <typename Problem, typename QZResult>
  decltype(auto) qz_eigen_pairs( Problem const& problem, QZResult const& qz_result ) {
    typedef typename std::tuple_element<0,QZResult>::type  qz_type ;
    typedef typename qz_type::value_type                   value_type ;
    typedef decltype(std::abs(value_type()))               real_value_type ;
    typedef std::complex<real_value_type>                  complex_value_type ;

    // Get eigenvalues and eigenvectors of linearization stored in qz_pair.
    auto const& qz_pair = std::get<0>(qz_result) ;

    glas2::shared_vector< complex_value_type > eigen_values( qz_pair.size() ) ;
    glas2::matrix< complex_value_type > eigen_vectors_lin( qz_pair.size(), qz_pair.size() ) ;
    glas2::matrix< value_type > A_copy( qz_pair.size(), qz_pair.size() ) ; A_copy = qz_pair.A_ ;
    glas2::matrix< value_type > B_copy( qz_pair.size(), qz_pair.size() ) ; B_copy = qz_pair.B_ ;
    std::cout << "A = " << A_copy << std::endl ;
    std::cout << "B = " << B_copy << std::endl ;

    CORK::lapack::eig( A_copy, B_copy, eigen_vectors_lin, eigen_values ) ;
    std::cout << eigen_values << std::endl ;

    // Select first block for eigenevctors of NEP:
    glas2::shared_matrix< complex_value_type > eigen_vectors( problem.size(), qz_pair.size() ) ;
    eigen_vectors = eigen_vectors_lin( glas2::range(0,eigen_vectors.num_rows()), glas2::all() ) ;

    auto const& nep_eig = matrix_valued_function::select_nep_from_pair< complex_value_type, std::false_type >( CORK::deref(problem) ) ;

    glas2::shared_vector< real_value_type > residual_norms( eigen_values.size() ) ;
    eigen_residual_norms( nep_eig, eigen_values, eigen_vectors, residual_norms ) ;

    return std::make_tuple( eigen_values, eigen_vectors, residual_norms ) ;
 } // cork_eigen_pairs()


} } // namespace CORK::eigs


#endif
