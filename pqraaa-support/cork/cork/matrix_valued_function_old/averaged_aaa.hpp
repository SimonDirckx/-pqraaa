//  (C) Copyright Karl Meerbergen 2018.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_matrix_valued_function_averaged_aaa_hpp
#define cork_matrix_valued_function_averaged_aaa_hpp

#include <cork/matrix_valued_function/matrix_polynomial.hpp>
#include <cork/matrix_valued_function/difference.hpp>
#include <cork/approximation/aaa.hpp>
#include <cork/approximation/sv_aaa_approximation.hpp>
#include <cork/options/value_of.hpp>
#include <cork/options/test_approximation_error.hpp>
#include <cork/basis/barycentric_rational.hpp>
#include <cork/basis/union_of_bases.hpp>
#include <cork/coefficient_matrices/combined.hpp>

namespace CORK { namespace matrix_valued_function {
  
  template <typename T, typename MatrixValuedFunction, typename Domain>
  decltype (auto) averaged_AAA ( MatrixValuedFunction const& nep, Domain const& domain ) {
    return averaged_AAA<T>( nep, domain, std::tuple<>() ) ;
  } // averaged_AAA()

  template <typename T, typename MatrixValuedFunction, typename Domain, typename Options>
  decltype(auto) averaged_AAA( MatrixValuedFunction const& nep, Domain const& domain, Options const& options ) {
    typedef T value_type ;

    glas2::vector<value_type> u( nep.size() ) ;
    glas2::vector<value_type> v( nep.size() ) ;
    glas2::vector<value_type> z( nep.size() ) ;
    randomize(u) ; u /= norm_2(u) ;
    randomize(v) ; v /= norm_2(v) ;
    auto fun = [&]( auto s ) { fill(z,0.0) ; nep.multiply_add( s, u, z ) ; return inner_prod(v, z) ; } ;
    auto aaa_approx = approximation::AAA<T>( fun, domain, options ) ;

    // Make AAA_approximation of the entire approximation.
    approximation::SV_AAA_approximation<typename decltype(aaa_approx)::value_type, typename decltype(aaa_approx)::value_type> sv_aaa_approx( nep.basis().basis_2().num_terms(), aaa_approx.n() ) ;
    glas2::vector<value_type> temp( nep.basis().basis_2().num_terms() ) ;
    for (int i=0; i<aaa_approx.n(); ++i) {
      nep.basis().basis_2().evaluate( aaa_approx.nodes()(i), temp ) ;
      sv_aaa_approx.add_node( aaa_approx.weights()(i), aaa_approx.nodes()(i), temp ) ;
    }

    // Make state_space basis
    glas2::shared_vector<value_type> nodes( copy( sv_aaa_approx.nodes() ) ) ;
    glas2::shared_vector<value_type> weights( copy( sv_aaa_approx.weights() ) ) ;
    basis::barycentric_rational< decltype(weights), decltype(nodes) > state_space_basis( weights, nodes ) ;

    // Make combinations of coefficient matrices
    // The first (1,1) block is the identity matrix corresponding to the polynomial part
    // The combinations in the (2,2) block are equal to the function values of AAA multiplied with the weights, for each nonlinear functions, i.e., for each row of 'combinations'
    typedef glas2::shared_matrix< value_type > matrix_type ;
    matrix_type combinations( nep.coefficient_matrices().num_matrices(), nep.basis().basis_1().num_terms()+nodes.size() ) ;
    fill( combinations, 0.0 ) ;

    glas2::range basis_range(0,nep.basis().basis_1().num_terms()) ;
    // (1,1) block
    combinations( basis_range, basis_range ) = glas2::identity_matrix< value_type >( basis_range.size(), basis_range.size() ) ;

    // (2,2) block: the weights are taken into the linearization, not the matrices.
    combinations( glas2::range_from_end(nep.basis().basis_1().num_terms(),0), glas2::range_from_end(nep.basis().basis_1().num_terms(),0) ) = transpose( sv_aaa_approx.coefficients() ) ;
/*    for (int i=nep.basis().num_terms(); i<combinations.num_rows(); ++i) {
      combinations( i, glas2::range_from_end(nep.basis().num_terms(),0) ) = combinations( i, glas2::range_from_end(nep.basis().num_terms(),0) ) * aaa_approx.weights() ;
    }*/

    // Make the coefficient matrices for the AAA representation.
    CORK::coefficient_matrices::combined< typename MatrixValuedFunction::coefficient_matrices_type const&, matrix_type > combined( nep.coefficient_matrices(), combinations ) ;

    // Make Rational eigenvalue problem, using a union of bases
    auto uni_basis = basis::make_union_of_bases( nep.basis().basis_1(), state_space_basis ) ;
    matrix_polynomial< decltype(uni_basis), decltype(combined) > pep( uni_basis, combined ) ;

    if (options::value_of<options::test_approximation_error>(options)) {
       std::cout << "Degree of the matrix polynomial " << pep.basis().num_terms() << std::endl ;
       std::cout << "Error on the AAA approximation " << CORK::matrix_valued_function::difference( pep, nep, domain, options ) << std::endl ;
    }

    return pep ;
  } // averaged_AAA()

} } // namespace CORK::matrix_valued_function

#endif
