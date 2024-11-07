//  (C) Copyright Karl Meerbergen 2022.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_matrix_valued_function_user_defined_hpp
#define cork_matrix_valued_function_user_defined_hpp

#include <cork/linear_solver/user_defined.hpp>
#include <cork/linear_solver/linear_solver.hpp>
#include <cork/coefficient_matrices/matrices_by_functions.hpp>
#include <cork/matrix_valued_function/linear_solver.hpp>
#include <cork/matrix_valued_function/nonlinear_matrix.hpp>
#include <cork/utility/ref.hpp>

namespace CORK { namespace matrix_valued_function {

  /*
  template <typename Basis, typename CoefficientMatrices>
  class user_defined
  : public nonlinear_matrix< Basis, CoefficientMatrices >
  {
    public:
      typedef typename CORK::deref_type< Basis >::type               basis_type ;
      typedef typename CORK::deref_type< CoefficientMatrices >::type coefficient_matrices_type ;

    public:
      user_defined( Basis const& basis, CoefficientMatrices const& matrices )
      : nonlinear_matrix< Basis, CoefficientMatrices >( basis, matrices )
      {}

    public:
      template <typename T>
      using linear_solver_type = CORK::linear_solver::user_defined< T, coefficient_matrices_type const* > ;

      template <typename T>
      auto linear_solver() const {
        return linear_solver_type<T>( &CORK::deref( this->matrices_ ) ) ;
      }
  } ; // class user_defined
*/

  template <typename ValueType, typename MatrixPolynomial, typename MatrixPolynomialType>
  struct linear_solver_struct< ValueType, MatrixPolynomial, MatrixPolynomialType, typename std::enable_if< coefficient_matrices::is_matrices_by_functions<typename MatrixPolynomialType::coefficient_matrices_type>::value >::type >  {
/*    typedef MatrixPolynomialType matrix_polynomial_type ;
    // Pass shared data as a reference.
    typedef typename matrix_polynomial_type::template value_type_for<ValueType> value_type ;
    typedef typename CORK::linear_solver::linear_solver_traits< value_type, std::reference_wrapper< typename matrix_polynomial_type::coefficient_matrices_type const> >::type solver_type ;
    typedef matrix_valued_function::linear_solver< std::reference_wrapper< typename matrix_polynomial_type::basis_type const>
                                                 , std::reference_wrapper< typename matrix_polynomial_type::coefficient_matrices_type const>
                                                 , solver_type
                                                 > type ;

    static type apply( MatrixPolynomial const& mp ) {
      typedef CORK::linear_solver::user_defined< T, coefficient_matrices_type const* > ;
      linear_solver_type solver( CORK::linear_solver::make_linear_solver<typename solver_type::value_type>( std::cref( CORK::deref(mp).coefficient_matrices() ) ) ) ;
      return type( std::cref(CORK::deref(mp).basis()), std::cref(CORK::deref(mp).coefficient_matrices()), solver ) ;
    }*/

    typedef typename MatrixPolynomialType::coefficient_matrices_type coefficient_matrices_type ;
    typedef CORK::linear_solver::user_defined< ValueType, typename MatrixPolynomialType::coefficient_matrices_data_type > type ;

    static type apply( MatrixPolynomial const& mp ) {
      return type( CORK::deref(mp).coefficient_matrices_data() ) ;
    }
  } ; // struct linear_solver_struct


} } // namespace CORK::matrix_valued_function

#endif
