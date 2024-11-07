//  (C) Copyright Karl Meerbergen 2022.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_matrix_valued_function_linear_solver_hpp
#define cork_matrix_valued_function_linear_solver_hpp

#include <cork/options/value_of.hpp>
#include <cork/options/relative_shift_modifier.hpp>
#include <cork/options/absolute_shift_modifier.hpp>
#include <cork/exception/linear_solver_failure.hpp>
#include <cork/linear_solver/linear_solver.hpp>
#include <cork/basis4cork/basis4cork.hpp>
#include <cork/utility/ref.hpp>
#include <glas2/matrix.hpp>
#include <glas2/vector.hpp>
#include <type_traits>
#include <cassert>

namespace CORK { namespace matrix_valued_function {

  template <typename Basis, typename CoefficientMatrices, typename CoefficientMatricesLinearSolver>
  class linear_solver
  {
    public:
      typedef typename CORK::deref_type< CoefficientMatricesLinearSolver >::type linear_solver_type ;
      typedef typename CORK::deref_type< CoefficientMatrices >::type             coefficient_matrices_type ;
      typedef typename CORK::deref_type< Basis >::type                           sequence_of_functions_type ;
      typedef typename linear_solver_type::size_type                             size_type ;

    public:
      typedef typename linear_solver_type::value_type value_type ;

    public:
      explicit linear_solver( Basis const& basis, CoefficientMatrices const& coefficient_matrices, CoefficientMatricesLinearSolver solver )
      : basis_( basis )
      , coefficient_matrices_( coefficient_matrices )
      , solver_( solver )
      {}

    private:
      size_type num_terms() const { return CORK::deref(basis_).num_terms() ; }

    public:
      size_type size() const { return CORK::deref(solver_).size() ; }

    public:
      // Old code, deprecated
      void prepare_solve( value_type const& shift ) {
        glas2::vector< value_type > coefs( num_terms() ) ;
        CORK::deref(basis_).evaluate( shift, coefs.pass_ref() ) ;
        CORK::deref(solver_).prepare_solve( [&]( auto& matrix ) { CORK::deref(coefficient_matrices_).accumulate( coefs, matrix ) ; } ) ;
      } // prepare_solve()

      void prepare_solve( basis4CORK::basis4CORK<typename std::decay<Basis>::type> const& basis4CORK ) {
        glas2::vector< value_type > coefs( num_terms() ) ;
        basis4CORK.evaluate( coefs.pass_ref() ) ;
        CORK::deref(solver_).prepare_solve( coefs ) ;
      } // prepare_solve()

      template <typename V>
      void solve( V v ) const {
        static_assert( std::is_same< typename std::decay<V>::type::value_type, value_type >::value, "V must have the same value_type as solver" ) ;
        CORK::deref(solver_).solve( v ) ;
      } // solve()

    public:
      // For 2020 code
      template <typename V, typename Options>
      void solve( value_type& shift, V v, bool is_new_shift, Options const& options ) {
        try {
          prepare_solve( shift ) ;
        } catch (exception::linear_solver_failure& e ) {
          // Set new shift and retry
          value_type new_shift = shift + options::value_of<options::relative_shift_modifier<value_type>>(options) * shift + options::value_of<options::absolute_shift_modifier<value_type>>(options) ;
          shift = new_shift ;
          prepare_solve( new_shift ) ;
        }

        solve( v ) ;
      }

    private:
      Basis                           basis_ ;
      CoefficientMatrices             coefficient_matrices_ ;
      CoefficientMatricesLinearSolver solver_ ;
  } ; // linear_solver

/*
  template <typename Basis, typename CoefficientMatricesLinearSolver>
  decltype (auto) make_linear_solver( Basis basis, CoefficientMatricesLinearSolver solver ) {
    return linear_solver<Basis,CoefficientMatricesLinearSolver>( basis, solver ) ;
  }*/

  template <typename ValueType, typename MatrixPolynomial, typename MatrixPolynomialType, typename EnableIf=void>
  struct linear_solver_struct {
    typedef MatrixPolynomialType matrix_polynomial_type ;
    // Pass shared data as a reference.
    typedef typename matrix_polynomial_type::template value_type_for<ValueType> value_type ;
    typedef typename CORK::linear_solver::linear_solver_traits< value_type, typename matrix_polynomial_type::coefficient_matrices_type >::type solver_type ;
    //typedef matrix_valued_function::linear_solver< std::reference_wrapper< typename matrix_polynomial_type::basis_type const>
    //                                             , std::reference_wrapper< typename matrix_polynomial_type::coefficient_matrices_type const>
    //                                             , solver_type
    //                                             > type ;
    typedef matrix_valued_function::linear_solver< typename matrix_polynomial_type::basis_data_type
                                                 , typename matrix_polynomial_type::coefficient_matrices_data_type
                                                 , solver_type
                                                 > type ;

    static type apply( MatrixPolynomial const& mp ) {
      //std::cout << "Do not select this one" << std::endl ;
      //std::cout << typeid(CoefficientMatrices).name() << std::endl ;
      solver_type solver( CORK::linear_solver::make_linear_solver<typename solver_type::value_type>( CORK::deref(mp).coefficient_matrices() ) ) ;
      return type( CORK::deref(mp).basis_data(), CORK::deref(mp).coefficient_matrices_data(), solver ) ;
    }
  } ; // struct linear_solver_struct

  template <typename ValueType, typename MatrixPolynomial>
  struct linear_solver_type
  : linear_solver_struct< ValueType,MatrixPolynomial,typename CORK::deref_type<MatrixPolynomial>::type>
  { } ;

  template <typename ValueType, typename MatrixPolynomial>
  decltype (auto) make_linear_solver( MatrixPolynomial mp ) {
    return linear_solver_struct<ValueType,MatrixPolynomial,typename CORK::deref_type<MatrixPolynomial>::type>::apply( mp ) ;
  }

} } // namespace CORK::matrix_valued_function

#endif
