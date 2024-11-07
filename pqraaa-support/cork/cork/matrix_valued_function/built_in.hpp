//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_matrix_valued_function_built_in_hpp
#define cork_matrix_valued_function_built_in_hpp

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
#include <typeinfo>

namespace CORK { namespace matrix_valued_function {

  template <typename Basis, typename CoefficientMatrices, typename CoefficientMatricesLinearSolver>
  class built_in
  {
    public:
      typedef typename CORK::deref_type< CoefficientMatricesLinearSolver >::type linear_solver_type ;
      typedef typename CORK::deref_type< CoefficientMatrices >::type             coefficient_matrices_type ;
      typedef typename CORK::deref_type< Basis >::type                           sequence_of_functions_type ;
      typedef typename linear_solver_type::size_type                             size_type ;

    public:
      typedef typename linear_solver_type::value_type value_type ;

    public:
      explicit built_in( Basis const& basis, CoefficientMatrices const& coefficient_matrices )
      : basis_( basis )
      , coefficient_matrices_( coefficient_matrices )
      , solver_( CORK::linear_solver::make_linear_solver<typename linear_solver_type::value_type>( std::cref( coefficient_matrices_ ) ) )
      {}

    private:
      size_type num_terms() const { return CORK::deref(basis_).num_terms() ; }

    public:
      size_type size() const {
        assert( CORK::deref(coefficient_matrices_).num_rows()==CORK::deref(coefficient_matrices_).num_columns() ) ;
        return CORK::deref(coefficient_matrices_).num_rows() ;
      }

    private:
      void prepare_solve( value_type const& shift ) {
        glas2::vector< value_type > coefs( num_terms() ) ;
        CORK::deref(basis_).evaluate( shift, coefs.pass_ref() ) ;
        solver_.prepare_solve( [&]( auto& matrix ) { CORK::deref(coefficient_matrices_).accumulate( coefs, matrix ) ; } ) ;
      } // prepare_solve()

      /*void prepare_solve( basis4CORK::basis4CORK<typename CORK::deref_type<Basis>::type> const& basis4CORK ) {
        glas2::vector< value_type > coefs( num_terms() ) ;
        basis4CORK.evaluate( coefs.pass_ref() ) ;
        solver_.prepare_solve( coefs ) ;
      } // prepare_solve()*/

      template <typename V>
      void solve( V v ) const {
        static_assert( std::is_same< typename std::decay<V>::type::value_type, value_type >::value, "V must have the same value_type as solver" ) ;
        solver_.solve( v ) ;
      } // solve()

    public:
      template <typename V, typename Options>
      void solve( value_type& shift, size_type clone, V v, bool is_new_shift, Options const& options ) const {
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

    public:
      typedef linear_solver<std::reference_wrapper<sequence_of_functions_type const>, std::reference_wrapper<coefficient_matrices_type const>, linear_solver_type > clone_type ;
      clone_type clone() const {
        return clone_type( std::cref(basis_), std::cref(coefficient_matrices_), CORK::linear_solver::make_linear_solver<value_type>( coefficient_matrices_ ) ) ;
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
    typedef typename CORK::linear_solver::linear_solver_traits< value_type, std::reference_wrapper< typename matrix_polynomial_type::coefficient_matrices_type const> >::type solver_type ;
    typedef matrix_valued_function::linear_solver< std::reference_wrapper< typename matrix_polynomial_type::basis_type const>
                                                 , std::reference_wrapper< typename matrix_polynomial_type::coefficient_matrices_type const>
                                                 , solver_type
                                                 > type ;

    static type apply( MatrixPolynomial const& mp ) {
      //std::cout << "Do not select this one" << std::endl ;
      //std::cout << typeid(CoefficientMatrices).name() << std::endl ;
      solver_type solver( CORK::linear_solver::make_linear_solver<typename solver_type::value_type>( std::cref( CORK::deref(mp).coefficient_matrices() ) ) ) ;
      return type( std::cref(CORK::deref(mp).basis()), std::cref(CORK::deref(mp).coefficient_matrices()), solver ) ;
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
