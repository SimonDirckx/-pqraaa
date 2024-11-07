//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_matrix_valued_function_matrix_polynomial_with_shadow_hpp
#define cork_matrix_valued_function_matrix_polynomial_with_shadow_hpp

#include <cork/matrix_valued_function/linear_solver.hpp>
#include <cork/utility/pass_value.hpp>
#include <cork/utility/pass_lvalue_reference.hpp>
#include <cork/utility/value_type.hpp>
#include <cork/utility/value_type_for.hpp>
#include <glas2/vector.hpp>
#include <cassert>
#include <type_traits>

namespace CORK { namespace matrix_valued_function {

  template <typename MatrixPolynomial, typename Shadow>
  class matrix_polynomial_with_shadow
  {
    public:
      typedef typename std::decay< Shadow >::type                        shadow_type ;
      typedef typename std::decay< MatrixPolynomial >::type              matrix_polynomial_type ;
      typedef typename matrix_polynomial_type::basis_type                basis_type ;
      typedef typename matrix_polynomial_type::coefficient_matrices_type coefficient_matrices_type ;

    public:
      template <typename ValueType>
      using value_type_for = typename shadow_type::template value_type_for<ValueType> ;

      //typedef typename shadow_type::native_value_type native_value_type ;

    public:
      matrix_polynomial_with_shadow( MatrixPolynomial polynomial, Shadow shadow )
      : polynomial_( polynomial )
      , shadow_( shadow )
      {
        assert( polynomial_.size()==shadow_.size() ) ;
      }

    public:
      auto size() const {
        return polynomial_.size() ;
      }

      auto num_terms() const {
        //std::cout << "matrix_polynomial_with_shadow: num_terms() is called" << std::endl ;
        return polynomial_.num_terms() ;
      } // num_terms()

      basis_type const& basis() const { return polynomial_.basis(); }
      coefficient_matrices_type const& coefficient_matrices() const { return polynomial_.coefficient_matrices() ; }

    public:
      shadow_type const& shadow() const {return shadow_ ;}

    public:
      template <typename Shift, typename X, typename W>
      void multiply_add( Shift const& shift, X const& x, W w ) const {
        shadow_.multiply_add( shift, x, w ) ;
      } // multiply_add()

    private:
      MatrixPolynomial polynomial_ ;
      Shadow           shadow_ ;
  } ; // class matrix_polynomial

  template <typename MatrixPolynomial, typename Shadow>
  decltype (auto) make_matrix_polynomial_with_shadow_lvalue( MatrixPolynomial const& poly, Shadow const& shadow ) {
    return matrix_polynomial_with_shadow< MatrixPolynomial const&, Shadow const&>( poly, shadow ) ;
  }

  template <typename MatrixPolynomial, typename Shadow>
  decltype (auto) make_matrix_polynomial_with_shadow( MatrixPolynomial const& poly, Shadow const& shadow ) {
    return matrix_polynomial_with_shadow< MatrixPolynomial, Shadow >( poly, shadow ) ;
  }


  template <typename ValueType, typename MatrixPolynomial, typename Shadow>
  struct linear_solver_struct< ValueType, matrix_polynomial_with_shadow< MatrixPolynomial, Shadow> >
  : linear_solver_struct< ValueType, typename std::decay<Shadow>::type >
  {
    static auto apply( matrix_polynomial_with_shadow< MatrixPolynomial, Shadow> const& mp ) {
      return make_linear_solver<ValueType>( mp.shadow() ) ;
    }
  } ; // struct linear_solver_struct


} } // namespace CORK::matrix_valued_function


namespace CORK {
      template <typename ValueType, typename MatrixPolynomial, typename Shadow>
      struct value_type_for< ValueType, CORK::matrix_valued_function::matrix_polynomial_with_shadow< MatrixPolynomial, Shadow > >
      : value_type_for< ValueType, Shadow >
      {} ;

      template <typename MatrixPolynomial, typename Shadow>
      struct value_type< CORK::matrix_valued_function::matrix_polynomial_with_shadow< MatrixPolynomial, Shadow > >
      : value_type< Shadow >
      {} ;
} // namespace CORK

#endif
