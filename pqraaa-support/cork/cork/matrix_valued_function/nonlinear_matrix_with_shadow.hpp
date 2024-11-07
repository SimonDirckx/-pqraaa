//  (C) Copyright Karl Meerbergen 2022.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_matrix_valued_function_nonlinear_matrix_with_shadow_hpp
#define cork_matrix_valued_function_nonlinear_matrix_with_shadow_hpp

#include <cork/coefficient_matrices/coefficient_matrices.hpp>
#include <cork/coefficient_matrices/coefficient_matrices4cork.hpp>
#include <cork/linear_solver/linear_solver.hpp>
#include <cork/utility/ref.hpp>
#include <glas2/vector.hpp>
#include <cassert>
#include <type_traits>

namespace CORK { namespace matrix_valued_function {

  //!
  //! Basis: type of polynomial basis, can be evaluated in any value_type in principle, unless the basis does not support this
  //! CoefficientMatrices: matrix coefficients, can have any value_type
  //! T: is the type of the matrix polynomial. Usually this is determined by the type of the matrices, but sometimes, the basis may have a specific
  //!    type so that T should be different. This allows the user to set the type manually if it is different from the default.
  //!
  template <typename NonlinearMatrix, typename Shadow>
  class nonlinear_matrix_with_shadow
  {
    public:
      typedef typename CORK::deref_type< NonlinearMatrix >::type               nonlinear_matrix_type ;
      typedef typename CORK::deref_type< Shadow >::type                        shadow_type ;
      typedef typename nonlinear_matrix_type::basis_type                       basis_type ;
      typedef typename nonlinear_matrix_type::coefficient_matrices_type        coefficient_matrices_type ;

    public:
      template <typename ValueType>
      using value_type_for = typename nonlinear_matrix_type::template value_type_for<ValueType> ;

      // For testing the shift
      template <typename ValueType>
      using has_value_type_for = typename nonlinear_matrix_type::template has_value_type_for<ValueType> ;

    public:
      nonlinear_matrix_with_shadow( NonlinearMatrix polynomial, Shadow shadow )
      : matrix_( polynomial )
      , shadow_( shadow )
      {
        assert( CORK::deref(matrix_).size()==CORK::deref(shadow_).size() ) ;
      }

    public:
      auto size() const {
        return CORK::deref(matrix_).size() ;
      }

      auto num_terms() const {
        //std::cout << "matrix_matrix_with_shadow: num_terms() is called" << std::endl ;
        return CORK::deref(matrix_).num_terms() ;
      } // num_terms()

      basis_type const& basis() const { return CORK::deref(matrix_).basis(); }
      coefficient_matrices_type const& coefficient_matrices() const { return CORK::deref(matrix_).coefficient_matrices() ; }

    public:
      shadow_type const& shadow() const {return CORK::deref(shadow_) ;}

    public:
      template <typename T>
      using linear_solver_type = typename shadow_type::template linear_solver_type<T> ;

      template <typename T>
      auto linear_solver() const {
        return CORK::deref(shadow_).template linear_solver<T>() ;
      }

    public:
      template <typename Shift, typename X, typename W>
      void multiply_add( Shift const& shift, X const& x, W w ) const {
        CORK::deref(shadow_).multiply_add( shift, x, w ) ;
      } // multiply_add()

    private:
      NonlinearMatrix matrix_ ;
      Shadow           shadow_ ;
  } ; // class nonlinear_matrix_with_shadow

} } // namespace CORK::matrix_valued_function

#endif
