//  (C) Copyright Karl Meerbergen 2022.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_matrix_valued_function_nonlinear_matrix_hpp
#define cork_matrix_valued_function_nonlinear_matrix_hpp

#include <cork/basis/concept.hpp>
#include <cork/coefficient_matrices/coefficient_matrices.hpp>
#include <cork/coefficient_matrices/coefficient_matrices4cork.hpp>
#include <cork/linear_solver/linear_solver.hpp>
#include <cork/matrix_valued_function/linear_solver.hpp>
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
#ifdef CORK_USE_CONCEPTS
  template <Concept::Basis Basis, typename CoefficientMatrices>
#else
  template <typename Basis, typename CoefficientMatrices>
#endif
  class nonlinear_matrix
  {
    public:
      typedef Basis                                                  basis_data_type ;
      typedef CoefficientMatrices                                    coefficient_matrices_data_type ;
      typedef typename CORK::deref_type< Basis >::type               basis_type ;
      typedef typename CORK::deref_type< CoefficientMatrices >::type coefficient_matrices_type ;

    public:
      template <typename ValueType>
      using value_type_for = typename std::common_type< typename coefficient_matrices_type::template value_type_for<ValueType>
                                                      , typename basis_type::template value_type_for<ValueType>
                                                      , typename coefficient_matrices_type::template value_type_for< typename basis_type::template value_type_for<ValueType> >
                                                      >::type ;

      // For testing the shift
      template <typename ValueType>
      using has_value_type_for = std::bool_constant< basis_type::template has_value_type_for<ValueType>::value
                                                   && coefficient_matrices_type::template has_value_type_for< ValueType >::value
                                                   && coefficient_matrices_type::template has_value_type_for< typename basis_type::template value_type_for<ValueType> >::value
                                                   > ;

    public:
      nonlinear_matrix( Basis const& basis, CoefficientMatrices const& matrices )
      : basis_( basis )
      , matrices_( matrices )
      {
        assert( (long int)(CORK::deref(matrices_).num_matrices()) == (long int)(CORK::deref(basis_).num_terms()) ) ;
      } // nonlinear_matrix( functions, matrices )

    public:
      auto size() const {
        assert( CORK::deref(matrices_).num_rows()==CORK::deref(matrices_).num_columns() ) ;
        return CORK::deref(matrices_).num_rows() ;
      }

      auto num_terms() const {
        return CORK::deref(basis_).num_terms() ;
      } // num_terms()

      basis_type const& basis() const { return CORK::deref(basis_); }
      coefficient_matrices_type const& coefficient_matrices() const { return CORK::deref(matrices_) ; }

      Basis const& basis_data() const { return basis_; }
      CoefficientMatrices const& coefficient_matrices_data() const { return matrices_ ; }

    public:
      template <typename T>
      using linear_solver_type = typename matrix_valued_function::linear_solver_type< T, nonlinear_matrix >::type ;
      //using linear_solver_type = matrix_valued_function::linear_solver< basis_type const*, coefficient_matrices_type const*, typename CORK::linear_solver::linear_solver_traits<T, coefficient_matrices_type>::type > ;

      template <typename T>
      auto linear_solver() const {
        return CORK::matrix_valued_function::make_linear_solver<T>( *this ) ;
      }

    public:
      template <typename Shift, typename X, typename W>
      void multiply_add_slow( Shift const& shift, X const& x, W w ) const {
        glas2::vector< typename basis_type::template value_type_for<Shift> > coefs( basis_.num_terms() ) ;
        basis_.evaluate( shift, coefs ) ;

        // Can lose efficiency. To be fixed.
#ifndef NDEBUG
// can be that the same multiplication is done multiple times
std::cerr << "nonlinear_matrix::multiply_add_slow() : Can lose efficiency here" << std::endl ;
#endif
        glas2::vector<typename W::value_type> xx(w.size()) ;
        for (typename decltype(coefs)::size_type i=0; i<coefs.size(); ++i) {
  //        xx = coefs(i) * x ;
  //        matrices_.multiply_add( i, xx, w ) ;
          matrices_.multiply_add( i, coefs(i)*x, w ) ;
        }
      } // multiply_add_slow()

      template <typename Shift, typename X, typename W>
      void multiply_add( Shift const& shift, X const& x, W w ) const {
        glas2::vector< typename basis_type::template value_type_for<Shift> > coefs( CORK::deref(basis_).num_terms() ) ;
        basis_.evaluate( shift, coefs ) ;

        coefficient_matrices::coefficient_matrices4CORK<typename W::value_type, coefficient_matrices_type> coef_mat4CORK( CORK::deref(matrices_) ) ;
        coef_mat4CORK.initialize_schedule( CORK::deref(matrices_) ) ;

        glas2::vector<typename W::value_type> xx(w.size()) ;
        for (typename decltype(coefs)::size_type i=0; i<coefs.size(); ++i) {
          coef_mat4CORK.schedule( i, coefs(i)*x ) ;
        }
        coef_mat4CORK.apply_scheduled( CORK::deref(matrices_), w ) ;
      } // multiply_add()

    protected:
      Basis               basis_ ;
      CoefficientMatrices matrices_ ;
  } ; // class nonlinear_matrix

  template <typename Basis, typename CoefficientMatrices>
  std::ostream& operator<<( std::ostream& s, nonlinear_matrix<Basis,CoefficientMatrices> const& nep ) {
    s << "nonlinear_matrix{"<< nep.basis() << "," << nep.coefficient_matrices() <<"}" ;
    return s ;
  }

/*  template <typename Basis, typename CoefficientMatrices>
  decltype (auto) make_nonlinear_matrix( Basis const& poly, CoefficientMatrices const& mat ) {
    return nonlinear_matrix< decltype(remove_pass(poly)), decltype(remove_pass(mat)) >( remove_pass(poly), remove_pass(mat) ) ;
  } // make_nonlinear_matrix()
*/

} } // namespace CORK::matrix_valued_function

#endif
