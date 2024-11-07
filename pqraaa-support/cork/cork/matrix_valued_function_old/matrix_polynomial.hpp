//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_matrix_valued_function_matrix_polynomial_hpp
#define cork_matrix_valued_function_matrix_polynomial_hpp

#include <cork/coefficient_matrices/coefficient_matrices.hpp>
#include <cork/coefficient_matrices/coefficient_matrices4cork.hpp>
#include <cork/utility/pass_value.hpp>
#include <cork/utility/pass_lvalue_reference.hpp>
#include <cork/utility/value_type.hpp>
#include <cork/utility/value_type_for.hpp>
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
  template <typename Basis, typename CoefficientMatrices>
  class matrix_polynomial
  {
    public:
      typedef typename std::decay< Basis >::type               basis_type ;
      typedef typename std::decay< CoefficientMatrices >::type coefficient_matrices_type ;

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
      matrix_polynomial( Basis basis, CoefficientMatrices matrices )
      : basis_( basis )
      , matrices_( matrices )
      {
        assert( (long int)(matrices_.num_matrices()) == (long int)(basis_.num_terms()) ) ;
      } // matrix_polynomial( functions, matrices )

    public:
      auto size() const {
        assert( matrices_.num_rows()==matrices_.num_columns() ) ;
        return matrices_.num_rows() ;
      }

      auto num_terms() const {
        return basis_.num_terms() ;
      } // num_terms()

      basis_type const& basis() const { return basis_; }
      coefficient_matrices_type const& coefficient_matrices() const { return matrices_ ; }

    public:
      template <typename Shift, typename X, typename W>
      void multiply_add_slow( Shift const& shift, X const& x, W w ) const {
        glas2::vector< typename basis_type::template value_type_for<Shift> > coefs( basis_.num_terms() ) ;
        basis_.evaluate( shift, coefs ) ;

        // Can lose efficiency. To be fixed.
#ifndef NDEBUG
// can be that the same multiplication is done multiple times
std::cerr << "matrix_polynomial::multiply_add_slow() : Can lose efficiency here" << std::endl ;
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
        glas2::vector< typename basis_type::template value_type_for<Shift> > coefs( basis_.num_terms() ) ;
        basis_.evaluate( shift, coefs ) ;

        coefficient_matrices::coefficient_matrices4CORK<typename W::value_type, CoefficientMatrices> coef_mat4CORK( matrices_ ) ;
        coef_mat4CORK.initialize_schedule( matrices_ ) ;

        glas2::vector<typename W::value_type> xx(w.size()) ;
        for (typename decltype(coefs)::size_type i=0; i<coefs.size(); ++i) {
          coef_mat4CORK.schedule( i, coefs(i)*x ) ;
        }
        coef_mat4CORK.apply_scheduled( matrices_, w ) ;
      } // multiply_add()

    private:
      Basis               basis_ ;
      CoefficientMatrices matrices_ ;
  } ; // class matrix_polynomial

  template <typename Basis, typename CoefficientMatrices>
  decltype (auto) make_matrix_polynomial( Basis const& poly, CoefficientMatrices const& mat ) {
    return matrix_polynomial< decltype(remove_pass(poly)), decltype(remove_pass(mat)) >( remove_pass(poly), remove_pass(mat) ) ;
  } // make_matrix_polynomial()


} } // namespace CORK::matrix_valued_function

#endif
