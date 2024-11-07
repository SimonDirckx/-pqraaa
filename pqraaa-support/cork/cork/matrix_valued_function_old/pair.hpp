//  (C) Copyright Karl Meerbergen 2021.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_matrix_valued_function_pair_hpp
#define cork_matrix_valued_function_pair_hpp

#include <cork/utility/pass_value.hpp>
#include <cork/utility/pass_lvalue_reference.hpp>
#include <cork/utility/value_type.hpp>
#include <cork/utility/value_type_for.hpp>
#include <glas2/vector.hpp>
#include <cassert>
#include <type_traits>

namespace CORK { namespace matrix_valued_function {

  template <typename NEP_R, typename NEP_C>
  class pair
  {
    public:
      template <typename ValueType>
      using value_type_for = typename std::conditional< NEP_R::has_value_type_for<ValueType>::value
                                                      , typename NEP_R::value_type_for<ValueType>::type
                                                      , typename NEP_C::value_type_for<ValueType>::type
                                                      >::type ;

      // For testing the shift
      template <typename ValueType>
      using has_value_type_for = std::bool_constant< NEP_R::has_value_type_for<ValueType>::value || NEP_C::has_value_type_for<ValueType>::value > ;

    public:
      pair( NEP_R const& nep_r, NEP_C const& nep_c )
      : nep_r_( nep_r )
      , nep_c_( nep_c )
      {}

    public:
      auto size() const {
        assert( neps_c_.size()==neps_r_.size() ) ;
        return neps_c_.size() ;
      }

      auto num_terms() const {
        assert( neps_c_.num_terms()==neps_r_.num_terms() ) ;
        return nep_c_.num_terms() ;
      } // num_terms()

      basis_type const& basis() const { return nep_c_.basis() ; }
      auto coefficient_matrices() const {
        return coefficient_matrices:matrices_by_functions_pair< typename NEP_C::coefficient_matrices_type, typename NEP_R::coefficient_matrices_type >(
                                                        nep_c_.coefficient_matrices(), nep_r_.coefficient_matrices() ) ;
      }

    private:
      template <typename Shift, typename X, typename W>
      void multiply_add( Shift const& shift, X const& x, W w, std::false_type ) const {
        nep_c_.multiply_add( shift, x, w ) ;
      }

      template <typename Shift, typename X, typename W>
      void multiply_add( Shift const& shift, X const& x, W w, std::true_type ) const {
        nep_r_.multiply_add( shift, x, w ) ;
      }

    public:
      template <typename Shift, typename X, typename W>
      void multiply_add( Shift const& shift, X const& x, W w ) const {
        multiply_add( shift, x, w, typename NEP_R::has_value_type_for<Shift>() ) ;
      } // multiply_add()

    private:
      NEP_R const&        nep_r_ ;
      NEP_C const&        nep_c_ ;
  } ; // class pair


} } // namespace CORK::matrix_valued_function

#endif
