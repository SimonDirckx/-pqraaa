//  (C) Copyright Karl Meerbergen 2021.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_matrix_valued_function_select_nep_from_pair_hpp
#define cork_matrix_valued_function_select_nep_from_pair_hpp

#include <cork/matrix_valued_function/empty.hpp>
#include <cork/matrix_valued_function/nonlinear_matrix.hpp>
#include <cork/utility/ref.hpp>
#include <utility> // For pair

namespace CORK { namespace matrix_valued_function {

  template <typename T, typename AllowFail, typename Basis, typename CoefficientMatrices>
  typename std::enable_if< nonlinear_matrix<Basis,CoefficientMatrices>::template has_value_type_for<T>::value
                         , nonlinear_matrix<Basis,CoefficientMatrices> const&
                         >:: type select_nep_from_pair( nonlinear_matrix<Basis,CoefficientMatrices> const& nep ) {
    return nep ;
  }

  template <typename T, typename AllowFail, typename Basis, typename CoefficientMatrices>
  typename std::enable_if< !nonlinear_matrix<Basis,CoefficientMatrices>::template has_value_type_for<T>::value
                         , decltype(empty<T>().matrix_polynomial())
                         >:: type select_nep_from_pair( nonlinear_matrix<Basis,CoefficientMatrices> const& nep ) {
    static_assert( AllowFail::value, "CORK::value_type of shifts is not compatible with nonlinear matrix" ) ;
    return empty<T>().matrix_polynomial() ;
  }


  template <typename T, typename AllowFail, typename NEP1, typename NEP2, typename EnableIf=void>
  struct select_nep_from_pair_traits {
  } ;

  template <typename T, typename AllowFail, typename NEP1, typename NEP2>
  struct select_nep_from_pair_traits< T, AllowFail, NEP1, NEP2, typename std::enable_if< NEP1::template has_value_type_for<T>::value >::type > {
    static auto const& apply( NEP1 const& nep1, NEP2 const& nep2 ) {
      return nep1 ;
    }
  } ;

  template <typename T, typename AllowFail, typename NEP1, typename NEP2>
  struct select_nep_from_pair_traits< T, AllowFail, NEP1, NEP2, typename std::enable_if< !NEP1::template has_value_type_for<T>::value >::type > {
    static auto const& apply( NEP1 const& nep1, NEP2 const& nep2 ) {
      static_assert( NEP2::template has_value_type_for<T>::value || AllowFail::value, "CORK::value_type of shifts is not compatible with nonlinear matrix" ) ;
      return nep2 ;
    }
  } ;

/*  template <typename T, typename NEP1, typename NEP2>
  struct select_nep_from_pair_traits< T, NEP1, NEP2, typename std::enable_if< !NEP1::template has_value_type_for<T>::value
                                                                           && !NEP2::template has_value_type_for<T>::value
                                                                           >::type > {
    static auto const& apply( NEP1 const& nep1, NEP2 const& nep2 ) {
      return empty<T>().matrix_polynomial() ;
    }
  } ;
*/
  template <typename T, typename AllowFail, typename NEP1, typename NEP2>
  auto const& select_nep_from_pair( std::pair<NEP1,NEP2> const& pair ) {
    return select_nep_from_pair_traits<T,AllowFail,NEP1,NEP2>::apply( CORK::deref(pair.first), CORK::deref(pair.second) ) ;
  }


} } // namespace CORK::matrix_valued_function

#endif
