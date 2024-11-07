//  (C) Copyright Karl Meerbergen 2021.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_matrix_valued_function_use_shadowed_nonlinear_matrix_hpp
#define cork_matrix_valued_function_use_shadowed_nonlinear_matrix_hpp

#include <cork/matrix_valued_function/matrix_polynomial_with_shadow.hpp>
#include <cork/options/use_shadowed_nonlinear_matrix.hpp>
#include <cork/options/has_option.hpp>

namespace CORK { namespace matrix_valued_function {

  namespace detail {
    template <int Options>
    struct use_shadowed_nonlinear_matrix {} ;

    template <>
    struct use_shadowed_nonlinear_matrix< true > {
      template <typename NEP, typename Shadow>
      static auto apply( NEP const& nep, Shadow const& shadow ) {
        return make_matrix_polynomial_with_shadow( nep, shadow ) ;
      }
    } ;

    template <>
    struct use_shadowed_nonlinear_matrix< false > {
      template <typename NEP, typename Shadow>
      static auto apply( NEP const& nep, Shadow const& shadow ) {
        return nep ;
      }
    } ;
  } // namespace detail

  template <typename Options, typename NEP, typename Shadow>
  auto use_shadowed_nonlinear_matrix( Options const& options, NEP const& nep, Shadow const& shadow ) {
    return detail::use_shadowed_nonlinear_matrix< options::has_option<options::use_shadowed_nonlinear_matrix,Options>::value >::apply( nep, shadow ) ;
  }

} } // namespace CORK::matrix_valued_function

#endif
