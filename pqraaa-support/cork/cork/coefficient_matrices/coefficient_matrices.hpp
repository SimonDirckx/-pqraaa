//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompany_glasing file 
//  LICENSE_1_0.txt)

#ifndef cork_coefficient_matrices_coefficient_matrices_hpp
#define cork_coefficient_matrices_coefficient_matrices_hpp

#include <cork/coefficient_matrices/any_glas.hpp>
#include <cassert>
#include <type_traits>

namespace CORK { namespace coefficient_matrices {

  template <typename Sequence, typename EnableIf=void>
  struct coefficient_matrices_traits
  {
    typedef any_glas< Sequence > type ;

    static type apply( Sequence const& seq ) {
      return type( seq ) ;
    }
  } ;

  template <typename Sequence>
  decltype (auto) make_coefficient_matrices( Sequence const& seq ) {
    return coefficient_matrices_traits< Sequence >::apply( seq ) ;
  }

  template <typename Sequence>
  decltype (auto) make_coefficient_matrices_lvalue_reference( Sequence const& seq ) {
    return coefficient_matrices_traits< Sequence const& >::apply( seq ) ;
  }

} } // namespace CORK::coefficient_matrices

#endif
