//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_bsp_backend_default_matrix_matrix_multiply_hpp
#define glas2_bsp_backend_default_matrix_matrix_multiply_hpp

#include <glas2/bsp/backend/default/default_backend.hpp>
#include <glas2/matrix/algorithm/fill.hpp>
#include <glas2/bsp/matrix/concept/bsp_matrix.hpp>
#include <type_traits>
#include <cassert>

namespace glas2 { namespace bsp {

  template <typename To, typename M1, typename M2>
  To assign( current_backend(), To to, glas2::bsp::matrix_multiply<M1,M2> const& that ) {
    fill( to, 0.0 ) ;
    plus_assign( current_backend(), to, that ) ;
    return to ;
  }

  template <typename To, typename M1, typename M2>
  To plus_assign( current_backend(), To to, glas2::bsp::matrix_multiply<M1,M2> const& that ) {
    return to ;
  }

} } // namespace glas2::bsp

#endif
