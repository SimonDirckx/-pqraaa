//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2015.
//
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_krylov_rhs_iteration_hpp
#define cork_krylov_rhs_iteration_hpp

#include <cork/krylov/options_type.hpp>

namespace CORK {

  template <typename TwoLevel, typename Strategy, typename Op>
  void rks_steps( Op const& op, Shifts const& shifts, int& step_number, Continuation const& continuation, QQ Q, UU U, HH H, KK K ) {
    assert( step_number+shifts.size()<=K.num_columns() ) ;
    assert( step_number+shifts.size()<=H.num_columns() ) ;
    assert( step_number+shifts.size()<K.numrows() ) ;
    assert( step_number+shifts.size()<H.numrows() ) ;
    assert( step_number+shifts.size()<U.num_columns() ) ;

    for (size_type i=0; i<shifts.size(); ++i) {
    }
  } // rks_step()
   
} // namespace CORK


#endif
