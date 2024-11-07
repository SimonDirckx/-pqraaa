//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_matrix_algorithm_solve_hpp
#define glas2_matrix_algorithm_solve_hpp

#include <glas2/vector/type/range.hpp>
#include <glas2/vector/algorithm/ops.hpp>
#include <glas2/vector/algorithm/ops_assign.hpp>
#include <glas2/matrix/algorithm/upper.hpp>
#include <glas2/matrix/algorithm/lower.hpp>
#include <cassert>

namespace glas2 {

  template <typename M, typename V>
  void solve( upper_view<M> const& matrix, V v ) {
    assert( matrix.num_rows()==matrix.num_columns() ) ;
    assert( v.size()==matrix.num_columns() ) ;

    auto m = matrix.matrix() ;

    for ( int i=m.num_columns()-1; i>=0; --i ) {
      v(i) /= m(i,i) ;
      v( range(0,i) ) -= v(i) * matrix_selection< decltype(m), range, int >::apply( m, range(0,i), i ) ;
    }
  }

  template <typename M, typename V>
  void solve( lower_view<M> const& matrix, V v ) {
    assert( matrix.num_rows()==matrix.num_columns() ) ;
    assert( v.size()==matrix.num_columns() ) ;

    auto m = matrix.matrix() ;

    for ( int i=0; i<m.num_columns(); ++i ) {
      v(i) /= m(i,i) ;
      v( range(i+1,v.size()) ) -= v(i) * matrix_selection< decltype(m), range, int >::apply( m, range(i+1,v.size()), i ) ;
    }
  }

} // namespace glas2

#endif
