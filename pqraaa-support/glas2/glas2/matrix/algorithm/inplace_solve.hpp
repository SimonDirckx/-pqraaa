//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_matrix_algorithm_inplace_solve_hpp
#define glas2_matrix_algorithm_inplace_solve_hpp

#include <glas2/vector/type/range.hpp>
#include <glas2/vector/algorithm/ops.hpp>
#include <glas2/vector/algorithm/ops_assign.hpp>
#include <glas2/matrix/algorithm/upper.hpp>
#include <glas2/matrix/algorithm/lower.hpp>
#include <cassert>

namespace glas2 {

  template <typename M, typename V>
  typename std::enable_if< is<Vector,V>::value >::type inplace_solve( upper_view<M> const& matrix, V v ) {
    assert( matrix.num_rows()==matrix.num_columns() ) ;
    assert( v.size()==matrix.num_columns() ) ;

    auto m = matrix.matrix() ;

    for ( int i=m.num_columns()-1; i>=0; --i ) {
      v(i) /= m(i,i) ;
      v( range(0,i) ) -= v(i) * matrix_selection< decltype(m), range, int >::apply( m, range(0,i), i ) ;
    }
  }

  template <typename M, typename V>
  typename std::enable_if< is<Vector,V>::value >::type inplace_solve( lower_view<M> const& matrix, V v ) {
    assert( matrix.num_rows()==matrix.num_columns() ) ;
    assert( v.size()==matrix.num_columns() ) ;

    auto m = matrix.matrix() ;

    for ( int i=0; i<m.num_columns(); ++i ) {
      v(i) /= m(i,i) ;
      v( range(i+1,v.size()) ) -= v(i) * matrix_selection< decltype(m), range, int >::apply( m, range(i+1,v.size()), i ) ;
    }
  }


  template <typename M, typename V>
  typename std::enable_if< is<DenseMatrix,V>::value >::type inplace_solve( upper_view<M> const& matrix, V v ) {
    assert( matrix.num_rows()==matrix.num_columns() ) ;
    assert( v.num_rows()==matrix.num_columns() ) ;

    auto m = matrix.matrix() ;

    for ( int i=m.num_columns()-1; i>=0; --i ) {
      v(i,glas2::all()) /= m(i,i) ;
      v( range(0,i),all() ) -= outer_prod( matrix_selection< decltype(m), range, int >::apply( m, range(0,i), i ), v( i,all() ) ) ;
    }
  }

  template <typename M, typename V>
  typename std::enable_if< is<DenseMatrix,V>::value >::type inplace_solve( V v, lower_view<M> const& matrix ) {
    assert( matrix.num_rows()==matrix.num_columns() ) ;
    assert( v.num_columns()==matrix.num_columns() ) ;

    auto m = matrix.matrix() ;

    for ( int i=m.num_columns()-1; i>=0; --i ) {
      v(all(),i) /= m(i,i) ;
      v(all(), range(0,i)) -= outer_prod( v( all(),i ), m( i, range(0,i) ) ) ;
    }
  }

  template <typename M, typename V>
  typename std::enable_if< is<DenseMatrix,V>::value >::type inplace_solve( lower_view<M> const& matrix, V v ) {
    assert( matrix.num_rows()==matrix.num_columns() ) ;
    assert( v.num_rows()==matrix.num_columns() ) ;

    auto m = matrix.matrix() ;

    for ( int i=0; i<m.num_columns(); ++i ) {
      v(i,glas2::all()) /= m(i,i) ;
      if (i+1<m.num_columns()) v( range_from_end(i+1,0),all() ) -= outer_prod( matrix_selection< decltype(m), range_from_end, int >::apply( m, range_from_end(i+1,0), i ), v(i,all()) ) ;
    }
  }

  template <typename M, typename V>
  typename std::enable_if< is<DenseMatrix,V>::value >::type inplace_solve( V v, upper_view<M> const& matrix ) {
    assert( matrix.num_rows()==matrix.num_columns() ) ;
    assert( v.num_columns()==matrix.num_columns() ) ;

    auto m = matrix.matrix() ;

    for ( int i=0; i<m.num_columns(); ++i ) {
      v(glas2::all(),i) /= m(i,i) ;
      if (i+1<m.num_columns()) v(all(), range_from_end(i+1,0)) -= outer_prod( v(all(),i), m( i, range_from_end(i+1,0) ) ) ;
      //if (i+1<m.num_columns()) v(all(), range_from_end(i+1,0)) -= outer_prod( v(all(),i), matrix_selection< decltype(m), int, range_from_end >::apply( m, i, range_from_end(i+1,0) ) ) ;
    }
  }

} // namespace glas2

#endif
