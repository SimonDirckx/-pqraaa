//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_matrix_algorithm_inplace_multiply_hpp
#define glas2_matrix_algorithm_inplace_multiply_hpp

#include <glas2/vector/type/range.hpp>
#include <glas2/vector/algorithm/ops.hpp>
#include <glas2/vector/algorithm/ops_assign.hpp>
#include <glas2/matrix/algorithm/upper.hpp>
#include <glas2/matrix/algorithm/lower.hpp>
#include <cassert>

namespace glas2 {

  template <typename M, typename V>
  typename std::enable_if< is<Vector,V>::value >::type inplace_multiply( upper_view<M> const& matrix, V v ) {
    assert( matrix.num_rows()==matrix.num_columns() ) ;
    assert( v.size()==matrix.num_columns() ) ;

    auto m = matrix.matrix() ;

    for ( int i=0; i<m.num_columns(); ++i ) {
      v(i) = inner_prod( m( i, range_from_end(i,0) ), v(range_from_end(i,0) ) ) ;
    }
  }

  template <typename M, typename V>
  typename std::enable_if< is<Vector,V>::value >::type inplace_multiply( lower_view<M> const& matrix, V v ) {
    assert( matrix.num_rows()==matrix.num_columns() ) ;
    assert( v.size()==matrix.num_columns() ) ;

    auto m = matrix.matrix() ;

    for ( int i=m.num_columns()-1; i>=0; --i ) {
      v(i) = inner_prod( v(range(0,i+1)), m( i, range(0,i+1) ) ) ;
    }
  }


  template <typename M, typename V>
  typename std::enable_if< is<DenseMatrix,V>::value >::type inplace_multiply( upper_view<M> const& matrix, V v ) {
    assert( matrix.num_rows()==matrix.num_columns() ) ;
    assert( v.num_rows()==matrix.num_columns() ) ;

    auto m = matrix.matrix() ;

    for ( int i=0; i<m.num_columns(); ++i ) {
      v(i,all()) = multiply( transpose(v(range_from_end(i,0),all())), m( i, range_from_end(i,0) ) ) ;
    }
  }

  template <typename M, typename V>
  typename std::enable_if< is<DenseMatrix,V>::value >::type inplace_multiply( V v, lower_view<M> const& matrix ) {
    assert( matrix.num_rows()==matrix.num_columns() ) ;
    assert( v.num_columns()==matrix.num_columns() ) ;

    auto m = matrix.matrix() ;

    for ( int i=0; i<m.num_columns(); ++i ) {
      v(all(),i) = multiply( v(all(),range_from_end(i,0)), m( range_from_end(i,0), i ) ) ;
    }
  }

  template <typename M, typename V>
  typename std::enable_if< is<DenseMatrix,V>::value >::type inplace_multiply( lower_view<M> const& matrix, V v ) {
    assert( matrix.num_rows()==matrix.num_columns() ) ;
    assert( v.num_rows()==matrix.num_columns() ) ;

    auto m = matrix.matrix() ;

    for ( int i=m.num_columns()-1; i>=0; --i ) {
      v(i,all()) = multiply( transpose(v(range(0,i+1),all())), m( i, range(0,i+1) ) ) ;
    }
  }

  template <typename M, typename V>
  typename std::enable_if< is<DenseMatrix,V>::value >::type inplace_multiply( V v, upper_view<M> const& matrix ) {
    assert( matrix.num_rows()==matrix.num_columns() ) ;
    assert( v.num_columns()==matrix.num_columns() ) ;

    auto m = matrix.matrix() ;

    for ( int i=m.num_columns()-1; i>=0; --i ) {
      v(all(),i) = multiply( v(all(),range(0,i+1)), m( range(0,i+1), i ) ) ;
    }
  }

} // namespace glas2

#endif
