#include <glas3/array/type/empty_array.hpp>
#include <glas3/array/dense_array/container/dense_scalar.hpp>
#include <glas3/array/dense_array/container/dense_vector.hpp>
#include <glas3/array/dense_array/container/dense_matrix.hpp>
#include <glas3/array/dense_array/container/dense_array.hpp>

#include <glas3/array/dense_array/algorithm/ops.hpp>
#include <glas3/array/dense_array/algorithm/permute.hpp>
#include <glas3/array/dense_array/algorithm/ttt.hpp>
#include <glas3/array/algorithm/shape.hpp>

#include <iostream>

int main() {

 	std::cout << '\n';

 	std::ptrdiff_t i ;

 	std::cout << "shape operation" << '\n' ;
 	glas3::empty_array e ;
 	glas3::dense_scalar<int> s = 1 ;
 	glas3::dense_vector<int> v = {1, 2, 3, 4, 5, 6} ;
 	glas3::dense_matrix<int> m({1, 2, 3, 4, 5, 6}, {2, 3}) ;
 	glas3::dense_array<int> a({1, 2, 3, 4, 5, 6}, {2, 3, 1}) ;
    auto b = ( a * ( v + a ) - v ) / v ;
    auto p = glas3::permute( a, {2, 0, 1} ) ;
    auto t = glas3::ttt( v, a ) ;

    for ( i = 0; i < e.shape().size(); ++i ) {
    	if ( e.shape()[i] != glas3::shape(e)[i] ) return 1 ;
    }

    for ( i = 0; i < s.shape().size(); ++i ) {
    	if ( s.shape()[i] != glas3::shape(s)[i] ) return 1 ;
    }

    for ( i = 0; i < v.shape().size(); ++i ) {
    	if ( v.shape()[i] != glas3::shape(v)[i] ) return 1 ;
    }

    for ( i = 0; i < m.shape().size(); ++i ) {
    	if ( m.shape()[i] != glas3::shape(m)[i] ) return 1 ;
    }

    for ( i = 0; i < a.shape().size(); ++i ) {
    	if ( a.shape()[i] != glas3::shape(a)[i] ) return 1 ;
    }

    for ( i = 0; i < b.shape().size(); ++i ) {
    	if ( b.shape()[i] != glas3::shape(b)[i] ) return 1 ;
    }

    for ( i = 0; i < p.shape().size(); ++i ) {
    	if ( p.shape()[i] != glas3::shape(p)[i] ) return 1 ;
    }

    for ( i = 0; i < t.shape().size(); ++i ) {
    	if ( t.shape()[i] != glas3::shape(t)[i] ) return 1 ;
    }

    std::cout << '\n' ;

	return 0;

}
