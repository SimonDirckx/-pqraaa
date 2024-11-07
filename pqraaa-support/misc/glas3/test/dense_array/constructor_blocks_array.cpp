#include <glas3/array/dense_array/container/dense_array.hpp>
#include <glas3/array/dense_array/type/blocks_array.hpp>

#include <glas3/array/dense_array/algorithm/iostream.hpp>

#include <iostream>
#include <vector>

int main() {

 	std::cout << '\n';

 	std::ptrdiff_t i, j ;

    glas3::dense_array<double> a ( {0, 1, 2, 3, 4, 5}, {2, 3} ) ;
    glas3::dense_array<double> b ( {7, 8, 9, 10, 11, 12}, {3, 2} ) ;
    glas3::dense_array<double> c ( {13, 14, 15, 16}, {1, 4} ) ;
    glas3::blocks_array<double> diag( {a, b, c}, {1, 0} ) ;
    if ( diag.size() != 6 * 9 || diag.ndof() != 6 + 6 + 4 || diag.shape().size() != 2 || diag.shape()[0] != 6 || diag.shape()[1] != 9 || diag.shape().shape().size() != 1 ) return 1 ;
    a[0] -= 7 ;
    b[3] -= 28 ;
    c[1] += 28 ;
    glas3::dense_array<double> d ( 0, {6, 9} ) ;
    d({{0, 1}, {0, 1, 2}}) = a ;
    d({{2, 3, 4}, {3, 4}}) = b ;
    d({5, {5, 6, 7, 8}}) = c ;
	for ( j = 0; j < 9; ++j ) {
		for ( i = 0; i < 6; ++i ) {
    		if ( diag[i + 6 * j] != d[i + 6 * j] || diag({i, j}) != d({i, j}) ) return 1 ;
    	}
    }

	std::cout << "glas3::blocks_array<double> diag( {a, b, c}, {1, 0} ): " << diag << '\n' ;

	glas3::dense_array<double> e ( {0, 1, 2, 3, 4, 5}, {2, 3} ) ;
	glas3::dense_array<double> f ( {7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18}, {4, 3} ) ;
	glas3::dense_array<double> g ( {13, 14, 15}, {1, 3} ) ;
    glas3::blocks_array<double, 2> cat( {e, f, g}, {0} ) ;
    if ( cat.size() != 7 * 3 || cat.ndof() != 6 + 12 + 3 || cat.shape().size() != 2 || cat.shape()[0] != 7 || cat.shape()[1] != 3 || cat.shape().shape().size() != 1 ) return 1 ;
    e[0] -= 7 ;
    f[3] -= 28 ;
    g[1] += 28 ;
    glas3::dense_array<double> h ( 0, {7, 3} ) ;
    h({{0, 1}, {0, 1, 2}}) = e ;
    h({{2, 3, 4, 5}, {0, 1, 2}}) = f ;
    h({6, {0, 1, 2}}) = g ;
	for ( j = 0; j < 3; ++j ) {
		for ( i = 0; i < 7; ++i ) {
    		if ( cat[i + 7 * j] != h[i + 7 * j] || cat({i, j}) != h({i, j}) ) return 1 ;
    	}
    }

	std::cout << "glas3::blocks_array<double, 2> cat( {e, f, g}, {0} ): " << cat << '\n' ;

	std::cout << '\n' ;

	return 0;
}
