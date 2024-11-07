#include <glas3/array/dense_array/container/dense_scalar.hpp>
#include <glas3/array/dense_array/container/dense_vector.hpp>
#include <glas3/array/dense_array/container/dense_array.hpp>

#include <glas3/array/dense_array/algorithm/tile.hpp>

#include <iostream>

int main() {

 	std::cout << '\n';

 	std::ptrdiff_t i, j ;

	std::cout << "tile operation (shallow copy)" << '\n' ;
	std::cout << "auto tile_a = glas3::tile( a({5, 7, 11, 13}, {2, 2}), {2, 3}) ) ;" << '\n' ;
 	glas3::dense_array<int> a({5, 7, 11, 13}, {2, 2}) ;
    auto tile_a = glas3::tile( a, {2, 3} ) ;
    if ( tile_a.size() != 4*2*3 || tile_a.shape().size() != 2 || tile_a.shape()[0] != 2*2 || tile_a.shape()[1] != 2*3 || tile_a.shape().shape().size() != 1 ) return 1 ;
    a[0] -= 1 ;
    for ( i = 0; i < 2*2; ++i ) {
    	for ( j = 0; j < 2*3; ++j ) {
    		if ( tile_a[i + 2*2 * j] != a[(i%2) + 2 * (j%2)] || tile_a({i, j}) != a({i%2, j%2}) ) return 1 ;
    	}
    }
    for ( i = 2*2; i > 0 ; --i ) {
    	for ( j = 2*3; j > 0; --j ) {
    		if ( tile_a[(i-1) + 2*2 * (j-1)] != a[((i-1)%2) + 2 * ((j-1)%2)] || tile_a({i-1, j-1}) != a({(i-1)%2, (j-1)%2}) ) return 1 ;
    	}
    }

	std::cout << '\n' ;

	return 0;

}
