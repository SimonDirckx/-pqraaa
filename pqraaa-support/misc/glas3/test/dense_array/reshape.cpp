#include <glas3/array/dense_array/container/dense_scalar.hpp>
#include <glas3/array/dense_array/container/dense_array.hpp>

#include <glas3/array/dense_array/algorithm/reshape.hpp>

#include <iostream>

int main() {

 	std::cout << '\n';

 	std::ptrdiff_t i, j, k ;

	std::cout << "reshape operation (shallow copy)" << '\n' ;
	std::cout << "auto reshape_a1 = glas3::reshape( glas3::dense_array<int>({1, 2, 3, 4, 5, 6}, {2, 3}), {1, 6, 1} ) ;" << '\n' ;
	std::cout << "auto reshape_a2 = glas3::reshape( glas3::dense_array<int>({1, 2, 3, 4, 5, 6}, {2, 3}), {6} ) ;" << '\n' ;
 	glas3::dense_array<int> a({1, 2, 3, 4, 5, 6}, {2, 3}) ;
 	auto reshape_a1 = glas3::reshape( a, {1, 3, 2} ) ;
 	auto reshape_a2 = glas3::reshape( a, {6} ) ;
    if ( reshape_a1.size() != 6 || reshape_a1.shape().size() != 3 || reshape_a1.shape()[0] != 1 || reshape_a1.shape()[1] != 3 || reshape_a1.shape()[2] != 2 || reshape_a1.shape().shape().size() != 1 ) return 1 ;
    if ( reshape_a2.size() != 6 || reshape_a2.shape().size() != 1 || reshape_a2.shape()[0] != 6 || reshape_a2.shape().shape().size() != 1 ) return 1 ;
    a[0] -= 1 ;
    for ( i = 0; i < 1; ++i ) {
    	for ( j = 0; j < 3; ++j ) {
    		for ( k = 0; k < 2; ++k ) {
    			if ( reshape_a1[i + 1 * j + 1 * 3 * k] != a[i + 1 * j + 1 * 3 * k] || reshape_a1({i, j, k}) != a[i + 1 * j + 1 * 3 * k] ) return 1 ;
    		}
    	}
    }
    for ( i = 0; i < 6; ++i ) {
		if ( reshape_a2[i] != a[i] || reshape_a2({i}) != a[i] ) return 1 ;
    }

	std::cout << "reshape operation (shallow copy)" << '\n' ;
	std::cout << "auto reshape_s = glas3::reshape( glas3::dense_scalar<int>(7), {1, 1, 1} ) ;" << '\n' ;
	glas3::dense_scalar<int> s = 7 ;
 	auto reshape_s = glas3::reshape( s, {1, 1, 1} ) ;
 	if ( reshape_s.size() != 1 || reshape_s.shape().size() != 3 || reshape_s.shape()[0] != 1 || reshape_s.shape()[1] != 1 || reshape_s.shape()[2] != 1 || reshape_s.shape().shape().size() != 1 ) return 1 ;
    s[0] -= 1 ;
    for ( i = 0; i < 1; ++i ) {
    	for ( j = 0; j < 1; ++j ) {
    		for ( k = 0; k < 1; ++k ) {
    			if ( reshape_s[i + 1 * j + 1 * 1 * k] != s[i + 1 * j + 1 * 1 * k] || reshape_s({i, j, k}) != s[i + 1 * j + 1 * 1 * k] ) return 1 ;
    		}
    	}
    }

    std::cout << '\n' ;

	return 0;

}
