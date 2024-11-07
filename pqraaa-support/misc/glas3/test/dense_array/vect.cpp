#include <glas3/array/dense_array/container/dense_scalar.hpp>
#include <glas3/array/dense_array/container/dense_array.hpp>

#include <glas3/array/dense_array/algorithm/vect.hpp>

#include <iostream>

int main() {

 	std::cout << '\n';

 	std::ptrdiff_t i, j ;

	std::cout << "vect operation (shallow copy)" << '\n' ;
	std::cout << "auto vect_a1 = glas3::vect( glas3::dense_array<int>({1, 2, 3, 4, 5, 6}, {2, 3}) ) ;" << '\n' ;
 	glas3::dense_array<int> a1({1, 2, 3, 4, 5, 6}, {2, 3}) ;
 	auto vect_a1 = glas3::vect( a1 ) ;
    if ( vect_a1.size() != 6 || vect_a1.shape().size() != 1 || vect_a1.shape()[0] != 6 || vect_a1.shape().shape().size() != 1 ) return 1 ;
    a1[0] -= 1 ;
    for ( i = 0; i < 2; ++i ) {
    	for ( j = 0; j < 3; ++j ) {
    		if ( vect_a1[i + 2 * j] != a1[i + 2 * j] || vect_a1({i + 2 * j}) != a1({i, j}) ) return 1 ;
    	}
    }

	std::cout << "vect operation (shallow copy)" << '\n' ;
	std::cout << "auto vect_s = glas3::vect( glas3::dense_scalar<int>(7) ) ;" << '\n' ;
	glas3::dense_scalar<int> s = 7 ;
 	auto vect_s = glas3::vect( s ) ;
 	if ( vect_s.size() != 1 || vect_s.shape().size() != 1 || vect_s.shape()[0] != 1 || vect_s.shape().shape().size() != 1 ) return 1 ;
    s[0] -= 1 ;
    if ( vect_s[0] != s[0] ) return 1 ;

    std::cout << '\n' ;

	return 0;

}
