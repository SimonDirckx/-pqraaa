#include <glas3/array/dense_array/container/dense_scalar.hpp>
#include <glas3/array/dense_array/container/dense_array.hpp>

#include <glas3/array/dense_array/algorithm/squeeze.hpp>

#include <iostream>

int main() {

 	std::cout << '\n';

 	std::ptrdiff_t i, j ;

	std::cout << "squeeze operation (shallow copy)" << '\n' ;
	std::cout << "auto squeeze_a = glas3::squeeze( glas3::dense_array<int>({1, 2, 3, 4, 5, 6}, {1, 2, 1, 3, 1, 1}) ) ;" << '\n' ;
 	glas3::dense_array<int> a({1, 2, 3, 4, 5, 6}, {1, 2, 1, 3, 1, 1}) ;
 	auto squeeze_a = glas3::squeeze( a ) ;
    if ( squeeze_a.size() != 6 || squeeze_a.shape().size() != 2 || squeeze_a.shape()[0] != 2 || squeeze_a.shape()[1] != 3 || squeeze_a.shape().shape().size() != 1 ) return 1 ;
    a[0] -= 1 ;
    for ( i = 0; i < 2; ++i ) {
    	for ( j = 0; j < 3; ++j ) {
    		if ( squeeze_a[i + 2 * j] != a[i + 2 * j] || squeeze_a({i, j}) != a({0, i, 0, j, 0, 0}) ) return 1 ;
    	}
    }

	std::cout << "squeeze operation (shallow copy)" << '\n' ;
	std::cout << "auto squeeze_s = glas3::squeeze( glas3::dense_scalar<int>(7) ) ;" << '\n' ;
	glas3::dense_scalar<int> s = 7 ;
 	auto squeeze_s = glas3::squeeze( s ) ;
 	if ( squeeze_s.size() != 1 || squeeze_s.shape().size() != 0 || squeeze_s.shape().shape().size() != 1 ) return 1 ;
    s[0] -= 1 ;
    if ( squeeze_s[0] != s[0] ) return 1 ;

    std::cout << '\n' ;

	return 0;

}
